/*****************************************************************************\
 *  agent.c - parallel background communication functions. This is where  
 *	logic could be placed for broadcast communications.
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-5 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>, et. al.
 *  Derived from pdsh written by Jim Garlick <garlick1@llnl.gov>
 *  UCRL-CODE-217948.
 *  
 *  This file is part of SLURM, a resource management program.
 *  For details, see <http://www.llnl.gov/linux/slurm/>.
 *  
 *  SLURM is free software; you can redistribute it and/or modify it under
 *  the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *  
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
 *****************************************************************************
 *  Theory of operation:
 *
 *  The functions below permit slurm to initiate parallel tasks as a 
 *  detached thread and let the functions below make sure the work happens. 
 *  For example, when a job's time limit is to be changed slurmctld needs 
 *  to notify the slurmd on every node to which the job was allocated.  
 *  We don't want to hang slurmctld's primary function (the job update RPC)  
 *  to perform this work, so it just initiates an agent to perform the work.  
 *  The agent is passed all details required to perform the work, so it will 
 *  be possible to execute the agent as an pthread, process, or even a daemon 
 *  on some other computer.
 *
 *  The main agent thread creates a separate thread for each node to be
 *  communicated with up to AGENT_THREAD_COUNT. A special watchdog thread  
 *  sends SIGLARM to any threads that have been active (in DSH_ACTIVE state)  
 *  for more than COMMAND_TIMEOUT seconds. 
 *  The agent responds to slurmctld via a function call or an RPC as required.
 *  For example, informing slurmctld that some node is not responding.
 *
 *  All the state for each thread is maintained in thd_t struct, which is 
 *  used by the watchdog thread as well as the communication threads.
\*****************************************************************************/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <pthread.h>
#include <pwd.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdlib.h>

#include "src/common/list.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/node_select.h"
#include "src/common/xsignal.h"
#include "src/common/xassert.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_protocol_interface.h"
#include "src/common/forward.h"
#include "src/slurmctld/agent.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/ping_nodes.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/state_save.h"
#include "src/slurmctld/srun_comm.h"

#define MAX_RETRIES	10

typedef enum {
	DSH_NEW,        /* Request not yet started */
	DSH_ACTIVE,     /* Request in progress */
	DSH_DONE,       /* Request completed normally */
	DSH_NO_RESP,    /* Request timed out */
	DSH_FAILED      /* Request resulted in error */
} state_t;

typedef struct thd_complete {
	bool work_done; 	/* assume all threads complete */
	int fail_cnt;   	/* assume no threads failures */
	int no_resp_cnt;	/* assume all threads respond */
	int retry_cnt;  	/* assume no required retries */
	int max_delay;
	time_t now;
} thd_complete_t;

typedef struct thd {
	pthread_t thread;		/* thread ID */
	pthread_attr_t attr;		/* thread attributes */
	state_t state;			/* thread state */
	time_t start_time;		/* start time */
	time_t end_time;		/* end time or delta time 
					 * upon termination */
	struct sockaddr_in slurm_addr;	/* network address */
	forward_t forward;	        /* structure holding info for all
					   forwarding info
					*/	
	char node_name[MAX_SLURM_NAME];	/* node's name */
	List ret_list;
} thd_t;

typedef struct agent_info {
	pthread_mutex_t thread_mutex;	/* agent specific mutex */
	pthread_cond_t thread_cond;	/* agent specific condition */
	uint32_t thread_count;		/* number of threads records */
	uint32_t threads_active;	/* currently active threads */
	uint16_t retry;			/* if set, keep trying */
	thd_t *thread_struct;		/* thread structures */
	bool get_reply;			/* flag if reply expected */
	slurm_msg_type_t msg_type;	/* RPC to be issued */
	void **msg_args_pptr;		/* RPC data to be used */
} agent_info_t;

typedef struct task_info {
	pthread_mutex_t *thread_mutex_ptr; /* pointer to agent specific 
					    * mutex */
	pthread_cond_t *thread_cond_ptr;/* pointer to agent specific
					 * condition */
	uint32_t *threads_active_ptr;	/* currently active thread ptr */
	thd_t *thread_struct_ptr;	/* thread structures ptr */
	bool get_reply;			/* flag if reply expected */
	slurm_msg_type_t msg_type;	/* RPC to be issued */
	void *msg_args_ptr;		/* ptr to RPC data to be used */
} task_info_t;

typedef struct queued_request {
	agent_arg_t* agent_arg_ptr;	/* The queued request */
	time_t       last_attempt;	/* Time of last xmit attempt */
} queued_request_t;

typedef struct mail_info {
	char *user_name;
	char *message;
} mail_info_t;

static void _alarm_handler(int dummy);
static inline void _comm_err(char *node_name);
static void _list_delete_retry(void *retry_entry);
static agent_info_t *_make_agent_info(agent_arg_t *agent_arg_ptr);
static task_info_t *_make_task_data(agent_info_t *agent_info_ptr, int inx);
static void _notify_slurmctld_jobs(agent_info_t *agent_ptr);
static void _notify_slurmctld_nodes(agent_info_t *agent_ptr, 
		int no_resp_cnt, int retry_cnt);
static void _purge_agent_args(agent_arg_t *agent_arg_ptr);
static void _queue_agent_retry(agent_info_t * agent_info_ptr, int count);
static void _slurmctld_free_job_launch_msg(batch_job_launch_msg_t * msg);
static void _spawn_retry_agent(agent_arg_t * agent_arg_ptr);
static void *_thread_per_group_rpc(void *args);
static int   _valid_agent_arg(agent_arg_t *agent_arg_ptr);
static void *_wdog(void *args);

static mail_info_t *_mail_alloc(void);
static void  _mail_free(void *arg);
static void  _mail_proc(mail_info_t *mi);
static char *_mail_type_str(uint16_t mail_type);

static pthread_mutex_t retry_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mail_mutex  = PTHREAD_MUTEX_INITIALIZER;
static List retry_list = NULL;		/* agent_arg_t list for retry */
static List mail_list = NULL;		/* pending e-mail requests */

static pthread_mutex_t agent_cnt_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  agent_cnt_cond  = PTHREAD_COND_INITIALIZER;
static int agent_cnt = 0;

static bool run_scheduler = false;

/*
 * agent - party responsible for transmitting an common RPC in parallel 
 *	across a set of nodes. Use agent_queue_request() if immediate 
 *	execution is not essential.
 * IN pointer to agent_arg_t, which is xfree'd (including slurm_addr, 
 *	node_names and msg_args) upon completion if AGENT_IS_THREAD is set
 * RET always NULL (function format just for use as pthread)
 */
void *agent(void *args)
{
	int i, rc, retries = 0;
	pthread_attr_t attr_wdog;
	pthread_t thread_wdog;
	agent_arg_t *agent_arg_ptr = args;
	agent_info_t *agent_info_ptr = NULL;
	thd_t *thread_ptr;
	task_info_t *task_specific_ptr;

	//info("I am here and agent_cnt is %d of %d",agent_cnt, MAX_AGENT_CNT);
	slurm_mutex_lock(&agent_cnt_mutex);
	while (1) {
		if (agent_cnt < MAX_AGENT_CNT) {
			agent_cnt++;
			break;
		} else {	/* wait for state change and retry */
			pthread_cond_wait(&agent_cnt_cond, &agent_cnt_mutex);
		}
	}
	slurm_mutex_unlock(&agent_cnt_mutex);

	/* basic argument value tests */
	if (_valid_agent_arg(agent_arg_ptr))
		goto cleanup;

	xsignal(SIGALRM, _alarm_handler);

	/* initialize the agent data structures */
	agent_info_ptr = _make_agent_info(agent_arg_ptr);
	thread_ptr = agent_info_ptr->thread_struct;

	/* start the watchdog thread */
	slurm_attr_init(&attr_wdog);
	if (pthread_attr_setdetachstate
	    (&attr_wdog, PTHREAD_CREATE_JOINABLE))
		error("pthread_attr_setdetachstate error %m");
	while (pthread_create(&thread_wdog, &attr_wdog, _wdog,
				(void *) agent_info_ptr)) {
		error("pthread_create error %m");
		if (++retries > MAX_RETRIES)
			fatal("Can't create pthread");
		sleep(1);	/* sleep and again */
	}
	slurm_attr_destroy(&attr_wdog);
#if 	AGENT_THREAD_COUNT < 1
	fatal("AGENT_THREAD_COUNT value is invalid");
#endif
	debug2("got %d threads to send out",agent_info_ptr->thread_count);
	/* start all the other threads (up to AGENT_THREAD_COUNT active) */
	for (i = 0; i < agent_info_ptr->thread_count; i++) {

		/* wait until "room" for another thread */
		slurm_mutex_lock(&agent_info_ptr->thread_mutex);
		while (agent_info_ptr->threads_active >=
		       AGENT_THREAD_COUNT) {
			pthread_cond_wait(&agent_info_ptr->thread_cond,
					  &agent_info_ptr->thread_mutex);
		}

		/* create thread specific data, NOTE: freed from 
		 *      _thread_per_group_rpc() */
		task_specific_ptr = _make_task_data(agent_info_ptr, i);

		slurm_attr_init(&thread_ptr[i].attr);
		if (pthread_attr_setdetachstate(&thread_ptr[i].attr,
						PTHREAD_CREATE_DETACHED))
			error("pthread_attr_setdetachstate error %m");
		while ((rc = pthread_create(&thread_ptr[i].thread,
					    &thread_ptr[i].attr,
					    _thread_per_group_rpc,
					    (void *) task_specific_ptr))) {
			error("pthread_create error %m");
			if (agent_info_ptr->threads_active)
				pthread_cond_wait(&agent_info_ptr->
						  thread_cond,
						  &agent_info_ptr->
						  thread_mutex);
			else {
				slurm_mutex_unlock(&agent_info_ptr->
						     thread_mutex);
				sleep(1);
				slurm_mutex_lock(&agent_info_ptr->
						   thread_mutex);
			}
		}
		slurm_attr_destroy(&thread_ptr[i].attr);
		agent_info_ptr->threads_active++;
		slurm_mutex_unlock(&agent_info_ptr->thread_mutex);
	}
		
	/* wait for termination of remaining threads */
	pthread_join(thread_wdog, NULL);
	slurm_mutex_lock(&agent_info_ptr->thread_mutex);
	while (agent_info_ptr->threads_active != 0) {
		pthread_cond_wait(&agent_info_ptr->thread_cond,
				&agent_info_ptr->thread_mutex);
	}
	slurm_mutex_unlock(&agent_info_ptr->thread_mutex);

      cleanup:
#if AGENT_IS_THREAD
	_purge_agent_args(agent_arg_ptr);
#endif

	if (agent_info_ptr) {
		xfree(agent_info_ptr->thread_struct);
		xfree(agent_info_ptr);
	}
	slurm_mutex_lock(&agent_cnt_mutex);
	if (agent_cnt > 0)
		agent_cnt--;
	else
		error("agent_cnt underflow");
	if (agent_cnt < MAX_AGENT_CNT)
		agent_retry(RPC_RETRY_INTERVAL);
	slurm_mutex_unlock(&agent_cnt_mutex);
	pthread_cond_broadcast(&agent_cnt_cond);

	return NULL;
}

/* Basic validity test of agent argument */
static int _valid_agent_arg(agent_arg_t *agent_arg_ptr)
{
	xassert(agent_arg_ptr);
	xassert(agent_arg_ptr->slurm_addr);
	xassert(agent_arg_ptr->node_names);
	xassert((agent_arg_ptr->msg_type == SRUN_PING) ||
		(agent_arg_ptr->msg_type == SRUN_TIMEOUT) || 
		(agent_arg_ptr->msg_type == SRUN_NODE_FAIL) || 
		(agent_arg_ptr->msg_type == REQUEST_SIGNAL_JOB) ||
		(agent_arg_ptr->msg_type == REQUEST_TERMINATE_JOB) ||
		(agent_arg_ptr->msg_type == REQUEST_KILL_TIMELIMIT) || 
		(agent_arg_ptr->msg_type == REQUEST_UPDATE_JOB_TIME) ||
		(agent_arg_ptr->msg_type == REQUEST_SIGNAL_TASKS) ||
		(agent_arg_ptr->msg_type == REQUEST_TERMINATE_TASKS) || 
		(agent_arg_ptr->msg_type == REQUEST_PING) || 
		(agent_arg_ptr->msg_type == REQUEST_BATCH_JOB_LAUNCH) || 
		(agent_arg_ptr->msg_type == REQUEST_SHUTDOWN) || 
		(agent_arg_ptr->msg_type == REQUEST_SUSPEND) || 
		(agent_arg_ptr->msg_type == REQUEST_RECONFIGURE) ||
	        (agent_arg_ptr->msg_type == RESPONSE_RESOURCE_ALLOCATION) ||
		(agent_arg_ptr->msg_type == REQUEST_NODE_REGISTRATION_STATUS));

	if (agent_arg_ptr->node_count == 0)
		return SLURM_FAILURE;	/* no messages to be sent */
	return SLURM_SUCCESS;
}

static agent_info_t *_make_agent_info(agent_arg_t *agent_arg_ptr)
{
	int i;
	agent_info_t *agent_info_ptr;
	thd_t *thread_ptr;
	int *span = set_span(agent_arg_ptr->node_count, 0);
	int thr_count = 0;
       	forward_t forward;

	agent_info_ptr = xmalloc(sizeof(agent_info_t));
	slurm_mutex_init(&agent_info_ptr->thread_mutex);
	if (pthread_cond_init(&agent_info_ptr->thread_cond, NULL))
		fatal("pthread_cond_init error %m");
	agent_info_ptr->thread_count   = agent_arg_ptr->node_count;
	agent_info_ptr->retry          = agent_arg_ptr->retry;
	agent_info_ptr->threads_active = 0;
	thread_ptr = xmalloc(agent_info_ptr->thread_count * sizeof(thd_t));
	agent_info_ptr->thread_struct  = thread_ptr;
	agent_info_ptr->msg_type       = agent_arg_ptr->msg_type;
	agent_info_ptr->msg_args_pptr  = &agent_arg_ptr->msg_args;
	
	if ((agent_arg_ptr->msg_type != REQUEST_SHUTDOWN) &&
	    (agent_arg_ptr->msg_type != REQUEST_RECONFIGURE))
		agent_info_ptr->get_reply = true;

	forward.cnt = agent_info_ptr->thread_count;
	forward.name = agent_arg_ptr->node_names;
	forward.addr = agent_arg_ptr->slurm_addr;
	forward.node_id = NULL;
	forward.timeout = SLURM_MESSAGE_TIMEOUT_MSEC_STATIC;

	for (i = 0; i < agent_info_ptr->thread_count; i++) {
		thread_ptr[thr_count].state      = DSH_NEW;
		thread_ptr[thr_count].slurm_addr = 
			agent_arg_ptr->slurm_addr[i];
		strncpy(thread_ptr[thr_count].node_name,
			&agent_arg_ptr->node_names[i * MAX_SLURM_NAME],
			MAX_SLURM_NAME);

		forward_set(&thread_ptr[thr_count].forward,
			    span[thr_count],
			    &i,
			    &forward);

		thr_count++;		       
	}
	xfree(span);
	agent_info_ptr->thread_count = thr_count;
	return agent_info_ptr;
}

static task_info_t *_make_task_data(agent_info_t *agent_info_ptr, int inx)
{
	task_info_t *task_info_ptr;
	task_info_ptr = xmalloc(sizeof(task_info_t));

	task_info_ptr->thread_mutex_ptr  = &agent_info_ptr->thread_mutex;
	task_info_ptr->thread_cond_ptr   = &agent_info_ptr->thread_cond;
	task_info_ptr->threads_active_ptr= &agent_info_ptr->threads_active;
	task_info_ptr->thread_struct_ptr = &agent_info_ptr->thread_struct[inx];
	task_info_ptr->get_reply         = agent_info_ptr->get_reply;
	task_info_ptr->msg_type          = agent_info_ptr->msg_type;
	task_info_ptr->msg_args_ptr      = *agent_info_ptr->msg_args_pptr;

	return task_info_ptr;
}

static void _update_wdog_state(thd_t *thread_ptr, 
			       state_t *state, 
			       thd_complete_t *thd_comp)
{
	switch(*state) {
	case DSH_ACTIVE:
		thd_comp->work_done = false;
		if (thread_ptr->end_time <= thd_comp->now) {
			debug3("agent thread %lu timed out\n", 
			       (unsigned long) 
			       thread_ptr->thread);
			if (pthread_kill(thread_ptr->thread, SIGALRM) == ESRCH)
				*state = DSH_NO_RESP;
		}
		break;
	case DSH_NEW:
		thd_comp->work_done = false;
		break;
	case DSH_DONE:
		if (thd_comp->max_delay < 
		    (int)thread_ptr->end_time)
			thd_comp->max_delay = (int)thread_ptr->end_time;
		break;
	case DSH_NO_RESP:
		thd_comp->no_resp_cnt++;
		thd_comp->retry_cnt++;
		break;
	case DSH_FAILED:
		thd_comp->fail_cnt++;
		break;
	}
}

/* 
 * _wdog - Watchdog thread. Send SIGALRM to threads which have been active 
 *	for too long. 
 * IN args - pointer to agent_info_t with info on threads to watch
 * Sleep between polls with exponential times (from 0.125 to 1.0 second) 
 */
static void *_wdog(void *args)
{
	bool srun_agent = false;
	int i;
	agent_info_t *agent_ptr = (agent_info_t *) args;
	thd_t *thread_ptr = agent_ptr->thread_struct;
	unsigned long usec = 1250000;
	ListIterator itr;
	ret_types_t *ret_type = NULL;
	thd_complete_t thd_comp;


	if ( (agent_ptr->msg_type == SRUN_PING) ||
	     (agent_ptr->msg_type == SRUN_TIMEOUT) ||
	     (agent_ptr->msg_type == RESPONSE_RESOURCE_ALLOCATION) ||
	     (agent_ptr->msg_type == SRUN_NODE_FAIL) )
		srun_agent = true;

	thd_comp.max_delay = 0;
		
	while (1) {
		thd_comp.work_done   = true;/* assume all threads complete */
		thd_comp.fail_cnt    = 0;   /* assume no threads failures */
		thd_comp.no_resp_cnt = 0;   /* assume all threads respond */
		thd_comp.retry_cnt   = 0;   /* assume no required retries */
		thd_comp.now         = time(NULL);
		
		usleep(usec);
		usec = MIN((usec * 2), 1000000);		

		slurm_mutex_lock(&agent_ptr->thread_mutex);
		for (i = 0; i < agent_ptr->thread_count; i++) {
			if(!thread_ptr[i].ret_list) {
				_update_wdog_state(&thread_ptr[i],
						   &thread_ptr[i].state,
						   &thd_comp);
			} else {
				itr = list_iterator_create(
					thread_ptr[i].ret_list);
				while((ret_type = list_next(itr)) != NULL) {
					_update_wdog_state(
						&thread_ptr[i],
						(state_t *)&ret_type->msg_rc,
						&thd_comp);
				}
				list_iterator_destroy(itr);
			}
		}
		if (thd_comp.work_done)
			break;
		
		slurm_mutex_unlock(&agent_ptr->thread_mutex);
	}

	if (srun_agent) {
	        _notify_slurmctld_jobs(agent_ptr);
	} else {
	        _notify_slurmctld_nodes(agent_ptr, 
					thd_comp.no_resp_cnt, 
					thd_comp.retry_cnt);
	}

	for (i = 0; i < agent_ptr->thread_count; i++) {
		if (thread_ptr[i].ret_list)
			list_destroy(thread_ptr[i].ret_list);
	}

	if (thd_comp.max_delay)
		debug2("agent maximum delay %d seconds", thd_comp.max_delay);
	
	slurm_mutex_unlock(&agent_ptr->thread_mutex);
	return (void *) NULL;
}

static void _notify_slurmctld_jobs(agent_info_t *agent_ptr)
{
#if AGENT_IS_THREAD
	/* Locks: Write job */
	slurmctld_lock_t job_write_lock =
	    { NO_LOCK, WRITE_LOCK, NO_LOCK, NO_LOCK };
	uint32_t job_id = 0, step_id = 0;
	thd_t *thread_ptr = agent_ptr->thread_struct;

	if        (agent_ptr->msg_type == SRUN_PING) {
		srun_ping_msg_t *msg = *agent_ptr->msg_args_pptr;
		job_id  = msg->job_id;
		step_id = msg->step_id;
	} else if (agent_ptr->msg_type == SRUN_TIMEOUT) {
		srun_timeout_msg_t *msg = *agent_ptr->msg_args_pptr;
		job_id  = msg->job_id;
		step_id = msg->step_id;
	} else if (agent_ptr->msg_type == SRUN_NODE_FAIL) {
		srun_node_fail_msg_t *msg = *agent_ptr->msg_args_pptr;
		job_id  = msg->job_id;
		step_id = msg->step_id;
	} else if (agent_ptr->msg_type == RESPONSE_RESOURCE_ALLOCATION) {
		resource_allocation_response_msg_t *msg =
			*agent_ptr->msg_args_pptr;
		job_id  = msg->job_id;
		step_id = NO_VAL;
	} else {
		error("_notify_slurmctld_jobs invalid msg_type %u",
			agent_ptr->msg_type);
		return;
	}
	lock_slurmctld(job_write_lock);
	if  (thread_ptr[0].state == DSH_DONE) {
		srun_response(job_id, step_id);
	}
	unlock_slurmctld(job_write_lock);
#else
	fatal("Code development needed here if agent is not thread");
#endif
}

static void _notify_slurmctld_nodes(agent_info_t *agent_ptr, 
				    int no_resp_cnt, int retry_cnt)
{
	ListIterator itr = NULL;
	ListIterator data_itr;
	ret_types_t *ret_type = NULL;
	ret_data_info_t *ret_data_info = NULL;
	state_t state;
	int is_ret_list = 1;

#if AGENT_IS_THREAD
	/* Locks: Read config, write job, write node */
	slurmctld_lock_t node_write_lock =
	    { READ_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };
#endif
	thd_t *thread_ptr = agent_ptr->thread_struct;
	int i;

	/* Notify slurmctld of non-responding nodes */
	if (no_resp_cnt) {
#if AGENT_IS_THREAD
		/* Update node table data for non-responding nodes */
		lock_slurmctld(node_write_lock);
		if (agent_ptr->msg_type == REQUEST_BATCH_JOB_LAUNCH) {
			/* Requeue the request */
			batch_job_launch_msg_t *launch_msg_ptr = 
					*agent_ptr->msg_args_pptr;
			uint32_t job_id = launch_msg_ptr->job_id;
			job_complete(job_id, 0, true, 0);
		}
		unlock_slurmctld(node_write_lock);
#else
		fatal("Code development needed here if agent is not thread");
#endif
	}
	if (retry_cnt && agent_ptr->retry)
		_queue_agent_retry(agent_ptr, retry_cnt);

	/* Update last_response on responding nodes */
#if AGENT_IS_THREAD
	lock_slurmctld(node_write_lock);
	for (i = 0; i < agent_ptr->thread_count; i++) {
		if(!thread_ptr[i].ret_list) {
			state = thread_ptr[i].state;
			is_ret_list = 0;
			goto switch_on_state;
		} 
		is_ret_list = 1;
		
		itr = list_iterator_create(thread_ptr[i].ret_list);
		while((ret_type = list_next(itr)) != NULL) {
			state = ret_type->msg_rc;
		switch_on_state:
			switch(state) {
			case DSH_NO_RESP:
				if(!is_ret_list) {
					node_not_resp(thread_ptr[i].node_name,
						      thread_ptr[i].
						      start_time);
					break;
				}
				data_itr = list_iterator_create(
					ret_type->ret_data_list);
				while((ret_data_info = list_next(data_itr)) 
				      != NULL) 
					node_not_resp(ret_data_info->node_name,
						      thread_ptr[i].
						      start_time);
				list_iterator_destroy(data_itr);
				break;
			case DSH_FAILED:
				if(!is_ret_list) {
					set_node_down(thread_ptr[i].node_name, 
						      "Prolog/epilog failure");
					break;
				}
				data_itr = list_iterator_create(
					ret_type->ret_data_list);
				while((ret_data_info = list_next(data_itr)) 
				      != NULL)
					set_node_down(ret_data_info->node_name,
						      "Prolog/epilog failure");
				list_iterator_destroy(data_itr);
				break;
			case DSH_DONE:
				if(!is_ret_list) {
					node_did_resp(thread_ptr[i].node_name);
					break;
				}
				data_itr = list_iterator_create(
					ret_type->ret_data_list);
				while((ret_data_info = list_next(data_itr)) 
				      != NULL)
					node_did_resp(
						ret_data_info->node_name);
				list_iterator_destroy(data_itr);
				break;
			default:
				if(!is_ret_list) {
					error("unknown state returned for %s",
					      thread_ptr[i].node_name);
					break;
				}
				data_itr = list_iterator_create(
					ret_type->ret_data_list);
				while((ret_data_info = list_next(data_itr)) 
				      != NULL)
					error("unknown state returned for %s",
					      ret_data_info->node_name);
				list_iterator_destroy(data_itr);
				break;
			}
			if(!is_ret_list)
				goto finished;
		}
		list_iterator_destroy(itr);
	}
finished:
	unlock_slurmctld(node_write_lock);
	if (run_scheduler) {
		run_scheduler = false;
		/* below functions all have their own locking */
		if (schedule())	{
			schedule_job_save();
			schedule_node_save();
		}
	}
	if ((agent_ptr->msg_type == REQUEST_PING) ||
	    (agent_ptr->msg_type == REQUEST_NODE_REGISTRATION_STATUS))
		ping_end();
#else
	fatal("Code development needed here if agent is not thread");
#endif
}

/* Report a communications error for specified node */
static inline void _comm_err(char *node_name)
{
#if AGENT_IS_THREAD
	if (is_node_resp (node_name))
#endif
		error("agent/send_recv_msg: %s: %m", node_name);
}

/*
 * _thread_per_group_rpc - thread to issue an RPC for a group of nodes
 *                         sending message out to one and forwarding it to
 *                         others if necessary.
 * IN/OUT args - pointer to task_info_t, xfree'd on completion
 */
static void *_thread_per_group_rpc(void *args)
{
	int rc = SLURM_SUCCESS;
	slurm_msg_t msg;
	task_info_t *task_ptr = (task_info_t *) args;
	/* we cache some pointers from task_info_t because we need 
	 * to xfree args before being finished with their use. xfree 
	 * is required for timely termination of this pthread because 
	 * xfree could lock it at the end, preventing a timely
	 * thread_exit */
	pthread_mutex_t *thread_mutex_ptr   = task_ptr->thread_mutex_ptr;
	pthread_cond_t  *thread_cond_ptr    = task_ptr->thread_cond_ptr;
	uint32_t        *threads_active_ptr = task_ptr->threads_active_ptr;
	thd_t           *thread_ptr         = task_ptr->thread_struct_ptr;
	state_t thread_state = DSH_NO_RESP;
	slurm_msg_type_t msg_type = task_ptr->msg_type;
	bool is_kill_msg, srun_agent;
	List ret_list = NULL;
	ListIterator itr;
	ListIterator data_itr;
	ret_types_t *ret_type = NULL;
	ret_data_info_t *ret_data_info = NULL;
	int found = 0;

#if AGENT_IS_THREAD
	/* Locks: Write job, write node */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };
	/* Locks: Read job */
	slurmctld_lock_t job_read_lock = {
		NO_LOCK, READ_LOCK, NO_LOCK, NO_LOCK };
#endif
	xassert(args != NULL);

	is_kill_msg = (	(msg_type == REQUEST_KILL_TIMELIMIT) ||
			(msg_type == REQUEST_TERMINATE_JOB) );
	srun_agent = (	(msg_type == SRUN_PING)    ||
			(msg_type == SRUN_TIMEOUT) ||
			(msg_type == RESPONSE_RESOURCE_ALLOCATION) ||
			(msg_type == SRUN_NODE_FAIL) );

	thread_ptr->start_time = time(NULL);

	/* don't try to communicate with defunct job */
#if AGENT_IS_THREAD
	if (srun_agent) {
		uint32_t          job_id   = 0;
		enum job_states    state   = JOB_END;
		struct job_record *job_ptr = NULL;

		if (msg_type == SRUN_PING) {
			srun_ping_msg_t *msg = task_ptr->msg_args_ptr;
			job_id  = msg->job_id;
		} else if (msg_type == SRUN_TIMEOUT) {
			srun_timeout_msg_t *msg = task_ptr->msg_args_ptr;
			job_id  = msg->job_id;
		} else if (msg_type == SRUN_NODE_FAIL) {
			srun_node_fail_msg_t *msg = task_ptr->msg_args_ptr;
			job_id  = msg->job_id;
		} else if (msg_type == RESPONSE_RESOURCE_ALLOCATION) {
			resource_allocation_response_msg_t *msg = 
				task_ptr->msg_args_ptr;
			job_id  = msg->job_id;
		}
		lock_slurmctld(job_read_lock);
		if (job_id)
			job_ptr = find_job_record(job_id);
		if (job_ptr)
			state = job_ptr->job_state;	
		unlock_slurmctld(job_read_lock);
		if ((state == JOB_RUNNING) ||
		    ((state & JOB_COMPLETING) 
		     && (msg_type == SRUN_NODE_FAIL))) {
			; /* proceed with the communication */
		} else {	
			thread_state = DSH_DONE;
			goto cleanup;
		}
	}
#endif

	slurm_mutex_lock(thread_mutex_ptr);
	thread_ptr->state = DSH_ACTIVE;
	slurm_mutex_unlock(thread_mutex_ptr);

	/* send request message */
	msg.address  = thread_ptr->slurm_addr;
	msg.msg_type = msg_type;
	msg.data     = task_ptr->msg_args_ptr;
	msg.forward  = thread_ptr->forward;
	msg.ret_list = NULL;
	msg.orig_addr.sin_addr.s_addr = 0;
	msg.srun_node_id = 0;
	msg.forward_struct_init = 0;

	//info("forwarding to %d",msg.forward.cnt);
	thread_ptr->end_time = thread_ptr->start_time + COMMAND_TIMEOUT;
	if (task_ptr->get_reply) {
		if ((ret_list = slurm_send_recv_rc_msg(&msg, 
						       msg.forward.timeout)) 
		    == NULL) {
			if (!srun_agent)
				_comm_err(thread_ptr->node_name);
			goto cleanup;
		}
	} else {
		if (slurm_send_only_node_msg(&msg) < 0) {
			if (!srun_agent)
				_comm_err(thread_ptr->node_name);
		} else
			thread_state = DSH_DONE;
		goto cleanup;
	}

	//info("got %d states back from the send", list_count(ret_list));
	found = 0;
	itr = list_iterator_create(ret_list);		
	while((ret_type = list_next(itr)) != NULL) {
		data_itr = list_iterator_create(ret_type->ret_data_list);
		while((ret_data_info = list_next(data_itr)) != NULL) {
			rc = ret_type->msg_rc;
			if(!found 
			   && !strcmp(ret_data_info->node_name,"localhost")) {
			  //info("got localhost");
				xfree(ret_data_info->node_name);
				ret_data_info->node_name = 
					xstrdup(thread_ptr->node_name);
				found = 1;
			}
			/* info("response for %s rc = %d", */
/* 			     ret_data_info->node_name, */
/* 			     ret_type->msg_rc); */
			if(rc == SLURM_ERROR) {
				errno = ret_type->err;
				continue;
			}
#if AGENT_IS_THREAD
			
			/* SPECIAL CASE: Mark node as IDLE if job already 
			   complete */
			if (is_kill_msg && 
			    (rc == ESLURMD_KILL_JOB_ALREADY_COMPLETE)) {
				kill_job_msg_t *kill_job;
				kill_job = (kill_job_msg_t *) 
					task_ptr->msg_args_ptr;
				rc = SLURM_SUCCESS;
				lock_slurmctld(job_write_lock);
				if (job_epilog_complete(kill_job->job_id, 
							ret_data_info->
							node_name, 
							rc))
					run_scheduler = true;
				unlock_slurmctld(job_write_lock);
			}
			
			/* SPECIAL CASE: Kill non-startable batch job */
			if ((msg_type == REQUEST_BATCH_JOB_LAUNCH)
			    && rc) {
				batch_job_launch_msg_t *launch_msg_ptr = 
					task_ptr->msg_args_ptr;
				uint32_t job_id = launch_msg_ptr->job_id;
				info("Killing non-startable batch job %u: %s", 
				     job_id, slurm_strerror(rc));
				thread_state = DSH_DONE;
				lock_slurmctld(job_write_lock);
				job_complete(job_id, 0, false, 1);
				unlock_slurmctld(job_write_lock);
				continue;
				//goto cleanup;
			}
#endif
		}
		list_iterator_destroy(data_itr);
	
		if (((msg_type == REQUEST_SIGNAL_TASKS) 
		     ||   (msg_type == REQUEST_TERMINATE_TASKS)) 
		    && (rc == ESRCH)) {
			/* process is already dead, not a real error */
			rc = SLURM_SUCCESS;
		}
		
		switch (rc) {
		case SLURM_SUCCESS:
			/*debug3("agent processed RPC to node %s", 
			  ret_data_info->node_name); */
			thread_state = DSH_DONE;
			break;
		case ESLURMD_EPILOG_FAILED:
			data_itr = 
				list_iterator_create(ret_type->ret_data_list);
			while((ret_data_info = list_next(data_itr)) != NULL) 
				error("Epilog failure on host %s, "
				      "setting DOWN", 
				      ret_data_info->node_name);
			list_iterator_destroy(data_itr);
	
			thread_state = DSH_FAILED;
			break;
		case ESLURMD_PROLOG_FAILED:
			data_itr = 
				list_iterator_create(ret_type->ret_data_list);
			while((ret_data_info = list_next(data_itr)) != NULL) 
				error("Prolog failure on host %s, "
				      "setting DOWN", 
				      ret_data_info->node_name);
			list_iterator_destroy(data_itr);
			thread_state = DSH_FAILED;
			break;
		case ESLURM_INVALID_JOB_ID:  
			/* Not indicative of a real error */
		case ESLURMD_JOB_NOTRUNNING: 
			/* Not indicative of a real error */
			data_itr = 
				list_iterator_create(ret_type->ret_data_list);
			while((ret_data_info = list_next(data_itr)) != NULL) 
				debug2("agent processed RPC to node %s: %s",
				       ret_data_info->node_name,
				       slurm_strerror(rc));
			list_iterator_destroy(data_itr);
	
			thread_state = DSH_DONE;
			break;
		default:
			data_itr = 
				list_iterator_create(ret_type->ret_data_list);
			while((ret_data_info = list_next(data_itr)) != NULL) 
				if (!srun_agent) {
					errno = ret_type->err;
					_comm_err(ret_data_info->node_name);
				}
			list_iterator_destroy(data_itr);
			if (srun_agent)
				thread_state = DSH_FAILED;
			else	 /* Not serious error, don't DRAIN node */
				thread_state = DSH_DONE;
		}	
		ret_type->msg_rc = thread_state;
	}
	list_iterator_destroy(itr);

cleanup:
	xfree(args);
	destroy_forward(&thread_ptr->forward);
	slurm_mutex_lock(thread_mutex_ptr);
	thread_ptr->ret_list = ret_list;
	thread_ptr->state = thread_state;
	thread_ptr->end_time = (time_t) difftime(time(NULL), 
						 thread_ptr->start_time);
	/* Signal completion so another thread can replace us */
	(*threads_active_ptr)--;
	slurm_mutex_unlock(thread_mutex_ptr);
	pthread_cond_signal(thread_cond_ptr);
	return (void *) NULL;
}

/*
 * SIGALRM handler.  We are really interested in interrupting hung communictions
 * and causing them to return EINTR. Multiple interupts might be required.
 */
static void _alarm_handler(int dummy)
{
	xsignal(SIGALRM, _alarm_handler);
}

/*
 * _queue_agent_retry - Queue any failed RPCs for later replay
 * IN agent_info_ptr - pointer to info on completed agent requests
 * IN count - number of agent requests which failed, count to requeue
 */
static void _queue_agent_retry(agent_info_t * agent_info_ptr, int count)
{
	agent_arg_t *agent_arg_ptr;
	queued_request_t *queued_req_ptr = NULL;
	thd_t *thread_ptr = agent_info_ptr->thread_struct;
	int i, j;

	if (count == 0)
		return;

	/* build agent argument with just the RPCs to retry */
	agent_arg_ptr = xmalloc(sizeof(agent_arg_t));
	agent_arg_ptr->node_count = count;
	agent_arg_ptr->retry = 1;
	agent_arg_ptr->slurm_addr = xmalloc(sizeof(slurm_addr) * count);
	agent_arg_ptr->node_names = 
		xmalloc(sizeof(char) * MAX_SLURM_NAME * count);
	agent_arg_ptr->msg_type = agent_info_ptr->msg_type;
	agent_arg_ptr->msg_args = *(agent_info_ptr->msg_args_pptr);
	*(agent_info_ptr->msg_args_pptr) = NULL;

	j = 0;
	for (i = 0; i < agent_info_ptr->thread_count; i++) {
		if (thread_ptr[i].state != DSH_NO_RESP)
			continue;
		agent_arg_ptr->slurm_addr[j] = thread_ptr[i].slurm_addr;
		strncpy(&agent_arg_ptr->node_names[j * MAX_SLURM_NAME],
			thread_ptr[i].node_name, MAX_SLURM_NAME);
		if ((++j) == count)
			break;
	}
	if (count != j) {
		error("agent: Retry count (%d) != actual count (%d)", 
			count, j);
		agent_arg_ptr->node_count = j;
	}
	debug2("Queue RPC msg_type=%u, nodes=%d for retry", 
	       agent_arg_ptr->msg_type, j);

	/* add the requeust to a list */
	queued_req_ptr = xmalloc(sizeof(queued_request_t));
	queued_req_ptr->agent_arg_ptr = agent_arg_ptr;
	queued_req_ptr->last_attempt  = time(NULL);
	slurm_mutex_lock(&retry_mutex);
	if (retry_list == NULL) {
		retry_list = list_create(&_list_delete_retry);
		if (retry_list == NULL)
			fatal("list_create failed");
	}
	if (list_append(retry_list, (void *) queued_req_ptr) == 0)
		fatal("list_append failed");
	slurm_mutex_unlock(&retry_mutex);
}

/*
 * _list_delete_retry - delete an entry from the retry list, 
 *	see common/list.h for documentation
 */
static void _list_delete_retry(void *retry_entry)
{
	queued_request_t *queued_req_ptr;

	if (! retry_entry)
		return;

	queued_req_ptr = (queued_request_t *) retry_entry;
	_purge_agent_args(queued_req_ptr->agent_arg_ptr);
	xfree(queued_req_ptr);
}


/*
 * agent_retry - Agent for retrying pending RPCs. One pending request is 
 *	issued if it has been pending for at least min_wait seconds
 * IN min_wait - Minimum wait time between re-issue of a pending RPC
 * RET count of queued requests remaining (zero if none are old enough 
 * to re-issue)
 */
extern int agent_retry (int min_wait)
{
	int list_size = 0;
	time_t now = time(NULL);
	queued_request_t *queued_req_ptr = NULL;

	if (retry_list)
		list_size = list_count(retry_list);
	if (agent_cnt >= MAX_AGENT_CNT)		/* too much work already */
		return list_size;

	slurm_mutex_lock(&retry_mutex);
	if (retry_list) {
		double age = 0;
		queued_req_ptr = (queued_request_t *) list_peek(retry_list);
		if (queued_req_ptr) {
			age = difftime(now, queued_req_ptr->last_attempt);
			if (age > min_wait) {
				queued_req_ptr = (queued_request_t *) 
					list_pop(retry_list);
			} else { /* too new */
				queued_req_ptr = NULL;
				list_size = 0;
			}
		}
	}
	slurm_mutex_unlock(&retry_mutex);

	if (queued_req_ptr) {
		agent_arg_t *agent_arg_ptr = queued_req_ptr->agent_arg_ptr;
		xfree(queued_req_ptr);
		if (agent_arg_ptr)
			_spawn_retry_agent(agent_arg_ptr);
		else
			error("agent_retry found record with no agent_args");
	} else if (mail_list) {
		mail_info_t *mi;
		slurm_mutex_lock(&mail_mutex);
		mi = (mail_info_t *) list_dequeue(mail_list);
		slurm_mutex_unlock(&mail_mutex);
		if (mi)
			_mail_proc(mi);
	}

	return list_size;
}

/*
 * agent_queue_request - put a new request on the queue for execution or
 * 	execute now if not too busy
 * IN agent_arg_ptr - the request to enqueue
 */
void agent_queue_request(agent_arg_t *agent_arg_ptr)
{
	queued_request_t *queued_req_ptr = NULL;

	if (agent_cnt < MAX_AGENT_CNT) {	/* execute now */
		pthread_attr_t attr_agent;
		pthread_t thread_agent;
		int rc;
		slurm_attr_init(&attr_agent);
		if (pthread_attr_setdetachstate
				(&attr_agent, PTHREAD_CREATE_DETACHED))
			error("pthread_attr_setdetachstate error %m");
		rc = pthread_create(&thread_agent, &attr_agent,
				    agent, (void *) agent_arg_ptr);
		slurm_attr_destroy(&attr_agent);
		if (rc == 0)
			return;
	}

	queued_req_ptr = xmalloc(sizeof(queued_request_t));
	queued_req_ptr->agent_arg_ptr = agent_arg_ptr;
/*	queued_req_ptr->last_attempt  = 0; Implicit */

	slurm_mutex_lock(&retry_mutex);
	if (retry_list == NULL) {
		retry_list = list_create(&_list_delete_retry);
		if (retry_list == NULL)
			fatal("list_create failed");
	}
	list_prepend(retry_list, (void *)queued_req_ptr);
	slurm_mutex_unlock(&retry_mutex);
}

/* _spawn_retry_agent - pthread_create an agent for the given task */
static void _spawn_retry_agent(agent_arg_t * agent_arg_ptr)
{
	int retries = 0;
	pthread_attr_t attr_agent;
	pthread_t thread_agent;

	if (agent_arg_ptr == NULL)
		return;

	debug2("Spawning RPC agent for msg_type %u", 
	       agent_arg_ptr->msg_type);
	slurm_attr_init(&attr_agent);
	if (pthread_attr_setdetachstate(&attr_agent,
					PTHREAD_CREATE_DETACHED))
		error("pthread_attr_setdetachstate error %m");
	while (pthread_create(&thread_agent, &attr_agent,
			agent, (void *) agent_arg_ptr)) {
		error("pthread_create error %m");
		if (++retries > MAX_RETRIES)
			fatal("Can't create pthread");
		sleep(1);	/* sleep and try again */
	}
	slurm_attr_destroy(&attr_agent);
}

/* _slurmctld_free_job_launch_msg is a variant of slurm_free_job_launch_msg
 *	because all environment variables currently loaded in one xmalloc 
 *	buffer (see get_job_env()), which is different from how slurmd 
 *	assembles the data from a message
 */
static void _slurmctld_free_job_launch_msg(batch_job_launch_msg_t * msg)
{
	if (msg) {
		if (msg->environment) {
			xfree(msg->environment[0]);
			xfree(msg->environment);
		}
		slurm_free_job_launch_msg(msg);
	}
}

/* agent_purge - purge all pending RPC requests */
void agent_purge(void)
{
	if (retry_list == NULL)
		return;

	slurm_mutex_lock(&retry_mutex);
	list_destroy(retry_list);
	retry_list = NULL;
	if (mail_list)
		list_destroy(mail_list);
	mail_list = NULL;
	slurm_mutex_unlock(&retry_mutex);
}

static void _purge_agent_args(agent_arg_t *agent_arg_ptr)
{
	if (agent_arg_ptr == NULL)
		return;

	xfree(agent_arg_ptr->slurm_addr);
	xfree(agent_arg_ptr->node_names);
	if (agent_arg_ptr->msg_args) {
		if (agent_arg_ptr->msg_type == REQUEST_BATCH_JOB_LAUNCH)
			_slurmctld_free_job_launch_msg(agent_arg_ptr->msg_args);
		else if (agent_arg_ptr->msg_type == 
				RESPONSE_RESOURCE_ALLOCATION)
			slurm_free_resource_allocation_response_msg(
					agent_arg_ptr->msg_args);
		else if ((agent_arg_ptr->msg_type == REQUEST_SIGNAL_JOB)
		||       (agent_arg_ptr->msg_type == REQUEST_TERMINATE_JOB)
		||       (agent_arg_ptr->msg_type == REQUEST_KILL_TIMELIMIT))
			slurm_free_kill_job_msg(agent_arg_ptr->msg_args);
		else
			xfree(agent_arg_ptr->msg_args);
	}
	xfree(agent_arg_ptr);
}

static mail_info_t *_mail_alloc(void)
{
	return xmalloc(sizeof(mail_info_t));
}

static void _mail_free(void *arg)
{
	mail_info_t *mi = (mail_info_t *) arg;

	if (mi) {
		xfree(mi->user_name);
		xfree(mi->message);
		xfree(mi);
	}
}

/* process an email request and free the record */
static void _mail_proc(mail_info_t *mi)
{
	pid_t pid;
	int pfd[2];

	if (pipe(pfd) == -1) {
		error("pipe(): %m");
		goto fini;
	}

	pid = fork();
	if (pid < 0) {
		error("fork(): %m");
		(void) close(pfd[0]);
		(void) close(pfd[1]);
	} else if (pid == 0) { /* child */
		(void) close(0);
		dup(pfd[0]);
		(void) close(pfd[0]);
		(void) close(pfd[1]);
		(void) close(1);
		(void) close(2);
		execle("/bin/mail", "mail", mi->user_name, 
			"-s", mi->message, NULL, NULL);
		exit(1);
	} else {	/* parent */
		(void) close(pfd[0]);
		(void) close(pfd[1]);
		waitpid(pid, NULL, 0);
	}
fini:	_mail_free(mi);
	return;
}

static char *_mail_type_str(uint16_t mail_type)
{
	if (mail_type == MAIL_JOB_BEGIN)
		return "Began";
	if (mail_type == MAIL_JOB_END)
		return "Ended";
	if (mail_type == MAIL_JOB_FAIL)
		return "Failed";
	return "unknown";
}

/*
 * mail_job_info - Send e-mail notice of job state change
 * IN job_ptr - job identification
 * IN state_type - job transition type, see MAIL_JOB in slurm.h
 */
extern void mail_job_info (struct job_record *job_ptr, uint16_t mail_type)
{
	mail_info_t *mi = _mail_alloc();

	if (!job_ptr->mail_user) {
		struct passwd *pw;
		pw = getpwuid((uid_t) job_ptr->user_id);
		if (pw && pw->pw_name)
			mi->user_name = xstrdup(pw->pw_name);
		else {
			error("getpwuid(%u): %m", job_ptr->user_id);
			_mail_free(mi);
			return;
		}
	} else
		mi->user_name = xstrdup(job_ptr->mail_user);

	mi->message = xmalloc(sizeof(char)*128);
	sprintf(mi->message, "SLURM Job_id=%u Name=%.24s %s",
		job_ptr->job_id, job_ptr->name, 
		_mail_type_str(mail_type));

	info ("msg to %s: %s", mi->user_name, mi->message);

	slurm_mutex_lock(&mail_mutex);
	if (!mail_list) {
		mail_list = list_create(&_mail_free);
		if (!mail_list)
			fatal("list_create failed");
	}
	if (!list_enqueue(mail_list, (void *) mi))
		fatal("list_enqueue failed");
	slurm_mutex_unlock(&mail_mutex);
	return;
}

