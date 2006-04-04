/*****************************************************************************\
 *  step_mgr.c - manage the job step information of slurm
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>, et. al.
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
\*****************************************************************************/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <unistd.h>

#include <slurm/slurm_errno.h>

#include "src/common/bitstring.h"
#include "src/common/checkpoint.h"
#include "src/common/dist_tasks.h"
#include "src/common/slurm_protocol_interface.h"
#include "src/common/switch.h"
#include "src/common/xstring.h"

#include "src/slurmctld/agent.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/node_scheduler.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/jobacct.h"

#define MAX_RETRIES 10

static int _job_step_ckpt_error(struct step_record *step_ptr, slurm_fd conn_fd);
static void _pack_ctld_job_step_info(struct step_record *step, Buf buffer);
static bitstr_t * _pick_step_nodes (struct job_record  *job_ptr, 
				    job_step_create_request_msg_t *step_spec );
/* 
 * create_step_record - create an empty step_record for the specified job.
 * IN job_ptr - pointer to job table entry to have step record added
 * RET a pointer to the record or NULL if error
 * NOTE: allocates memory that should be xfreed with delete_step_record
 */
struct step_record * 
create_step_record (struct job_record *job_ptr) 
{
	struct step_record *step_ptr;

	xassert(job_ptr);
	step_ptr = (struct step_record *) xmalloc(sizeof (struct step_record));

	last_job_update = time(NULL);
	step_ptr->job_ptr = job_ptr; 
	step_ptr->step_id = (job_ptr->next_step_id)++;
	step_ptr->start_time = time ( NULL ) ;

	if (list_append (job_ptr->step_list, step_ptr) == NULL)
		fatal ("create_step_record: unable to allocate memory");

	return step_ptr;
}


/* 
 * delete_all_step_records - delete all step record for specified job_ptr
 * IN job_ptr - pointer to job table entry to have step record added
 */
void 
delete_all_step_records (struct job_record *job_ptr) 
{
	ListIterator step_iterator;
	struct step_record *step_ptr;

	xassert(job_ptr);
	step_iterator = list_iterator_create (job_ptr->step_list);

	last_job_update = time(NULL);
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		list_remove (step_iterator);
		if (step_ptr->switch_job) {
			switch_g_job_step_complete(step_ptr->switch_job,
				step_ptr->step_node_list);
			switch_free_jobinfo(step_ptr->switch_job);
		}
		checkpoint_free_jobinfo(step_ptr->check_job);
		xfree(step_ptr->host);
		xfree(step_ptr->name);
		xfree(step_ptr->step_node_list);
		FREE_NULL_BITMAP(step_ptr->step_node_bitmap);
		FREE_NULL_BITMAP(step_ptr->exit_node_bitmap);
		if (step_ptr->network)
			xfree(step_ptr->network);
		xfree(step_ptr);
	}		

	list_iterator_destroy (step_iterator);
}


/* 
 * delete_step_record - delete record for job step for specified job_ptr 
 *	and step_id
 * IN job_ptr - pointer to job table entry to have step record removed
 * IN step_id - id of the desired job step
 * RET 0 on success, errno otherwise
 */
int 
delete_step_record (struct job_record *job_ptr, uint32_t step_id) 
{
	ListIterator step_iterator;
	struct step_record *step_ptr;
	int error_code;

	xassert(job_ptr);
	error_code = ENOENT;
	step_iterator = list_iterator_create (job_ptr->step_list);		

	last_job_update = time(NULL);
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		if (step_ptr->step_id == step_id) {
			list_remove (step_iterator);
/* FIXME: If job step record is preserved after completion, 
 * the switch_g_job_step_complete() must be called upon completion 
 * and not upon record purging. Presently both events occur 
 * simultaneously. */
			if (step_ptr->switch_job) {
				switch_g_job_step_complete(step_ptr->switch_job, 
					step_ptr->step_node_list);
				switch_free_jobinfo (step_ptr->switch_job);
			}
			checkpoint_free_jobinfo (step_ptr->check_job);
			xfree(step_ptr->host);
			xfree(step_ptr->name);
			xfree(step_ptr->step_node_list);
			FREE_NULL_BITMAP(step_ptr->step_node_bitmap);
			FREE_NULL_BITMAP(step_ptr->exit_node_bitmap);
			if (step_ptr->network)
				xfree(step_ptr->network);
			xfree(step_ptr);
			error_code = 0;
			break;
		}
	}		

	list_iterator_destroy (step_iterator);
	return error_code;
}


/*
 * dump_step_desc - dump the incoming step initiate request message
 * IN step_spec - job step request specification from RPC
 */
void
dump_step_desc(job_step_create_request_msg_t *step_spec)
{
	if (step_spec == NULL) 
		return;

	debug3("StepDesc: user_id=%u job_id=%u node_count=%u, cpu_count=%u", 
		step_spec->user_id, step_spec->job_id, 
		step_spec->node_count, step_spec->cpu_count);
	debug3("   num_tasks=%u relative=%u task_dist=%u node_list=%s", 
		step_spec->num_tasks, step_spec->relative, 
		step_spec->task_dist, step_spec->node_list);
	debug3("   host=%s port=%u name=%s network=%s", 
		step_spec->host, step_spec->port, step_spec->name,
		step_spec->network);
}


/* 
 * find_step_record - return a pointer to the step record with the given 
 *	job_id and step_id
 * IN job_ptr - pointer to job table entry to have step record added
 * IN step_id - id of the desired job step or NO_VAL for first one
 * RET pointer to the job step's record, NULL on error
 */
struct step_record *
find_step_record(struct job_record *job_ptr, uint16_t step_id) 
{
	ListIterator step_iterator;
	struct step_record *step_ptr;

	if (job_ptr == NULL)
		return NULL;

	step_iterator = list_iterator_create (job_ptr->step_list);
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		if ((step_ptr->step_id == step_id)
		||  ((uint16_t) step_id == (uint16_t) NO_VAL)) {
			break;
		}
	}		
	list_iterator_destroy (step_iterator);

	return step_ptr;
}


/* 
 * job_step_signal - signal the specified job step
 * IN job_id - id of the job to be cancelled
 * IN step_id - id of the job step to be cancelled
 * IN signal - user id of user issuing the RPC
 * IN uid - user id of user issuing the RPC
 * RET 0 on success, otherwise ESLURM error code 
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 */
int job_step_signal(uint32_t job_id, uint32_t step_id, 
		    uint16_t signal, uid_t uid)
{
	struct job_record *job_ptr;
	struct step_record *step_ptr;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		error("job_step_cancel: invalid job id %u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}

	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;
	if (job_ptr->job_state != JOB_RUNNING) {
		verbose("job_step_signal: step %u.%u can not be sent signal "
			"%u from state=%s", job_id, step_id, signal,
			job_state_string(job_ptr->job_state));
		return ESLURM_TRANSITION_STATE_NO_UPDATE;
	}

	if ((job_ptr->user_id != uid) && (uid != 0) && (uid != getuid())) {
		error("Security violation, JOB_CANCEL RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}

	step_ptr = find_step_record(job_ptr, step_id);
	if (step_ptr == NULL) {
		info("job_step_cancel step %u.%u not found",
		     job_id, step_id);
		return ESLURM_INVALID_JOB_ID;
	}

	signal_step_tasks(step_ptr, signal);
	return SLURM_SUCCESS;
}

/*
 * signal_step_tasks - send specific signal to specific job step
 * IN step_ptr - step record pointer
 * IN signal - signal to send
 */
void signal_step_tasks(struct step_record *step_ptr, uint16_t signal)
{
	int i;
	kill_tasks_msg_t *kill_tasks_msg;
	agent_arg_t *agent_args;
	int buf_rec_size = 0;

	xassert(step_ptr);
	agent_args = xmalloc(sizeof(agent_arg_t));
	if (signal == SIGKILL)
		agent_args->msg_type = REQUEST_TERMINATE_TASKS;
	else
		agent_args->msg_type = REQUEST_SIGNAL_TASKS;
	agent_args->retry = 1;
	kill_tasks_msg = xmalloc(sizeof(kill_tasks_msg_t));
	kill_tasks_msg->job_id      = step_ptr->job_ptr->job_id;
	kill_tasks_msg->job_step_id = step_ptr->step_id;
	kill_tasks_msg->signal      = signal;

	for (i = 0; i < node_record_count; i++) {
		if (bit_test(step_ptr->step_node_bitmap, i) == 0)
			continue;
		if ((agent_args->node_count + 1) > buf_rec_size) {
			buf_rec_size += 128;
			xrealloc((agent_args->slurm_addr),
				 (sizeof(struct sockaddr_in) *
				  buf_rec_size));
			xrealloc((agent_args->node_names),
				 (MAX_SLURM_NAME * buf_rec_size));
		}
		agent_args->slurm_addr[agent_args->node_count] =
		    node_record_table_ptr[i].slurm_addr;
		strncpy(&agent_args->
			node_names[MAX_SLURM_NAME * agent_args->node_count],
			node_record_table_ptr[i].name, MAX_SLURM_NAME);
		agent_args->node_count++;
#ifdef HAVE_FRONT_END		/* Operate only on front-end */
		break;
#endif
	}

	if (agent_args->node_count == 0) {
		xfree(kill_tasks_msg);
		xfree(agent_args);
		return;
	}

	agent_args->msg_args = kill_tasks_msg;
	agent_queue_request(agent_args);
	return;
}


/* 
 * job_step_complete - note normal completion the specified job step
 * IN job_id - id of the job to be completed
 * IN step_id - id of the job step to be completed
 * IN uid - user id of user issuing the RPC
 * IN requeue - job should be run again if possible
 * IN job_return_code - job's return code, if set then set state to JOB_FAILED
 * RET 0 on success, otherwise ESLURM error code 
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 */
int job_step_complete(uint32_t job_id, uint32_t step_id, uid_t uid,
		      bool requeue, uint32_t job_return_code)
{
	struct job_record *job_ptr;
	struct step_record *step_ptr;
	int error_code;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		info("job_step_complete: invalid job id %u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}
	
	step_ptr = find_step_record(job_ptr, step_id);
	if (step_ptr == NULL)
		return ESLURM_INVALID_JOB_ID;
	else 
		jobacct_step_complete(step_ptr);
	
	if ((job_ptr->kill_on_step_done)
	&&  (list_count(job_ptr->step_list) <= 1)
	&&  (!IS_JOB_FINISHED(job_ptr)))
		return job_complete(job_id, uid, requeue, job_return_code);

	if ((job_ptr->user_id != uid) && (uid != 0) && (uid != getuid())) {
		error("Security violation, JOB_COMPLETE RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}

	last_job_update = time(NULL);
	error_code = delete_step_record(job_ptr, step_id);
	if (error_code == ENOENT) {
		info("job_step_complete step %u.%u not found", job_id,
		     step_id);
		return ESLURM_ALREADY_DONE;
	}
	return SLURM_SUCCESS;
}

/*
 * aggregate_step_data - given a step_record, aggregate all process info
 * IN step - pointer to step record
 * IN rusage - rusage struct
 * IN psize - psize
 * IN vsize - vsize 
 * RET NONE
 */
void aggregate_step_data(struct step_record *step,
			 struct rusage rusage, int psize, int vsize)
{
	step->max_psize = MAX(step->max_psize, psize);
	step->max_vsize = MAX(step->max_vsize, vsize);
	/* sum up all rusage stuff */
	step->rusage.ru_utime.tv_sec	+= rusage.ru_utime.tv_sec;
	step->rusage.ru_utime.tv_usec	+= rusage.ru_utime.tv_usec;
	while (step->rusage.ru_utime.tv_usec >= 1E6) {
		step->rusage.ru_utime.tv_sec++;
		step->rusage.ru_utime.tv_usec -= 1E6;
	}
	step->rusage.ru_stime.tv_sec	+= rusage.ru_stime.tv_sec;
	step->rusage.ru_stime.tv_usec	+= rusage.ru_stime.tv_usec;
	while (step->rusage.ru_stime.tv_usec >= 1E6) {
		step->rusage.ru_stime.tv_sec++;
		step->rusage.ru_stime.tv_usec -= 1E6;
	}
	step->rusage.ru_maxrss		+= rusage.ru_maxrss;
	step->rusage.ru_ixrss		+= rusage.ru_ixrss;
	step->rusage.ru_idrss		+= rusage.ru_idrss;
	step->rusage.ru_isrss		+= rusage.ru_isrss;
	step->rusage.ru_minflt		+= rusage.ru_minflt;
	step->rusage.ru_majflt		+= rusage.ru_majflt;
	step->rusage.ru_nswap		+= rusage.ru_nswap;
	step->rusage.ru_inblock		+= rusage.ru_inblock;
	step->rusage.ru_oublock		+= rusage.ru_oublock;
	step->rusage.ru_msgsnd		+= rusage.ru_msgsnd;
	step->rusage.ru_msgrcv		+= rusage.ru_msgrcv;
	step->rusage.ru_nsignals	+= rusage.ru_nsignals;
	step->rusage.ru_nvcsw		+= rusage.ru_nvcsw;
	step->rusage.ru_nivcsw		+= rusage.ru_nivcsw;
}

/* 
 * _pick_step_nodes - select nodes for a job step that satify its requirements
 *	we satify the super-set of constraints.
 * IN job_ptr - pointer to job to have new step started
 * IN step_spec - job step specification
 * global: node_record_table_ptr - pointer to global node table
 * NOTE: returns all of a job's nodes if step_spec->node_count == INFINITE
 * NOTE: returned bitmap must be freed by the caller using bit_free()
 */
static bitstr_t *
_pick_step_nodes (struct job_record  *job_ptr, 
		job_step_create_request_msg_t *step_spec )
{

	bitstr_t *nodes_avail = NULL, *nodes_idle = NULL;
	bitstr_t *nodes_picked = NULL, *node_tmp = NULL;
	int error_code, nodes_picked_cnt = 0, cpus_picked_cnt, i;

	if (job_ptr->node_bitmap == NULL)
		return NULL;
	
	nodes_avail = bit_copy (job_ptr->node_bitmap);
	if (nodes_avail == NULL)
		fatal("bit_copy malloc failure");
	bit_and (nodes_avail, avail_node_bitmap);

	if ( step_spec->node_count == INFINITE)	/* use all nodes */
		return nodes_avail;
try_again:
	if (step_spec->node_list) {
		error_code = node_name2bitmap (step_spec->node_list, false, 
						&nodes_picked);
		if (error_code) {
			info ("_pick_step_nodes: invalid node list %s", 
				step_spec->node_list);
			goto cleanup;
		}
		if (bit_super_set (nodes_picked, job_ptr->node_bitmap) == 0) {
			info ("_pick_step_nodes: requested nodes %s not part "
				"of job %u", 
				step_spec->node_list, job_ptr->job_id);
			goto cleanup;
		}
		if(step_spec->task_dist == SLURM_DIST_ARBITRARY) {
			if (!strcmp(slurmctld_conf.switch_type,
				    "switch/elan")) {
				error("Can't do an ARBITRARY task layout with "
				      "switch type elan. Switching DIST type "
				      "to BLOCK");
				xfree(step_spec->node_list);
				step_spec->task_dist == SLURM_DIST_BLOCK;
				FREE_NULL_BITMAP(nodes_picked);
				goto try_again;
			}
			FREE_NULL_BITMAP(nodes_avail);
			return nodes_picked;
		}
	} else if (step_spec->relative) {
		/* Remove first (step_spec->relative) nodes from  
		 * available list */
		bitstr_t *relative_nodes = NULL;
		relative_nodes = 
			bit_pick_cnt (nodes_avail, step_spec->relative);
		if (relative_nodes == NULL) {
			info ("_pick_step_nodes: "
			      "Invalid relative value (%u) for job %u",
			      step_spec->relative, job_ptr->job_id);
			goto cleanup;
		}
		bit_not (relative_nodes);
		bit_and (nodes_avail, relative_nodes);
		bit_free (relative_nodes);
	} else {
		ListIterator step_iterator;
		struct step_record *step_p;
		nodes_picked = bit_alloc (bit_size (nodes_avail) );
		nodes_idle = bit_alloc (bit_size (nodes_avail) );
		if ((nodes_picked == NULL) || (nodes_idle == NULL))
			fatal("bit_alloc malloc failure");
		step_iterator = list_iterator_create (job_ptr->step_list);
		while ((step_p = (struct step_record *) list_next (step_iterator)))
			bit_or (nodes_idle, step_p->step_node_bitmap); 
		list_iterator_destroy (step_iterator);
		bit_not(nodes_idle);
		bit_and(nodes_idle, nodes_avail);
	}

	/* if user specifies step needs a specific processor count and 
	 * all nodes have the same processor count, just translate this to
	 * a node count */
	if (step_spec->cpu_count && (job_ptr->num_cpu_groups == 1)) {
		i = (step_spec->cpu_count + (job_ptr->cpus_per_node[0] - 1) ) 
				/ job_ptr->cpus_per_node[0];
		step_spec->node_count = (i > step_spec->node_count) ? 
						i : step_spec->node_count ;
		step_spec->cpu_count = 0;
	}

	if (step_spec->node_count) {
		nodes_picked_cnt = bit_set_count(nodes_picked);
		if (nodes_idle 
		&&  (step_spec->node_count > nodes_picked_cnt)) {
			node_tmp = bit_pick_cnt(nodes_idle,
					(step_spec->node_count -
					nodes_picked_cnt));
			if (node_tmp == NULL)
				goto cleanup;
			bit_or  (nodes_picked, node_tmp);
			bit_not (node_tmp);
			bit_and (nodes_idle, node_tmp);
			bit_and (nodes_avail, node_tmp);
			bit_free (node_tmp);
			node_tmp = NULL;
			nodes_picked_cnt = step_spec->node_count;
		}
		if (step_spec->node_count > nodes_picked_cnt) {
			node_tmp = bit_pick_cnt(nodes_avail, 
					(step_spec->node_count - 
					nodes_picked_cnt));
			if (node_tmp == NULL)
				goto cleanup;
			bit_or  (nodes_picked, node_tmp);
			bit_not (node_tmp);
			bit_and (nodes_avail, node_tmp);
			bit_free (node_tmp);
			node_tmp = NULL;
			nodes_picked_cnt = step_spec->node_count;
		}
	}

	if (step_spec->cpu_count) {
		cpus_picked_cnt = count_cpus(nodes_picked);
		if (nodes_idle
		&&  (step_spec->cpu_count > cpus_picked_cnt)) {
			int first_bit, last_bit;
			first_bit = bit_ffs(nodes_idle);
			last_bit  = bit_fls(nodes_idle);
			for (i = first_bit; i <= last_bit; i++) {
				if (bit_test (nodes_idle, i) != 1)
					continue;
				bit_set (nodes_picked, i);
				bit_clear (nodes_avail, i);
				/* bit_clear (nodes_idle, i);	unused */
				cpus_picked_cnt +=
					node_record_table_ptr[i].cpus;
				if (cpus_picked_cnt >= step_spec->cpu_count)
					break;
			}
			if (step_spec->cpu_count > cpus_picked_cnt)
				goto cleanup;
		}
		if (step_spec->cpu_count > cpus_picked_cnt) {
			int first_bit, last_bit;
			first_bit = bit_ffs(nodes_avail);
			last_bit  = bit_fls(nodes_avail);
 			for (i = first_bit; i <= last_bit; i++) {
				if (bit_test (nodes_avail, i) != 1)
					continue;
				bit_set (nodes_picked, i);
				cpus_picked_cnt += 
					node_record_table_ptr[i].cpus;
				if (cpus_picked_cnt >= step_spec->cpu_count)
					break;
			}
			if (step_spec->cpu_count > cpus_picked_cnt)
				goto cleanup;
		}
	}

	FREE_NULL_BITMAP(nodes_avail);
	FREE_NULL_BITMAP(nodes_idle);
	return nodes_picked;

cleanup:
	FREE_NULL_BITMAP(nodes_avail);
	FREE_NULL_BITMAP(nodes_idle);
	FREE_NULL_BITMAP(nodes_picked);
	return NULL;
}


/*
 * step_create - creates a step_record in step_specs->job_id, sets up the
 *	according to the step_specs.
 * IN step_specs - job step specifications
 * OUT new_step_record - pointer to the new step_record (NULL on error)
 * IN kill_job_when_step_done - if set kill the job on step completion
 * IN batch_step - if set then step is a batch script
 * RET - 0 or error code
 * NOTE: don't free the returned step_record because that is managed through
 * 	the job.
 */
extern int
step_create(job_step_create_request_msg_t *step_specs, 
	    struct step_record** new_step_record,
	    bool kill_job_when_step_done, bool batch_step)
{
	struct step_record *step_ptr;
	struct job_record  *job_ptr;
	bitstr_t *nodeset;
	int node_count;
	time_t now = time(NULL);

	*new_step_record = NULL;
	job_ptr = find_job_record (step_specs->job_id);
	if (job_ptr == NULL)
		return ESLURM_INVALID_JOB_ID ;

	if ((step_specs->user_id != job_ptr->user_id) &&
	    (step_specs->user_id != 0))
		return ESLURM_ACCESS_DENIED ;

	if (IS_JOB_PENDING(job_ptr))
		return ESLURM_INVALID_JOB_ID ;

	if (IS_JOB_FINISHED(job_ptr) || 
	    (job_ptr->end_time <= time(NULL)))
		return ESLURM_ALREADY_DONE;

	if ((step_specs->task_dist != SLURM_DIST_CYCLIC) &&
	    (step_specs->task_dist != SLURM_DIST_BLOCK) &&
	    (step_specs->task_dist != SLURM_DIST_ARBITRARY))
		return ESLURM_BAD_DIST;

	if (job_ptr->kill_on_step_done)
		/* Don't start more steps, job already being cancelled */
		return ESLURM_ALREADY_DONE;
	job_ptr->kill_on_step_done = kill_job_when_step_done;

	job_ptr->time_last_active = now;
	nodeset = _pick_step_nodes (job_ptr, step_specs);
	if (nodeset == NULL)
		return ESLURM_REQUESTED_NODE_CONFIG_UNAVAILABLE ;
	node_count = bit_set_count(nodeset);
	
	if (step_specs->num_tasks == NO_VAL) {
		if (step_specs->cpu_count != NO_VAL)
			step_specs->num_tasks = step_specs->cpu_count;
		else
			step_specs->num_tasks = node_count;
	}
	
	if ((step_specs->num_tasks < 1)
	|| (step_specs->num_tasks > (node_count*MAX_TASKS_PER_NODE)))
		return ESLURM_BAD_TASK_COUNT;

	step_ptr = create_step_record (job_ptr);
	if (step_ptr == NULL)
		fatal ("create_step_record failed with no memory");

	/* set the step_record values */
	/* Here is where the node list is set for the job */
	if(step_specs->node_list)
		step_ptr->step_node_list = xstrdup(step_specs->node_list);
	else
		step_ptr->step_node_list = bitmap2node_name(nodeset);
	xfree(step_specs->node_list);
	step_specs->node_list = bitmap2node_name(nodeset);
	step_ptr->step_node_bitmap = nodeset;
	step_ptr->cyclic_alloc = 
		(uint16_t) (step_specs->task_dist == SLURM_DIST_CYCLIC);
	step_ptr->num_tasks = step_specs->num_tasks;
	step_ptr->num_cpus = step_specs->cpu_count;
	step_ptr->time_last_active = now;
	step_ptr->port = step_specs->port;
	step_ptr->host = xstrdup(step_specs->host);
	step_ptr->batch_step = batch_step;
	step_ptr->exit_code = NO_VAL;

	/* step's name and network default to job's values if not 
	 * specified in the step specification */
	if (step_specs->name && step_specs->name[0])
		step_ptr->name = xstrdup(step_specs->name);
	else
		step_ptr->name = xstrdup(job_ptr->name);
	if (step_specs->network && step_specs->network[0])
		step_ptr->network = xstrdup(step_specs->network);
	else
		step_ptr->network = xstrdup(job_ptr->network);

	/* a batch script does not need switch info */
	if (!batch_step) {
		uint32_t *tasks_per_node;
		tasks_per_node = distribute_tasks(job_ptr->nodes,
						  job_ptr->num_cpu_groups,
						  job_ptr->cpus_per_node,
						  job_ptr->cpu_count_reps,
						  step_ptr->step_node_list,
						  step_ptr->num_tasks);

		if (switch_alloc_jobinfo (&step_ptr->switch_job) < 0)
			fatal ("step_create: switch_alloc_jobinfo error");
		
		if (switch_build_jobinfo(step_ptr->switch_job, 
					 step_ptr->step_node_list,
					 tasks_per_node, 
					 step_ptr->cyclic_alloc,
					 step_ptr->network) < 0) {
			error("switch_build_jobinfo: %m");
			xfree(tasks_per_node);
			delete_step_record (job_ptr, step_ptr->step_id);
			return ESLURM_INTERCONNECT_FAILURE;
		}
		xfree(tasks_per_node);
	}
	if (checkpoint_alloc_jobinfo (&step_ptr->check_job) < 0)
		fatal ("step_create: checkpoint_alloc_jobinfo error");

	*new_step_record = step_ptr;
	jobacct_step_start(step_ptr);
	return SLURM_SUCCESS;
}

/* Pack the data for a specific job step record
 * IN step - pointer to a job step record
 * IN/OUT buffer - location to store data, pointers automatically advanced
 */
static void _pack_ctld_job_step_info(struct step_record *step, Buf buffer)
{
	pack_job_step_info_members(step->job_ptr->job_id,
				   step->step_id,
				   step->job_ptr->user_id,
				   step->num_tasks,
				   step->start_time,
				   step->job_ptr->partition,
				   step->step_node_list, 
				   step->name, step->network, buffer);
}

/* 
 * pack_ctld_job_step_info_response_msg - packs job step info
 * IN job_id - specific id or zero for all
 * IN step_id - specific id or zero for all
 * IN uid - user issuing request
 * IN show_flags - job step filtering options
 * OUT buffer - location to store data, pointers automatically advanced 
 * RET - 0 or error code
 * NOTE: MUST free_buf buffer
 */
extern int pack_ctld_job_step_info_response_msg(uint32_t job_id, 
			uint32_t step_id, uid_t uid, 
			uint16_t show_flags, Buf buffer)
{
	ListIterator job_iterator;
	ListIterator step_iterator;
	int error_code = 0;
	uint32_t steps_packed = 0, tmp_offset;
	struct step_record *step_ptr;
	struct job_record *job_ptr;
	time_t now = time(NULL);

	pack_time(now, buffer);
	pack32(steps_packed, buffer);	/* steps_packed placeholder */

	part_filter_set(uid);
	if (job_id == 0) {
		/* Return all steps for all jobs */
		job_iterator = list_iterator_create(job_list);
		while ((job_ptr = 
				(struct job_record *) 
				list_next(job_iterator))) {
			if (((show_flags & SHOW_ALL) == 0) && 
			    (job_ptr->part_ptr) && 
			    (job_ptr->part_ptr->hidden))
				continue;
			step_iterator =
			    list_iterator_create(job_ptr->step_list);
			while ((step_ptr =
					(struct step_record *)
					list_next(step_iterator))) {
				_pack_ctld_job_step_info(step_ptr, buffer);
				steps_packed++;
			}
			list_iterator_destroy(step_iterator);
		}
		list_iterator_destroy(job_iterator);

	} else if (step_id == 0) {
		/* Return all steps for specific job_id */
		job_ptr = find_job_record(job_id);
		if (((show_flags & SHOW_ALL) == 0) && 
		    (job_ptr->part_ptr) && 
		    (job_ptr->part_ptr->hidden))
			job_ptr = NULL;
		if (job_ptr) {
			step_iterator = 
				list_iterator_create(job_ptr->step_list);
			while ((step_ptr =
					(struct step_record *)
					list_next(step_iterator))) {
				_pack_ctld_job_step_info(step_ptr, buffer);
				steps_packed++;
			}
			list_iterator_destroy(step_iterator);
		} else
			error_code = ESLURM_INVALID_JOB_ID;
	} else {
		/* Return data for specific job_id.step_id */
		job_ptr = find_job_record(job_id);
		if (((show_flags & SHOW_ALL) == 0) && 
		    (job_ptr->part_ptr) && 
		    (job_ptr->part_ptr->hidden))
			job_ptr = NULL;
		step_ptr = find_step_record(job_ptr, step_id);
		if (step_ptr == NULL)
			error_code = ESLURM_INVALID_JOB_ID;
		else {
			_pack_ctld_job_step_info(step_ptr, buffer);
			steps_packed++;
		}
	}
	part_filter_clear();

	/* put the real record count in the message body header */
	tmp_offset = get_buf_offset(buffer);
	set_buf_offset(buffer, 0);
	pack_time(now, buffer);
	pack32(steps_packed, buffer);
	set_buf_offset(buffer, tmp_offset);

	return error_code;
}

/* 
 * step_on_node - determine if the specified job has any job steps allocated to 
 * 	the specified node 
 * IN job_ptr - pointer to an active job record
 * IN node_ptr - pointer to a node record
 * RET true of job has step on the node, false otherwise 
 */
bool step_on_node(struct job_record  *job_ptr, struct node_record *node_ptr)
{
	ListIterator step_iterator;
	struct step_record *step_ptr;
	bool found = false;
	int bit_position;

	if ((job_ptr == NULL) || (node_ptr == NULL))
		return false;

	bit_position = node_ptr - node_record_table_ptr;
	step_iterator = list_iterator_create (job_ptr->step_list);	
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		if (bit_test(step_ptr->step_node_bitmap, bit_position)) {
			found = true;
			break;
		}
	}		

	list_iterator_destroy (step_iterator);
	return found;
}

/*
 * job_step_checkpoint - perform some checkpoint operation
 * IN ckpt_ptr - checkpoint request message 
 * IN uid - user id of the user issuing the RPC
 * IN conn_fd - file descriptor on which to send reply
 * RET 0 on success, otherwise ESLURM error code
 */
extern int job_step_checkpoint(checkpoint_msg_t *ckpt_ptr,
		uid_t uid, slurm_fd conn_fd)
{
	int rc = SLURM_SUCCESS;
	struct job_record *job_ptr;
	struct step_record *step_ptr;
	checkpoint_resp_msg_t resp_data;
	slurm_msg_t resp_msg;

	resp_msg.forward.cnt = 0;
	resp_msg.ret_list = NULL;

	/* find the job */
	job_ptr = find_job_record (ckpt_ptr->job_id);
	if (job_ptr == NULL) {
		rc = ESLURM_INVALID_JOB_ID;
		goto reply;
	}
	if ((uid != job_ptr->user_id) && (uid != 0)) {
		rc = ESLURM_ACCESS_DENIED ;
		goto reply;
	}
	if (job_ptr->job_state == JOB_PENDING) {
		rc = ESLURM_JOB_PENDING;
		goto reply;
	} else if (job_ptr->job_state == JOB_SUSPENDED) {
		/* job can't get cycles for checkpoint 
		 * if it is already suspended */
		rc = ESLURM_DISABLED;
		goto reply;
	} else if (job_ptr->job_state != JOB_RUNNING) {
		rc = ESLURM_ALREADY_DONE;
		goto reply;
	}

	bzero((void *)&resp_data, sizeof(checkpoint_resp_msg_t));
	/* find the individual job step */
	if (ckpt_ptr->step_id != NO_VAL) {
		step_ptr = find_step_record(job_ptr, ckpt_ptr->step_id);
		if (step_ptr == NULL) {
			rc = ESLURM_INVALID_JOB_ID;
			goto reply;
		} else {
			rc = checkpoint_op(ckpt_ptr->op, ckpt_ptr->data, 
				(void *)step_ptr, &resp_data.event_time, 
				&resp_data.error_code, &resp_data.error_msg);
			last_job_update = time(NULL);
		}
	}

	/* operate on all of a job's steps */
	else {
		int update_rc = -2;
		bool error_reply = false;
		ListIterator step_iterator;

		step_iterator = list_iterator_create (job_ptr->step_list);
		while ((step_ptr = (struct step_record *) 
					list_next (step_iterator))) {
			update_rc = checkpoint_op(ckpt_ptr->op, 
						  ckpt_ptr->data,
						  (void *)step_ptr,
						  &resp_data.event_time,
						  &resp_data.error_code,
						  &resp_data.error_msg);
			rc = MAX(rc, update_rc);
		}
		if (update_rc != -2)	/* some work done */
			last_job_update = time(NULL);
		list_iterator_destroy (step_iterator);
	}

    reply:
	if ((rc == SLURM_SUCCESS) &&
	    ((ckpt_ptr->op == CHECK_ABLE) || (ckpt_ptr->op == CHECK_ERROR))) {
		resp_msg.msg_type = RESPONSE_CHECKPOINT;
		resp_msg.data = &resp_data;
		 (void) slurm_send_node_msg(conn_fd, &resp_msg);
	} else {
		return_code_msg_t rc_msg;
		rc_msg.return_code = rc;
		resp_msg.msg_type  = RESPONSE_SLURM_RC;
		resp_msg.data      = &rc_msg;
		(void) slurm_send_node_msg(conn_fd, &resp_msg);
	}
	return rc;
}

/*
 * job_step_checkpoint_comp - note job step checkpoint completion
 * IN ckpt_ptr - checkpoint complete status message
 * IN uid - user id of the user issuing the RPC
 * IN conn_fd - file descriptor on which to send reply
 * RET 0 on success, otherwise ESLURM error code
 */
extern int job_step_checkpoint_comp(checkpoint_comp_msg_t *ckpt_ptr,
		uid_t uid, slurm_fd conn_fd)
{
	int rc = SLURM_SUCCESS;
	struct job_record *job_ptr;
	struct step_record *step_ptr;
	slurm_msg_t resp_msg;
	return_code_msg_t rc_msg;
	
	resp_msg.forward.cnt = 0;
	resp_msg.ret_list = NULL;

	/* find the job */
	job_ptr = find_job_record (ckpt_ptr->job_id);
	if (job_ptr == NULL) {
		rc = ESLURM_INVALID_JOB_ID;
		goto reply;
	}
	if ((uid != job_ptr->user_id) && (uid != 0)) {
		rc = ESLURM_ACCESS_DENIED;
		goto reply;
	}
	if (job_ptr->job_state == JOB_PENDING) {
		rc = ESLURM_JOB_PENDING;
		goto reply;
	} else if ((job_ptr->job_state != JOB_RUNNING)
	&&         (job_ptr->job_state != JOB_SUSPENDED)) {
		rc = ESLURM_ALREADY_DONE;
		goto reply;
	}
 
	step_ptr = find_step_record(job_ptr, ckpt_ptr->step_id);
	if (step_ptr == NULL) {
		rc = ESLURM_INVALID_JOB_ID;
		goto reply;
	} else {
		rc = checkpoint_comp((void *)step_ptr, ckpt_ptr->begin_time, 
			ckpt_ptr->error_code, ckpt_ptr->error_msg);
		last_job_update = time(NULL);
	}

    reply:
	rc_msg.return_code = rc;
	resp_msg.msg_type  = RESPONSE_SLURM_RC;
	resp_msg.data      = &rc_msg;
	(void) slurm_send_node_msg(conn_fd, &resp_msg);
	return rc;
}

/*
 * step_partial_comp - Note the completion of a job step on at least
 *	some of its nodes
 * IN req     - step_completion_msg RPC from slurmstepd
 * OUT rem    - count of nodes for which responses are still pending
 * OUT max_rc - highest return code for any step thus far
 * RET 0 on success, otherwise ESLURM error code
 */
extern int step_partial_comp(step_complete_msg_t *req, int *rem, 
		int *max_rc)
{
	struct job_record *job_ptr;
	struct step_record *step_ptr;
	int nodes;

	/* find the job, step, and validate input */
	job_ptr = find_job_record (req->job_id);
	if (job_ptr == NULL)
		return ESLURM_INVALID_JOB_ID;
	if (job_ptr->job_state == JOB_PENDING)
		return ESLURM_JOB_PENDING;
	step_ptr = find_step_record(job_ptr, req->job_step_id);
	if (step_ptr == NULL)
		return ESLURM_INVALID_JOB_ID;
	if (req->range_last < req->range_first)
		return EINVAL;

	aggregate_step_data(step_ptr, 
			    req->rusage, req->max_psize, req->max_vsize);

	if (step_ptr->exit_code == NO_VAL) {
		/* initialize the node bitmap for exited nodes */
		nodes = bit_set_count(step_ptr->step_node_bitmap);
		if (req->range_last >= nodes)	/* range is zero origin */
			return EINVAL;
		xassert(step_ptr->exit_node_bitmap == NULL);
		step_ptr->exit_node_bitmap = bit_alloc(nodes);
		if (step_ptr->exit_node_bitmap == NULL)
			fatal("bit_alloc: %m");
		step_ptr->exit_code = req->step_rc;
	} else {
		xassert(step_ptr->exit_node_bitmap);
		nodes = _bitstr_bits(step_ptr->exit_node_bitmap);
		if (req->range_last >= nodes)	/* range is zero origin */
			return EINVAL;
		step_ptr->exit_code = MAX(step_ptr->exit_code, req->step_rc);
	}

	bit_nset(step_ptr->exit_node_bitmap, req->range_first,
		req->range_last);
	if (rem)
		*rem = bit_clear_count(step_ptr->exit_node_bitmap);
	if (max_rc)
		*max_rc = step_ptr->exit_code;

	return SLURM_SUCCESS;
}

