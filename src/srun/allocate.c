/*****************************************************************************\
 * src/srun/allocate.c - srun functions for managing node allocations
 * $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <mgrondona@llnl.gov>.
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
 *  In addition, as a special exception, the copyright holders give permission 
 *  to link the code of portions of this program with the OpenSSL library under 
 *  certain conditions as described in each individual source file, and 
 *  distribute linked combinations including the two. You must obey the GNU 
 *  General Public License in all respects for all of the code used other than 
 *  OpenSSL. If you modify file(s) with this exception, you may extend this 
 *  exception to your version of the file(s), but you are not obligated to do 
 *  so. If you do not wish to do so, delete this exception statement from your
 *  version.  If you delete this exception statement from all source files in 
 *  the program, then also delete it here.
 *  
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA.
\*****************************************************************************/

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdlib.h>
#include <sys/poll.h>


#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/slurm_auth.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"
#include "src/common/xsignal.h"
#include "src/common/xstring.h"
#include "src/common/forward.h"

#include "src/srun/allocate.h"
#include "src/srun/msg.h"
#include "src/srun/opt.h"
#include "src/srun/attach.h"

#define MAX_ALLOC_WAIT 60	/* seconds */
#define MIN_ALLOC_WAIT  5	/* seconds */
#define MAX_RETRIES    10

/*
 * Static Prototypes
 */
static int   _accept_msg_connection(slurm_fd slurmctld_fd,
		resource_allocation_response_msg_t **resp);
static int   _handle_msg(slurm_msg_t *msg, \
		resource_allocation_response_msg_t **resp);
static int   _wait_for_alloc_rpc(int sleep_time,
		resource_allocation_response_msg_t **resp);
static void  _wait_for_resources(resource_allocation_response_msg_t **resp);
static bool  _retry();
static void  _intr_handler(int signo);

static job_step_create_request_msg_t * _step_req_create(srun_job_t *j);

static sig_atomic_t destroy_job = 0;
static srun_job_t *allocate_job = NULL;

int
allocate_test(void)
{
	int rc;
	job_desc_msg_t *j = job_desc_msg_create_from_opts (NULL);
	if(!j)
		return SLURM_ERROR;
	
	rc = slurm_job_will_run(j);
	job_desc_msg_destroy(j);
	return rc;
}

resource_allocation_response_msg_t *
allocate_nodes(void)
{
	int rc = 0;
	static int sigarray[] = { SIGQUIT, SIGINT, SIGTERM, 0 };
	SigFunc *oquitf, *ointf, *otermf;
	sigset_t oset;
	resource_allocation_response_msg_t *resp = NULL;
	job_desc_msg_t *j = job_desc_msg_create_from_opts (NULL);
	if(!j)
		return NULL;
	
	oquitf = xsignal(SIGQUIT, _intr_handler);
	ointf  = xsignal(SIGINT,  _intr_handler);
	otermf = xsignal(SIGTERM, _intr_handler);

	xsignal_save_mask(&oset);
	xsignal_unblock(sigarray);

	/* Do not re-use existing job id when submitting new job
	 * from within a running job */
	if (getenv("SLURM_JOBID") != NULL) {
		info("WARNING: Creating SLURM job allocation from within "
			"another allocation");
		info("WARNING: You are attempting to initiate a second job");
		if (!opt.jobid_set)	/* Let slurmctld set jobid */
			j->job_id = NO_VAL;
	}
	
	while ((rc = slurm_allocate_resources(j, &resp) < 0) && _retry()) {
		if (destroy_job)
			goto done;
	} 

	if(!resp)
		goto done;
	
	if ((rc == 0) && (resp->node_list == NULL)) {
		if (resp->error_code)
			verbose("Warning: %s", slurm_strerror(resp->error_code));
		_wait_for_resources(&resp);
	}
	/* For diagnosing a node problem, administrators need to sometimes
	 * run a job on N nodes one of which must be the node believed to 
	 * have a problem (e.g. "srun -N4 -w bad_node diagnostic"). The 
	 * below logic prevents this from working and necessiates the 
	 * admin identify four specific nodes to use for the above test
	 * instead of just the one bad node. Otherwise only the one 
	 * bad node is used in the job's allocation. */
	if(resp->node_list && j->req_nodes) {
		xfree(resp->node_list);
		resp->node_list = xstrdup(j->req_nodes);
	}
		
    done:
	xsignal_set_mask(&oset);
	xsignal(SIGINT,  ointf);
	xsignal(SIGTERM, otermf);
	xsignal(SIGQUIT, oquitf);

	job_desc_msg_destroy(j);

	return resp;
}

/* 
 * Returns jobid if SLURM_JOBID was set in the user's environment
 *  or if --jobid option was given, else returns 0
 */
uint32_t
jobid_from_env(void)
{
	if (opt.jobid != NO_VAL)
		return ((uint32_t) opt.jobid);
	else 
		return (0);
}

static void
_wait_for_resources(resource_allocation_response_msg_t **resp)
{
	resource_allocation_response_msg_t *r = *resp;
	int sleep_time = MIN_ALLOC_WAIT;
	int job_id = r->job_id;

	if (!opt.quiet)
		info ("job %u queued and waiting for resources", r->job_id);

	slurm_free_resource_allocation_response_msg(r);

	/* Keep polling until the job is allocated resources */
	while (_wait_for_alloc_rpc(sleep_time, resp) <= 0) {

		if (slurm_allocation_lookup_lite(job_id, resp) >= 0)
			break;

		if (slurm_get_errno() == ESLURM_JOB_PENDING) 
			debug3 ("Still waiting for allocation");
		else 
			fatal ("Unable to confirm allocation for job %u: %m", 
			       job_id);

		if (destroy_job) {
			verbose("cancelling job %u", job_id);
			slurm_complete_job(job_id, 0);
			debugger_launch_failure(allocate_job);
			exit(0);
		}

		if (sleep_time < MAX_ALLOC_WAIT)
			sleep_time++;
	}
	if (!opt.quiet)
		info ("job %u has been allocated resources", (*resp)->job_id);
}

/* Wait up to sleep_time for RPC from slurmctld indicating resource allocation
 * has occured.
 * IN sleep_time: delay in seconds
 * OUT resp: resource allocation response message
 * RET 1 if resp is filled in, 0 otherwise */
static int
_wait_for_alloc_rpc(int sleep_time, resource_allocation_response_msg_t **resp)
{
	struct pollfd fds[1];
	slurm_fd slurmctld_fd;

	if ((slurmctld_fd = slurmctld_msg_init()) < 0) {
		sleep (sleep_time);
		return (0);
	}

	fds[0].fd = slurmctld_fd;
	fds[0].events = POLLIN;

	while (poll (fds, 1, (sleep_time * 1000)) < 0) {
		switch (errno) {
			case EAGAIN:
			case EINTR:
				return (-1);
			case ENOMEM:
			case EINVAL:
			case EFAULT:
				fatal("poll: %m");
			default:
				error("poll: %m. Continuing...");
		}
	}

	if (fds[0].revents & POLLIN)
		return (_accept_msg_connection(slurmctld_fd, resp));

	return (0);
}

/* Accept RPC from slurmctld and process it.
 * IN slurmctld_fd: file descriptor for slurmctld communications
 * OUT resp: resource allocation response message
 * RET 1 if resp is filled in, 0 otherwise */
static int 
_accept_msg_connection(slurm_fd slurmctld_fd, 
		resource_allocation_response_msg_t **resp)
{
	slurm_fd     fd;
	slurm_msg_t *msg = NULL;
	slurm_addr   cli_addr;
	char         host[256];
	uint16_t     port;
	int          rc = 0;
	List ret_list;

	fd = slurm_accept_msg_conn(slurmctld_fd, &cli_addr);
	if (fd < 0) {
		error("Unable to accept connection: %m");
		return rc;
	}

	slurm_get_addr(&cli_addr, &port, host, sizeof(host));
	debug2("got message connection from %s:%d", host, port);

	msg = xmalloc(sizeof(slurm_msg_t));
	slurm_init_slurm_msg(msg, NULL);
	msg->conn_fd = fd;
		
  again:
	ret_list = slurm_receive_msg(fd, msg, 0);

	if(!ret_list || errno != SLURM_SUCCESS) {
		if (errno == EINTR) {
			goto again;
		}
		if(ret_list)
			list_destroy(ret_list);
			
		error("_accept_msg_connection[%s]: %m", host);
		slurm_free_msg(msg);
		return SLURM_ERROR;
	}
	if(list_count(ret_list)>0) {
		error("_accept_msg_connection: "
		      "got %d from receive, expecting 0",
		      list_count(ret_list));
	}
	msg->ret_list = ret_list;
	
	
	rc = _handle_msg(msg, resp); /* handle_msg frees msg */
	slurm_free_msg(msg);
		
	slurm_close_accepted_conn(fd);
	return rc;
}

/* process RPC from slurmctld
 * IN msg: message recieved
 * OUT resp: resource allocation response message
 * RET 1 if resp is filled in, 0 otherwise */
static int
_handle_msg(slurm_msg_t *msg, resource_allocation_response_msg_t **resp)
{
	uid_t req_uid   = g_slurm_auth_get_uid(msg->auth_cred);
	uid_t uid       = getuid();
	uid_t slurm_uid = (uid_t) slurm_get_slurm_user_id();
	int rc = 0;
	srun_timeout_msg_t *to;

	if ((req_uid != slurm_uid) && (req_uid != 0) && (req_uid != uid)) {
		error ("Security violation, slurm message from uid %u",
			(unsigned int) req_uid);
		return 0;
	}

	switch (msg->msg_type) {
		case SRUN_PING:
			debug3("slurmctld ping received");
			slurm_send_rc_msg(msg, SLURM_SUCCESS);
			slurm_free_srun_ping_msg(msg->data);
			break;
		case SRUN_JOB_COMPLETE:
			debug3("job complete received");
			/* FIXME: do something here */
			slurm_free_srun_job_complete_msg(msg->data);	
			break;
		case RESPONSE_RESOURCE_ALLOCATION:
			debug2("resource allocation response received");
			slurm_send_rc_msg(msg, SLURM_SUCCESS);
			*resp = msg->data;
			rc = 1;
			break;
		case SRUN_TIMEOUT:
			debug2("timeout received");
			to = msg->data;
			timeout_handler(to->timeout);
			slurm_free_srun_timeout_msg(msg->data);
			break;
		default:
			error("received spurious message type: %d\n",
				 msg->msg_type);
	}
	return rc;
}

static bool
_retry()
{
	static int  retries = 0;
	static char *msg = "Slurm controller not responding, "
		           "sleeping and retrying.";

	if (errno == ESLURM_ERROR_ON_DESC_TO_RECORD_COPY) {
		if (retries == 0)
			error (msg);
		else if (retries < MAX_RETRIES)
			debug (msg);
		else
			return false;
		sleep (++retries);
	} else {
		error("Unable to allocate resources: %m");
		return false;
	}

	return true;
}

/*
 * SIGINT handler while waiting for resources to become available.
 */
static void
_intr_handler(int signo)
{
	destroy_job = 1;
}


/*
 * Create job description structure based off srun options
 * (see opt.h)
 */
job_desc_msg_t *
job_desc_msg_create_from_opts (char *script)
{
	extern char **environ;
	job_desc_msg_t *j = xmalloc(sizeof(*j));
	char buf[8192];
	hostlist_t hl = NULL;
	
	slurm_init_job_desc_msg(j);
	
	j->contiguous     = opt.contiguous;
	j->features       = opt.constraints;
	j->immediate      = opt.immediate;
	j->name           = opt.job_name;
	j->req_nodes      = xstrdup(opt.nodelist);
	if (j->req_nodes == NULL) {
		char *nodelist = NULL;
		char *hostfile = getenv("SLURM_HOSTFILE");
		
		if (hostfile != NULL) {
			nodelist = slurm_read_hostfile(hostfile, opt.nprocs);
			if (nodelist == NULL) {
				error("Failure getting NodeNames from "
				      "hostfile");
				/* FIXME - need to fail somehow */
			} else {
				debug("loading nodes from hostfile %s",
				      hostfile);
				j->req_nodes = xstrdup(nodelist);
				free(nodelist);
				opt.distribution = SLURM_DIST_ARBITRARY;
			}
		}
	}
	/* simplify the job allocation nodelist, 
	  not laying out tasks until step */
	if(j->req_nodes) {
		hl = hostlist_create(j->req_nodes);
		hostlist_ranged_string(hl, sizeof(buf), buf);
		xfree(opt.nodelist);
		opt.nodelist = xstrdup(buf);
		hostlist_uniq(hl);
		hostlist_ranged_string(hl, sizeof(buf), buf);
		hostlist_destroy(hl);

		xfree(j->req_nodes);
		j->req_nodes = xstrdup(buf);
	}
	
	if(opt.distribution == SLURM_DIST_ARBITRARY
	   && !j->req_nodes) {
		error("With Arbitrary distribution you need to "
		      "specify a nodelist or hostfile with the -w option");
		return NULL;
	}
	j->exc_nodes      = opt.exc_nodes;
	j->partition      = opt.partition;
	j->min_nodes      = opt.min_nodes;
	j->min_sockets    = opt.min_sockets_per_node;
	j->min_cores      = opt.min_cores_per_socket;
	j->min_threads    = opt.min_threads_per_core;
	j->user_id        = opt.uid;
	j->dependency     = opt.dependency;
	if (opt.nice)
		j->nice   = NICE_OFFSET + opt.nice;
	j->task_dist      = opt.distribution;
	j->plane_size     = opt.plane_size;
	j->group_id       = opt.gid;
	j->mail_type      = opt.mail_type;
	if (opt.mail_user)
		j->mail_user = xstrdup(opt.mail_user);
	if (opt.begin)
		j->begin_time = opt.begin;
	if (opt.network)
		j->network = xstrdup(opt.network);
	if (opt.account)
		j->account = xstrdup(opt.account);

	if (opt.hold)
		j->priority     = 0;
	if (opt.jobid != NO_VAL)
		j->job_id	= opt.jobid;
#if SYSTEM_DIMENSIONS
	if (opt.geometry[0] > 0) {
		int i;
		for (i=0; i<SYSTEM_DIMENSIONS; i++)
			j->geometry[i] = opt.geometry[i];
	}
#endif

	if (opt.conn_type != -1)
		j->conn_type = opt.conn_type;
			
	if (opt.no_rotate)
		j->rotate = 0;

	if (opt.max_nodes)
		j->max_nodes    = opt.max_nodes;
	if (opt.max_sockets_per_node)
		j->max_sockets  = opt.max_sockets_per_node;
	if (opt.max_cores_per_socket)
		j->max_cores    = opt.max_cores_per_socket;
	if (opt.max_threads_per_core)
		j->max_threads  = opt.max_threads_per_core;

	if (opt.job_min_cpus > -1)
		j->job_min_procs    = opt.job_min_cpus;
	if (opt.job_min_sockets > -1)
		j->job_min_sockets  = opt.job_min_sockets;
	if (opt.job_min_cores > -1)
		j->job_min_cores    = opt.job_min_cores;
	if (opt.job_min_threads > -1)
		j->job_min_threads  = opt.job_min_threads;
	if (opt.job_min_memory > -1)
		j->job_min_memory   = opt.job_min_memory;
	if (opt.job_min_tmp_disk > -1)
		j->job_min_tmp_disk = opt.job_min_tmp_disk;
	if (opt.overcommit) {
		j->num_procs    = opt.min_nodes;
		j->overcommit	= opt.overcommit;
	} else
		j->num_procs    = opt.nprocs * opt.cpus_per_task;
	if (opt.nprocs_set)
		j->num_tasks    = opt.nprocs;

	if (opt.cpus_set)
		j->cpus_per_task = opt.cpus_per_task;

	if (opt.no_kill)
		j->kill_on_node_fail   = 0;
	if (opt.time_limit > -1)
		j->time_limit          = opt.time_limit;
	j->shared = opt.shared;

	/* srun uses the same listening port for the allocation response
	 * message as all other messages */
	j->alloc_resp_port = slurmctld_comm_addr.port;
	j->other_port = slurmctld_comm_addr.port;
	if (slurmctld_comm_addr.hostname) {
		j->alloc_resp_hostname = xstrdup(slurmctld_comm_addr.hostname);
		j->other_hostname = xstrdup(slurmctld_comm_addr.hostname);
	} else {
		j->alloc_resp_hostname = NULL;
		j->other_hostname = NULL;
	}

	if (script) {
		/*
		 * If script is set then we are building a request for
		 *  a batch job
		 */
		xassert (opt.batch);

		j->environment = environ;
		j->env_size = envcount (environ);
		j->script = script;
		j->argv = remote_argv;
		j->argc = remote_argc;
		j->err  = opt.efname;
		j->in   = opt.ifname;
		j->out  = opt.ofname;
		j->work_dir = opt.cwd;
		j->no_requeue = opt.no_requeue;
	}

	return (j);
}

void
job_desc_msg_destroy(job_desc_msg_t *j)
{
	if (j) {
		xfree(j->account);
		xfree(j->alloc_resp_hostname);
		xfree(j->other_hostname);
		xfree(j);
	}
}

static job_step_create_request_msg_t *
_step_req_create(srun_job_t *j)
{
	job_step_create_request_msg_t *r = xmalloc(sizeof(*r));
	r->job_id     = j->jobid;
	r->user_id    = opt.uid;

	r->node_count = j->nhosts;
	/* info("send %d or %d? sending %d", opt.max_nodes, */
/* 		     j->nhosts, r->node_count); */
	if(r->node_count > j->nhosts) {
		error("Asking for more nodes that allocated");
		return NULL;
	}
	r->cpu_count  = opt.overcommit ? r->node_count
		                       : (opt.nprocs*opt.cpus_per_task);
	r->num_tasks  = opt.nprocs;
	r->node_list  = xstrdup(opt.nodelist);
	r->network    = xstrdup(opt.network);
	r->name       = xstrdup(opt.job_name);
	r->relative   = (uint16_t)opt.relative;
	r->overcommit = opt.overcommit ? 1 : 0;
	debug("requesting job %d, user %d, nodes %d (%s)", 
	      r->job_id, r->user_id, r->node_count, r->node_list);
	debug("cpus %d, tasks %d, name %s, relative %d", 
	      r->cpu_count, r->num_tasks, r->name, r->relative);
	
	switch (opt.distribution) {
	case SLURM_DIST_BLOCK:
		r->task_dist = SLURM_DIST_BLOCK;
		break;
	case SLURM_DIST_ARBITRARY:
		r->task_dist = SLURM_DIST_ARBITRARY;
		break;
	case SLURM_DIST_CYCLIC:
		r->task_dist = SLURM_DIST_CYCLIC;
		break;
	case SLURM_DIST_CYCLIC_CYCLIC:
		r->task_dist = SLURM_DIST_CYCLIC_CYCLIC;
		break;
	case SLURM_DIST_CYCLIC_BLOCK:
		r->task_dist = SLURM_DIST_CYCLIC_BLOCK;
		break;
	case SLURM_DIST_BLOCK_CYCLIC:
		r->task_dist = SLURM_DIST_BLOCK_CYCLIC;
		break;
	case SLURM_DIST_BLOCK_BLOCK:
		r->task_dist = SLURM_DIST_BLOCK_BLOCK;
		break;
	case SLURM_DIST_PLANE:
		r->task_dist = SLURM_DIST_PLANE;
		r->plane_size = opt.plane_size;
		break;
	default:
		r->task_dist = (opt.nprocs <= r->node_count) 
			? SLURM_DIST_CYCLIC : SLURM_DIST_BLOCK;
		break;

	}
	/* make sure we set the env correctly */
	opt.distribution = r->task_dist;
	
	if (slurmctld_comm_addr.port) {
		r->host = xstrdup(slurmctld_comm_addr.hostname);
		r->port = slurmctld_comm_addr.port;
	}

	return(r);
}

int
create_job_step(srun_job_t *job)
{
	job_step_create_request_msg_t  *req  = NULL;
	job_step_create_response_msg_t *resp = NULL;
	
	if (!(req = _step_req_create(job))) {
		error ("Unable to allocate step request message");
		return -1;
	}

	if ((slurm_job_step_create(req, &resp) < 0) || (resp == NULL)) {
		error ("Unable to create job step: %m");
		return -1;
	}
	
	job->stepid  = resp->job_step_id;
	job->step_layout = resp->step_layout;
	job->cred    = resp->cred;
	job->switch_job = resp->switch_job;
		
	if(!job->step_layout) {
		error("step_layout not returned");
		return -1;
	}
	
	/*
	 * Recreate filenames which may depend upon step id
	 */
	job_update_io_fnames(job);

	slurm_free_job_step_create_request_msg(req);
	
	return 0;
}

void 
set_allocate_job(srun_job_t *job) 
{
	allocate_job = job;
	return;
}
