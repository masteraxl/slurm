/*****************************************************************************\
 *  step_launch.c - launch a parallel job step
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher J. Morrone <morrone2@llnl.gov>
 *  UCRL-CODE-226842.
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

#include <errno.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <netinet/in.h>
#include <sys/param.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <netdb.h> /* for gethostbyname */

#include <slurm/slurm.h>

#include "src/common/hostlist.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_protocol_defs.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/eio.h"
#include "src/common/net.h"
#include "src/common/fd.h"
#include "src/common/slurm_auth.h"
#include "src/common/forward.h"
#include "src/common/plugstack.h"
#include "src/common/slurm_cred.h"
#include "src/common/mpi.h"

#include "src/api/step_launch.h"
#include "src/api/step_ctx.h"
#include "src/api/pmi_server.h"

extern char **environ;

/**********************************************************************
 * General declarations for step launch code
 **********************************************************************/
static int _launch_tasks(slurm_step_ctx_t *ctx,
			 launch_tasks_request_msg_t *launch_msg,
			 uint32_t timeout);
static char *_lookup_cwd(void);

/**********************************************************************
 * Message handler declarations
 **********************************************************************/
static uid_t  slurm_uid;
static int _msg_thr_create(struct step_launch_state *sls, int num_nodes);
static void _handle_msg(struct step_launch_state *sls, slurm_msg_t *msg);
static bool _message_socket_readable(eio_obj_t *obj);
static int _message_socket_accept(eio_obj_t *obj, List objs);

static struct io_operations message_socket_ops = {
	readable:	&_message_socket_readable,
	handle_read:	&_message_socket_accept
};


/**********************************************************************
 * API functions
 **********************************************************************/

/* 
 * slurm_step_launch_params_t_init - initialize a user-allocated
 *      slurm_job_step_launch_t structure with default values.
 *	default values.  This function will NOT allocate any new memory.
 * IN ptr - pointer to a structure allocated by the user.
 *      The structure will be intialized.
 */
void slurm_step_launch_params_t_init (slurm_step_launch_params_t *ptr)
{
	static slurm_step_io_fds_t fds = SLURM_STEP_IO_FDS_INITIALIZER;

	/* Set all values to zero (in other words, "NULL" for pointers) */
	memset(ptr, 0, sizeof(slurm_step_launch_params_t));

	ptr->buffered_stdio = true;
	memcpy(&ptr->local_fds, &fds, sizeof(fds));
	ptr->gid = getgid();
}

/*
 * slurm_step_launch - launch a parallel job step
 * IN ctx - job step context generated by slurm_step_ctx_create
 * IN callbacks - Identify functions to be called when various events occur
 * RET SLURM_SUCCESS or SLURM_ERROR (with errno set)
 */
int slurm_step_launch (slurm_step_ctx_t *ctx,
		       const slurm_step_launch_params_t *params,
		       const slurm_step_launch_callbacks_t *callbacks)
{
	launch_tasks_request_msg_t launch;
	int i;
	char **env = NULL;
	char **mpi_env = NULL;
	int rc = SLURM_SUCCESS;

	debug("Entering slurm_step_launch");
	memset(&launch, 0, sizeof(launch));

	if (ctx == NULL || ctx->magic != STEP_CTX_MAGIC) {
		error("Not a valid slurm_step_ctx_t!");

		slurm_seterrno(EINVAL);
		return SLURM_ERROR;
	}

	/* Initialize the callback pointers */
	if (callbacks != NULL) {
		/* copy the user specified callback pointers */
		memcpy(&(ctx->launch_state->callback), callbacks,
		       sizeof(slurm_step_launch_callbacks_t));
	} else {
		/* set all callbacks to NULL */
		memset(&(ctx->launch_state->callback), 0,
		       sizeof(slurm_step_launch_callbacks_t));
	}

	if (mpi_hook_client_init(params->mpi_plugin_name) == SLURM_ERROR) {
		slurm_seterrno(SLURM_MPI_PLUGIN_NAME_INVALID);
		return SLURM_ERROR;
	}
	/* Now, hack the step_layout struct if the following it true.
	   This looks like an ugly hack to support LAM/MPI's lamboot. */
	if (mpi_hook_client_single_task_per_node()) {
		for (i = 0; i < ctx->step_resp->step_layout->node_cnt; i++)
			ctx->step_resp->step_layout->tasks[i] = 1;
	}
	if ((ctx->launch_state->mpi_state =
	     mpi_hook_client_prelaunch(ctx->launch_state->mpi_info, &mpi_env))
	    == NULL) {
		slurm_seterrno(SLURM_MPI_PLUGIN_PRELAUNCH_SETUP_FAILED);
		return SLURM_ERROR;
	}

	/* Create message receiving sockets and handler thread */
	_msg_thr_create(ctx->launch_state,
			ctx->step_resp->step_layout->node_cnt);

	/* Start tasks on compute nodes */
	launch.job_id = ctx->step_req->job_id;
	launch.uid = ctx->step_req->user_id;
	launch.gid = params->gid;
	launch.argc = params->argc;
	launch.argv = params->argv;
	launch.cred = ctx->step_resp->cred;
	launch.job_step_id = ctx->step_resp->job_step_id;
	if (params->env == NULL) {
		/* if the user didn't specify an environment, grab the
		   environment of the running process */
		env_array_merge(&env, (const char **)environ);
	} else {
		env_array_merge(&env, (const char **)params->env);
	}
	{
		/* FIXME - hostname and IP need to be user settable */
		char *launcher_hostname = xshort_hostname();
		struct hostent *ent = gethostbyname(launcher_hostname);

		env_array_for_step(&env,
				   ctx->step_resp,
				   launcher_hostname,
				   ctx->launch_state->resp_port[0],
				   ent->h_addr_list[0]);
		xfree(launcher_hostname);
	}
	env_array_merge(&env, (const char **)mpi_env);
	env_array_free(mpi_env);

	launch.envc = envcount(env);
	launch.env = env;
	if (params->cwd != NULL) {
		launch.cwd = xstrdup(params->cwd);
	} else {
		launch.cwd = _lookup_cwd();
	}
	launch.nnodes = ctx->step_resp->step_layout->node_cnt;
	launch.nprocs = ctx->step_resp->step_layout->task_cnt;
	launch.slurmd_debug = params->slurmd_debug;
	launch.switch_job = ctx->step_resp->switch_job;
	launch.task_prolog = params->task_prolog;
	launch.task_epilog = params->task_epilog;
	launch.cpu_bind_type = params->cpu_bind_type;
	launch.cpu_bind = params->cpu_bind;
	launch.mem_bind_type = params->mem_bind_type;
	launch.mem_bind = params->mem_bind;
	launch.multi_prog = params->multi_prog ? 1 : 0;
	launch.cpus_per_task	= params->cpus_per_task;
	launch.ntasks_per_node	= params->ntasks_per_node;
	launch.ntasks_per_socket= params->ntasks_per_socket;
	launch.ntasks_per_core	= params->ntasks_per_core;
	launch.task_dist	= params->task_dist;
	launch.plane_size	= params->plane_size;
	launch.options = job_options_create();
	launch.complete_nodelist = 
		xstrdup(ctx->step_resp->step_layout->node_list);
	spank_set_remote_options (launch.options);
	launch.task_flags = 0;
	if (params->parallel_debug)
		launch.task_flags |= TASK_PARALLEL_DEBUG;

	launch.tasks_to_launch = ctx->step_resp->step_layout->tasks;
	launch.cpus_allocated  = ctx->step_resp->step_layout->tasks;
	launch.global_task_ids = ctx->step_resp->step_layout->tids;
	
	launch.user_managed_io = params->user_managed_io ? 1 : 0;
	ctx->launch_state->user_managed_io = params->user_managed_io;
	if (!ctx->launch_state->user_managed_io) {
		launch.ofname = params->remote_output_filename;
		launch.efname = params->remote_error_filename;
		launch.ifname = params->remote_input_filename;
		launch.buffered_stdio = params->buffered_stdio ? 1 : 0;
		ctx->launch_state->io.normal =
			client_io_handler_create(params->local_fds,
						 ctx->step_req->num_tasks,
						 ctx->step_req->node_count,
						 ctx->step_resp->cred,
						 params->labelio);
		if (ctx->launch_state->io.normal == NULL) {
			rc = SLURM_ERROR;
			goto fail1;
		}
		if (client_io_handler_start(ctx->launch_state->io.normal) 
		    != SLURM_SUCCESS) {
			rc = SLURM_ERROR;
			goto fail1;
		}
		launch.num_io_port = ctx->launch_state->io.normal->num_listen;
		launch.io_port = xmalloc(sizeof(uint16_t)*launch.num_io_port);
		for (i = 0; i < launch.num_io_port; i++) {
			launch.io_port[i] =
				ctx->launch_state->io.normal->listenport[i];
		}
	} else { /* user_managed_io is true */
		/* initialize user_managed_io_t */
		ctx->launch_state->io.user =
			(user_managed_io_t *)xmalloc(sizeof(user_managed_io_t));
		ctx->launch_state->io.user->connected = 0;
		ctx->launch_state->io.user->sockets =
			(int *)xmalloc(sizeof(int)*ctx->step_req->num_tasks);
	}

	launch.num_resp_port = ctx->launch_state->num_resp_port;
	launch.resp_port = xmalloc(sizeof(uint16_t) * launch.num_resp_port);
	for (i = 0; i < launch.num_resp_port; i++) {
		launch.resp_port[i] = ctx->launch_state->resp_port[i];
	}

	_launch_tasks(ctx, &launch, params->msg_timeout);

	/* clean up */
	xfree(launch.resp_port);
	if (!ctx->launch_state->user_managed_io) {
		xfree(launch.io_port);
	}
	goto done;
fail1:

done:
	xfree(launch.complete_nodelist);
	xfree(launch.cwd);
	env_array_free(env);
	job_options_destroy(launch.options);
	return rc;
}

/*
 * Block until all tasks have started.
 */
int slurm_step_launch_wait_start(slurm_step_ctx_t *ctx)
{
	struct step_launch_state *sls = ctx->launch_state;
	/* Wait for all tasks to start */
	pthread_mutex_lock(&sls->lock);
	while (bit_set_count(sls->tasks_started) < sls->tasks_requested) {
		if (sls->abort) {
			if (!sls->abort_action_taken) {
				slurm_kill_job_step(ctx->job_id,
						    ctx->step_resp->job_step_id,
						    SIGKILL);
				sls->abort_action_taken = true;
			}
			pthread_mutex_unlock(&sls->lock);
			return SLURM_ERROR;
		}
		pthread_cond_wait(&sls->cond, &sls->lock);
	}

	if (sls->user_managed_io) {
		while(sls->io.user->connected < sls->tasks_requested) {
			if (sls->abort) {
				if (!sls->abort_action_taken) {
					slurm_kill_job_step(
						ctx->job_id,
						ctx->step_resp->job_step_id,
						SIGKILL);
					sls->abort_action_taken = true;
				}
				pthread_mutex_unlock(&sls->lock);
				return SLURM_ERROR;
			}
			pthread_cond_wait(&sls->cond, &sls->lock);
		}
	}

	pthread_mutex_unlock(&sls->lock);
	return SLURM_SUCCESS;
}

/*
 * Block until all tasks have finished (or failed to start altogether).
 */
void slurm_step_launch_wait_finish(slurm_step_ctx_t *ctx)
{
	struct step_launch_state *sls = ctx->launch_state;
	struct timespec ts = {0, 0};
	bool time_set = false;
	int errnum;

	/* Wait for all tasks to complete */
	pthread_mutex_lock(&sls->lock);
	while (bit_set_count(sls->tasks_exited) < sls->tasks_requested) {
		if (!sls->abort) {
			pthread_cond_wait(&sls->cond, &sls->lock);
		} else {
			if (!sls->abort_action_taken) {
				slurm_kill_job_step(ctx->job_id,
						    ctx->step_resp->job_step_id,
						    SIGKILL);
				sls->abort_action_taken = true;
			}
			if (!time_set) {
				/* Only set the time once, because we only
				 * want to wait 10 seconds, no matter how many
				 * times the condition variable is signalled.
				 */
				ts.tv_sec = time(NULL) + 10;
				time_set = true;
				/* FIXME - should this be a callback? */
				info("Job step aborted: Waiting up to "
				     "10 seconds for job step to finish.");
			}

			errnum = pthread_cond_timedwait(&sls->cond,
							&sls->lock, &ts);
			if (errnum == ETIMEDOUT) {
				error("Timed out waiting for job step to "
				      "complete");
				/* 
				 * Send kill again, in case steps were still
				 * launching the first time.
				 * FIXME - eventually the slurmd should
				 *   be made smart enough to really ensure
				 *   that a killed step never starts.
				 */
				slurm_kill_job_step(
					ctx->job_id,
					ctx->step_resp->job_step_id,
					SIGKILL);
				if (!sls->user_managed_io)
					client_io_handler_abort(sls->io.normal);
				break;
			} else if (errnum != 0) {
				error("Error waiting on condition in"
				      " slurm_step_launch_wait_finish: %m");
				if (!sls->user_managed_io)
					client_io_handler_abort(sls->io.normal);
				break;
			}
		}
	}
	
	/* Then shutdown the message handler thread */
	eio_signal_shutdown(sls->msg_handle);
	pthread_join(sls->msg_thread, NULL);
	eio_handle_destroy(sls->msg_handle);

	/* Then wait for the IO thread to finish */
	if (!sls->user_managed_io) {
		client_io_handler_finish(sls->io.normal);
		client_io_handler_destroy(sls->io.normal);
	}

	mpi_hook_client_fini(sls->mpi_state);

	pthread_mutex_unlock(&sls->lock);
}

/*
 * Abort an in-progress launch, or terminate the fully launched job step.
 *
 * Can be called from a signal handler.
 */
void slurm_step_launch_abort(slurm_step_ctx_t *ctx)
{
	struct step_launch_state *sls = ctx->launch_state;

	sls->abort = true;
	pthread_cond_signal(&sls->cond);
}

/* 
 * Forward a signal to all those nodes with running tasks 
 */
void slurm_step_launch_fwd_signal(slurm_step_ctx_t *ctx, int signo)
{
	int node_id, j, active, num_tasks;
	slurm_msg_t req;
	kill_tasks_msg_t msg;
	hostlist_t hl;
	char *name = NULL;
	char buf[8192];
	List ret_list = NULL;
	ListIterator itr;
	ret_data_info_t *ret_data_info = NULL;
	int rc = SLURM_SUCCESS;
	struct step_launch_state *sls = ctx->launch_state;
	
	debug2("forward signal %d to job", signo);
	
	/* common to all tasks */
	msg.job_id      = ctx->job_id;
	msg.job_step_id = ctx->step_resp->job_step_id;
	msg.signal      = (uint32_t) signo;
	
	pthread_mutex_lock(&sls->lock);
	
	hl = hostlist_create("");
	for (node_id = 0;
	     node_id < ctx->step_resp->step_layout->node_cnt;
	     node_id++) {
		active = 0;		
		num_tasks = sls->layout->tasks[node_id];
		for (j = 0; j < num_tasks; j++) {
			if(bit_test(sls->tasks_started,
				    sls->layout->tids[node_id][j]) &&
			   !bit_test(sls->tasks_exited,
				     sls->layout->tids[node_id][j])) {
				/* this one has active tasks */
				active = 1;
				break;
			}
		}
		
		if (!active)
			continue;
		
		name = nodelist_nth_host(sls->layout->node_list, node_id);
		hostlist_push(hl, name);
		free(name);
	}

	pthread_mutex_unlock(&sls->lock);
	
	if(!hostlist_count(hl)) {
		hostlist_destroy(hl);
		goto nothing_left;
	}
	hostlist_ranged_string(hl, sizeof(buf), buf);
	hostlist_destroy(hl);
	name = xstrdup(buf);
	
	slurm_msg_t_init(&req);	
	req.msg_type = REQUEST_SIGNAL_TASKS;
	req.data     = &msg;
	
	debug3("sending signal to host %s", name);
	
	if (!(ret_list = slurm_send_recv_msgs(name, &req, 0))) { 
		error("fwd_signal: slurm_send_recv_msgs really failed bad");
		xfree(name);
		return;
	}
	xfree(name);
	itr = list_iterator_create(ret_list);		
	while((ret_data_info = list_next(itr))) {
		rc = slurm_get_return_code(ret_data_info->type, 
					   ret_data_info->data);
		/*
		 *  Report error unless it is "Invalid job id" which 
		 *    probably just means the tasks exited in the meanwhile.
		 */
		if ((rc != 0) && (rc != ESLURM_INVALID_JOB_ID)
		    &&  (rc != ESLURMD_JOB_NOTRUNNING) && (rc != ESRCH)) {
			error("%s: signal: %s", 
			      ret_data_info->node_name, 
			      slurm_strerror(rc));
		}
	}
	list_iterator_destroy(itr);
	list_destroy(ret_list);
nothing_left:
	debug2("All tasks have been signalled");
	
}

/**********************************************************************
 * Functions used by step_ctx code, but not exported throught the API
 **********************************************************************/
/*
 * Create a launch state structure for a specified step context, "ctx".
 */
struct step_launch_state *step_launch_state_create(slurm_step_ctx_t *ctx)
{
	struct step_launch_state *sls;
	slurm_step_layout_t *layout = ctx->step_resp->step_layout;

	sls = xmalloc(sizeof(struct step_launch_state));
	if (sls != NULL) {
		sls->slurmctld_socket_fd = -1;
		sls->tasks_requested = layout->task_cnt;
		sls->tasks_started = bit_alloc(layout->task_cnt);
		sls->tasks_exited = bit_alloc(layout->task_cnt);
		sls->layout = layout;
		sls->resp_port = NULL;
		sls->abort = false;
		sls->abort_action_taken = false;
		sls->mpi_info->jobid = ctx->step_req->job_id;
		sls->mpi_info->stepid = ctx->step_resp->job_step_id;
		sls->mpi_info->step_layout = layout;
		sls->mpi_state = NULL;
		pthread_mutex_init(&sls->lock, NULL);
		pthread_cond_init(&sls->cond, NULL);
	}
	return sls;
}

/*
 * Free the memory associated with the a launch state structure.
 */
void step_launch_state_destroy(struct step_launch_state *sls)
{
	/* First undo anything created in step_launch_state_create() */
	pthread_mutex_destroy(&sls->lock);
	pthread_cond_destroy(&sls->cond);
	bit_free(sls->tasks_started);
	bit_free(sls->tasks_exited);

	/* Now clean up anything created by slurm_step_launch() */
	if (sls->resp_port != NULL) {
		xfree(sls->resp_port);
	}
}

/**********************************************************************
 * Message handler functions
 **********************************************************************/
static void *_msg_thr_internal(void *arg)
{
	struct step_launch_state *sls = (struct step_launch_state *)arg;

	eio_handle_mainloop(sls->msg_handle);

	return NULL;
}

static inline int
_estimate_nports(int nclients, int cli_per_port)
{
	div_t d;
	d = div(nclients, cli_per_port);
	return d.rem > 0 ? d.quot + 1 : d.quot;
}

static int _msg_thr_create(struct step_launch_state *sls, int num_nodes)
{
	int sock = -1;
	short port = -1;
	eio_obj_t *obj;
	int i;

	debug("Entering _msg_thr_create()");
	slurm_uid = (uid_t) slurm_get_slurm_user_id();

	sls->msg_handle = eio_handle_create();
	sls->num_resp_port = _estimate_nports(num_nodes, 48);
	sls->resp_port = xmalloc(sizeof(uint16_t) * sls->num_resp_port);
	for (i = 0; i < sls->num_resp_port; i++) {
		if (net_stream_listen(&sock, &port) < 0) {
			error("unable to intialize step launch listening socket: %m");
			return SLURM_ERROR;
		}
		sls->resp_port[i] = port;
		obj = eio_obj_create(sock, &message_socket_ops, (void *)sls);
		eio_new_initial_obj(sls->msg_handle, obj);
	}
	/* finally, add the listening port that we told the slurmctld about
	   eariler in the step context creation phase */
	if (sls->slurmctld_socket_fd > -1) {
		obj = eio_obj_create(sls->slurmctld_socket_fd,
				     &message_socket_ops, (void *)sls);
		eio_new_initial_obj(sls->msg_handle, obj);
	}

	if (pthread_create(&sls->msg_thread, NULL,
			   _msg_thr_internal, (void *)sls) != 0) {
		error("pthread_create of message thread: %m");
		return SLURM_ERROR;
	}
	return SLURM_SUCCESS;
}

static bool _message_socket_readable(eio_obj_t *obj)
{
	debug3("Called _message_socket_readable");
	if (obj->shutdown == true) {
		if (obj->fd != -1) {
			debug2("  false, shutdown");
			close(obj->fd);
			obj->fd = -1;
			/*_wait_for_connections();*/
		} else {
			debug2("  false");
		}
		return false;
	}
	return true;
}

static int _message_socket_accept(eio_obj_t *obj, List objs)
{
	struct step_launch_state *sls = (struct step_launch_state *)obj->arg;

	int fd;
	unsigned char *uc;
	short port;
	struct sockaddr_un addr;
	slurm_msg_t *msg = NULL;
	int len = sizeof(addr);
	int timeout = 0;	/* slurm default value */
	int rc = 0;
	
	debug3("Called _msg_socket_accept");

	while ((fd = accept(obj->fd, (struct sockaddr *)&addr,
			    (socklen_t *)&len)) < 0) {
		if (errno == EINTR)
			continue;
		if (errno == EAGAIN
		    || errno == ECONNABORTED
		    || errno == EWOULDBLOCK) {
			return SLURM_SUCCESS;
		}
		error("Error on msg accept socket: %m");
		obj->shutdown = true;
		return SLURM_SUCCESS;
	}

	fd_set_close_on_exec(fd);
	fd_set_blocking(fd);

	/* Should not call slurm_get_addr() because the IP may not be
	   in /etc/hosts. */
	uc = (unsigned char *)&((struct sockaddr_in *)&addr)->sin_addr.s_addr;
	port = ((struct sockaddr_in *)&addr)->sin_port;
	debug2("got message connection from %u.%u.%u.%u:%hu",
	       uc[0], uc[1], uc[2], uc[3], ntohs(port));
	fflush(stdout);

	msg = xmalloc(sizeof(slurm_msg_t));
	slurm_msg_t_init(msg);

	/* multiple jobs (easily induced via no_alloc) and highly
	 * parallel jobs using PMI sometimes result in slow message 
	 * responses and timeouts. Raise the default timeout for srun. */
	timeout = slurm_get_msg_timeout() * 8000;
again:
	if((rc = slurm_receive_msg(fd, msg, timeout)) != 0) {
		if (errno == EINTR) {
			goto again;
		}
		error("slurm_receive_msg[%u.%u.%u.%u]: %m",
		      uc[0],uc[1],uc[2],uc[3]);
		goto cleanup;
	}

	_handle_msg(sls, msg); /* handle_msg frees msg */
cleanup:
	if ((msg->conn_fd >= 0) && slurm_close_accepted_conn(msg->conn_fd) < 0)
		error ("close(%d): %m", msg->conn_fd);
	slurm_free_msg(msg);

	return SLURM_SUCCESS;
}

static void
_launch_handler(struct step_launch_state *sls, slurm_msg_t *resp)
{
	launch_tasks_response_msg_t *msg = resp->data;
	int i;

	pthread_mutex_lock(&sls->lock);

	for (i = 0; i < msg->count_of_pids; i++) {
		bit_set(sls->tasks_started, msg->task_ids[i]);
	}

	if (sls->callback.task_start != NULL)
		(sls->callback.task_start)(msg);

	pthread_cond_signal(&sls->cond);
	pthread_mutex_unlock(&sls->lock);

}

static void 
_exit_handler(struct step_launch_state *sls, slurm_msg_t *exit_msg)
{
	task_exit_msg_t *msg = (task_exit_msg_t *) exit_msg->data;
	int i;

	pthread_mutex_lock(&sls->lock);

	for (i = 0; i < msg->num_tasks; i++) {
		debug("task %d done", msg->task_id_list[i]);
		bit_set(sls->tasks_exited, msg->task_id_list[i]);
	}

	if (sls->callback.task_finish != NULL)
		(sls->callback.task_finish)(msg);

	pthread_cond_signal(&sls->cond);
	pthread_mutex_unlock(&sls->lock);
}

/*
 * Take the list of node names of down nodes and convert into an
 * array of nodeids for the step.  The nodeid array is passed to
 * client_io_handler_downnodes to notify the IO handler to expect no
 * further IO from that node.
 */
static void
_node_fail_handler(struct step_launch_state *sls, slurm_msg_t *fail_msg)
{
	srun_node_fail_msg_t *nf = fail_msg->data;
	hostset_t fail_nodes, all_nodes;
	hostlist_iterator_t fail_itr;
	char *node;
	int num_node_ids;
	int *node_ids;
	int i, j;
	int node_id, num_tasks;

	fail_nodes = hostset_create(nf->nodelist);
	fail_itr = hostset_iterator_create(fail_nodes);
	num_node_ids = hostset_count(fail_nodes);
	node_ids = xmalloc(sizeof(int) * num_node_ids);

	pthread_mutex_lock(&sls->lock);
	all_nodes = hostset_create(sls->layout->node_list);
	/* find the index number of each down node */
	for (i = 0; i < num_node_ids; i++) {
		node = hostlist_next(fail_itr);
		node_id = node_ids[i] = hostset_find(all_nodes, node);
		free(node);

		/* find all of the task that should run on this node and
		 * mark them as having started and exited.  If they haven't
		 * started yet, they never will, and likewise for exiting.
		 */
		num_tasks = sls->layout->tasks[node_id];
		for (j = 0; j < num_tasks; j++) {
			debug2("marking task %d done on failed node %d",
			       sls->layout->tids[node_id][j], node_id);
			bit_set(sls->tasks_started,
				sls->layout->tids[node_id][j]);
			bit_set(sls->tasks_exited,
				sls->layout->tids[node_id][j]);
		}
	}

	if (!sls->user_managed_io) {
		client_io_handler_downnodes(sls->io.normal, node_ids,
					    num_node_ids);
	}
	pthread_cond_signal(&sls->cond);
	pthread_mutex_unlock(&sls->lock);

	xfree(node_ids);
	hostlist_iterator_destroy(fail_itr);
	hostset_destroy(fail_nodes);
	hostset_destroy(all_nodes);
}

/*
 * The TCP connection that was used to send the task_spawn_io_msg_t message
 * will be used as the user managed IO stream.  The remote end of the TCP stream
 * will be connected to the stdin, stdout, and stderr of the task.  The
 * local end of the stream is stored in the user_managed_io_t structure, and
 * is left to the user to manage (the user can retrieve the array of
 * socket descriptors using slurm_step_ctx_get()).
 *
 * To allow the message TCP stream to be reused for spawn IO traffic we
 * set the slurm_msg_t's conn_fd to -1 to avoid having the caller close the
 * TCP stream.
 */
static void
_task_user_managed_io_handler(struct step_launch_state *sls,
			      slurm_msg_t *user_io_msg)
{
	task_user_managed_io_msg_t *msg =
		(task_user_managed_io_msg_t *) user_io_msg->data;

	pthread_mutex_lock(&sls->lock);

	debug("task %d user managed io stream established", msg->task_id);
	/* sanity check */
	if (msg->task_id >= sls->tasks_requested) {
		error("_task_user_managed_io_handler:"
		      " bad task ID %u (of %d tasks)",
		      msg->task_id, sls->tasks_requested);
	}

	sls->io.user->connected++;
	fd_set_blocking(user_io_msg->conn_fd);
	sls->io.user->sockets[msg->task_id] = user_io_msg->conn_fd;

	/* prevent the caller from closing the user managed IO stream */
	user_io_msg->conn_fd = -1;

	pthread_cond_signal(&sls->cond);
	pthread_mutex_unlock(&sls->lock);
}

/*
 * Identify the incoming message and call the appropriate handler function.
 */
static void
_handle_msg(struct step_launch_state *sls, slurm_msg_t *msg)
{
	uid_t req_uid = g_slurm_auth_get_uid(msg->auth_cred);
	uid_t uid = getuid();
	int rc;
	
	if ((req_uid != slurm_uid) && (req_uid != 0) && (req_uid != uid)) {
		error ("Security violation, slurm message from uid %u", 
		       (unsigned int) req_uid);
 		return;
	}

	switch (msg->msg_type) {
	case RESPONSE_LAUNCH_TASKS:
		debug2("received task launch");
		_launch_handler(sls, msg);
		slurm_free_launch_tasks_response_msg(msg->data);
		break;
	case MESSAGE_TASK_EXIT:
		debug2("received task exit");
		_exit_handler(sls, msg);
		slurm_free_task_exit_msg(msg->data);
		break;
	case SRUN_NODE_FAIL:
		debug2("received srun node fail");
		_node_fail_handler(sls, msg);
		slurm_free_srun_node_fail_msg(msg->data);
		break;
	case SRUN_TIMEOUT:
		debug2("received job step timeout message");
		/* FIXME - does nothing yet */
		slurm_free_srun_timeout_msg(msg->data);
		break;
	case SRUN_JOB_COMPLETE:
		debug2("received job step complete message");
		/* FIXME - does nothing yet */
		slurm_free_srun_job_complete_msg(msg->data);
		break;
	case PMI_KVS_PUT_REQ:
		debug2("PMI_KVS_PUT_REQ received");
		rc = pmi_kvs_put((struct kvs_comm_set *) msg->data);
		slurm_send_rc_msg(msg, rc);
		break;
	case PMI_KVS_GET_REQ:
		debug2("PMI_KVS_GET_REQ received");
		rc = pmi_kvs_get((kvs_get_msg_t *) msg->data);
		slurm_send_rc_msg(msg, rc);
		slurm_free_get_kvs_msg((kvs_get_msg_t *) msg->data);
		break;
	case TASK_USER_MANAGED_IO_STREAM:
		debug2("TASK_USER_MANAGED_IO_STREAM");
		_task_user_managed_io_handler(sls, msg);
		break;
	default:
		error("received spurious message type: %d",
		      msg->msg_type);
		break;
	}
	return;
}

/**********************************************************************
 * Task launch functions
 **********************************************************************/
static int _launch_tasks(slurm_step_ctx_t *ctx,
			 launch_tasks_request_msg_t *launch_msg,
			 uint32_t timeout)
{
	slurm_msg_t msg;
	List ret_list = NULL;
	ListIterator ret_itr;
	ret_data_info_t *ret_data = NULL;
	int rc = SLURM_SUCCESS;

	debug("Entering _launch_tasks");
	slurm_msg_t_init(&msg);
	msg.msg_type = REQUEST_LAUNCH_TASKS;
	msg.data = launch_msg;
	
	if(!(ret_list = slurm_send_recv_msgs(
		     ctx->step_resp->step_layout->node_list,
		     &msg, timeout))) {
		error("slurm_send_recv_msgs failed miserably: %m");
		return SLURM_ERROR;
	}
	ret_itr = list_iterator_create(ret_list);
	while ((ret_data = list_next(ret_itr))) {
		rc = slurm_get_return_code(ret_data->type, 
					   ret_data->data);
		debug("launch returned msg_rc=%d err=%d type=%d",
		      rc, ret_data->err, ret_data->type);
		if (rc != SLURM_SUCCESS) {
			errno = ret_data->err;
			error("Task launch failed on node %s: %m",
			      ret_data->node_name);
		} else {
#if 0 /* only for debugging, might want to make this a callback */
			errno = ret_data->err;
			info("Launch success on node %s",
			     ret_data->node_name);
#endif
		}
	}
	list_iterator_destroy(ret_itr);
	list_destroy(ret_list);
	return SLURM_SUCCESS;
}

/* returns an xmalloc cwd string, or NULL if lookup failed. */
static char *_lookup_cwd(void)
{
	char buf[PATH_MAX];

	if (getcwd(buf, PATH_MAX) != NULL) {
		return xstrdup(buf);
	} else {
		return NULL;
	}
}
