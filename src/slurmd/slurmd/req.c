/*****************************************************************************\
 *  src/slurmd/slurmd/req.c - slurmd request handling
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <mgrondona@llnl.gov>.
 *  UCRL-CODE-2002-040.
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
#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <fcntl.h>
#include <pthread.h>
#include <sched.h>
#include <signal.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/param.h>		/* MAXPATHLEN */
#include <sys/poll.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "src/common/hostlist.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/node_select.c"
#include "src/common/slurm_auth.h"
#include "src/common/slurm_cred.h"
#include "src/common/slurm_jobacct.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_protocol_interface.h"
#include "src/common/xstring.h"
#include "src/common/xmalloc.h"
#include "src/common/list.h"
#include "src/common/util-net.h"

#include "src/slurmd/slurmd/slurmd.h"
#include "src/slurmd/common/proctrack.h"
#include "src/slurmd/common/slurmstepd_init.h"
#include "src/slurmd/common/stepd_api.h"

#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN	64
#endif


static int  _abort_job(uint32_t job_id);
static char ** _build_env(uint32_t jobid, uid_t uid, char *bgl_part_id);
static bool _slurm_authorized_user(uid_t uid);
static bool _job_still_running(uint32_t job_id);
static int  _kill_all_active_steps(void *auth_cred, uint32_t jobid,
				   int sig, bool batch);
static void _rpc_launch_tasks(slurm_msg_t *, slurm_addr *);
static void _rpc_spawn_task(slurm_msg_t *, slurm_addr *);
static void _rpc_batch_job(slurm_msg_t *, slurm_addr *);
static void _rpc_signal_tasks(slurm_msg_t *, slurm_addr *);
static void _rpc_terminate_tasks(slurm_msg_t *, slurm_addr *);
static void _rpc_timelimit(slurm_msg_t *, slurm_addr *);
static void _rpc_reattach_tasks(slurm_msg_t *, slurm_addr *);
static void _rpc_signal_job(slurm_msg_t *, slurm_addr *);
static void _rpc_terminate_job(slurm_msg_t *, slurm_addr *);
static void _rpc_update_time(slurm_msg_t *, slurm_addr *);
static void _rpc_shutdown(slurm_msg_t *msg, slurm_addr *cli_addr);
static void _rpc_reconfig(slurm_msg_t *msg, slurm_addr *cli_addr);
static void _rpc_pid2jid(slurm_msg_t *msg, slurm_addr *);
static int  _rpc_ping(slurm_msg_t *, slurm_addr *);
static int  _run_prolog(uint32_t jobid, uid_t uid, char *bgl_part_id);
static int  _run_epilog(uint32_t jobid, uid_t uid, char *bgl_part_id);

static bool _pause_for_job_completion(void *auth_cred, uint32_t jobid,
				      int maxtime);
static int _waiter_init (uint32_t jobid);
static int _waiter_complete (uint32_t jobid);

static bool _steps_completed_now(uint32_t jobid);
static void _wait_state_completed(uint32_t jobid, int max_delay);

/*
 *  List of threads waiting for jobs to complete
 */
static List waiters;

static pthread_mutex_t launch_mutex = PTHREAD_MUTEX_INITIALIZER;

void
slurmd_req(slurm_msg_t *msg, slurm_addr *cli)
{
	int rc;

	switch(msg->msg_type) {
	case REQUEST_BATCH_JOB_LAUNCH:
		/* Mutex locking moved into _rpc_batch_job() due to 
		 * very slow prolog on Blue Gene system. Only batch 
		 * jobs are supported on Blue Gene (no job steps). */
		_rpc_batch_job(msg, cli);
		slurm_free_job_launch_msg(msg->data);
		break;
	case REQUEST_LAUNCH_TASKS:
		slurm_mutex_lock(&launch_mutex);
		_rpc_launch_tasks(msg, cli);
		slurm_free_launch_tasks_request_msg(msg->data);
		slurm_mutex_unlock(&launch_mutex);
		break;
	case REQUEST_SPAWN_TASK:
		slurm_mutex_lock(&launch_mutex);
		_rpc_spawn_task(msg, cli);
		slurm_free_spawn_task_request_msg(msg->data);
		slurm_mutex_unlock(&launch_mutex);
		break;
	case REQUEST_SIGNAL_TASKS:
		debug2("Processing RPC: REQUEST_SIGNAL_TASKS");
		_rpc_signal_tasks(msg, cli);
		slurm_free_kill_tasks_msg(msg->data);
		break;
	case REQUEST_TERMINATE_TASKS:
		debug2("Processing RPC: REQUEST_TERMINATE_TASKS");
		_rpc_terminate_tasks(msg, cli);
		slurm_free_kill_tasks_msg(msg->data);
		break;
	case REQUEST_KILL_TIMELIMIT:
		debug2("Processing RPC: REQUEST_KILL_TIMELIMIT");
		_rpc_timelimit(msg, cli);
		slurm_free_timelimit_msg(msg->data);
		break; 
	case REQUEST_REATTACH_TASKS:
		_rpc_reattach_tasks(msg, cli);
		slurm_free_reattach_tasks_request_msg(msg->data);
		break;
	case REQUEST_SIGNAL_JOB:
		debug2("Processing RPC: REQUEST_SIGNAL_JOB");
		_rpc_signal_job(msg, cli);
		slurm_free_kill_job_msg(msg->data);
		break;
	case REQUEST_TERMINATE_JOB:
		debug2("Processing RPC: REQUEST_TERMINATE_JOB");
		_rpc_terminate_job(msg, cli);
		slurm_free_kill_job_msg(msg->data);
		break;
	case REQUEST_UPDATE_JOB_TIME:
		_rpc_update_time(msg, cli);
		slurm_free_update_job_time_msg(msg->data);
		break;
	case REQUEST_SHUTDOWN:
		_rpc_shutdown(msg, cli);
		slurm_free_shutdown_msg(msg->data);
		break;
	case REQUEST_RECONFIGURE:
		_rpc_reconfig(msg, cli);
		/* No body to free */
		break;
	case REQUEST_NODE_REGISTRATION_STATUS:
		/* Treat as ping (for slurmctld agent, just return SUCCESS) */
		rc = _rpc_ping(msg, cli);
		/* No body to free */
		/* Then initiate a separate node registration */
		if (rc == SLURM_SUCCESS)
			send_registration_msg(SLURM_SUCCESS, true);
		break;
	case REQUEST_PING:
		_rpc_ping(msg, cli);
		/* No body to free */
		break;
	case REQUEST_JOB_ID:
		_rpc_pid2jid(msg, cli);
		slurm_free_job_id_request_msg(msg->data);
		break;
	case MESSAGE_JOBACCT_DATA:
		{
			int rc=SLURM_SUCCESS;
			debug3("jobacct(%i) received jobacct message",
					getpid());
			slurm_send_rc_msg(msg,rc); /* ACK the message */
			debug3("jobacct(%i) sent jobacct rc=%d message",
					getpid(), rc);
			rc=g_slurm_jobacct_process_message(msg);
			debug3("jobacct(%i) slurm_jobacct_process_message "
					"rc=%d",
					getpid(), rc);
			slurm_free_jobacct_msg(msg->data);
		}
		break;
	default:
		error("slurmd_req: invalid request msg type %d\n",
		      msg->msg_type);
		slurm_send_rc_msg(msg, EINVAL);
		break;
	}
	return;
}

/*
 *  Need to close all open fds
 */
static void
_close_fds(void)
{
	int i;
	int maxfd = 1024;
	for (i = 4; i < maxfd; i++) {
		close(i);
	}
}

static int
_send_slurmstepd_init(int fd, slurmd_step_type_t type, void *req, 
		      slurm_addr *cli, slurm_addr *self)
{
	int rc;
	int len = 0;
	Buf buffer;
	slurm_msg_t *msg = NULL;

	/* send type over to slurmstepd */
	safe_write(fd, &type, sizeof(int));

	/* send conf over to slurmstepd */
	buffer = init_buf(0);
	pack_slurmd_conf_lite(conf, buffer);
	len = get_buf_offset(buffer);
	safe_write(fd, &len, sizeof(int));
	safe_write(fd, get_buf_data(buffer), len);
	free_buf(buffer);

	/* send cli over to slurmstepd */
	buffer = init_buf(0);
	slurm_pack_slurm_addr(cli, buffer);
	len = get_buf_offset(buffer);
	safe_write(fd, &len, sizeof(int));
	safe_write(fd, get_buf_data(buffer), len);
	free_buf(buffer);

	/* send self over to slurmstepd */
	if(self) {
		buffer = init_buf(0);
		slurm_pack_slurm_addr(self, buffer);
		len = get_buf_offset(buffer);
		safe_write(fd, &len, sizeof(int));
		safe_write(fd, get_buf_data(buffer), len);
		free_buf(buffer);
	} else {
		len = 0;
		safe_write(fd, &len, sizeof(int));
	}

	/* send req over to slurmstepd */
	msg = xmalloc(sizeof(slurm_msg_t));
	switch(type) {
	case LAUNCH_BATCH_JOB:
		msg->msg_type = REQUEST_BATCH_JOB_LAUNCH;
		break;
	case LAUNCH_TASKS:
		msg->msg_type = REQUEST_LAUNCH_TASKS;
		break;
	case SPAWN_TASKS:
		msg->msg_type = REQUEST_SPAWN_TASK;
		break;
	default:
		error("Was sent a task I didn't understand");
		break;
	}
	buffer = init_buf(0);
	msg->data = req;
	pack_msg(msg, buffer);
	len = get_buf_offset(buffer);
	safe_write(fd, &len, sizeof(int));
	safe_write(fd, get_buf_data(buffer), len);
	free_buf(buffer);
	slurm_free_msg(msg);

	return 0;

rwfail:
	error("_send_slurmstepd_init failed");
	return -1;
}


/*
 * Fork and exec the slurmstepd, then send the slurmstepd its
 * initialization data.  Then wait for slurmstepd to send an "ok"
 * message before returning.  When the "ok" message is received,
 * the slurmstepd has created and begun listening on its unix
 * domain socket.
 *
 * Note that this code forks twice and it is the grandchild that
 * becomes the slurmstepd process, so the slurmstepd's parent process
 * will be init, not slurmd.
 */
static int
_forkexec_slurmstepd(slurmd_step_type_t type, void *req, 
		     slurm_addr *cli, slurm_addr *self)
{
	pid_t pid;
	int to_stepd[2] = {-1, -1};
	int to_slurmd[2] = {-1, -1};

	if (pipe(to_stepd) < 0 || pipe(to_slurmd) < 0) {
		error("_forkexec_slurmstepd pipe failed: %m");
		return -1;
	}

	if ((pid = fork()) < 0) {
		error("fork_slurmd: fork: %m");
		close(to_stepd[0]);
		close(to_stepd[1]);
		close(to_slurmd[0]);
		close(to_slurmd[1]);
		return -1;
	} else if (pid > 0) {
		int ok;
		int rc = 0;
		/*
		 * Parent sends initialization data to the slurmstepd
		 * over the to_stepd pipe, and waits for an "ok"
		 * reply on the to_slurmd pipe.
		 */
		if (close(to_stepd[0]) < 0)
			error("Unable to close read to_stepd in parent: %m");
		if (close(to_slurmd[1]) < 0)
			error("Unable to close write to_slurmd in parent: %m");

		if ((rc = _send_slurmstepd_init(to_stepd[1], type,
						req, cli, self)) < 0) {
			rc = -1;
		}
		if (read(to_slurmd[0], &ok, sizeof(int)) != sizeof(int)) {
			error("Error reading \"ok\" message "
			      " from slurmstepd: %m");
			rc = -2;
		}

		/* Reap child */
		if (waitpid(pid, NULL, 0) < 0)
			error("Unable to reap slurmd child process");
		if (close(to_stepd[1]) < 0)
			error("close write to_stepd in parent: %m");
		if (close(to_slurmd[0]) < 0)
			error("close read to_slurmd in parent: %m");

		return rc;
	} else {
		char **argv = NULL;
		/*
		 * Child forks and exits
		 */
		if (setsid() < 0)
			error("fork_slurmd: setsid: %m");
		if ((pid = fork()) < 0)
			error("fork_slurmd: Unable to fork grandchild: %m");
		else if (pid > 0) { /* child */
			exit(0);
		}

		/*
		 * Grandchild exec's the slurmstepd
		 */
		slurm_shutdown_msg_engine(conf->lfd);
		
		if (close(to_stepd[1]) < 0)
			error("close write to_stepd in grandchild: %m");
		if (close(to_slurmd[0]) < 0)
			error("close read to_slurmd in parent: %m");
		if (dup2(to_stepd[0], STDIN_FILENO) == -1) {
			error("dup2 over STDIN_FILENO: %m");
			exit(1);
		}
		if (dup2(to_slurmd[1], STDOUT_FILENO) == -1) {
			error("dup2 over STDOUT_FILENO: %m");
			exit(1);
		}
		argv = xmalloc(2 * sizeof(char *));
		argv[0] = SLURMD_STEP_PATH;
		argv[1] = NULL;
		execvp(argv[0], argv);

		fatal("exec of slurmstepd failed: %m");
		exit(2);
	}
}

static int
_check_job_credential(slurm_cred_t cred, uint32_t jobid, 
		      uint32_t stepid, uid_t uid, int tasks_to_launch)
{
	slurm_cred_arg_t arg;
	hostset_t        hset    = NULL;
	bool             user_ok = _slurm_authorized_user(uid); 
	int              host_index = -1;
	int              rc;

	/*
	 * First call slurm_cred_verify() so that all valid
	 * credentials are checked
	 */
	if (((rc = slurm_cred_verify(conf->vctx, cred, &arg)) < 0) && !user_ok)
		return SLURM_ERROR;

	/*
	 * If uid is the slurm user id or root, do not bother
	 * performing validity check of the credential
	 */
	if (user_ok) {
		if (rc >= 0)
			xfree(arg.hostlist);
		return SLURM_SUCCESS;
	}

	if ((arg.jobid != jobid) || (arg.stepid != stepid)) {
		error("job credential for %d.%d, expected %d.%d",
		      arg.jobid, arg.stepid, jobid, stepid); 
		goto fail;
	}

	if (arg.uid != uid) {
		error("job credential created for uid %ld, expected %ld",
		      (long) arg.uid, (long) uid);
		goto fail;
	}

	/*
	 * Check that credential is valid for this host
	 */
	if (!(hset = hostset_create(arg.hostlist))) {
		error("Unable to parse credential hostlist: `%s'", 
		      arg.hostlist);
		goto fail;
	}

	if (!hostset_within(hset, conf->node_name)) {
		error("job credential invald for this host [%d.%d %ld %s]",
		      arg.jobid, arg.stepid, (long) arg.uid, arg.hostlist);
		goto fail;
	}

        if ((arg.ntask_cnt > 0) && (tasks_to_launch > 0)) {

                host_index = hostset_index(hset, conf->node_name, jobid);

                /* Left in here for debugging purposes
                if(host_index >= 0)
                  debug3(" cons_res %u ntask_cnt %d task[%d] = %d = task_to_launch %d host %s ", 
                         arg.jobid, arg.ntask_cnt, host_index, arg.ntask[host_index], 
                         tasks_to_launch, conf->node_name);
		*/

                if (host_index < 0) { 
                        error("job cr credential invalid host_index %d for job %d",
                              host_index, arg.jobid);
                        goto fail; 
                }
                
                if (!(arg.ntask[host_index] == tasks_to_launch)) {
                        error("job cr credential (%d != %d) invalid for this host [%d.%d %ld %s]",
                              arg.ntask[host_index], tasks_to_launch, arg.jobid, arg.stepid, 
                              (long) arg.uid, arg.hostlist);
                        goto fail;
                }
        }

	hostset_destroy(hset);
	xfree(arg.hostlist);
	arg.ntask_cnt = 0;
	xfree(arg.ntask);

	return SLURM_SUCCESS;

    fail:
	if (hset)
		hostset_destroy(hset);
	xfree(arg.hostlist);
        arg.ntask_cnt = 0;
        xfree(arg.ntask);
	slurm_seterrno_ret(ESLURMD_INVALID_JOB_CREDENTIAL);
}


static void 
_rpc_launch_tasks(slurm_msg_t *msg, slurm_addr *cli)
{
	int      errnum = 0;
	uint16_t port;
	char     host[MAXHOSTNAMELEN];
	uid_t    req_uid;
	launch_tasks_request_msg_t *req = msg->data;
	uint32_t jobid  = req->job_id;
	uint32_t stepid = req->job_step_id;
	bool     super_user = false, run_prolog = false;
	slurm_addr self;
	socklen_t adlen;

	req_uid = g_slurm_auth_get_uid(msg->cred);

	super_user = _slurm_authorized_user(req_uid);

	if ((super_user == false) && (req_uid != req->uid)) {
		error("launch task request from uid %u",
		      (unsigned int) req_uid);
		errnum = ESLURM_USER_ID_MISSING;	/* or invalid user */
		goto done;
	}

	slurm_get_ip_str(cli, &port, host, sizeof(host));
	info("launch task %u.%u request from %u.%u@%s", req->job_id, 
	     req->job_step_id, req->uid, req->gid, host);

#ifndef HAVE_FRONT_END
	if (!slurm_cred_jobid_cached(conf->vctx, req->job_id)) 
		run_prolog = true;
#endif

	if (_check_job_credential(req->cred, jobid, stepid, req_uid,
			req->tasks_to_launch) < 0) {
		errnum = errno;
		error("Invalid job credential from %ld@%s: %m", 
		      (long) req_uid, host);
		goto done;
	}
	if (slurm_cred_revoked(conf->vctx, jobid)) {
		info("Job credential revoked for %u", jobid);
		errnum = ESLURMD_CREDENTIAL_REVOKED;
		goto done;
	}

	/* xassert(slurm_cred_jobid_cached(conf->vctx, req->job_id));*/

	/* Run job prolog if necessary */
	if (run_prolog && (_run_prolog(req->job_id, req->uid, NULL) != 0)) {
		error("[job %u] prolog failed", req->job_id);
		errnum = ESLURMD_PROLOG_FAILED;
		goto done;
	}

	adlen = sizeof(self);
	_slurm_getsockname(msg->conn_fd, (struct sockaddr *)&self, &adlen);
	errnum = _forkexec_slurmstepd(LAUNCH_TASKS, (void *)req, cli, &self);

    done:
	if (slurm_send_rc_msg(msg, errnum) < 0) {

		error("launch_tasks: unable to send return code: %m");

		/*
		 * Rewind credential so that srun may perform retry
		 */
		slurm_cred_rewind(conf->vctx, req->cred); /* ignore errors */

	} else if (errnum == SLURM_SUCCESS)
		save_cred_state(conf->vctx);

	/*
	 *  If job prolog failed, indicate failure to slurmctld
	 */
	if (errnum == ESLURMD_PROLOG_FAILED)
		send_registration_msg(errnum, false);	
}


static void 
_rpc_spawn_task(slurm_msg_t *msg, slurm_addr *cli)
{
	int      errnum = 0;
	uint16_t port;
	char     host[MAXHOSTNAMELEN];
	uid_t    req_uid;
	spawn_task_request_msg_t *req = msg->data;
	uint32_t jobid  = req->job_id;
	uint32_t stepid = req->job_step_id;
	bool     super_user = false, run_prolog = false;
	slurm_addr self;
	socklen_t adlen;
        int spawn_tasks_to_launch = -1;

	req_uid = g_slurm_auth_get_uid(msg->cred);

	super_user = _slurm_authorized_user(req_uid);

	if ((super_user == false) && (req_uid != req->uid)) {
		error("spawn task request from uid %u",
		      (unsigned int) req_uid);
		errnum = ESLURM_USER_ID_MISSING;	/* or invalid user */
		goto done;
	}

	slurm_get_ip_str(cli, &port, host, sizeof(host));
	info("spawn task %u.%u request from %u@%s", req->job_id, 
	     req->job_step_id, req->uid, host);

#ifndef HAVE_FRONT_END
	if (!slurm_cred_jobid_cached(conf->vctx, req->job_id)) 
		run_prolog = true;
#endif

	if (_check_job_credential(req->cred, jobid, stepid, req_uid, 
			spawn_tasks_to_launch) < 0) {
		errnum = errno;
		error("Invalid job credential from %ld@%s: %m", 
		      (long) req_uid, host);
		goto done;
	}
	if (slurm_cred_revoked(conf->vctx, jobid)) {
		info("Job credential revoked for %u", jobid);
		errnum = ESLURMD_CREDENTIAL_REVOKED;
		goto done;
	}

	/* xassert(slurm_cred_jobid_cached(conf->vctx, req->job_id));*/

	/* Run job prolog if necessary */
	if (run_prolog && (_run_prolog(req->job_id, req->uid, NULL) != 0)) {
		error("[job %u] prolog failed", req->job_id);
		errnum = ESLURMD_PROLOG_FAILED;
		goto done;
	}

	adlen = sizeof(self);
	_slurm_getsockname(msg->conn_fd, (struct sockaddr *)&self, &adlen);
	errnum = _forkexec_slurmstepd(SPAWN_TASKS, (void *)req, cli, &self);

    done:
	if (slurm_send_rc_msg(msg, errnum) < 0) {

		error("spawn_task: unable to send return code: %m");

		/*
		 * Rewind credential so that srun may perform retry
		 */
		slurm_cred_rewind(conf->vctx, req->cred); /* ignore errors */

	} else if (errnum == SLURM_SUCCESS)
		save_cred_state(conf->vctx);

	/*
	 *  If job prolog failed, indicate failure to slurmctld
	 */
	if (errnum == ESLURMD_PROLOG_FAILED)
		send_registration_msg(errnum, false);
}

static void
_prolog_error(batch_job_launch_msg_t *req, int rc)
{
	char *err_name_ptr, err_name[128], path_name[MAXPATHLEN];
	int fd;

	if (req->err)
		err_name_ptr = req->err;
	else {
		snprintf(err_name, sizeof(err_name), "slurm-%u.err", req->job_id);
		err_name_ptr = err_name;
	}
	if (err_name_ptr[0] == '/')
		snprintf(path_name, MAXPATHLEN, "%s", err_name_ptr);
	else if (req->work_dir)
		snprintf(path_name, MAXPATHLEN, "%s/%s", req->work_dir, err_name_ptr);
	else
		snprintf(path_name, MAXPATHLEN, "/%s", err_name_ptr);

	if ((fd = open(path_name, (O_CREAT|O_APPEND|O_WRONLY), 0644)) == -1) {
		error("Unable to open %s: %s", path_name, slurm_strerror(errno));
		return;
	}
	snprintf(err_name, 128, "Error running slurm prolog: %d\n", WEXITSTATUS(rc));
	write(fd, err_name, strlen(err_name));
	fchown(fd, (uid_t) req->uid, (gid_t) req->gid);
	close(fd);
}

static void
_rpc_batch_job(slurm_msg_t *msg, slurm_addr *cli)
{
	batch_job_launch_msg_t *req = (batch_job_launch_msg_t *)msg->data;
	bool     first_job_run = true;
	int      rc = SLURM_SUCCESS;
	uid_t    req_uid = g_slurm_auth_get_uid(msg->cred);
	char    *bgl_part_id = NULL;
	bool	replied = false;

	if (!_slurm_authorized_user(req_uid)) {
		error("Security violation, batch launch RPC from uid %u",
		      (unsigned int) req_uid);
		rc = ESLURM_USER_ID_MISSING;	/* or bad in this case */
		goto done;
	} 

	if (req->step_id != NO_VAL && req->step_id != 0)
		first_job_run = false;
	
	/*
	 * Insert jobid into credential context to denote that
	 * we've now "seen" an instance of the job
	 */
	if (first_job_run) {
		slurm_cred_insert_jobid(conf->vctx, req->job_id);

		/* 
	 	 * Run job prolog on this node
	 	 */
		select_g_get_jobinfo(req->select_jobinfo, SELECT_DATA_PART_ID, 
			&bgl_part_id);

#ifdef HAVE_BGL
		/* BlueGene prolog waits for partition boot and is very slow.
		 * Just reply now and send a separate kill job request if the 
		 * prolog or launch fail. */
		slurm_send_rc_msg(msg, rc);
		replied = true;
#endif

		rc = _run_prolog(req->job_id, req->uid, bgl_part_id);
		xfree(bgl_part_id);
		if (rc != 0) {
			error("[job %u] prolog failed", req->job_id);
			_prolog_error(req, rc);
			rc = ESLURMD_PROLOG_FAILED;
			goto done;
		}
	}

	/* Since job could have been killed while the prolog was 
	 * running (especially on BlueGene, which can wait  minutes
	 * for partition booting). Test if the credential has since
	 * been revoked and exit as needed. */
	if (slurm_cred_revoked(conf->vctx, req->job_id)) {
		info("Job %u already killed, do not launch tasks",  
			req->job_id);
		rc = ESLURMD_CREDENTIAL_REVOKED;     /* job already ran */
		goto done;
	}

	slurm_mutex_lock(&launch_mutex);
	if (req->step_id == NO_VAL)
		info("Launching batch job %u for UID %d",
			req->job_id, req->uid);
	else
		info("Launching batch job %u.%u for UID %d",
			req->job_id, req->step_id, req->uid);
	rc = _forkexec_slurmstepd(LAUNCH_BATCH_JOB, (void *)req, cli, NULL);
	slurm_mutex_unlock(&launch_mutex);

    done:
	if (!replied)
		slurm_send_rc_msg(msg, rc);
	else if (rc != 0) {
		/* prolog or job launch failure, 
		 * tell slurmctld that the job failed */
		(void) _abort_job(req->job_id);
	}
}

static int
_abort_job(uint32_t job_id)
{
	complete_job_step_msg_t  resp;
	slurm_msg_t resp_msg;

	resp.job_id       = job_id;
	resp.job_step_id  = NO_VAL;
	resp.job_rc       = 1;
	resp.slurm_rc     = 0;
	resp.node_name    = NULL;	/* unused */
	resp_msg.msg_type = REQUEST_COMPLETE_JOB_STEP;
	resp_msg.data     = &resp;
	return slurm_send_only_controller_msg(&resp_msg);
}

static void
_rpc_reconfig(slurm_msg_t *msg, slurm_addr *cli_addr)
{
	uid_t req_uid = g_slurm_auth_get_uid(msg->cred);

	if (!_slurm_authorized_user(req_uid))
		error("Security violation, reconfig RPC from uid %u",
		      (unsigned int) req_uid);
	else
		kill(conf->pid, SIGHUP);

	/* Never return a message, slurmctld does not expect one */
}

static void
_rpc_shutdown(slurm_msg_t *msg, slurm_addr *cli_addr)
{
	uid_t req_uid = g_slurm_auth_get_uid(msg->cred);

	if (!_slurm_authorized_user(req_uid))
		error("Security violation, shutdown RPC from uid %u",
		      (unsigned int) req_uid);
	else {
		if (kill(conf->pid, SIGTERM) != 0)
			error("kill(%u,SIGTERM): %m", conf->pid);
	}

	/* Never return a message, slurmctld does not expect one */
}

static int
_rpc_ping(slurm_msg_t *msg, slurm_addr *cli_addr)
{
	int        rc = SLURM_SUCCESS;
	uid_t req_uid = g_slurm_auth_get_uid(msg->cred);

	if (!_slurm_authorized_user(req_uid)) {
		error("Security violation, ping RPC from uid %u",
		      (unsigned int) req_uid);
		rc = ESLURM_USER_ID_MISSING;	/* or bad in this case */
	}

	/* Return result. If the reply can't be sent this indicates that
	 * 1. The network is broken OR
	 * 2. slurmctld has died    OR
	 * 3. slurmd was paged out due to full memory
	 * If the reply request fails, we send an registration message to 
	 * slurmctld in hopes of avoiding having the node set DOWN due to
	 * slurmd paging and not being able to respond in a timely fashion. */
	if (slurm_send_rc_msg(msg, rc) < 0) {
		error("Error responding to ping: %m");
		send_registration_msg(SLURM_SUCCESS, false);
	}
	return rc;
}

static void
_rpc_signal_tasks(slurm_msg_t *msg, slurm_addr *cli_addr)
{
	int               rc = SLURM_SUCCESS;
	uid_t             req_uid;
	kill_tasks_msg_t *req = (kill_tasks_msg_t *) msg->data;
	step_loc_t loc;

	/*
	 * Use slurmstepd API to request that the slurmdstepd handle
	 * the signalling.
	 */
	loc.directory = conf->spooldir;
	loc.nodename = "nodename";
	loc.jobid = req->job_id;
	loc.stepid = req->job_step_id;

#ifdef HAVE_AIX
#  ifdef SIGMIGRATE
#    ifdef SIGSOUND
	/* SIGMIGRATE and SIGSOUND are used to initiate job checkpoint on AIX.
	 * These signals are not sent to the entire process group, but just a
	 * single process, namely the PMD. */
	if (req->signal == SIGMIGRATE || req->signal == SIGSOUND) {
		rc = stepd_signal_task_local(loc, msg->cred, req->signal, 0);
		goto done;
	}
#    endif
#  endif
#endif

	rc = stepd_signal(loc, msg->cred, req->signal);
	
  done:
	slurm_send_rc_msg(msg, rc);
}

static void
_rpc_terminate_tasks(slurm_msg_t *msg, slurm_addr *cli_addr)
{
	kill_tasks_msg_t *req = (kill_tasks_msg_t *) msg->data;
	step_loc_t        loc;
	int               rc = SLURM_SUCCESS;

	debug3("Entering _rpc_terminate_tasks");
	/*
	 * Use slurmstepd API to request that the slurmdstepd handle
	 * the signalling.
	 */
	loc.directory = conf->spooldir;
	loc.nodename = "nodename";
	loc.jobid = req->job_id;
	loc.stepid = req->job_step_id;
	rc = stepd_signal_container(loc, msg->cred, req->signal);

  done:
	slurm_send_rc_msg(msg, rc);
}

/* 
 *  For the specified job_id: reply to slurmctld, 
 *   sleep(configured kill_wait), then send SIGKILL 
 *  FIXME! - Perhaps we should send SIGXCPU first?
 */
static void
_rpc_timelimit(slurm_msg_t *msg, slurm_addr *cli_addr)
{
	uid_t           uid = g_slurm_auth_get_uid(msg->cred);
	kill_job_msg_t *req = msg->data;
	int             nsteps;

	if (!_slurm_authorized_user(uid)) {
		error ("Security violation: rpc_timelimit req from uid %ld", 
		       (long) uid);
		slurm_send_rc_msg(msg, ESLURM_USER_ID_MISSING);
		return;
	}

	/*
	 *  Indicate to slurmctld that we've received the message
	 */
	slurm_send_rc_msg(msg, SLURM_SUCCESS);
	slurm_close_accepted_conn(msg->conn_fd);
	msg->conn_fd = -1;

	nsteps = _kill_all_active_steps(msg->cred, req->job_id, SIGTERM, false);
	verbose( "Job %u: timeout: sent SIGTERM to %d active steps", 
	         req->job_id, nsteps );

	/* Revoke credential, send SIGKILL, run epilog, etc. */
	_rpc_terminate_job(msg, cli_addr); 
}

static void  _rpc_pid2jid(slurm_msg_t *msg, slurm_addr *cli)
{
	job_id_request_msg_t *req = (job_id_request_msg_t *) msg->data;
	slurm_msg_t           resp_msg;
	job_id_response_msg_t resp;
	bool         found = false; 
	List         steps;
	ListIterator i;
	step_loc_t *stepd;

	steps = stepd_available(conf->spooldir, "nodename");
	i = list_iterator_create(steps);
	while (stepd = list_next(i)) {
		if (stepd_pid_in_container(*stepd, req->job_pid)
		    || req->job_pid == stepd_daemon_pid(*stepd)) {
			resp.job_id = stepd->jobid;
			found = true;
			break;
		}
	}
	list_iterator_destroy(i);
	list_destroy(steps);

	if (found) {
		debug3("_rpc_pid2jid: pid(%u) found in %u",
		       req->job_pid, resp.job_id);
		resp_msg.address      = msg->address;
		resp_msg.msg_type     = RESPONSE_JOB_ID;
		resp_msg.data         = &resp;
		slurm_send_node_msg(msg->conn_fd, &resp_msg);
	} else {
		/* We could possibly scan the proc table and figure 
		 * out which job this pid belongs to, but for now 
		 * we only handle the job's top level pid */
		debug3("_rpc_pid2jid: pid(%u) not found", req->job_pid);
		slurm_send_rc_msg(msg, ESLURM_INVALID_JOB_ID);
	}
}


static void 
_rpc_reattach_tasks(slurm_msg_t *msg, slurm_addr *cli)
{
	reattach_tasks_request_msg_t  *req = msg->data;
	reattach_tasks_response_msg_t  resp;
	slurm_msg_t                    resp_msg;
	int          rc   = SLURM_SUCCESS;
	uint16_t     port = 0;
	char         host[MAXHOSTNAMELEN];
	int          i;
	slurm_addr   ioaddr;
	step_loc_t loc;

	memset(&resp, 0, sizeof(reattach_tasks_response_msg_t));
	slurm_get_ip_str(cli, &port, host, sizeof(host));

	/* 
	 * Set response addr by resp_port and client address
	 */
	memcpy(&resp_msg.address, cli, sizeof(slurm_addr));
	slurm_set_addr(&resp_msg.address, req->resp_port, NULL); 

	/* 
	 * Set IO address and by io_port and client address
	 */
	memcpy(&ioaddr, cli, sizeof(slurm_addr));
	slurm_set_addr(&ioaddr, req->io_port, NULL);

	/*
	 * Use slurmstepd API to request that the slurmdstepd handle
	 * the attach.
	 */
	loc.directory = conf->spooldir;
	loc.nodename = "nodename";
	loc.jobid = req->job_id;
	loc.stepid = req->job_step_id;
	rc = stepd_attach(loc, &ioaddr, &resp_msg.address,
			  msg->cred, req->cred, &resp);
	if (rc != SLURM_SUCCESS) {
		debug2("stepd_attach call failed");
		goto done;
	}

    done:
	debug2("update step addrs rc = %d", rc);
	resp_msg.data         = &resp;
	resp_msg.msg_type     = RESPONSE_REATTACH_TASKS;
	resp.node_name        = conf->node_name;
	resp.srun_node_id     = req->srun_node_id;
	resp.return_code      = rc;

	slurm_send_only_node_msg(&resp_msg);

	xfree(resp.gtids);
	xfree(resp.local_pids);

}

/*
 * _kill_all_active_steps - signals all steps of a job
 * jobid IN - id of job to signal
 * sig   IN - signal to send
 * batch IN - if true signal batch script, otherwise skip it
 * RET count of signaled job steps (plus batch script, if applicable)
 */
static int
_kill_all_active_steps(void *auth_cred, uint32_t jobid, int sig, bool batch)
{
	List steps;
	ListIterator i;
	step_loc_t *stepd;
	int step_cnt  = 0;  

	steps = stepd_available(conf->spooldir, "nodename");
	i = list_iterator_create(steps);
	while (stepd = list_next(i)) {
		if (stepd->jobid != jobid) {
			/* multiple jobs expected on shared nodes */
			debug3("Step from other job: jobid=%d (this jobid=%d)",
			      stepd->jobid, jobid);
			continue;
		}

		if ((stepd->stepid == NO_VAL) && (!batch))
			continue;

		step_cnt++;

		debug2("container signal %d to job %u.%u",
		       sig, jobid, stepd->stepid);
		if (stepd_signal_container(*stepd, auth_cred, sig) < 0)
			error("kill jid %d: %m", jobid);
	}
	list_iterator_destroy(i);
	list_destroy(steps);
	if (step_cnt == 0)
		debug2("No steps in jobid %d to send signal %d", jobid, sig);
	return step_cnt;
}

static bool
_job_still_running(uint32_t job_id)
{
	bool         retval = false;
	List         steps;
	ListIterator i;
	step_loc_t  *s     = NULL;   

	steps = stepd_available(conf->spooldir, "nodename");
	i = list_iterator_create(steps);
	while ((s = list_next(i))) {
		if (s->jobid == job_id
		    && stepd_state(*s) != SLURMSTEPD_NOT_RUNNING) {
			retval = true;
			break;
		}
	}
	list_iterator_destroy(i);
	list_destroy(steps);

	return retval;
}

/*
 * Wait until all job steps are in SLURMD_JOB_COMPLETE state.
 * This indicates that interconnect_postfini has completed and 
 * freed the switch windows (as needed only for Federation switch).
 */
static void
_wait_state_completed(uint32_t jobid, int max_delay)
{
	char *switch_type = slurm_get_switch_type();
	int i;

	if (strcmp(switch_type, "switch/federation")) {
		xfree(switch_type);
		return;
	}
	xfree(switch_type);

	for (i=0; i<max_delay; i++) {
		if (_steps_completed_now(jobid))
			break;
		sleep(1);
	}
	if (i >= max_delay)
		error("timed out waiting for job %u to complete", jobid);
}

static bool
_steps_completed_now(uint32_t jobid)
{
	List steps;
	ListIterator i;
	step_loc_t *stepd;
	bool rc = true;

	steps = stepd_available(conf->spooldir, "nodename");
	i = list_iterator_create(steps);
	while (stepd = list_next(i)) {
		if (stepd->jobid == jobid
		    && stepd_state(*stepd) != SLURMSTEPD_NOT_RUNNING) {
			rc = false;
			break;
		}
	}
	list_iterator_destroy(i);
	list_destroy(steps);

	return rc;
}

/*
 *  Send epilog complete message to currently active comtroller.
 *   Returns SLURM_SUCCESS if message sent successfully,
 *           SLURM_FAILURE if epilog complete message fails to be sent.
 */
static int
_epilog_complete(uint32_t jobid, int rc)
{
	int                    ret = SLURM_SUCCESS;
	slurm_msg_t            msg;
	epilog_complete_msg_t  req;

	_wait_state_completed(jobid, 5);

	req.job_id      = jobid;
	req.return_code = rc;
	req.node_name   = conf->node_name;
	if (switch_g_alloc_node_info(&req.switch_nodeinfo))
		error("switch_g_alloc_node_info: %m");
	if (switch_g_build_node_info(req.switch_nodeinfo))
		error("switch_g_build_node_info: %m");

	msg.msg_type    = MESSAGE_EPILOG_COMPLETE;
	msg.data        = &req;

	if (slurm_send_only_controller_msg(&msg) < 0) {
		error("Unable to send epilog complete message: %m");
		ret = SLURM_ERROR;
	} else
		debug ("Job %u: sent epilog complete msg: rc = %d", jobid, rc);

	switch_g_free_node_info(&req.switch_nodeinfo);
	return ret;
}


/* FIXME - kill_job_msg_t doesn't have a signal number in the payload
 *         so it won't suffice for this RPC.  Fortunately, no one calls
 *         this RPC.
 */
static void 
_rpc_signal_job(slurm_msg_t *msg, slurm_addr *cli)
{
	int             rc     = SLURM_SUCCESS;
	kill_job_msg_t *req    = msg->data;
	uid_t           uid    = g_slurm_auth_get_uid(msg->cred);
	int             nsteps = 0;
	int		delay;
	char           *bgl_part_id = NULL;

	error("_rpc_signal_job not yet implemented");
	/* 
	 * check that requesting user ID is the SLURM UID
	 */
	if (!_slurm_authorized_user(uid)) {
		error("Security violation: kill_job(%ld) from uid %ld",
		      req->job_id, (long) uid);
		if (msg->conn_fd >= 0)
			slurm_send_rc_msg(msg, ESLURM_USER_ID_MISSING);
		return;
	} 

	/*nsteps = _kill_all_active_steps(req->job_id, SIGTERM, true);*/

	/*
	 *  At this point, if connection still open, we send controller
	 *   a "success" reply to indicate that we've recvd the msg.
	 */
	if (msg->conn_fd >= 0) {
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
		if (slurm_close_accepted_conn(msg->conn_fd) < 0)
			error ("_rpc_signal_job: close(%d): %m", msg->conn_fd);
		msg->conn_fd = -1;
	}
}

static void 
_rpc_terminate_job(slurm_msg_t *msg, slurm_addr *cli)
{
	int             rc     = SLURM_SUCCESS;
	kill_job_msg_t *req    = msg->data;
	uid_t           uid    = g_slurm_auth_get_uid(msg->cred);
	int             nsteps = 0;
	int		delay;
	char           *bgl_part_id = NULL;

	/* 
	 * check that requesting user ID is the SLURM UID
	 */
	if (!_slurm_authorized_user(uid)) {
		error("Security violation: kill_job(%ld) from uid %ld",
		      req->job_id, (long) uid);
		if (msg->conn_fd >= 0)
			slurm_send_rc_msg(msg, ESLURM_USER_ID_MISSING);
		return;
	} 

	/*
	 *  Initialize a "waiter" thread for this jobid. If another
	 *   thread is already waiting on termination of this job, 
	 *   _waiter_init() will return SLURM_ERROR. In this case, just 
	 *   notify slurmctld that we recvd the message successfully,
	 *   then exit this thread.
	 */
	if (_waiter_init(req->job_id) == SLURM_ERROR) {
		if (msg->conn_fd >= 0)
			slurm_send_rc_msg (msg, SLURM_SUCCESS);
		return;
	}


	/*
	 * "revoke" all future credentials for this jobid
	 */
	if (slurm_cred_revoke(conf->vctx, req->job_id) < 0) {
		debug("revoking cred for job %u: %m", req->job_id);
	} else {
		save_cred_state(conf->vctx);
		debug("credential for job %u revoked", req->job_id);
	}

	/*
	 * Tasks might be stopped (possibly by a debugger)
	 * so send SIGCONT first.
	 */
	_kill_all_active_steps(msg->cred, req->job_id, SIGCONT, true);

	nsteps = _kill_all_active_steps(msg->cred, req->job_id, SIGTERM, true);

	/*
	 *  If there are currently no active job steps and no
	 *    configured epilog to run, bypass asynchronous reply and
	 *    notify slurmctld that we have already completed this
	 *    request. We need to send current switch state on AIX
	 *    systems, so this bypass can not be used.
	 */
#ifndef HAVE_AIX
	if ((nsteps == 0) && !conf->epilog) {
		if (msg->conn_fd >= 0)
			slurm_send_rc_msg(msg, 
				ESLURMD_KILL_JOB_ALREADY_COMPLETE);
		slurm_cred_begin_expiration(conf->vctx, req->job_id);
		_waiter_complete(req->job_id);
		return;
	}
#endif

	/*
	 *  At this point, if connection still open, we send controller
	 *   a "success" reply to indicate that we've recvd the msg.
	 */
	if (msg->conn_fd >= 0) {
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
		if (slurm_close_accepted_conn(msg->conn_fd) < 0)
			error ("rpc_kill_job: close(%d): %m", msg->conn_fd);
		msg->conn_fd = -1;
		
	}

	/*
	 *  Check for corpses
	 */
	delay = MAX(conf->cf.kill_wait, 5);
	if ( !_pause_for_job_completion (msg->cred, req->job_id, delay)
	     && _kill_all_active_steps(msg->cred, req->job_id, SIGKILL, true) ) {
		/*
		 *  Block until all user processes are complete.
		 */
		_pause_for_job_completion (msg->cred, req->job_id, 0);
	}

	/*
	 *  Begin expiration period for cached information about job.
	 *   If expiration period has already begun, then do not run
	 *   the epilog again, as that script has already been executed 
	 *   for this job.
	 */
	if (slurm_cred_begin_expiration(conf->vctx, req->job_id) < 0) {
		debug("Not running epilog for jobid %d: %m", req->job_id);
		goto done;
	}

	save_cred_state(conf->vctx);
	select_g_get_jobinfo(req->select_jobinfo, SELECT_DATA_PART_ID,
		&bgl_part_id);
	rc = _run_epilog(req->job_id, req->job_uid, bgl_part_id);
	xfree(bgl_part_id);
	if (rc != 0) {
		error ("[job %u] epilog failed", req->job_id);
		rc = ESLURMD_EPILOG_FAILED;
	} else
		debug("completed epilog for jobid %u", req->job_id);
	
    done:
	_epilog_complete(req->job_id, rc);
	_waiter_complete(req->job_id);
}

/*
 *  Returns true if "uid" is a "slurm authorized user" - i.e. uid == 0
 *   or uid == slurm user id at this time.
 */
static bool
_slurm_authorized_user(uid_t uid)
{
	return ((uid == (uid_t) 0) || (uid == conf->slurm_user_id));
}


struct waiter {
	uint32_t jobid;
	pthread_t thd;
};


static struct waiter *
_waiter_create(uint32_t jobid)
{
	struct waiter *wp = xmalloc(sizeof(struct waiter));

	wp->jobid = jobid;
	wp->thd   = pthread_self();

	return wp;
}

static int _find_waiter(struct waiter *w, uint32_t *jp)
{
	return (w->jobid == *jp);
}

static void _waiter_destroy(struct waiter *wp)
{
	xfree(wp);
}

static int _waiter_init (uint32_t jobid)
{
	if (!waiters)
		waiters = list_create((ListDelF) _waiter_destroy);

	/* 
	 *  Exit this thread if another thread is waiting on job
	 */
	if (list_find_first (waiters, (ListFindF) _find_waiter, &jobid))
		return SLURM_ERROR;
	else 
		list_append(waiters, _waiter_create(jobid));

	return (SLURM_SUCCESS);
}

static int _waiter_complete (uint32_t jobid)
{
	return (list_delete_all (waiters, (ListFindF) _find_waiter, &jobid));
}

/*
 *  Like _wait_for_procs(), but only wait for up to max_time seconds
 *  if max_time == 0, send SIGKILL to tasks repeatedly
 *    
 *  Returns true if all job 
 */
static bool
_pause_for_job_completion (void *auth_cred, uint32_t job_id, int max_time)
{
	int sec = 0, rc = 0;

	while ( ((sec++ < max_time) || (max_time == 0))
	      && (rc = _job_still_running (job_id))) {
		if ((max_time == 0) && (sec > 1))
			_kill_all_active_steps(auth_cred, job_id,
					       SIGKILL, true);
		sleep (1);
	}
	/* 
	 * Return true if job is NOT running
	 */
	return (!rc);
}

static void 
_rpc_update_time(slurm_msg_t *msg, slurm_addr *cli)
{
	int   rc      = SLURM_SUCCESS;
	uid_t req_uid = g_slurm_auth_get_uid(msg->cred);
	job_time_msg_t *req = (job_time_msg_t *) msg->data;

	if ((req_uid != conf->slurm_user_id) && (req_uid != 0)) {
		rc = ESLURM_USER_ID_MISSING;
		error("Security violation, uid %u can't update time limit",
		      (unsigned int) req_uid);
		goto done;
	} 

/* FIXME */
/* 	if (shm_update_job_timelimit(req->job_id, req->expiration_time) < 0) { */
/* 		error("updating lifetime for job %u: %m", req->job_id); */
/* 		rc = ESLURM_INVALID_JOB_ID; */
/* 	} else */
/* 		debug("reset job %u lifetime", req->job_id); */

    done:
	slurm_send_rc_msg(msg, rc);
}

/* NOTE: xfree returned value */
static char **
_build_env(uint32_t jobid, uid_t uid, char *bgl_part_id)
{
	char **env = xmalloc(sizeof(char *));
	env[0]  = NULL;
	setenvf(&env, "SLURM_JOBID", "%u", jobid);
	setenvf(&env, "SLURM_UID",   "%u", uid);
	if (bgl_part_id) {
		setenvf(&env, "MPIRUN_PARTITION",
			"%s", bgl_part_id);
	}
	return env;
}

static int 
_run_prolog(uint32_t jobid, uid_t uid, char *bgl_part_id)
{
	int error_code;
	char *my_prolog;
	char **my_env = _build_env(jobid, uid, bgl_part_id);

	slurm_mutex_lock(&conf->config_mutex);
	my_prolog = xstrdup(conf->prolog);
	slurm_mutex_unlock(&conf->config_mutex);

	error_code = run_script("prolog", my_prolog, jobid, uid, 
			-1, my_env);
	xfree(my_prolog);
	xfree(my_env);

	return error_code;
}

static int 
_run_epilog(uint32_t jobid, uid_t uid, char *bgl_part_id)
{
	int error_code;
	char *my_epilog;
	char **my_env = _build_env(jobid, uid, bgl_part_id);

	slurm_mutex_lock(&conf->config_mutex);
	my_epilog = xstrdup(conf->epilog);
	slurm_mutex_unlock(&conf->config_mutex);

	error_code = run_script("epilog", my_epilog, jobid, uid, 
			-1, my_env);
	xfree(my_epilog);
	xfree(my_env);

	return error_code;
}
