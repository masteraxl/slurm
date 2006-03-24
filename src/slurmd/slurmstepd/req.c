/*****************************************************************************\
 *  src/slurmd/slurmstepd/req.c - slurmstepd domain socket request handling
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher Morrone <morrone2@llnl.gov>
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "src/common/xstring.h"
#include "src/common/xmalloc.h"
#include "src/common/fd.h"
#include "src/common/eio.h"
#include "src/common/slurm_auth.h"

#include "src/slurmd/slurmd/slurmd.h"
#include "src/slurmd/common/stepd_api.h"
#include "src/slurmd/slurmstepd/slurmstepd.h"
#include "src/slurmd/slurmstepd/slurmstepd_job.h"
#include "src/slurmd/slurmstepd/req.h"

static void *_handle_accept(void *arg);
static int _handle_request(int fd, slurmd_job_t *job, uid_t uid, gid_t gid);
static int _handle_state(int fd, slurmd_job_t *job);
static int _handle_info(int fd, slurmd_job_t *job);
static int _handle_signal_process_group(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_signal_task_local(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_signal_container(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_attach(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_pid_in_container(int fd, slurmd_job_t *job);
static int _handle_daemon_pid(int fd, slurmd_job_t *job);
static int _handle_suspend(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_resume(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_terminate(int fd, slurmd_job_t *job, uid_t uid);
static int _handle_completion(int fd, slurmd_job_t *job, uid_t uid);
static bool _msg_socket_readable(eio_obj_t *obj);
static int _msg_socket_accept(eio_obj_t *obj, List objs);

struct io_operations msg_socket_ops = {
	readable:	&_msg_socket_readable,
	handle_read:	&_msg_socket_accept
};

static char *socket_name;
static pthread_mutex_t suspend_mutex = PTHREAD_MUTEX_INITIALIZER;
static bool suspended = false;

struct request_params {
	int fd;
	slurmd_job_t *job;
};

/*
 *  Returns true if "uid" is a "slurm authorized user" - i.e. uid == 0
 *   or uid == slurm user id at this time.
 */
static bool
_slurm_authorized_user(uid_t uid)
{
	return ((uid == (uid_t) 0) || (uid == conf->slurm_user_id));
}

/*
 * Create a named unix domain listening socket.
 * (cf, Stevens APUE 1st ed., section 15.5.2)
 */
static int
_create_socket(const char *name)
{
	int fd;
	int len;
	struct sockaddr_un addr;

	/* create a unix domain stream socket */
	if ((fd = socket(AF_UNIX, SOCK_STREAM, 0)) < 0)
		return -1;
	fd_set_close_on_exec(fd);

	memset(&addr, 0, sizeof(addr));
	addr.sun_family = AF_UNIX;
	strcpy(addr.sun_path, name);
	len = strlen(addr.sun_path)+1 + sizeof(addr.sun_family);

	/* bind the name to the descriptor */
	if (bind(fd, (struct sockaddr *) &addr, len) < 0)
		return -2;

	if (listen(fd, 5) < 0)
		return -3;

	return fd;
}

static int
_domain_socket_create(const char *dir, const char *nodename,
		     uint32_t jobid, uint32_t stepid)
{
	int fd;
	char *name = NULL;
	struct stat stat_buf;

	/*
	 * Make sure that "dir" exists and is a directory.
	 */
	if (stat(dir, &stat_buf) < 0) {
		error("Domain socket directory %s: %m", dir);
		return -1;
	} else if (!S_ISDIR(stat_buf.st_mode)) {
		error("%s is not a directory", dir);
		return -1;
	}

	/*
	 * Now build the the name of socket, and create the socket.
	 */
	xstrfmtcat(name, "%s/%s_%u.%u", dir, nodename, jobid, stepid);

	/*
	 * First check to see if the named socket already exists.
	 */
	if (stat(name, &stat_buf) == 0) {
		error("Socket %s already exists", name);
		xfree(name);
		errno = ESLURMD_STEP_EXISTS;
		return -1;
	}

	fd = _create_socket(name);
	if (fd < 0)
		fatal("Could not create domain socket: %m");

	chmod(name, 0777);
	socket_name = name;
	return fd;
}

static void
_domain_socket_destroy(int fd)
{
	if (close(fd) < 0)
		error("Unable to close domain socket: %m");

	if (unlink(socket_name) == -1)
		error("Unable to unlink domain socket: %m");
}


static void *
_msg_thr_internal(void *job_arg)
{
	slurmd_job_t *job = (slurmd_job_t *) job_arg;

	debug("Message thread started pid = %lu", (unsigned long) getpid());
	eio_handle_mainloop(job->msg_handle);
	debug("Message thread exited");
}

int
msg_thr_create(slurmd_job_t *job)
{
	int fd;
	eio_obj_t *eio_obj;
	pthread_attr_t attr;

	errno = 0;
	fd = _domain_socket_create(conf->spooldir, conf->node_name,
				  job->jobid, job->stepid);
	if (fd == -1)
		return SLURM_ERROR;

	fd_set_nonblocking(fd);

	eio_obj = eio_obj_create(fd, &msg_socket_ops, (void *)job);
	job->msg_handle = eio_handle_create();
	eio_new_initial_obj(job->msg_handle, eio_obj);

	slurm_attr_init(&attr);
	if (pthread_create(&job->msgid, &attr,
			   &_msg_thr_internal, (void *)job) != 0) {
		error("pthread_create: %m");
		return SLURM_ERROR;
	}

	return SLURM_SUCCESS;
}

static bool 
_msg_socket_readable(eio_obj_t *obj)
{
	debug3("Called _msg_socket_readable");
	if (obj->shutdown == true) {
		if (obj->fd != -1) {
			debug2("  false, shutdown");
			_domain_socket_destroy(obj->fd);
			obj->fd = -1;
		} else {
			debug2("  false");
		}
		return false;
	}
	return true;
}

static int
_msg_socket_accept(eio_obj_t *obj, List objs)
{
	slurmd_job_t *job = (slurmd_job_t *)obj->arg;
	int fd;
	struct sockaddr_un addr;
	int len = sizeof(addr);
	struct request_params *param = NULL;
	pthread_attr_t attr;
	pthread_t id;

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

	slurm_attr_init(&attr);
	if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED) != 0) {
		close(fd);
		error("Unable to set detachstate on attr: %m");
		return SLURM_ERROR;
	}

	param = xmalloc(sizeof(struct request_params));
	param->fd = fd;
	param->job = job;
	if (pthread_create(&id, &attr, &_handle_accept, (void *)param) != 0) {
		error("stepd_api message engine pthread_create: %m");
		_handle_accept((void *)param);
	}
	param = NULL;

	debug3("Leaving _msg_socket_accept");
	return SLURM_SUCCESS;
}

static void *
_handle_accept(void *arg)
{
	/*struct request_params *param = (struct request_params *)arg;*/
	int fd = ((struct request_params *)arg)->fd;
	slurmd_job_t *job = ((struct request_params *)arg)->job;
	int req;
	int len;
	Buf buffer;
	void *auth_cred;
	int rc;
	uid_t uid;
	gid_t gid;

	debug3("Entering _handle_accept (new thread)");
	xfree(arg);

	safe_read(fd, &req, sizeof(int));
	if (req != REQUEST_CONNECT) {
		error("First message must be REQUEST_CONNECT");
		goto fail;
	}

	safe_read(fd, &len, sizeof(int));
	buffer = init_buf(len);
	safe_read(fd, get_buf_data(buffer), len);

	/* Unpack and verify the auth credential */
	auth_cred = g_slurm_auth_unpack(buffer);
	if (auth_cred == NULL) {
		error("Unpacking authentication credential: %s",
		      g_slurm_auth_errstr(g_slurm_auth_errno(NULL)));
		free_buf(buffer);
		goto fail;
	}
	rc = g_slurm_auth_verify(auth_cred, NULL, 2);
	if (rc != SLURM_SUCCESS) {
		error("Verifying authentication credential: %s",
		      g_slurm_auth_errstr(g_slurm_auth_errno(auth_cred)));
		(void) slurm_free_cred(auth_cred);
		free_buf(buffer);
		goto fail;
	}

	/* Get the uid & gid from the credential, then destroy it. */
	uid = g_slurm_auth_get_uid(auth_cred);
	gid = g_slurm_auth_get_gid(auth_cred);
	debug3("  Identity: uid=%d, gid=%d", uid, gid);
	slurm_free_cred(auth_cred);
	free_buf(buffer);

	rc = SLURM_SUCCESS;
	safe_write(fd, &rc, sizeof(int));
	
	while (1) {
		rc = _handle_request(fd, job, uid, gid);
		if (rc != SLURM_SUCCESS)
			break;
	}

	if (close(fd) == -1)
		error("Closing accepted fd: %m");
	debug3("Leaving  _handle_accept");
	return NULL;

fail:
	rc = SLURM_FAILURE;
	safe_write(fd, &rc, sizeof(int));
rwfail:
	if (close(fd) == -1)
		error("Closing accepted fd after error: %m");
	debug("Leaving  _handle_accept on an error");
	return NULL;
}


int
_handle_request(int fd, slurmd_job_t *job, uid_t uid, gid_t gid)
{
	int rc = 0;
	int req;

	debug3("Entering _handle_request");
	if ((rc = read(fd, &req, sizeof(int))) != sizeof(int)) {
		if (rc == 0) { /* EOF, normal */
			return -1;
		} else {
			debug3("Leaving _handle_request on read error");
			return SLURM_FAILURE;
		}
	}
	debug3("Got request");
	rc = SLURM_SUCCESS;
	switch (req) {
	case REQUEST_SIGNAL_PROCESS_GROUP:
		debug("Handling REQUEST_SIGNAL_PROCESS_GROUP");
		rc = _handle_signal_process_group(fd, job, uid);
		break;
	case REQUEST_SIGNAL_TASK_LOCAL:
		debug("Handling REQUEST_SIGNAL_TASK_LOCAL");
		rc = _handle_signal_task_local(fd, job, uid);
		break;
	case REQUEST_SIGNAL_TASK_GLOBAL:
		debug("Handling REQUEST_SIGNAL_TASK_LOCAL (not implemented)");
		break;
	case REQUEST_SIGNAL_CONTAINER:
		debug("Handling REQUEST_SIGNAL_CONTAINER");
		rc = _handle_signal_container(fd, job, uid);
		break;
	case REQUEST_STATE:
		debug("Handling REQUEST_STATE");
		rc = _handle_state(fd, job);
		break;
	case REQUEST_INFO:
		debug("Handling REQUEST_INFO");
		rc = _handle_info(fd, job);
		break;
	case REQUEST_ATTACH:
		debug("Handling REQUEST_ATTACH");
		rc = _handle_attach(fd, job, uid);
		break;
	case REQUEST_PID_IN_CONTAINER:
		debug("Handling REQUEST_PID_IN_CONTAINER");
		rc = _handle_pid_in_container(fd, job);
		break;
	case REQUEST_DAEMON_PID:
		debug("Handling REQUEST_DAEMON_PID");
		rc = _handle_daemon_pid(fd, job);
		break;
	case REQUEST_STEP_SUSPEND:
		debug("Handling REQUEST_STEP_SUSPEND");
		rc = _handle_suspend(fd, job, uid);
		break;
	case REQUEST_STEP_RESUME:
		debug("Handling REQUEST_STEP_RESUME");
		rc = _handle_resume(fd, job, uid);
		break;
	case REQUEST_STEP_TERMINATE:
		debug("Handling REQUEST_STEP_TERMINATE");
		rc = _handle_terminate(fd, job, uid);
		break;
	case REQUEST_STEP_COMPLETION:
		debug("Handling REQUEST_STEP_COMPLETION");
		rc = _handle_completion(fd, job, uid);
		break;
	default:
		error("Unrecognized request: %d", req);
		rc = SLURM_FAILURE;
		break;
	}

	debug3("Leaving  _handle_request: %s",
	       rc ? "SLURM_FAILURE" : "SLURM_SUCCESS");
	return rc;
}

static int
_handle_state(int fd, slurmd_job_t *job)
{
	safe_write(fd, &job->state, sizeof(slurmstepd_state_t));

	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_info(int fd, slurmd_job_t *job)
{
	safe_write(fd, &job->uid, sizeof(uid_t));
	safe_write(fd, &job->jobid, sizeof(uint32_t));
	safe_write(fd, &job->stepid, sizeof(uint32_t));

	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_signal_process_group(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int signal;

	debug("_handle_signal_process_group for job %u.%u",
	      job->jobid, job->stepid);

	safe_read(fd, &signal, sizeof(int));

	debug3("  uid = %d", uid);
	if (uid != job->uid && !_slurm_authorized_user(uid)) {
		debug("kill req from uid %ld for job %u.%u owned by uid %ld",
		      (long)uid, job->jobid, job->stepid, (long)job->uid);
		rc = EPERM;
		goto done;
	}

	/*
	 * Sanity checks
	 */
	if (job->pgid <= (pid_t)1) {
		debug ("step %u.%u invalid [jmgr_pid:%d pgid:%u]", 
                       job->jobid, job->stepid, job->jmgr_pid, job->pgid);
		rc = ESLURMD_JOB_NOTRUNNING;
		goto done;
	}

	/*
	 * Signal the process group
	 */
	pthread_mutex_lock(&suspend_mutex);
	if (suspended) {
		rc = ESLURMD_STEP_SUSPENDED;
		pthread_mutex_unlock(&suspend_mutex);
		goto done;
	}

	if (killpg(job->pgid, signal) == -1) {
		rc = -1;
		verbose("Error sending signal %d to %u.%u, pgid %d: %s", 
			signal, job->jobid, job->stepid, job->pgid,
			slurm_strerror(rc));
	} else {
		verbose("Sent signal %d to %u.%u, pgid %d", 
			signal, job->jobid, job->stepid, job->pgid);
	}
	pthread_mutex_unlock(&suspend_mutex);
	
done:
	/* Send the return code */
	safe_write(fd, &rc, sizeof(int));
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_signal_task_local(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int signal;
	int ltaskid; /* local task index */

	debug("_handle_signal_task_local for job %u.%u",
	      job->jobid, job->stepid);

	safe_read(fd, &signal, sizeof(int));
	safe_read(fd, &ltaskid, sizeof(int));

	debug3("  uid = %d", uid);
	if (uid != job->uid && !_slurm_authorized_user(uid)) {
		debug("kill req from uid %ld for job %u.%u owned by uid %ld",
		      (long)uid, job->jobid, job->stepid, (long)job->uid);
		rc = EPERM;
		goto done;
	}

	/*
	 * Sanity checks
	 */
	if (ltaskid < 0 || ltaskid >= job->ntasks) {
		debug("step %u.%u invalid local task id %d", 
		      job->jobid, job->stepid, ltaskid);
		rc = SLURM_ERROR;
		goto done;
	}
	if (!job->task
	    || !job->task[ltaskid]) {
		debug("step %u.%u no task info for task id %d",
		      job->jobid, job->stepid, ltaskid);
		rc = SLURM_ERROR;
		goto done;
	}
	if (job->task[ltaskid]->pid <= 1) {
		debug("step %u.%u invalid pid %d for task %d",
		      job->jobid, job->stepid,
		      job->task[ltaskid]->pid, ltaskid);
		rc = SLURM_ERROR;
		goto done;
	}

	/*
	 * Signal the task
	 */
	pthread_mutex_lock(&suspend_mutex);
	if (suspended) {
		rc = ESLURMD_STEP_SUSPENDED;
		pthread_mutex_unlock(&suspend_mutex);
		goto done;
	}

	if (kill(job->task[ltaskid]->pid, signal) == -1) {
		rc = -1;
		verbose("Error sending signal %d to %u.%u, pid %d: %s", 
			signal, job->jobid, job->stepid,
			job->task[ltaskid]->pid, slurm_strerror(rc));
	} else {
		verbose("Sent signal %d to %u.%u, pid %d", 
			signal, job->jobid, job->stepid,
			job->task[ltaskid]->pid);
	}
	pthread_mutex_unlock(&suspend_mutex);
	
done:
	/* Send the return code */
	safe_write(fd, &rc, sizeof(int));
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_signal_container(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int errnum = 0;
	int signal;

	debug("_handle_signal_container for job %u.%u",
	      job->jobid, job->stepid);

	safe_read(fd, &signal, sizeof(int));

	debug3("  uid = %d", uid);
	if (uid != job->uid && !_slurm_authorized_user(uid)) {
		debug("kill container req from uid %ld for job %u.%u "
		      "owned by uid %ld",
		      (long)uid, job->jobid, job->stepid, (long)job->uid);
		rc = -1;
		errnum = EPERM;
		goto done;
	}

	/*
	 * Sanity checks
	 */
	if (job->cont_id == 0) {
		debug ("step %u.%u invalid container [cont_id:%u]", 
			job->jobid, job->stepid, job->cont_id);
		rc = -1;
		errnum = ESLURMD_JOB_NOTRUNNING;
		goto done;
	}

	/*
	 * Signal the container
	 */
	pthread_mutex_lock(&suspend_mutex);
	if (suspended) {
		rc = -1;
		errnum = ESLURMD_STEP_SUSPENDED;
		pthread_mutex_unlock(&suspend_mutex);
		goto done;
	}

	if (slurm_container_signal(job->cont_id, signal) < 0) {
		rc = -1;
		errnum = errno;
		verbose("Error sending signal %d to %u.%u: %s", 
			signal, job->jobid, job->stepid, 
			slurm_strerror(rc));
	} else {
		verbose("Sent signal %d to %u.%u", 
			signal, job->jobid, job->stepid);
	}
	pthread_mutex_unlock(&suspend_mutex);

done:
	/* Send the return code and errnum */
	safe_write(fd, &rc, sizeof(int));
	safe_write(fd, &errnum, sizeof(int));
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_terminate(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int errnum = 0;

	debug("_handle_terminate for job %u.%u",
	      job->jobid, job->stepid);

	debug3("  uid = %d", uid);
	if (uid != job->uid && !_slurm_authorized_user(uid)) {
		debug("terminate req from uid %ld for job %u.%u "
		      "owned by uid %ld",
		      (long)uid, job->jobid, job->stepid, (long)job->uid);
		rc = -1;
		errnum = EPERM;
		goto done;
	}

	/*
	 * Sanity checks
	 */
	if (job->cont_id == 0) {
		debug ("step %u.%u invalid container [cont_id:%u]", 
			job->jobid, job->stepid, job->cont_id);
		rc = -1;
		errnum = ESLURMD_JOB_NOTRUNNING;
		goto done;
	}

	/*
	 * Signal the container with SIGKILL
	 */
	pthread_mutex_lock(&suspend_mutex);
	if (suspended) {
		debug("Terminating suspended job step %u.%u",
		      job->jobid, job->stepid);
	}

	if (slurm_container_signal(job->cont_id, SIGKILL) < 0) {
		rc = -1;
		errnum = errno;
		verbose("Error sending signal %d to %u.%u: %s", 
			signal, job->jobid, job->stepid, 
			slurm_strerror(rc));
	} else {
		verbose("Sent signal %d to %u.%u", 
			signal, job->jobid, job->stepid);
	}
	pthread_mutex_unlock(&suspend_mutex);

done:
	/* Send the return code and errnum */
	safe_write(fd, &rc, sizeof(int));
	safe_write(fd, &errnum, sizeof(int));
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_attach(int fd, slurmd_job_t *job, uid_t uid)
{
	srun_info_t *srun;
	int rc = SLURM_SUCCESS;
	int sig_len;

	debug("_handle_attach for job %u.%u", job->jobid, job->stepid);

	srun       = xmalloc(sizeof(srun_info_t));
	srun->key  = (srun_key_t *)xmalloc(SLURM_CRED_SIGLEN);

	debug("sizeof(srun_info_t) = %d, sizeof(slurm_addr) = %d",
	      sizeof(srun_info_t), sizeof(slurm_addr));
	safe_read(fd, &srun->ioaddr, sizeof(slurm_addr));
	safe_read(fd, &srun->resp_addr, sizeof(slurm_addr));
	safe_read(fd, srun->key, SLURM_CRED_SIGLEN);

	/*
	 * Check if jobstep is actually running.
	 */
	if (job->state != SLURMSTEPD_STEP_RUNNING) {
		rc = ESLURMD_JOB_NOTRUNNING;
		goto done;
	}

	/*
	 * At the moment, it only makes sense for the slurmd to make this
	 * call, so only _slurm_authorized_user is allowed.
	 */
	if (!_slurm_authorized_user(uid)) {
		error("uid %ld attempt to attach to job %u.%u owned by %ld",
		      (long) uid, job->jobid, job->stepid, (long)job->uid);
		rc = EPERM;
		goto done;
	}

	list_prepend(job->sruns, (void *) srun);
	rc = io_client_connect(srun, job);
	debug("  back from io_client_connect, rc = %d", rc);
done:
	/* Send the return code */
	safe_write(fd, &rc, sizeof(int));

	debug("  in _handle_attach rc = %d", rc);
	if (rc == SLURM_SUCCESS) {
		/* Send response info */
		uint32_t *pids, *gtids;
		int len, i;

		debug("  in _handle_attach sending response info");
		len = job->ntasks * sizeof(uint32_t);
		pids = xmalloc(len);
		gtids = xmalloc(len);
		
		if (job->task != NULL) {
			for (i = 0; i < job->ntasks; i++) {
				if (job->task[i] == NULL)
					continue;
				pids[i] = (uint32_t)job->task[i]->pid;
				gtids[i] = job->task[i]->gtid;
			}
		}

		safe_write(fd, &job->ntasks, sizeof(uint32_t));
		safe_write(fd, pids, len);
		safe_write(fd, gtids, len);
		xfree(pids);
		xfree(gtids);

		len = strlen(job->argv[0]) + 1; /* +1 to include the \0 */
		safe_write(fd, &len, sizeof(int));
		safe_write(fd, job->argv[0], len);
	}

	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_pid_in_container(int fd, slurmd_job_t *job)
{
	bool rc = false;
	pid_t pid;

	debug("_handle_pid_in_container for job %u.%u",
	      job->jobid, job->stepid);

	safe_read(fd, &pid, sizeof(pid_t));
	
	/*
	 * FIXME - we should add a new call in the proctrack API
	 *         that simply returns "true" if a pid is in the step
	 */
	if (job->cont_id == slurm_container_find(pid))
		rc = true;

	/* Send the return code */
	safe_write(fd, &rc, sizeof(bool));

	debug("Leaving _handle_pid_in_container");
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_daemon_pid(int fd, slurmd_job_t *job)
{
	safe_write(fd, &job->jmgr_pid, sizeof(pid_t));

	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_suspend(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int errnum = 0;

	debug("_handle_suspend for job %u.%u",
	      job->jobid, job->stepid);

	debug3("  uid = %d", uid);
	if (!_slurm_authorized_user(uid)) {
		debug("job step suspend request from uid %ld for job %u.%u ",
		      (long)uid, job->jobid, job->stepid);
		rc = -1;
		errnum = EPERM;
		goto done;
	}

	/*
	 * Signal the container
	 */
	pthread_mutex_lock(&suspend_mutex);
	if (suspended) {
		rc = -1;
		errnum = ESLURMD_STEP_SUSPENDED;
		pthread_mutex_unlock(&suspend_mutex);
		goto done;
	} else {
		if (slurm_container_signal(job->cont_id, SIGSTOP) < 0) {
			verbose("Error suspending %u.%u: %m", 
			        job->jobid, job->stepid);
		} else {
			verbose("Suspended %u.%u", job->jobid, job->stepid);
		}
		suspended = true;
	}
	pthread_mutex_unlock(&suspend_mutex);

done:
	/* Send the return code and errno */
	safe_write(fd, &rc, sizeof(int));
	safe_write(fd, &errnum, sizeof(int));
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_resume(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int errnum = 0;

	debug("_handle_resume for job %u.%u",
	      job->jobid, job->stepid);

	debug3("  uid = %d", uid);
	if (!_slurm_authorized_user(uid)) {
		debug("job step resume request from uid %ld for job %u.%u ",
		      (long)uid, job->jobid, job->stepid);
		rc = -1;
		errnum = EPERM;
		goto done;
	}

	/*
	 * Signal the container
	 */
	pthread_mutex_lock(&suspend_mutex);
	if (!suspended) {
		rc = -1;
		errnum = ESLURMD_STEP_NOTSUSPENDED;
		pthread_mutex_unlock(&suspend_mutex);
		goto done;
	} else {
		if (slurm_container_signal(job->cont_id, SIGCONT) < 0) {
			verbose("Error resuming %u.%u: %m", 
			        job->jobid, job->stepid);
		} else {
			verbose("Resumed %u.%u", job->jobid, job->stepid);
		}
		suspended = false;
	}
	pthread_mutex_unlock(&suspend_mutex);

done:
	/* Send the return code and errno */
	safe_write(fd, &rc, sizeof(int));
	safe_write(fd, &errnum, sizeof(int));
	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}

static int
_handle_completion(int fd, slurmd_job_t *job, uid_t uid)
{
	int rc = SLURM_SUCCESS;
	int errnum = 0;
	int first;
	int last;

	debug("_handle_completion for job %u.%u",
	      job->jobid, job->stepid);

	debug3("  uid = %d", uid);
	if (!_slurm_authorized_user(uid)) {
		debug("job step completion message from uid %ld for job %u.%u ",
		      (long)uid, job->jobid, job->stepid);
		rc = -1;
		errnum = EPERM;
		/* Send the return code and errno */
		safe_write(fd, &rc, sizeof(int));
		safe_write(fd, &errnum, sizeof(int));
		return SLURM_SUCCESS;
	}

	safe_read(fd, &first, sizeof(int));
	safe_read(fd, &last, sizeof(int));

	/*
	 * Record the completed nodes
	 */
	pthread_mutex_lock(&step_complete.lock);
	bit_nset(step_complete.bits,
		 first - (step_complete.rank+1),
		 last - (step_complete.rank+1));
	/* Send the return code and errno, we do this within the locked
	 * region to ensure that the stepd doesn't exit before we can
	 * perform this send. */
	safe_write(fd, &rc, sizeof(int));
	safe_write(fd, &errnum, sizeof(int));
	pthread_cond_signal(&step_complete.cond);
	pthread_mutex_unlock(&step_complete.lock);

	return SLURM_SUCCESS;
rwfail:
	return SLURM_FAILURE;
}
