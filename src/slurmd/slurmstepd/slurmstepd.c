/*****************************************************************************\
 *  src/slurmd/slurmstepd/slurmstepd.c - SLURM job-step manager.
 *  $Id: $
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <da@llnl.gov>.
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

#include <unistd.h>
#include <stdlib.h>
#include <signal.h>

#include "src/common/xmalloc.h"
#include "src/common/xsignal.h"

#include "src/slurmd/slurmd/slurmd.h"
#include "src/slurmd/common/slurmstepd_init.h"
#include "src/slurmd/common/stepd_api.h"
#include "src/slurmd/slurmstepd/slurmstepd.h"
#include "src/slurmd/slurmstepd/mgr.h"
#include "src/slurmd/slurmstepd/slurmstepd_job.h"

static int _init_from_slurmd(int sock, char **argv, slurm_addr **_cli,
			    slurm_addr **_self, slurm_msg_t **_msg);
static void _send_ok_to_slurmd(int sock);
static void _send_fail_to_slurmd(int sock);
static slurmd_job_t *_step_setup(slurm_addr *cli, slurm_addr *self,
				 slurm_msg_t *msg);
static void _step_cleanup(slurmd_job_t *job, slurm_msg_t *msg, int rc);

int slurmstepd_blocked_signals[] = {
	SIGPIPE, 0
};

int 
main (int argc, char *argv[])
{
	slurm_addr *cli;
	slurm_addr *self;
	slurm_msg_t *msg;
	slurmd_job_t *job;
	int rc;

	xsignal_block(slurmstepd_blocked_signals);
	conf = xmalloc(sizeof(*conf));
	conf->argv = &argv;
	conf->argc = &argc;
	conf->task_prolog = slurm_get_task_prolog();
	conf->task_epilog = slurm_get_task_epilog();
	init_setproctitle(argc, argv);
	if (slurm_proctrack_init() != SLURM_SUCCESS)
		return SLURM_FAILURE;

	_init_from_slurmd(STDIN_FILENO, argv, &cli, &self, &msg);
	close(STDIN_FILENO);

	job = _step_setup(cli, self, msg);

	/* sets job->msg_handle and job->msgid */
	if (msg_thr_create(job) == SLURM_ERROR) {
		_send_fail_to_slurmd(STDOUT_FILENO);
		return -1;
	}

	_send_ok_to_slurmd(STDOUT_FILENO);
	close(STDOUT_FILENO);

	rc = job_manager(job); /* blocks until step is complete */

	/* signal the message thread to shutdown, and wait for it */
	eio_signal_shutdown(job->msg_handle);
	pthread_join(job->msgid, NULL);

	_step_cleanup(job, msg, rc);

	xfree(cli);
	xfree(self);
	xfree(conf->hostname);
	xfree(conf->spooldir);
	xfree(conf->node_name);
	xfree(conf->logfile);
	xfree(conf->cf.job_acct_parameters);
	xfree(conf);
	info("done with job");
	return 0;
}

static void
_send_ok_to_slurmd(int sock)
{
	int ok = SLURM_SUCCESS;
	safe_write(sock, &ok, sizeof(int));
	return;
rwfail:
	error("Unable to send \"ok\" to slurmd");
}

static void
_send_fail_to_slurmd(int sock)
{
	int fail = SLURM_FAILURE;

	if (errno)
		fail = errno;
	safe_write(sock, &fail, sizeof(int));
	return;
rwfail:
	error("Unable to send \"fail\" to slurmd");
}

static int
_init_from_slurmd(int sock, char **argv,
		  slurm_addr **_cli, slurm_addr **_self, slurm_msg_t **_msg)
{
	char *incoming_buffer = NULL;
	Buf buffer;
	int step_type;
	int len;
	int rc;	
	char c;
	slurm_addr *cli = NULL;
	slurm_addr *self = NULL;
	slurm_msg_t *msg = NULL;

	/* receive job type from slurmd */
	safe_read(sock, &step_type, sizeof(int));
	debug3("step_type = %d", step_type);
	
	/* receive conf from slurmd */
	safe_read(sock, &len, sizeof(int));
	incoming_buffer = xmalloc(len);
	safe_read(sock, incoming_buffer, len);
	buffer = create_buf(incoming_buffer,len);
	if(unpack_slurmd_conf_lite_no_alloc(conf, buffer) == SLURM_ERROR) {
		fatal("slurmstepd: problem with unpack of slurmd_conf");
	}
	free_buf(buffer);
				
	debug2("debug level is %d.", conf->debug_level);
	conf->log_opts.stderr_level = conf->debug_level;
	conf->log_opts.logfile_level = conf->debug_level;
	conf->log_opts.syslog_level = conf->debug_level;
	/* forward the log options to slurmstepd */
	//log_alter(conf->log_opts, 0, NULL);
	/*
	 * If daemonizing, turn off stderr logging -- also, if
	 * logging to a file, turn off syslog.
	 *
	 * Otherwise, if remaining in foreground, turn off logging
	 * to syslog (but keep logfile level)
	 */
	if (conf->daemonize) {
		conf->log_opts.stderr_level = LOG_LEVEL_QUIET;
		if (conf->logfile)
			conf->log_opts.syslog_level = LOG_LEVEL_QUIET;
	} else 
		conf->log_opts.syslog_level  = LOG_LEVEL_QUIET;

	log_init(argv[0],conf->log_opts, LOG_DAEMON, conf->logfile);
	g_slurmd_jobacct_init(conf->cf.job_acct_parameters);
	switch_g_slurmd_step_init();

	/* receive cli from slurmd */
	safe_read(sock, &len, sizeof(int));
	incoming_buffer = xmalloc(sizeof(char) * len);
	safe_read(sock, incoming_buffer, len);
	buffer = create_buf(incoming_buffer,len);	
	cli = xmalloc(sizeof(slurm_addr));
	if(slurm_unpack_slurm_addr_no_alloc(cli, buffer) == SLURM_ERROR) {
		fatal("slurmstepd: problem with unpack of slurmd_conf");
	}
	free_buf(buffer);

	/* receive self from slurmd */
	safe_read(sock, &len, sizeof(int));
	if(len > 0) {
		/* receive packed self from main slurmd */
		incoming_buffer = xmalloc(sizeof(char) * len);
		safe_read(sock, incoming_buffer, len);
		buffer = create_buf(incoming_buffer,len);
		self = xmalloc(sizeof(slurm_addr));
		if(slurm_unpack_slurm_addr_no_alloc(self, buffer)
		   == SLURM_ERROR) {
			fatal("slurmstepd: problem with unpack of "
			      "slurmd_conf");
		}
		free_buf(buffer);
	}
		
	/* receive req from slurmd */
	safe_read(sock, &len, sizeof(int));
	incoming_buffer = xmalloc(sizeof(char) * len);
	safe_read(sock, incoming_buffer, len);
	buffer = create_buf(incoming_buffer,len);

	msg = xmalloc(sizeof(slurm_msg_t));
	switch(step_type) {
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
		fatal("Unrecognized launch/spawn RPC");
		break;
	}
	if(unpack_msg(msg, buffer) == SLURM_ERROR) 
		fatal("slurmstepd: we didn't unpack the request correctly");
	free_buf(buffer);

	*_cli = cli;
	*_self = self;
	*_msg = msg;

	return 1;

rwfail:
	fatal("Error reading initialization data from slurmd");
	exit(1);
}

static slurmd_job_t *
_step_setup(slurm_addr *cli, slurm_addr *self, slurm_msg_t *msg)
{
	slurmd_job_t *job;

	switch(msg->msg_type) {
	case REQUEST_BATCH_JOB_LAUNCH:
		debug2("setup for a batch_job");
		job = mgr_launch_batch_job_setup(msg->data, cli);
		break;
	case REQUEST_LAUNCH_TASKS:
		debug2("setup for a launch_task");
		job = mgr_launch_tasks_setup(msg->data, cli, self);
		break;
	case REQUEST_SPAWN_TASK:
		debug2("setup for a spawn_task");
		job = mgr_spawn_task_setup(msg->data, cli, self);
		break;
	default:
		fatal("handle_launch_message: Unrecognized launch/spawn RPC");
		break;
	}
	job->jmgr_pid = getpid();

	return job;
}

static void
_step_cleanup(slurmd_job_t *job, slurm_msg_t *msg, int rc)
{
	if (job->batch)
		mgr_launch_batch_job_cleanup(job, rc);
	else
		job_destroy(job);

	/* 
	 * The message cannot be freed until the jobstep is complete
	 * because the job struct has pointers into the msg, such
	 * as the switch jobinfo pointer.
	 */
	switch(msg->msg_type) {
	case REQUEST_BATCH_JOB_LAUNCH:
		slurm_free_job_launch_msg(msg->data);
		break;
	case REQUEST_LAUNCH_TASKS:
		slurm_free_launch_tasks_request_msg(msg->data);
		break;
	case REQUEST_SPAWN_TASK:
		slurm_free_spawn_task_request_msg(msg->data);
		break;
	default:
		fatal("handle_launch_message: Unrecognized launch/spawn RPC");
		break;
	}
	xfree(msg);
}
