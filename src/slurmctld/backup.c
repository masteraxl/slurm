/*****************************************************************************\
 *  backup.c - backup slurm controller
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette@llnl.gov>, Kevin Tew <tew1@llnl.gov>, et. al.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif				/* WITH_PTHREADS */

#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>

#include <slurm/slurm_errno.h>

#include "src/common/daemonize.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/slurm_auth.h"
#include "src/common/switch.h"
#include "src/common/xsignal.h"
#include "src/common/xstring.h"

#include "src/slurmctld/locks.h"
#include "src/slurmctld/read_config.h"
#include "src/slurmctld/slurmctld.h"

static int          _background_process_msg(slurm_msg_t * msg);
static void *       _background_rpc_mgr(void *no_data);
static void *       _background_signal_hand(void *no_data);
static int          _ping_controller(void);
inline static void  _update_cred_key(void);

/* Local variables */
static bool     dump_core = false;

/*
 * Static list of signals to block in this process
 * *Must be zero-terminated*
 */
static int backup_sigarray[] = {
	SIGINT,  SIGTERM, SIGCHLD, SIGUSR1,
	SIGUSR2, SIGTSTP, SIGXCPU, SIGQUIT,
	SIGPIPE, SIGALRM, SIGABRT, SIGHUP, 0
};

/* run_backup - this is the backup controller, it should run in standby 
 *	mode, assuming control when the primary controller stops responding */
void run_backup(void)
{
	time_t last_controller_response = time(NULL), last_ping = 0;
	pthread_attr_t thread_attr_sig, thread_attr_rpc;
	slurmctld_lock_t config_read_lock = { 
		READ_LOCK, NO_LOCK, NO_LOCK, NO_LOCK };

	info("slurmctld running in background mode");
	/* default: don't resume if shutdown */
	slurmctld_config.resume_backup = false;
	if (xsignal_block(backup_sigarray) < 0)
		error("Unable to block signals");

	/*
	 * create attached thread to process RPCs
	 */
	slurm_attr_init(&thread_attr_rpc);
	if (pthread_create(&slurmctld_config.thread_id_rpc, 
			&thread_attr_rpc, _background_rpc_mgr, NULL))
		fatal("pthread_create error %m");

	/*
	 * create attached thread for signal handling
	 */
	slurm_attr_init(&thread_attr_sig);
	if (pthread_create(&slurmctld_config.thread_id_sig,
			&thread_attr_sig, _background_signal_hand, NULL))
		fatal("pthread_create %m");

	sleep(5);       /* Give the primary slurmctld set-up time */
	/* repeatedly ping ControlMachine */
	while (slurmctld_config.shutdown_time == 0) {
		sleep(1);
		/* Lock of slurmctld_conf below not important */
		if (difftime(time(NULL), last_ping) <
		    slurmctld_conf.heartbeat_interval)
			continue;

		last_ping = time(NULL);
		if (_ping_controller() == 0)
			last_controller_response = time(NULL);
		else {
			uint32_t timeout;
			lock_slurmctld(config_read_lock);
			timeout = slurmctld_conf.slurmctld_timeout;
			unlock_slurmctld(config_read_lock);

			if (difftime(time(NULL), last_controller_response) >
					timeout)
				break;
		}
	}

	/* Since pidfile is created as user root (its owner is
	 *   changed to SlurmUser) SlurmUser may not be able to 
	 *   remove it, so this is not necessarily an error. 
	 * No longer need slurmctld_conf lock after above join. */
	if (unlink(slurmctld_conf.slurmctld_pidfile) < 0)
		verbose("Unable to remove pidfile '%s': %m",
			slurmctld_conf.slurmctld_pidfile);

	if (slurmctld_config.shutdown_time != 0) {
		info("BackupController terminating");
		pthread_join(slurmctld_config.thread_id_sig, NULL);
		log_fini();
		if (dump_core)
			abort();
		else
			exit(0);
	}

	lock_slurmctld(config_read_lock);
	error("ControlMachine %s not responding, "
		"BackupController %s taking over",
		slurmctld_conf.control_machine,
		slurmctld_conf.backup_controller);
	unlock_slurmctld(config_read_lock);

	pthread_kill(slurmctld_config.thread_id_sig, SIGTERM);
	pthread_join(slurmctld_config.thread_id_sig, NULL);
	pthread_join(slurmctld_config.thread_id_rpc, NULL);

	/* clear old state and read new state */
	job_fini();
	switch_clear();
	if (read_slurm_conf(2)) {	/* Recover all state */
		error("Unable to recover slurm state");
		abort();
	}
	slurmctld_config.shutdown_time = (time_t) 0;
	return;
}

/* _background_signal_hand - Process daemon-wide signals for the 
 *	backup controller */
static void *_background_signal_hand(void *no_data)
{
	int sig, rc;
	sigset_t set;
	/* Locks: Write configuration, job, node, and partition */
	slurmctld_lock_t config_write_lock = { 
		WRITE_LOCK, WRITE_LOCK, WRITE_LOCK, WRITE_LOCK };

	(void) pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
	(void) pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

	/* No need for slurmctld_conf lock yet */
	while ( (create_pidfile(slurmctld_conf.slurmctld_pidfile) < 0) && 
	        (errno == EAGAIN) ) {
		verbose("Retrying create_pidfile: %m");
		sleep(1);
	}

	while (slurmctld_config.shutdown_time == 0) {
		xsignal_sigset_create(backup_sigarray, &set);
		rc = sigwait(&set, &sig);
		if (rc == EINTR)
			continue;
		switch (sig) {
		case SIGINT:	/* kill -2  or <CTRL-C> */
		case SIGTERM:	/* kill -15 */
			info("Terminate signal (SIGINT or SIGTERM) received");
			slurmctld_config.shutdown_time = time(NULL);
			slurmctld_shutdown();
			return NULL;	/* Normal termination */
			break;
		case SIGHUP:    /* kill -1 */
			info("Reconfigure signal (SIGHUP) received");
			/*
			 * XXX - need to shut down the scheduler
			 * plugin, re-read the configuration, and then
			 * restart the (possibly new) plugin.
			 */
			lock_slurmctld(config_write_lock);
			rc = read_slurm_conf(0);
			if (rc)
				error("read_slurm_conf: %s",
					slurm_strerror(rc));
			else {
				/* Leave config lock set through this */
				_update_cred_key();
			}
			unlock_slurmctld(config_write_lock);
			break;
		case SIGABRT:   /* abort */
			info("SIGABRT received");
			slurmctld_config.shutdown_time = time(NULL);
			slurmctld_shutdown();
			dump_core = true;
			return NULL;    /* Normal termination */
			break;
		default:
			error("Invalid signal (%d) received", sig);
		}
	}
	return NULL;
}

/* Reset the job credential key based upon configuration parameters.
 * slurmctld_conf is locked on entry. */
static void _update_cred_key(void)
{	
	slurm_cred_ctx_key_update(slurmctld_config.cred_ctx, 
			slurmctld_conf.job_credential_private_key);
}

static void _sig_handler(int signal)
{
}

/* _background_rpc_mgr - Read and process incoming RPCs to the background 
 *	controller (that's us) */
static void *_background_rpc_mgr(void *no_data)
{
	slurm_fd newsockfd;
	slurm_fd sockfd;
	slurm_addr cli_addr;
	slurm_msg_t *msg = NULL;
	int error_code;
	/* Read configuration only */
	slurmctld_lock_t config_read_lock = { 
		READ_LOCK, NO_LOCK, NO_LOCK, NO_LOCK };
	int sigarray[] = {SIGUSR1, 0};

	(void) pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
	(void) pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
	debug3("_background_rpc_mgr pid = %lu", (unsigned long) getpid());

	/* initialize port for RPCs */
	lock_slurmctld(config_read_lock);
	if ((sockfd =
	     slurm_init_msg_engine_port(slurmctld_conf.slurmctld_port))
	    == SLURM_SOCKET_ERROR)
		fatal("slurm_init_msg_engine_port error %m");
	unlock_slurmctld(config_read_lock);

	/* Prepare to catch SIGUSR1 to interrupt accept().
	 * This signal is generated by the slurmctld signal
	 * handler thread upon receipt of SIGABRT, SIGINT, 
	 * or SIGTERM. That thread does all processing of  
	 * all signals. */ 
	xsignal(SIGUSR1, _sig_handler);
	xsignal_unblock(sigarray);

	/*
	 * Process incoming RPCs indefinitely
	 */
	while (slurmctld_config.shutdown_time == 0) {
		/* accept needed for stream implementation 
		 * is a no-op in message implementation that just passes 
		 * sockfd to newsockfd */
		if ((newsockfd = slurm_accept_msg_conn(sockfd, &cli_addr))
		    == SLURM_SOCKET_ERROR) {
			if (errno != EINTR)
				error("slurm_accept_msg_conn: %m");
			continue;
		}

		msg = xmalloc(sizeof(slurm_msg_t));
		msg->conn_fd = newsockfd;
		if (slurm_receive_msg(newsockfd, msg, 0) < 0)
			error("slurm_receive_msg: %m");
		else {
			error_code = _background_process_msg(msg);
			if ((error_code == SLURM_SUCCESS) &&
			    (msg->msg_type == REQUEST_SHUTDOWN_IMMEDIATE) &&
			    (slurmctld_config.shutdown_time == 0))
				 slurmctld_config.shutdown_time = time(NULL);
		}
		slurm_free_msg(msg);

		/* close should only be called when the socket 
		 * implementation is being used the following call will 
		 * be a no-op in a message/mongo implementation */
		slurm_close_accepted_conn(newsockfd);	/* close new socket */
	}

	debug3("_background_rpc_mgr shutting down");
	slurm_close_accepted_conn(sockfd);	/* close the main socket */
	pthread_exit((void *) 0);
	return NULL;
}

/* _background_process_msg - process an RPC to the backup_controller */
static int _background_process_msg(slurm_msg_t * msg)
{
	int error_code = SLURM_SUCCESS;

	if (msg->msg_type != REQUEST_PING) {
		bool super_user = false;
		uid_t uid = g_slurm_auth_get_uid(msg->cred);
		if ((uid == 0) || (uid == getuid()))
			super_user = true;

		if (super_user && 
		    (msg->msg_type == REQUEST_SHUTDOWN_IMMEDIATE)) {
			info("Performing RPC: REQUEST_SHUTDOWN_IMMEDIATE");
		} else if (super_user && 
			   (msg->msg_type == REQUEST_SHUTDOWN)) {
			info("Performing RPC: REQUEST_SHUTDOWN");
			pthread_kill(slurmctld_config.thread_id_sig, SIGTERM);
		} else if (super_user && 
			   (msg->msg_type == REQUEST_CONTROL)) {
			debug3("Ignoring RPC: REQUEST_CONTROL");
		} else {
			error("Invalid RPC received %d from uid %u", 
			      msg->msg_type, uid);
			error_code = ESLURM_IN_STANDBY_MODE;
		}
	}
	if (msg->msg_type != REQUEST_SHUTDOWN_IMMEDIATE)
		slurm_send_rc_msg(msg, error_code);
	return error_code;
}

/* Ping primary ControlMachine
 * RET 0 if no error */
static int _ping_controller(void)
{
	int rc;
	slurm_msg_t req;
	/* Locks: Read configuration */
	slurmctld_lock_t config_read_lock = { 
		READ_LOCK, NO_LOCK, NO_LOCK, NO_LOCK };

	/* 
	 *  Set address of controller to ping
	 */
	lock_slurmctld(config_read_lock);
	debug3("pinging slurmctld at %s", slurmctld_conf.control_addr);
	slurm_set_addr(&req.address, slurmctld_conf.slurmctld_port, 
	               slurmctld_conf.control_addr);
	unlock_slurmctld(config_read_lock);

	req.msg_type = REQUEST_PING;

	if (slurm_send_recv_rc_msg(&req, &rc, 0) < 0) {
		error("_ping_controller/slurm_send_node_msg error: %m");
		return SLURM_ERROR;
	}

	if (rc) {
		error("_ping_controller/response error %d", rc);
		return SLURM_PROTOCOL_ERROR;
	}

	return SLURM_PROTOCOL_SUCCESS;
}

