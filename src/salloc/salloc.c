/*****************************************************************************\
 *  salloc.c - Request a SLURM job allocation and
 *             launch a user-specified command.
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>

#include <slurm/slurm.h>

#include "src/common/xstring.h"
#include "src/common/xmalloc.h"
#include "src/common/xsignal.h"
#include "src/common/read_config.h"
#include "src/common/env.h"

#include "src/salloc/salloc.h"
#include "src/salloc/opt.h"

#define MAX_RETRIES 3

char **command_argv;
int command_argc;
pid_t command_pid = -1;
enum possible_allocation_states allocation_state = NOT_GRANTED;
pthread_mutex_t allocation_state_lock = PTHREAD_MUTEX_INITIALIZER;

static bool exit_flag = false;
static bool allocation_interrupted = false;
static uint32_t pending_job_id = 0;

static int fill_job_desc_from_opts(job_desc_msg_t *desc);
static void ring_terminal_bell(void);
static int fork_command(char **command);
static void _pending_callback(uint32_t job_id);
static void _ignore_signal(int signo);
static void _exit_on_signal(int signo);
static void _signal_while_allocating(int signo);
static void _job_complete_handler(srun_job_complete_msg_t *msg);
static void _timeout_handler(srun_timeout_msg_t *msg);
static void _user_msg_handler(srun_user_msg_t *msg);
static void _ping_handler(srun_ping_msg_t *msg);
static void _node_fail_handler(srun_node_fail_msg_t *msg);

int main(int argc, char *argv[])
{
	log_options_t logopt = LOG_OPTS_STDERR_ONLY;
	job_desc_msg_t desc;
	resource_allocation_response_msg_t *alloc;
	time_t before, after;
	allocation_msg_thread_t *msg_thr;
	char **env = NULL;
	int status = 0;
	int errnum = 0;
	int retries = 0;
	pid_t pid;
	pid_t rc_pid = 0;
	int rc = 0;
	static char *msg = "Slurm job queue full, sleeping and retrying.";
	slurm_allocation_callbacks_t callbacks;

	log_init(xbasename(argv[0]), logopt, 0, NULL);
	if (initialize_and_process_args(argc, argv) < 0) {
		fatal("salloc parameter parsing");
	}
	/* reinit log with new verbosity (if changed by command line) */
	if (opt.verbose || opt.quiet) {
		logopt.stderr_level += opt.verbose;
		logopt.stderr_level -= opt.quiet;
		logopt.prefix_level = 1;
		log_alter(logopt, 0, NULL);
	}

	/*
	 * Request a job allocation
	 */
	slurm_init_job_desc_msg(&desc);
	if (fill_job_desc_from_opts(&desc) == -1) {
		exit(1);
	}

	callbacks.ping = _ping_handler;
	callbacks.timeout = _timeout_handler;
	callbacks.job_complete = _job_complete_handler;
	callbacks.user_msg = _user_msg_handler;
	callbacks.node_fail = _node_fail_handler;
	/* create message thread to handle pings and such from slurmctld */
	msg_thr = slurm_allocation_msg_thr_create(&desc.other_port, &callbacks);

	xsignal(SIGHUP, _signal_while_allocating);
	xsignal(SIGINT, _signal_while_allocating);
	xsignal(SIGQUIT, _signal_while_allocating);
	xsignal(SIGPIPE, _signal_while_allocating);
	xsignal(SIGTERM, _signal_while_allocating);
	xsignal(SIGUSR1, _signal_while_allocating);
	xsignal(SIGUSR2, _signal_while_allocating);

	before = time(NULL);
	while ((alloc = slurm_allocate_resources_blocking(&desc, opt.max_wait,
					_pending_callback)) == NULL) {
		if ((errno != ESLURM_ERROR_ON_DESC_TO_RECORD_COPY) ||
		    (retries >= MAX_RETRIES))
			break;
		if (retries == 0)
			error(msg);
		else
			debug(msg);
		sleep (++retries);
	}

	if (alloc == NULL) {
		if (allocation_interrupted) {
			/* cancelled by signal */
		} else if (errno == EINTR) {
			error("Interrupted by signal."
			      "  Allocation request rescinded.");
		} else {
			error("Failed to allocate resources: %m");
		}
		slurm_allocation_msg_thr_destroy(msg_thr);
		exit(1);
	}
	after = time(NULL);

	xsignal(SIGHUP, _exit_on_signal);
	xsignal(SIGINT, _ignore_signal);
	xsignal(SIGQUIT, _ignore_signal);
	xsignal(SIGPIPE, _ignore_signal);
	xsignal(SIGTERM, _ignore_signal);
	xsignal(SIGUSR1, _ignore_signal);
	xsignal(SIGUSR2, _ignore_signal);

	/*
	 * Allocation granted!
	 */
	info("Granted job allocation %d", alloc->job_id);
	if (opt.bell == BELL_ALWAYS
	    || (opt.bell == BELL_AFTER_DELAY
		&& ((after - before) > DEFAULT_BELL_DELAY))) {
		ring_terminal_bell();
	}
	if (allocation_interrupted) {
		/* salloc process received a signal after
		 * slurm_allocate_resources_blocking returned with the
		 * allocation, but before the new signal handlers were
		 * registered.
		 */
		goto relinquish;
	}

	/*
	 * Run the user's command.
	 */
	env_array_for_job(&env, alloc);
	/* Add default task count for srun, if not already set */
	if (opt.nprocs_set)
		env_array_append_fmt(&env, "SLURM_NPROCS", "%d", opt.nprocs);
	if (opt.overcommit) {
		env_array_append_fmt(&env, "SLURM_OVERCOMMIT", "%d", 
			opt.overcommit);
	}
	if (opt.acctg_freq >= 0) {
		env_array_append_fmt(&env, "SLURM_ACCTG_FREQ", "%d",
			opt.acctg_freq);
	}
	if (opt.task_mem >= 0) {
		env_array_append_fmt(&env, "SLURM_TASK_MEM", "%d",
			opt.task_mem);
	}
	env_array_set_environment(env);
	env_array_free(env);
	pthread_mutex_lock(&allocation_state_lock);
	if (allocation_state == REVOKED) {
		error("Allocation was revoked before command could be run");
		pthread_mutex_unlock(&allocation_state_lock);
		return 1;
	} else {
		allocation_state = GRANTED;
		command_pid = pid = fork_command(command_argv);
	}
	pthread_mutex_unlock(&allocation_state_lock);

	/*
	 * Wait for command to exit, OR for waitpid to be interrupted by a
	 * signal.  Either way, we are going to release the allocation next.
	 */
	if (pid > 0) {
		while ((rc_pid = waitpid(pid, &status, 0)) == -1) {
			if (exit_flag)
				break;
			if (errno == EINTR)
				continue;
		}
		errnum = errno;
		if (rc_pid == -1 && errnum != EINTR)
			error("waitpid for %s failed: %m", command_argv[0]);
	}

	/*
	 * Relinquish the job allocation (if not already revoked).
	 */
relinquish:
	pthread_mutex_lock(&allocation_state_lock);
	if (allocation_state != REVOKED) {
		info("Relinquishing job allocation %d", alloc->job_id);
		if (slurm_complete_job(alloc->job_id, status)
		    != 0)
			error("Unable to clean up job allocation %d: %m",
			      alloc->job_id);
		else
			allocation_state = REVOKED;
	}
	pthread_mutex_unlock(&allocation_state_lock);

	slurm_free_resource_allocation_response_msg(alloc);
	slurm_allocation_msg_thr_destroy(msg_thr);

	/*
	 * Figure out what return code we should use.  If the user's command
	 * exitted normally, return the user's return code.
	 */
	rc = 1;
	if (rc_pid != -1) {
		
		if (WIFEXITED(status)) {
			rc = WEXITSTATUS(status);
		} else if (WIFSIGNALED(status)) {
			verbose("Command \"%s\" was terminated by signal %d",
				command_argv[0], WTERMSIG(status));
		}
	}

	return rc;
}


/* Returns 0 on success, -1 on failure */
static int fill_job_desc_from_opts(job_desc_msg_t *desc)
{
	desc->contiguous = opt.contiguous ? 1 : 0;
	desc->features = opt.constraints;
	desc->immediate = opt.immediate ? 1 : 0;
	desc->name = opt.job_name;
	desc->req_nodes = opt.nodelist;
	desc->exc_nodes = opt.exc_nodes;
	desc->partition = opt.partition;
	desc->min_nodes = opt.min_nodes;
	if (opt.max_nodes)
		desc->max_nodes = opt.max_nodes;
	desc->user_id = opt.uid;
	desc->group_id = opt.gid;
	if (opt.dependency)
		desc->dependency = xstrdup(opt.dependency);
	desc->task_dist  = opt.distribution;
	if (opt.plane_size != NO_VAL)
		desc->plane_size = opt.plane_size;
	if (opt.nice)
		desc->nice = NICE_OFFSET + opt.nice;
	desc->mail_type = opt.mail_type;
	if (opt.mail_user)
		desc->mail_user = xstrdup(opt.mail_user);
	if (opt.begin)
		desc->begin_time = opt.begin;
	if (opt.account)
		desc->account = xstrdup(opt.account);
	if (opt.comment)
		desc->comment = xstrdup(opt.comment);

	if (opt.hold)
		desc->priority     = 0;
#if SYSTEM_DIMENSIONS
	if (opt.geometry[0] > 0) {
		int i;
		for (i=0; i<SYSTEM_DIMENSIONS; i++)
			desc->geometry[i] = opt.geometry[i];
	}
#endif
	if (opt.conn_type != -1)
		desc->conn_type = opt.conn_type;
	if (opt.reboot)
		desc->reboot = 1;
	if (opt.no_rotate)
		desc->rotate = 0;
	if (opt.blrtsimage)
		desc->blrtsimage = xstrdup(opt.blrtsimage);
	if (opt.linuximage)
		desc->linuximage = xstrdup(opt.linuximage);
	if (opt.mloaderimage)
		desc->mloaderimage = xstrdup(opt.mloaderimage);
	if (opt.ramdiskimage)
		desc->ramdiskimage = xstrdup(opt.ramdiskimage);

	/* job constraints */
	if (opt.mincpus > -1)
		desc->job_min_procs = opt.mincpus;
	if (opt.minsockets > -1)
		desc->job_min_sockets = opt.minsockets;
	if (opt.mincores > -1)
		desc->job_min_cores = opt.mincores;
	if (opt.minthreads > -1)
		desc->job_min_threads = opt.minthreads;
	if (opt.realmem > -1)
		desc->job_min_memory = opt.realmem;
	if (opt.tmpdisk > -1)
		desc->job_min_tmp_disk = opt.tmpdisk;
	if (opt.overcommit) {
		desc->num_procs = opt.min_nodes;
		desc->overcommit = opt.overcommit;
	} else
		desc->num_procs = opt.nprocs * opt.cpus_per_task;
	if (opt.nprocs_set)
		desc->num_tasks = opt.nprocs;
	if (opt.cpus_set)
		desc->cpus_per_task = opt.cpus_per_task;
	if (opt.ntasks_per_node > -1)
		desc->ntasks_per_node = opt.ntasks_per_node;
	if (opt.ntasks_per_socket > -1)
		desc->ntasks_per_socket = opt.ntasks_per_socket;
	if (opt.ntasks_per_core > -1)
		desc->ntasks_per_core = opt.ntasks_per_core;

	/* node constraints */
	if (opt.min_sockets_per_node > -1)
		desc->min_sockets = opt.min_sockets_per_node;
	if (opt.max_sockets_per_node > -1)
		desc->max_sockets = opt.max_sockets_per_node;
	if (opt.min_cores_per_socket > -1)
		desc->min_cores = opt.min_cores_per_socket;
	if (opt.max_cores_per_socket > -1)
		desc->max_cores = opt.max_cores_per_socket;
	if (opt.min_threads_per_core > -1)
		desc->min_threads = opt.min_threads_per_core;
	if (opt.max_threads_per_core > -1)
		desc->max_threads = opt.max_threads_per_core;

	if (opt.no_kill)
		desc->kill_on_node_fail = 0;
	if (opt.time_limit  != NO_VAL)
		desc->time_limit = opt.time_limit;
	desc->shared = opt.shared;
	desc->job_id = opt.jobid;

	return 0;
}

static void ring_terminal_bell(void)
{
        if (isatty(STDOUT_FILENO)) {
                fprintf(stdout, "\a");
                fflush(stdout);
        }
}

/* returns the pid of the forked command, or <0 on error */
static pid_t fork_command(char **command)
{
	pid_t pid;

	pid = fork();
	if (pid < 0) {
		error("fork failed: %m");
	} else if (pid == 0) {
		/* child */
		execvp(command[0], command);

		/* should only get here if execvp failed */
		error("Unable to exec command \"%s\"", command[0]);
		exit(1);
	}
	/* parent returns */
	return pid;
}

static void _pending_callback(uint32_t job_id)
{
	info("Pending job allocation %u", job_id);
	pending_job_id = job_id;
}

static void _signal_while_allocating(int signo)
{
	allocation_interrupted = true;
	if (pending_job_id != 0) {
		slurm_complete_job(pending_job_id, 0);
	}
}

static void _ignore_signal(int signo)
{
	/* do nothing */
}

static void _exit_on_signal(int signo)
{
	exit_flag = true;
}

/* This typically signifies the job was cancelled by scancel */
static void _job_complete_handler(srun_job_complete_msg_t *comp)
{
	if (comp->step_id == NO_VAL) {
		pthread_mutex_lock(&allocation_state_lock);
		if (allocation_state != REVOKED) {
			/* If the allocation_state is already REVOKED, then
			 * no need to print this message.  We probably
			 * relinquished the allocation ourself.
			 */
			info("Job allocation %u has been revoked.",
			     comp->job_id);
		}
		if (allocation_state == GRANTED
		    && command_pid > -1
		    && opt.kill_command_signal_set) {
			verbose("Sending signal %d to command \"%s\", pid %d",
				opt.kill_command_signal,
				command_argv[0], command_pid);
			kill(command_pid, opt.kill_command_signal);
		}
		allocation_state = REVOKED;
		pthread_mutex_unlock(&allocation_state_lock);
	} else {
		verbose("Job step %u.%u is finished.",
			comp->job_id, comp->step_id);
	}
}

/*
 * Job has been notified of it's approaching time limit. 
 * Job will be killed shortly after timeout.
 * This RPC can arrive multiple times with the same or updated timeouts.
 * FIXME: We may want to signal the job or perform other action for this.
 * FIXME: How much lead time do we want for this message? Some jobs may 
 *	require tens of minutes to gracefully terminate.
 */
static void _timeout_handler(srun_timeout_msg_t *msg)
{
	static time_t last_timeout = 0;

	if (msg->timeout != last_timeout) {
		last_timeout = msg->timeout;
		verbose("Job allocation time limit to be reached at %s", 
			ctime(&msg->timeout));
	}
}

static void _user_msg_handler(srun_user_msg_t *msg)
{
	info("%s", msg->msg);
}

static void _ping_handler(srun_ping_msg_t *msg) 
{
	/* the api will respond so there really isn't anything to do
	   here */
}

static void _node_fail_handler(srun_node_fail_msg_t *msg)
{
	error("Node failure on %s", msg->nodelist);
}
