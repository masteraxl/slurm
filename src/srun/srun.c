/*****************************************************************************\
 *  srun.c - user interface to allocate resources, submit jobs, and execute 
 *	parallel jobs.
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona@llnl.gov>, et. al.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif

#ifdef HAVE_AIX
#  undef HAVE_UNSETENV
#  include <sys/checkpnt.h>
#endif
#ifndef HAVE_UNSETENV
#  include "src/common/unsetenv.h"
#endif

#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <ctype.h>
#include <fcntl.h>
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <grp.h>


#include "src/common/fd.h"
#include "src/common/log.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/switch.h"
#include "src/common/xmalloc.h"
#include "src/common/xsignal.h"
#include "src/common/xstring.h"
#include "src/common/net.h"
#include "src/common/mpi.h"
#include "src/common/slurm_rlimits_info.h"
#include "src/common/plugstack.h"

#include "src/srun/allocate.h"
#include "src/srun/srun_job.h"
#include "src/srun/opt.h"
#include "src/srun/sigstr.h"
#include "src/srun/debugger.h"
#include "src/srun/srun.h"
#include "src/srun/srun_pty.h"
#include "src/srun/multi_prog.h"
#include "src/api/pmi_server.h"

#define MAX_RETRIES 20
#define MAX_ENTRIES 50

#define	TYPE_NOT_TEXT	0
#define	TYPE_TEXT	1
#define	TYPE_SCRIPT	2

mpi_plugin_client_info_t mpi_job_info[1];
pid_t srun_ppid = 0;
static struct termios termdefaults;
int global_rc;
srun_job_t *job = NULL;

struct {
	bitstr_t *start_success;
	bitstr_t *start_failure;
	bitstr_t *finish_normal;
	bitstr_t *finish_abnormal;
} task_state;

/*
 * forward declaration of static funcs
 */
static void  _print_job_information(resource_allocation_response_msg_t *resp);
static void  _set_prio_process_env(void);
static int   _set_rlimit_env(void);
static int   _set_umask_env(void);
static char *_uint16_array_to_str(int count, const uint16_t *array);
static int   _become_user (void);
static void  _run_srun_prolog (srun_job_t *job);
static void  _run_srun_epilog (srun_job_t *job);
static int   _run_srun_script (srun_job_t *job, char *script);
static int   _change_rlimit_rss(void);
static int   _slurm_debug_env_val (void);
static int   _call_spank_local_user (srun_job_t *job);
static void  _set_stdio_fds(srun_job_t *job, slurm_step_io_fds_t *cio_fds);
static void  _define_symbols(void);
static void  _pty_restore(void);
static void  _step_opt_exclusive(void);
static void _task_start(launch_tasks_response_msg_t *msg);
static void _task_finish(task_exit_msg_t *msg);
static void _job_complete();
static void _task_state_struct_init(int num_tasks);
static void _task_state_struct_print(void);
static void _task_state_struct_free(void);
static void _mpir_init(int num_tasks);
static void _mpir_cleanup(void);
static void _mpir_set_executable_names(const char *executable_name);
static void _mpir_dump_proctable(void);
static void _handle_intr();
static void _handle_signal(int signo);
static int _setup_signals();

int srun(int ac, char **av)
{
	resource_allocation_response_msg_t *resp;
	env_t *env = xmalloc(sizeof(env_t));
	uint32_t job_id = 0;
	log_options_t logopt = LOG_OPTS_STDERR_ONLY;
	slurm_step_launch_params_t launch_params;
	slurm_step_launch_callbacks_t callbacks;
	
	env->stepid = -1;
	env->procid = -1;
	env->localid = -1;
	env->nodeid = -1;
	env->cli = NULL;
	env->env = NULL;

	logopt.stderr_level += _slurm_debug_env_val();
	log_init(xbasename(av[0]), logopt, 0, NULL);

/* 	xsignal(SIGQUIT, _ignore_signal); */
/* 	xsignal(SIGPIPE, _ignore_signal); */
/* 	xsignal(SIGUSR1, _ignore_signal); */
/* 	xsignal(SIGUSR2, _ignore_signal); */

	/* Initialize plugin stack, read options from plugins, etc.
	 */
	if (spank_init(NULL) < 0) {
		fatal("Plug-in initialization failed");
		_define_symbols();
	}

	/* Be sure to call spank_fini when srun exits.
	 */
	if (atexit((void (*) (void)) spank_fini) < 0)
		error("Failed to register atexit handler for plugins: %m");
		
	/* set default options, process commandline arguments, and
	 * verify some basic values
	 */
	if (initialize_and_process_args(ac, av) < 0) {
		error ("srun initialization failed");
		exit (1);
	}
	srun_ppid = getppid();
	
	/* reinit log with new verbosity (if changed by command line)
	 */
	if (_verbose || opt.quiet) {
		/* If log level is already increased, only increment the
		 *   level to the difference of _verbose an LOG_LEVEL_INFO
		 */
		if ((_verbose -= (logopt.stderr_level - LOG_LEVEL_INFO)) > 0)
			logopt.stderr_level += _verbose;
		logopt.stderr_level -= opt.quiet;
		logopt.prefix_level = 1;
		log_alter(logopt, 0, NULL);
	}

	(void) _set_rlimit_env();
	_set_prio_process_env();
	(void) _set_umask_env();
	
	/* now global "opt" should be filled in and available,
	 * create a job from opt
	 */
	if (opt.test_only) {
		int rc = allocate_test();
		if (rc) {
			slurm_perror("allocation failure");
			exit (1);
		}
		info("allocation success");
		exit (0);

	} else if (opt.no_alloc) {
		info("do not allocate resources");
		job = job_create_noalloc(); 
		if (create_job_step(job) < 0) {
			exit(1);
		}
	} else if ((resp = existing_allocation())) {
		job_id = resp->job_id;
		if (opt.alloc_nodelist == NULL)
                       opt.alloc_nodelist = xstrdup(resp->node_list);
		if (opt.exclusive)
			_step_opt_exclusive();

		job = job_step_create_allocation(resp);
		slurm_free_resource_allocation_response_msg(resp);

		if (!job || create_job_step(job) < 0)
			exit(1);
	} else {
		/* Combined job allocation and job step launch */
#ifdef HAVE_FRONT_END
		uid_t my_uid = getuid();
		if ((my_uid != 0)
		&&  (my_uid != slurm_get_slurm_user_id())) {
			error("srun task launch not supported on this system");
			exit(1);
		}
#endif
		if (opt.job_max_memory > 0) {		
			(void) _change_rlimit_rss();
		}
	
		if ( !(resp = allocate_nodes()) ) 
			exit(1);
		_print_job_information(resp);
		job = job_create_allocation(resp);
		opt.exclusive = false;	/* not applicable for this step */
		if (!job || create_job_step(job) < 0) {
			exit(1);
		}
		
		slurm_free_resource_allocation_response_msg(resp);
	}

	/*
	 *  Become --uid user
	 */
	if (_become_user () < 0)
		info ("Warning: Unable to assume uid=%lu\n", opt.uid);

	/*
	 *  Enhance environment for job
	 */
	env->nprocs = opt.nprocs;
	env->cpus_per_task = opt.cpus_per_task;
	if (opt.ntasks_per_node != NO_VAL)
		env->ntasks_per_node = opt.ntasks_per_node;
	if (opt.ntasks_per_socket != NO_VAL)
		env->ntasks_per_socket = opt.ntasks_per_socket;
	if (opt.ntasks_per_core != NO_VAL)
		env->ntasks_per_core = opt.ntasks_per_core;
	env->distribution = opt.distribution;
	if (opt.plane_size != NO_VAL)
		env->plane_size = opt.plane_size;
	env->cpu_bind_type = opt.cpu_bind_type;
	env->cpu_bind = opt.cpu_bind;
	env->mem_bind_type = opt.mem_bind_type;
	env->mem_bind = opt.mem_bind;
	env->overcommit = opt.overcommit;
	env->slurmd_debug = opt.slurmd_debug;
	env->labelio = opt.labelio;
	env->comm_port = slurmctld_comm_addr.port;
	env->comm_hostname = slurmctld_comm_addr.hostname;
	if(job) {
		uint16_t *tasks = NULL;
		slurm_step_ctx_get(job->step_ctx, SLURM_STEP_CTX_TASKS, 
				   &tasks);

		env->select_jobinfo = job->select_jobinfo;
		env->nhosts = job->nhosts;
		env->nodelist = job->nodelist;
		env->task_count = _uint16_array_to_str(
			job->nhosts, tasks);
		env->jobid = job->jobid;
		env->stepid = job->stepid;
	}
	if (opt.pty) {
		struct termios term;
		int fd = STDIN_FILENO;

		/* Save terminal settings for restore */
		tcgetattr(fd, &termdefaults); 
		tcgetattr(fd, &term);
		/* Set raw mode on local tty */
		cfmakeraw(&term);
		tcsetattr(fd, TCSANOW, &term);
		atexit(&_pty_restore);

		set_winsize(job);
		block_sigwinch();
		pty_thread_create(job);
		env->pty_port = job->pty_port;
		env->ws_col   = job->ws_col;
		env->ws_row   = job->ws_row;
	}
	setup_env(env);
	xfree(env->task_count);
	xfree(env);
	
	_task_state_struct_init(opt.nprocs);
	slurm_step_launch_params_t_init(&launch_params);
	launch_params.gid = opt.gid;
	launch_params.argc = opt.argc;
	launch_params.argv = opt.argv;
	launch_params.multi_prog = opt.multi_prog ? true : false;
	launch_params.cwd = opt.cwd;
	launch_params.slurmd_debug = opt.slurmd_debug;
	launch_params.buffered_stdio = !opt.unbuffered;
	launch_params.labelio = opt.labelio ? true : false;
	launch_params.remote_output_filename =fname_remote_string(job->ofname);
	launch_params.remote_input_filename = fname_remote_string(job->ifname);
	launch_params.remote_error_filename = fname_remote_string(job->efname);
	launch_params.task_prolog = opt.task_prolog;
	launch_params.task_epilog = opt.task_epilog;
	launch_params.cpu_bind = opt.cpu_bind;
	launch_params.cpu_bind_type = opt.cpu_bind_type;
	launch_params.mem_bind = opt.mem_bind;
	launch_params.mem_bind_type = opt.mem_bind_type;	
	launch_params.pty = opt.pty;
	launch_params.max_sockets     = opt.max_sockets_per_node;
	launch_params.max_cores       = opt.max_cores_per_socket;
	launch_params.max_threads     = opt.max_threads_per_core;
	launch_params.cpus_per_task = opt.cpus_per_task;
	launch_params.ntasks_per_node   = opt.ntasks_per_node;
	launch_params.ntasks_per_socket = opt.ntasks_per_socket;
	launch_params.ntasks_per_core   = opt.ntasks_per_core;

	/* job structure should now be filled in */
	_setup_signals();

	_set_stdio_fds(job, &launch_params.local_fds);

	if (MPIR_being_debugged) {
		launch_params.parallel_debug = true;
		pmi_server_max_threads(1);
	} else {
		launch_params.parallel_debug = false;
	}
	callbacks.task_start = _task_start;
	callbacks.task_finish = _task_finish;
	callbacks.job_complete = _job_complete;
	callbacks.timeout_handler = timeout_handler;

	_run_srun_prolog(job);

	_mpir_init(job->ctx_params.task_count);

	if (_call_spank_local_user (job) < 0) {
		error("Failure in local plugin stack");
		slurm_step_launch_abort(job->step_ctx);
		exit(1);
	}

	update_job_state(job, SRUN_JOB_LAUNCHING);
	if (slurm_step_launch(job->step_ctx, &launch_params, &callbacks)
	    != SLURM_SUCCESS) {
		error("Application launch failed: %m");
		goto cleanup;
	}

	update_job_state(job, SRUN_JOB_STARTING);
	if (slurm_step_launch_wait_start(job->step_ctx) == SLURM_SUCCESS) {
		update_job_state(job, SRUN_JOB_RUNNING);
		/* Only set up MPIR structures if the step launched
		   correctly. */
		if (opt.multi_prog)
			mpir_set_multi_name(job->ctx_params.task_count,
					    launch_params.argv[0]);
		else
			_mpir_set_executable_names(launch_params.argv[0]);
		MPIR_debug_state = MPIR_DEBUG_SPAWNED;
		MPIR_Breakpoint();
		if (opt.debugger_test)
			_mpir_dump_proctable();
	} else {
		info("Job step aborted before step completely launched.");
	}

	slurm_step_launch_wait_finish(job->step_ctx);

cleanup:
	_run_srun_epilog(job);
	slurm_step_ctx_destroy(job->step_ctx);
	_mpir_cleanup();
	_task_state_struct_free();
	log_fini();

	return global_rc;
}

static int _call_spank_local_user (srun_job_t *job)
{
	struct spank_launcher_job_info info[1];
	job_step_create_response_msg_t *step_resp;

	info->uid = opt.uid;
	info->gid = opt.gid;
	info->jobid = job->jobid;
	info->stepid = job->stepid;
	slurm_step_ctx_get(job->step_ctx, SLURM_STEP_CTX_RESP, &step_resp);
	info->step_layout = step_resp->step_layout;
	info->argc = opt.argc;
	info->argv = opt.argv;

	return spank_local_user(info);
}


static int _slurm_debug_env_val (void)
{
	long int level = 0;
	const char *val;

	if ((val = getenv ("SLURM_DEBUG"))) {
		char *p;
		if ((level = strtol (val, &p, 10)) < -LOG_LEVEL_INFO)
			level = -LOG_LEVEL_INFO;
		if (p && *p != '\0')
			level = 0;
	}
	return ((int) level);
}

/*
 * Return a string representation of an array of uint32_t elements.
 * Each value in the array is printed in decimal notation and elements
 * are seperated by a comma.  If sequential elements in the array
 * contain the same value, the value is written out just once followed
 * by "(xN)", where "N" is the number of times the value is repeated.
 *
 * Example:
 *   The array "1, 2, 1, 1, 1, 3, 2" becomes the string "1,2,1(x3),3,2"
 *
 * Returns an xmalloc'ed string.  Free with xfree().
 */
static char *_uint16_array_to_str(int array_len, const uint16_t *array)
{
	int i;
	int previous = 0;
	char *sep = ",";  /* seperator */
	char *str = xstrdup("");

	if(array == NULL)
		return str;

	for (i = 0; i < array_len; i++) {
		if ((i+1 < array_len)
		    && (array[i] == array[i+1])) {
				previous++;
				continue;
		}

		if (i == array_len-1) /* last time through loop */
			sep = "";
		if (previous > 0) {
			xstrfmtcat(str, "%u(x%u)%s",
				   array[i], previous+1, sep);
		} else {
			xstrfmtcat(str, "%u%s", array[i], sep);
		}
		previous = 0;
	}
	
	return str;
}

static void 
_print_job_information(resource_allocation_response_msg_t *resp)
{
	int i;
	char tmp_str[10], job_details[4096];

	sprintf(job_details, "jobid %d: nodes(%d):`%s', cpu counts: ", 
	        resp->job_id, resp->node_cnt, resp->node_list);

	for (i = 0; i < resp->num_cpu_groups; i++) {
		sprintf(tmp_str, ",%u(x%u)", resp->cpus_per_node[i], 
		        resp->cpu_count_reps[i]);
		if (i == 0)
			strcat(job_details, &tmp_str[1]);
		else if ((strlen(tmp_str) + strlen(job_details)) < 
		         sizeof(job_details))
			strcat(job_details, tmp_str);
		else
			break;
	}
	verbose("%s",job_details);
}

/* Set SLURM_UMASK environment variable with current state */
static int _set_umask_env(void)
{
	char mask_char[5];
	mode_t mask;

	if (getenv("SLURM_UMASK"))	/* use this value */
		return SLURM_SUCCESS;

	mask = (int)umask(0);
	umask(mask);

	sprintf(mask_char, "0%d%d%d", 
		((mask>>6)&07), ((mask>>3)&07), mask&07);
	if (setenvf(NULL, "SLURM_UMASK", "%s", mask_char) < 0) {
		error ("unable to set SLURM_UMASK in environment");
		return SLURM_FAILURE;
	}
	debug ("propagating UMASK=%s", mask_char); 
	return SLURM_SUCCESS;
}

/*
 * _set_prio_process_env
 *
 * Set the internal SLURM_PRIO_PROCESS environment variable to support
 * the propagation of the users nice value and the "PropagatePrioProcess"
 * config keyword.
 */
static void  _set_prio_process_env(void)
{
	int retval;

	errno = 0; /* needed to detect a real failure since prio can be -1 */

	if ((retval = getpriority (PRIO_PROCESS, 0)) == -1)  {
		if (errno) {
			error ("getpriority(PRIO_PROCESS): %m");
			return;
		}
	}

	if (setenvf (NULL, "SLURM_PRIO_PROCESS", "%d", retval) < 0) {
		error ("unable to set SLURM_PRIO_PROCESS in environment");
		return;
	}

	debug ("propagating SLURM_PRIO_PROCESS=%d", retval);
}

/* 
 *  Change SLURM_RLIMIT_RSS to the user specified value --job-mem 
 *  or opt.job_max_memory 
 */
static int _change_rlimit_rss(void)
{
	struct rlimit        rlim[1];
	long                 new_cur;
	int                  rc = SLURM_SUCCESS;
	
	if (getrlimit (RLIMIT_RSS, rlim) < 0)
		return (error ("getrlimit (RLIMIT_RSS): %m"));

	new_cur = opt.job_max_memory*1024; 
	if((new_cur > rlim->rlim_max) || (new_cur < 0))
		rlim->rlim_cur = rlim->rlim_max;
	else
		rlim->rlim_cur = new_cur;

	if (setenvf (NULL, "SLURM_RLIMIT_RSS", "%lu", rlim->rlim_cur) < 0)
		error ("unable to set %s in environment", "RSS");

	if (setrlimit (RLIMIT_RSS, rlim) < 0) 
		return (error ("Unable to change memoryuse: %m"));

	return rc;
}

/* Set SLURM_RLIMIT_* environment variables with current resource 
 * limit values, reset RLIMIT_NOFILE to maximum possible value */
static int _set_rlimit_env(void)
{
	int                  rc = SLURM_SUCCESS;
	struct rlimit        rlim[1];
	unsigned long        cur;
	char                 name[64], *format;
	slurm_rlimits_info_t *rli;

	for (rli = get_slurm_rlimits_info(); rli->name != NULL; rli++ ) {

		if (getrlimit (rli->resource, rlim) < 0) {
			error ("getrlimit (RLIMIT_%s): %m", rli->name);
			rc = SLURM_FAILURE;
			continue;
		}
		
		cur = (unsigned long) rlim->rlim_cur;
		snprintf(name, sizeof(name), "SLURM_RLIMIT_%s", rli->name);
		if (opt.propagate && rli->propagate_flag == PROPAGATE_RLIMITS)
			/*
			 * Prepend 'U' to indicate user requested propagate
			 */
			format = "U%lu";
		else
			format = "%lu";
		
		if (setenvf (NULL, name, format, cur) < 0) {
			error ("unable to set %s in environment", name);
			rc = SLURM_FAILURE;
			continue;
		}
		
		debug ("propagating RLIMIT_%s=%lu", rli->name, cur);
	}

	/* 
	 *  Now increase NOFILE to the max available for this srun
	 */
	if (getrlimit (RLIMIT_NOFILE, rlim) < 0)
	 	return (error ("getrlimit (RLIMIT_NOFILE): %m"));

	if (rlim->rlim_cur < rlim->rlim_max) {
		rlim->rlim_cur = rlim->rlim_max;
		if (setrlimit (RLIMIT_NOFILE, rlim) < 0) 
			return (error ("Unable to increase max no. files: %m"));
	}

	return rc;
}

static int _become_user (void)
{
	struct passwd *pwd = getpwuid (opt.uid);

	if (opt.uid == getuid ())
		return (0);

	if ((opt.egid != (gid_t) -1) && (setgid (opt.egid) < 0))
		return (error ("setgid: %m"));

	initgroups (pwd->pw_name, pwd->pw_gid); /* Ignore errors */

	if (setuid (opt.uid) < 0)
		return (error ("setuid: %m"));

	return (0);
}

static void _run_srun_prolog (srun_job_t *job)
{
	int rc;

	if (opt.prolog && strcasecmp(opt.prolog, "none") != 0) {
		rc = _run_srun_script(job, opt.prolog);
		debug("srun prolog rc = %d", rc);
	}
}

static void _run_srun_epilog (srun_job_t *job)
{
	int rc;

	if (opt.epilog && strcasecmp(opt.epilog, "none") != 0) {
		rc = _run_srun_script(job, opt.epilog);
		debug("srun epilog rc = %d", rc);
	}
}

static int _run_srun_script (srun_job_t *job, char *script)
{
	int status;
	pid_t cpid;
	int i;
	char **args = NULL;

	if (script == NULL || script[0] == '\0')
		return 0;

	if (access(script, R_OK | X_OK) < 0) {
		info("Access denied for %s: %m", script);
		return 0;
	}

	if ((cpid = fork()) < 0) {
		error ("run_srun_script: fork: %m");
		return -1;
	}
	if (cpid == 0) {

		/* set the scripts command line arguments to the arguments
		 * for the application, but shifted one higher
		 */
		args = xmalloc(sizeof(char *) * 1024);
		args[0] = script;
		for (i = 0; i < opt.argc; i++) {
			args[i+1] = opt.argv[i];
		}
		args[i+1] = NULL;
		execv(script, args);
		error("help! %m");
		exit(127);
	}

	do {
		if (waitpid(cpid, &status, 0) < 0) {
			if (errno == EINTR)
				continue;
			error("waidpid: %m");
			return 0;
		} else
			return status;
	} while(1);

	/* NOTREACHED */
}

static int
_is_local_file (fname_t *fname)
{
	if (fname->name == NULL)
		return 1;
	
	if (fname->taskid != -1)
		return 1;

	return ((fname->type != IO_PER_TASK) && (fname->type != IO_ONE));
}

static void
_set_stdio_fds(srun_job_t *job, slurm_step_io_fds_t *cio_fds)
{
	bool err_shares_out = false;

	/*
	 * create stdin file descriptor
	 */
	if (_is_local_file(job->ifname)) {
		if (job->ifname->name == NULL || job->ifname->taskid != -1) {
			cio_fds->in.fd = STDIN_FILENO;
		} else {
			cio_fds->in.fd = open(job->ifname->name, O_RDONLY);
			if (cio_fds->in.fd == -1)
				fatal("Could not open stdin file: %m");
		}
		if (job->ifname->type == IO_ONE) {
			job_step_create_response_msg_t *step_resp = NULL;
			
			slurm_step_ctx_get(job->step_ctx, SLURM_STEP_CTX_RESP,
					   &step_resp);
		
			cio_fds->in.taskid = job->ifname->taskid;
			cio_fds->in.nodeid = slurm_step_layout_host_id(
				step_resp->step_layout, job->ifname->taskid);
		}
	}

	/*
	 * create stdout file descriptor
	 */
	if (_is_local_file(job->ofname)) {
		if (job->ofname->name == NULL) {
			cio_fds->out.fd = STDOUT_FILENO;
		} else {
			cio_fds->out.fd = open(job->ofname->name,
					       O_CREAT|O_WRONLY|O_TRUNC, 0644);
			if (cio_fds->out.fd == -1)
				fatal("Could not open stdout file: %m");
		}
		if (job->ofname->name != NULL
		    && job->efname->name != NULL
		    && !strcmp(job->ofname->name, job->efname->name)) {
			err_shares_out = true;
		}
	}

	/*
	 * create seperate stderr file descriptor only if stderr is not sharing
	 * the stdout file descriptor
	 */
	if (err_shares_out) {
		debug3("stdout and stderr sharing a file");
		cio_fds->err.fd = cio_fds->out.fd;
		cio_fds->err.taskid = cio_fds->out.taskid;
	} else if (_is_local_file(job->efname)) {
		if (job->efname->name == NULL) {
			cio_fds->err.fd = STDERR_FILENO;
		} else {
			cio_fds->err.fd = open(job->efname->name,
					       O_CREAT|O_WRONLY|O_TRUNC, 0644);
			if (cio_fds->err.fd == -1)
				fatal("Could not open stderr file: %m");
		}
	}
}

/* Plugins must be able to resolve symbols.
 * Since srun statically links with src/api/libslurmhelper rather than 
 * dynamicaly linking with libslurm, we need to reference all needed 
 * symbols within srun. None of the functions below are actually 
 * used, but we need to load the symbols. */
static void _define_symbols(void)
{
	slurm_signal_job_step(0,0,0);	/* needed by mvapich and mpichgm */
}

static void _pty_restore(void)
{
	/* STDIN is probably closed by now */
	if (tcsetattr(STDOUT_FILENO, TCSANOW, &termdefaults) < 0)
		fprintf(stderr, "tcsetattr: %s\n", strerror(errno));
}

/* opt.exclusive is set, disable user task layout controls */
static void _step_opt_exclusive(void)
{
	if (!opt.nprocs_set)
		fatal("--nprocs must be set with --exclusive");
	if (opt.relative_set)
		fatal("--relative disabled, incompatible with --exclusive");
	if (opt.exc_nodes)
		fatal("--exclude is incompatible with --exclusive");
	if (opt.nodelist)
		fatal("--nodelist is incompatible with --exclusive");
}

static void
_task_start(launch_tasks_response_msg_t *msg)
{
	MPIR_PROCDESC *table;
	int taskid;
	int i;

	verbose("Node %s (%d), %d tasks started",
		msg->node_name, msg->srun_node_id, msg->count_of_pids);

	for (i = 0; i < msg->count_of_pids; i++) {
		taskid = msg->task_ids[i];
		table = &MPIR_proctable[taskid];
		table->host_name = xstrdup(msg->node_name);
		/* table->executable_name is set elsewhere */
		table->pid = msg->local_pids[i];

		if (msg->return_code == 0) {
			bit_set(task_state.start_success, taskid);
		} else {
			bit_set(task_state.start_failure, taskid);
		}
	}

}

static void
_terminate_job_step(slurm_step_ctx_t *step_ctx)
{
	uint32_t job_id, step_id;

	slurm_step_ctx_get(step_ctx, SLURM_STEP_CTX_JOBID, &job_id);
	slurm_step_ctx_get(step_ctx, SLURM_STEP_CTX_STEPID, &step_id);
	info("Terminating job step %u.%u", job_id, step_id);
	slurm_kill_job_step(job_id, step_id, SIGKILL);
}

static void
_handle_max_wait(int signo)
{
	info("First task exited %ds ago", opt.max_wait);
	_task_state_struct_print();
	_terminate_job_step(job->step_ctx);
}

static void
_task_finish(task_exit_msg_t *msg)
{
	static bool first_done = true;
	static bool first_error = true;
	int rc = 0;
	int i;

	verbose("%d tasks finished (rc=%u)",
		msg->num_tasks, msg->return_code);
	if (WIFEXITED(msg->return_code)) {
		rc = WEXITSTATUS(msg->return_code);
		if (rc != 0) {
			for (i = 0; i < msg->num_tasks; i++) {
				error("task %u exited with exit code %d",
				      msg->task_id_list[i], rc);
				bit_set(task_state.finish_abnormal,
					msg->task_id_list[i]);
			}
		} else {
			for (i = 0; i < msg->num_tasks; i++) {
				bit_set(task_state.finish_normal,
					msg->task_id_list[i]);
			}
		}
	} else if (WIFSIGNALED(msg->return_code)) {
		for (i = 0; i < msg->num_tasks; i++) {
			verbose("task %u killed by signal %d",
				msg->task_id_list[i],
				WTERMSIG(msg->return_code));
			bit_set(task_state.finish_abnormal,
				msg->task_id_list[i]);
		}
		rc = 1;
	}
	global_rc = MAX(global_rc, rc);

	if (first_error && rc > 0 && opt.kill_bad_exit) {
		first_error = false;
		_terminate_job_step(job->step_ctx);
	} else if (first_done && opt.max_wait > 0) {
		/* If these are the first tasks to finish we need to
		 * start a timer to kill off the job step if the other
		 * tasks don't finish within opt.max_wait seconds.
		 */
		first_done = false;
		debug2("First task has exited");
		xsignal(SIGALRM, _handle_max_wait);
		verbose("starting alarm of %d seconds", opt.max_wait);
		alarm(opt.max_wait);
	}
}

/* This typically signifies the job was cancelled by scancel */
static void
_job_complete()
{
	info("Force Terminated job");
}

static void
_task_state_struct_init(int num_tasks)
{
	task_state.start_success = bit_alloc(num_tasks);
	task_state.start_failure = bit_alloc(num_tasks);
	task_state.finish_normal = bit_alloc(num_tasks);
	task_state.finish_abnormal = bit_alloc(num_tasks);
}

/*
 * Tasks will most likely have bits set in multiple of the task_state
 * bit strings (e.g. a task can start normally and then later exit normally)
 * so we ensure that a task is only "seen" once.
 */
static void
_task_state_struct_print(void)
{
	bitstr_t *tmp, *seen, *not_seen;
	char buf[BUFSIZ];
	int len;

	len = bit_size(task_state.finish_abnormal); /* all the same length */
	tmp = bit_alloc(len);
	seen = bit_alloc(len);
	not_seen = bit_alloc(len);
	bit_not(not_seen);

	if (bit_set_count(task_state.finish_abnormal) > 0) {
		bit_copybits(tmp, task_state.finish_abnormal);
		bit_and(tmp, not_seen);
		bit_fmt(buf, BUFSIZ, tmp);
		info("task%s: exited abnormally", buf);
		bit_or(seen, tmp);
		bit_copybits(not_seen, seen);
		bit_not(not_seen);
	}

	if (bit_set_count(task_state.finish_normal) > 0) {
		bit_copybits(tmp, task_state.finish_normal);
		bit_and(tmp, not_seen);
		bit_fmt(buf, BUFSIZ, tmp);
		info("task%s: exited", buf);
		bit_or(seen, tmp);
		bit_copybits(not_seen, seen);
		bit_not(not_seen);
	}

	if (bit_set_count(task_state.start_failure) > 0) {
		bit_copybits(tmp, task_state.start_failure);
		bit_and(tmp, not_seen);
		bit_fmt(buf, BUFSIZ, tmp);
		info("task%s: failed to start", buf);
		bit_or(seen, tmp);
		bit_copybits(not_seen, seen);
		bit_not(not_seen);
	}

	if (bit_set_count(task_state.start_success) > 0) {
		bit_copybits(tmp, task_state.start_success);
		bit_and(tmp, not_seen);
		bit_fmt(buf, BUFSIZ, tmp);
		info("task%s: running", buf);
		bit_or(seen, tmp);
		bit_copybits(not_seen, seen);
		bit_not(not_seen);
	}
}

static void
_task_state_struct_free(void)
{
	bit_free(task_state.start_success);
	bit_free(task_state.start_failure);
	bit_free(task_state.finish_normal);
	bit_free(task_state.finish_abnormal);
}

/**********************************************************************
 * Functions for manipulating the MPIR_* global variables which
 * are accessed by parallel debuggers which trace slaunch.
 **********************************************************************/
static void
_mpir_init(int num_tasks)
{
	MPIR_proctable_size = num_tasks;
	MPIR_proctable = xmalloc(sizeof(MPIR_PROCDESC) * num_tasks);
	if (MPIR_proctable == NULL)
		fatal("Unable to initialize MPIR_proctable: %m");
}

static void
_mpir_cleanup()
{
	int i;

	for (i = 0; i < MPIR_proctable_size; i++) {
		xfree(MPIR_proctable[i].host_name);
		xfree(MPIR_proctable[i].executable_name);
	}
	xfree(MPIR_proctable);
}

static void
_mpir_set_executable_names(const char *executable_name)
{
	int i;

	for (i = 0; i < MPIR_proctable_size; i++) {
		MPIR_proctable[i].executable_name = xstrdup(executable_name);
		if (MPIR_proctable[i].executable_name == NULL)
			fatal("Unable to set MPI_proctable executable_name:"
			      " %m");
	}
}

static void
_mpir_dump_proctable()
{
	MPIR_PROCDESC *tv;
	int i;

	for (i = 0; i < MPIR_proctable_size; i++) {
		tv = &MPIR_proctable[i];
		if (!tv)
			break;
		info("task:%d, host:%s, pid:%d, executable:%s",
		     i, tv->host_name, tv->pid, tv->executable_name);
	}
}
	
static void _handle_intr()
{
	static time_t last_intr      = 0;
	static time_t last_intr_sent = 0;
	if (opt.quit_on_intr) {
		job_force_termination(job);
		slurm_step_launch_abort(job->step_ctx);
		return;
	}

	if (((time(NULL) - last_intr) > 1) && !opt.disable_status) {
		info("interrupt (one more within 1 sec to abort)");
		_task_state_struct_print();
		last_intr = time(NULL);
	} else  { /* second Ctrl-C in half as many seconds */
		update_job_state(job, SRUN_JOB_CANCELLED);
		/* terminate job */
		if (job->state < SRUN_JOB_FORCETERM) {
			if ((time(NULL) - last_intr_sent) < 1) {
				job_force_termination(job);
				slurm_step_launch_abort(job->step_ctx);
				return;
			}

			info("sending Ctrl-C to job");
			last_intr_sent = time(NULL);
			slurm_step_launch_fwd_signal(job->step_ctx, SIGINT);

		} else {
			job_force_termination(job);
			slurm_step_launch_abort(job->step_ctx);
		}
	}
}

static void _handle_signal(int signo)
{
	debug2("got signal %d", signo);

	switch (signo) {
	case SIGINT:
		_handle_intr();
		break;
	case SIGQUIT:
		info("Quit");
		/* continue with skurm_step_launch_abort */
	case SIGTERM:
	case SIGHUP:
		job_force_termination(job);
		slurm_step_launch_abort(job->step_ctx);
		break;
	/* case SIGTSTP: */
/* 		debug3("got SIGTSTP"); */
/* 		break; */
	case SIGCONT:
		debug3("got SIGCONT");
		break;
	default:
		slurm_step_launch_fwd_signal(job->step_ctx, signo);
		break;
	}
}

static int _setup_signals()
{
	int sigarray[] = {
		SIGINT,  SIGQUIT, /*SIGTSTP,*/ SIGCONT, SIGTERM,
		SIGALRM, SIGUSR1, SIGUSR2, SIGPIPE, 0
	};
	int rc = SLURM_SUCCESS, i=0, signo;

	xassert(job);
	xassert(job->step_ctx);

	while ((signo = sigarray[i++])) 
		xsignal(signo, _handle_signal);

	return rc;
}

