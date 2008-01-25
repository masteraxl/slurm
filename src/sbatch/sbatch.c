/*****************************************************************************\
 *  sbatch.c - Submit a SLURM batch script.
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

#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <slurm/slurm.h>

#include "src/common/env.h"
#include "src/common/plugstack.h"
#include "src/common/read_config.h"
#include "src/common/slurm_rlimits_info.h"
#include "src/common/xstring.h"
#include "src/common/xmalloc.h"

#include "src/sbatch/opt.h"

#define MAX_RETRIES 3

static void _call_spank_local_user(job_desc_msg_t desc,
				   submit_response_msg_t *resp);
static int   fill_job_desc_from_opts(job_desc_msg_t *desc);
static void *get_script_buffer(const char *filename, int *size);
static void  set_prio_process_env(void);
static int   set_umask_env(void);
static char *script_wrap(char *command_string);
static int  _set_rlimit_env(void);

int main(int argc, char *argv[])
{
	log_options_t logopt = LOG_OPTS_STDERR_ONLY;
	job_desc_msg_t desc;
	submit_response_msg_t *resp;
	char *script_name;
	void *script_body;
	int script_size = 0;
	int retries = 0;

	log_init(xbasename(argv[0]), logopt, 0, NULL);
	if (spank_init(NULL) < 0)
		fatal("Plug-in initialization failed");

	script_name = process_options_first_pass(argc, argv);
	/* reinit log with new verbosity (if changed by command line) */
	if (opt.verbose || opt.quiet) {
		logopt.stderr_level += opt.verbose;
		logopt.stderr_level -= opt.quiet;
		logopt.prefix_level = 1;
		log_alter(logopt, 0, NULL);
	}

	if (opt.wrap != NULL) {
		script_body = script_wrap(opt.wrap);
	} else {
		script_body = get_script_buffer(script_name, &script_size);
	}
	if (script_body == NULL)
		exit(1);

	if (process_options_second_pass((argc - opt.script_argc), argv,
					script_body, script_size) < 0) {
		fatal("sbatch parameter parsing");
	}

	(void) _set_rlimit_env();
	set_prio_process_env();
	set_umask_env();
	slurm_init_job_desc_msg(&desc);
	if (fill_job_desc_from_opts(&desc) == -1) {
		exit(2);
	}

	desc.script = (char *)script_body;

	while (slurm_submit_batch_job(&desc, &resp) < 0) {
		static char *msg = "Slurm job queue full, sleeping and retrying.";

		if ((errno != ESLURM_ERROR_ON_DESC_TO_RECORD_COPY) ||
		    (retries >= MAX_RETRIES)) {
			error("Batch job submission failed: %m");
			exit(3);
		}

		if (retries)
			debug(msg);
		else
			error(msg);
		sleep (++retries);
        }
	_call_spank_local_user(desc, resp);
	info("Submitted batch job %d", resp->job_id);
	xfree(desc.script);
	slurm_free_submit_response_response_msg(resp);
	spank_fini(NULL);
	return 0;
}

static void _call_spank_local_user(job_desc_msg_t desc,
				   submit_response_msg_t *resp)
{
	struct spank_launcher_job_info info[1];

	info->uid = desc.user_id;
	info->gid = desc.group_id;
	info->jobid = resp->job_id;
	info->stepid = SLURM_BATCH_SCRIPT;
	info->step_layout = NULL;
	info->argc = desc.argc;
	info->argv = desc.argv;

	if (spank_local_user(info) < 0)
		error("spank_local_user: %m");
}

/* Returns 0 on success, -1 on failure */
static int fill_job_desc_from_opts(job_desc_msg_t *desc)
{
	extern char **environ;

	if (opt.jobid_set)
		desc->job_id = opt.jobid;
	desc->contiguous = opt.contiguous ? 1 : 0;
	desc->features = opt.constraints;
	desc->immediate = opt.immediate;
	if (opt.job_name != NULL)
		desc->name = opt.job_name;
	else
		desc->name = xstrdup("sbatch");
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
	if (opt.conn_type != (uint16_t) NO_VAL)
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
	if (opt.time_limit != NO_VAL)
		desc->time_limit = opt.time_limit;
	desc->shared = opt.shared;

	desc->environment = NULL;
	if (opt.get_user_env_time >= 0) {
		struct passwd *pw = NULL;
		pw = getpwuid(opt.uid);
		if (pw != NULL) {
			desc->environment = env_array_user_default(pw->pw_name,
						opt.get_user_env_time,
						opt.get_user_env_mode);
			/* FIXME - should we abort if j->environment
			 * is NULL? */
		}
	}
	env_array_merge(&desc->environment, (const char **)environ);
	desc->env_size = envcount (desc->environment);
	desc->argv = opt.script_argv;
	desc->argc = opt.script_argc;
	desc->err  = opt.efname;
	desc->in   = opt.ifname;
	desc->out  = opt.ofname;
	desc->work_dir = opt.cwd;
	desc->no_requeue = opt.no_requeue;
	if (opt.open_mode)
		desc->open_mode = opt.open_mode;
	if (opt.acctg_freq >= 0)
		desc->acctg_freq = opt.acctg_freq;

	return 0;
}

/* Set SLURM_UMASK environment variable with current state */
static int set_umask_env(void)
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
 * set_prio_process_env
 *
 * Set the internal SLURM_PRIO_PROCESS environment variable to support
 * the propagation of the users nice value and the "PropagatePrioProcess"
 * config keyword.
 */
static void  set_prio_process_env(void)
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
 * Checks if the buffer starts with a shebang (#!).
 */
static bool has_shebang(const void *buf, int size)
{
	char *str = (char *)buf;

	if (size < 2)
		return false;

	if (str[0] != '#' || str[1] != '!')
		return false;

	return true;
}

/*
 * Checks if the buffer contains a NULL character (\0).
 */
static bool contains_null_char(const void *buf, int size)
{
	char *str = (char *)buf;
	int i;

	for (i = 0; i < size; i++) {
		if (str[i] == '\0')
			return true;
	}

	return false;
}

/*
 * If "filename" is NULL, the batch script is read from standard input.
 */
static void *get_script_buffer(const char *filename, int *size)
{
	int fd;
	char *buf = NULL;
	int buf_size = BUFSIZ;
	int buf_left;
	int script_size = 0;
	char *ptr = NULL;
	int tmp_size;

	/*
	 * First figure out whether we are reading from STDIN_FILENO
	 * or from a file.
	 */
	if (filename == NULL) {
		fd = STDIN_FILENO;
	} else {
		fd = open(filename, O_RDONLY);
		if (fd == -1) {
			error("Unable to open file %s", filename);
			goto fail;
		}
	}

	/*
	 * Then read in the script.
	 */
	buf = ptr = xmalloc(buf_size);
	buf_left = buf_size;
	while((tmp_size = read(fd, ptr, buf_left)) > 0) {
		buf_left -= tmp_size;
		script_size += tmp_size;
		if (buf_left == 0) {
			buf_size += BUFSIZ;
			xrealloc(buf, buf_size);
		}
		ptr = buf + script_size;
		buf_left = buf_size - script_size;
	}
	close(fd);

	/*
	 * Finally we perform some sanity tests on the script.
	 */
	if (script_size == 0) {
		error("Batch script is empty!");
		goto fail;
	} else if (script_size >= 0xffff) {
		error("Job script exceeds size supported by slurm");
		goto fail;
	} else if (xstring_is_whitespace(buf)) {
		error("Batch script contains only whitespace!");
		goto fail;
	} else if (!has_shebang(buf, script_size)) {
		error("This does not look like a batch script.  The first");
		error("line must start with #! followed by the path"
		      " to an interpreter.");
		error("For instance: #!/bin/sh");
		goto fail;
	} else if (contains_null_char(buf, script_size)) {
		error("The SLURM controller does not allow scripts that");
		error("contain a NULL character '\\0'.");
		goto fail;
	}

	*size = script_size;
	return buf;
fail:
	xfree(buf);
	*size = 0;
	return NULL;
}

/* Wrap a single command string in a simple shell script */
static char *script_wrap(char *command_string)
{
	char *script = NULL;

	xstrcat(script, "#!/bin/sh\n");
	xstrcat(script, "# This script was created by sbatch --wrap.\n\n");
	xstrcat(script, command_string);
	xstrcat(script, "\n");

	return script;
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

	/* Load default limits to be propagated from slurm.conf */
	slurm_conf_lock();
	slurm_conf_unlock();

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
