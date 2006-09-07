/*****************************************************************************\
 *  sbatch.c - Submit a SLURM batch script.
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher J. Morrone <morrone2@llnl.gov>
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <slurm/slurm.h>

#include "src/common/xstring.h"
#include "src/common/xmalloc.h"
#include "src/common/env.h"

#include "src/sbatch/opt.h"

static int fill_job_desc_from_opts(job_desc_msg_t *desc);
static void *get_script_buffer(const char *filename, int *size);
static int set_umask_env(void);

int main(int argc, char *argv[])
{
	log_options_t logopt = LOG_OPTS_STDERR_ONLY;
	job_desc_msg_t desc;
	submit_response_msg_t *resp;
	char *script_name;
	void *script_body;
	int script_size = 0;

	log_init(xbasename(argv[0]), logopt, 0, NULL);
	script_name = process_options_first_pass(argc, argv);
	/* reinit log with new verbosity (if changed by command line) */
	if (opt.verbose || opt.quiet) {
		logopt.stderr_level += opt.verbose;
		logopt.stderr_level -= opt.quiet;
		logopt.prefix_level = 1;
		log_alter(logopt, 0, NULL);
	}

	script_body = get_script_buffer(script_name, &script_size);
	if (script_body == NULL)
		exit(1);

	if (process_options_second_pass((argc - opt.script_argc), argv,
					script_body, script_size) < 0) {
		fatal("sbatch parameter parsing");
	}

	set_umask_env();
	slurm_init_job_desc_msg(&desc);
	if (fill_job_desc_from_opts(&desc) == -1) {
		exit(2);
	}

	desc.script = (char *)script_body;

	if (slurm_submit_batch_job(&desc, &resp) == -1) {
		error("Batch job submission failed: %m");
		exit(3);
	}

	info("Submitted batch job %d", resp->job_id);
	xfree(desc.script);
	slurm_free_submit_response_response_msg(resp);

	return 0;
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
	desc->dependency = opt.dependency;
	if (opt.nice)
		desc->nice = NICE_OFFSET + opt.nice;
	desc->mail_type = opt.mail_type;
	if (opt.mail_user)
		desc->mail_user = xstrdup(opt.mail_user);
	if (opt.begin)
		desc->begin_time = opt.begin;
	if (opt.account)
		desc->account = xstrdup(opt.account);

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
	if (opt.no_rotate)
		desc->rotate = 0;
	if (opt.mincpus > -1)
		desc->min_procs = opt.mincpus;
	if (opt.realmem > -1)
		desc->min_memory = opt.realmem;
	if (opt.tmpdisk > -1)
		desc->min_tmp_disk = opt.tmpdisk;
	desc->num_procs = opt.nprocs * opt.cpus_per_task;
	if (opt.nprocs_set)
		desc->num_tasks = opt.nprocs;
	if (opt.cpus_set)
		desc->cpus_per_task = opt.cpus_per_task;
	if (opt.no_kill)
		desc->kill_on_node_fail = 0;
	if (opt.time_limit > -1)
		desc->time_limit = opt.time_limit;
	desc->shared = opt.shared;

	desc->environment = environ;
	desc->env_size = envcount (environ);
	desc->argv = opt.script_argv;
	desc->argc = opt.script_argc;
	desc->err  = opt.efname;
	desc->in   = opt.ifname;
	desc->out  = opt.ofname;
	desc->work_dir = opt.cwd;
	desc->no_requeue = opt.no_requeue;

	return 0;
}

/* Set SLURM_UMASK environment variable with current state */
static int set_umask_env(void)
{
	char mask_char[5];
	mode_t mask = (int)umask(0);
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
