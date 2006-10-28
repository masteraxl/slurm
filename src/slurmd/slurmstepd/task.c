/*****************************************************************************\
 *  slurmd/slurmstepd/task.c - task launching functions for slurmstepd
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark A. Grondona <mgrondona@llnl.gov>.
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

#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <unistd.h>
#include <pwd.h>
#include <grp.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#if HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#if HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_AIX
#  include <sys/checkpnt.h>
#endif

#include <sys/resource.h>

#include <slurm/slurm_errno.h>

#include "src/common/env.h"
#include "src/common/fd.h"
#include "src/common/log.h"
#include "src/common/slurm_jobacct.h"
#include "src/common/switch.h"
#include "src/common/xsignal.h"
#include "src/common/xstring.h"
#include "src/common/mpi.h"
#include "src/common/xmalloc.h"
#include "src/common/plugstack.h"

#include "src/slurmd/slurmd/slurmd.h"
#include "src/slurmd/common/proctrack.h"
#include "src/slurmd/common/task_plugin.h"
#include "src/slurmd/slurmstepd/task.h"
#include "src/slurmd/slurmstepd/ulimits.h"
#include "src/slurmd/slurmstepd/io.h"
#include "src/slurmd/slurmstepd/pdebug.h"
#include "src/slurmd/slurmstepd/task_exec.h"

/*
 * Static prototype definitions.
 */
static void  _make_tmpdir(slurmd_job_t *job);
static int   _run_script(const char *name, const char *path, 
		slurmd_job_t *job);
static void  _update_env(char *buf, char ***env);

/* Search for "export NAME=value" records in buf and 
 * use them to add environment variables to env */
static void
_update_env(char *buf, char ***env)
{
	char *tmp_ptr, *name_ptr, *val_ptr, *buf_ptr = buf;

	while ((tmp_ptr = strstr(buf_ptr, "export"))) {
		buf_ptr += 6;
		while (isspace(buf_ptr[0]))
			buf_ptr++;
		if (buf_ptr[0] == '=')	/* mal-formed */
			continue;
		name_ptr = buf_ptr;	/* start of env var name */
		while ((buf_ptr[0] != '=') && (buf_ptr[0] != '\0'))
			buf_ptr++;
		if (buf_ptr[0] == '\0')	/* mal-formed */
			continue;
		buf_ptr[0] = '\0';	/* end of env var name */
		buf_ptr++;
		val_ptr = buf_ptr;	/* start of env var value */
		while ((!isspace(buf_ptr[0])) && (buf_ptr[0] != '\0'))
			buf_ptr++;
		if (isspace(buf_ptr[0])) {
			buf_ptr[0] = '\0';/* end of env var value */
			buf_ptr++;
		}
		debug("name:%s:val:%s:",name_ptr,val_ptr);
		if (setenvf(env, name_ptr, "%s", val_ptr))
			error("Unable to set %s environment variable", name_ptr);
	}		
}

/*
 * Run a task prolog script
 * name IN: class of program ("system prolog", "user prolog", etc.),
 *	if prefix is "user" then also set uid
 * path IN: pathname of program to run
 * job IN/OUT: pointer to associated job, can update job->env 
 *	if prolog
 * RET 0 on success, -1 on failure.
 */
static int
_run_script(const char *name, const char *path, slurmd_job_t *job)
{
	int status, rc, nread;
	pid_t cpid;
	int pfd[2];
	char buf[4096];

	xassert(job->env);
	if (path == NULL || path[0] == '\0')
		return 0;

	debug("[job %u] attempting to run %s [%s]", job->jobid, name, path);

	if (access(path, R_OK | X_OK) < 0) {
		debug("Not running %s [%s]: %m", name, path);
		return 0;
	}
	if (pipe(pfd) < 0) {
		error("executing %s: pipe: %m", name);
		return -1;
	}
	if ((cpid = fork()) < 0) {
		error("executing %s: fork: %m", name);
		return -1;
	}
	if (cpid == 0) {
		char *argv[2];

		argv[0] = xstrdup(path);
		argv[1] = NULL;
		close(1);
		dup(pfd[1]);
		close(2);
		close(0);
#ifdef SETPGRP_TWO_ARGS
		setpgrp(0, 0);
#else
		setpgrp();
#endif
		execve(path, argv, job->env);
		error("execve(): %m");
		exit(127);
	}

	close(pfd[1]);
	while ((nread = read(pfd[0], buf, sizeof(buf))) > 0) {
		buf[nread] = 0;
		//debug("read %d:%s:", nread, buf);
		_update_env(buf, &job->env);
	}

	close(pfd[0]);
	while (1) {
		rc = waitpid(cpid, &status, 0);
		if (rc < 0) {
			if (errno == EINTR)
				continue;
			error("waidpid: %m");
			return 0;
		} else  {
			killpg(cpid, SIGKILL);  /* kill children too */
			return status;
		}
	}

	/* NOTREACHED */
}

/*
 *  Current process is running as the user when this is called.
 */
void
exec_task(slurmd_job_t *job, int i, int waitfd)
{
	char c;
	int rc;
	slurmd_task_info_t *task = job->task[i];

	if (set_user_limits(job) < 0) {
		debug("Unable to set user limits");
		log_fini();
		exit(5);
	}

	if (i == 0)
		_make_tmpdir(job);

        /*
	 * Stall exec until all tasks have joined the same process group
	 */
        if ((rc = read (waitfd, &c, sizeof (c))) != 1) {
	        error ("_exec_task read failed, fd = %d, rc=%d: %m", waitfd, rc);
		log_fini();
		exit(1);
	}
	close(waitfd);

	job->envtp->jobid = job->jobid;
	job->envtp->stepid = job->stepid;
	job->envtp->nodeid = job->nodeid;
	job->envtp->cpus_on_node = job->cpus;
	job->envtp->env = job->env;
	
	job->envtp->procid = task->gtid;
	job->envtp->localid = task->id;
	job->envtp->task_pid = getpid();

	job->envtp->distribution = job->task_dist;
	job->envtp->plane_size   = job->plane_size;

	job->envtp->cpu_bind = xstrdup(job->cpu_bind);
	job->envtp->cpu_bind_type = job->cpu_bind_type;
	job->envtp->mem_bind = xstrdup(job->mem_bind);
	job->envtp->mem_bind_type = job->mem_bind_type;

	job->envtp->distribution = -1;
	setup_env(job->envtp);
	setenvf(&job->envtp->env, "SLURMD_NODENAME", "%s", conf->node_name);
	
	job->env = job->envtp->env;
	job->envtp->env = NULL;
	xfree(job->envtp->task_count);
	
	if (!job->batch) {
		if (interconnect_attach(job->switch_job, &job->env,
				job->nodeid, (uint32_t) i, job->nnodes,
				job->nprocs, task->gtid) < 0) {
			error("Unable to attach to interconnect: %m");
			log_fini();
			exit(1);
		}

		slurmd_mpi_init (job, task->gtid);
	
		pdebug_stop_current(job);
	}

	io_dup_stdio(task);

	/* task-specific pre-launch activities */

	if (spank_user_task (job, i) < 0) {
		error ("Failed to invoke task plugin stack\n");
		exit (1);
	}

	pre_launch(job);

	if (conf->task_prolog) {
		char *my_prolog;
		slurm_mutex_lock(&conf->config_mutex);
		my_prolog = xstrdup(conf->task_prolog);
		slurm_mutex_unlock(&conf->config_mutex);
		_run_script("slurm task_prolog", my_prolog, job);
		xfree(my_prolog);
	}
	if (job->task_prolog) {
		_run_script("user task_prolog", job->task_prolog, job); 
	}

	if (job->env == NULL) {
		debug("job->env is NULL");
		job->env = (char **)xmalloc(sizeof(char *));
		job->env[0] = (char *)NULL;
	}

	log_fini();
	execve(task->argv[0], task->argv, job->env);

	/* 
	 * error() and clean up if execve() returns:
	 */
	error("execve(): %s: %m", job->argv[0]); 
	printf("execve failed, %d %s\n", job->argc, job->argv[0]);
	exit(errno);
}

static void
_make_tmpdir(slurmd_job_t *job)
{
	char *tmpdir;

	if (!(tmpdir = getenvp(job->env, "TMPDIR")))
		return;

	if ((mkdir(tmpdir, 0700) < 0) && (errno != EEXIST))
		error ("Unable to create TMPDIR [%s]: %m", tmpdir);

	return;
}
