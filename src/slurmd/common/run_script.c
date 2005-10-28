/*****************************************************************************\
 * src/slurmd/common/run_script.c - code shared between slurmd and slurmstepd
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher Morrone <morrone2@llnl.gov>
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

#include <sys/wait.h>
#include <sys/errno.h>

#include "src/common/xmalloc.h"
#include "src/common/xassert.h"

#include "src/slurmd/common/run_script.h"


/*
 * Run a prolog or epilog script
 * name IN: class of program (prolog, epilog, etc.), 
 *	if prefix is "user" then also set uid
 * path IN: pathname of program to run
 * jobid, uid IN: info on associated job
 * max_wait IN: maximum time to wait in seconds, -1 for no limit
 * env IN: environment variables to use on exec, sets minimal environment
 *	if NULL
 * RET 0 on success, -1 on failure. 
 */
int
run_script(const char *name, const char *path, uint32_t jobid, uid_t uid, 
	   int max_wait, char **env)
{
	int status, rc, opt;
	pid_t cpid;

	xassert(env);
	if (path == NULL || path[0] == '\0')
		return 0;

	debug("[job %u] attempting to run %s [%s]", jobid, name, path);

	if (access(path, R_OK | X_OK) < 0) {
		debug("Not running %s [%s]: %m", name, path);
		return 0;
	}

	if ((cpid = fork()) < 0) {
		error ("executing %s: fork: %m", name);
		return -1;
	}
	if (cpid == 0) {
		char *argv[2];

		argv[0] = (char *)xstrdup(path);
		argv[1] = NULL;

		if (strncmp(name, "user", 4) == 0)
			setuid(uid);
		setpgrp();
		execve(path, argv, env);
		error("execve(): %m");
		exit(127);
	}

	if (max_wait < 0)
		opt = 0;
	else
		opt = WNOHANG;

	while (1) {
		rc = waitpid(cpid, &status, opt);
		if (rc < 0) {
			if (errno == EINTR)
				continue;
			error("waidpid: %m");
			return 0;
		} else if (rc == 0) {
			sleep(1);
			if ((--max_wait) == 0) {
				killpg(cpid, SIGKILL);
				opt = 0;
			}
		} else  {
			killpg(cpid, SIGKILL);	/* kill children too */
			return status;
		}
	}

	/* NOTREACHED */
}
