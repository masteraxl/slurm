/*****************************************************************************\
 *  daemonize.c - daemonization routine
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
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
#  include <config.h>
#endif

#include <fcntl.h>
#include <unistd.h>

#include "src/common/daemonize.h"
#include "src/common/fd.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/xassert.h"

/* closeall FDs >= a specified value */
static void
closeall(int fd)
{
	int fdlimit = sysconf(_SC_OPEN_MAX);

	while (fd < fdlimit) 
		close(fd++);
}

/* detach and go into background.
 * caller is responsible for umasks
 *
 * if nochdir == 0, will do a chdir to /
 * if noclose == 0, will close all FDs
 */
int
daemon(int nochdir, int noclose)
{
	switch (fork()) {
		case  0 : break;        /* child */
		case -1 : return -1;
		default : _exit(0);     /* exit parent */
	}

	if (setsid() < 0)
		return -1;

	switch (fork()) {
		case 0 : break;         /* child */
		case -1: return -1;
		default: _exit(0);      /* exit parent */
	}

	if(!nochdir && chdir("/") < 0) {
		error("chdir(/): %m");
		return -1;
	}

	/* Close all file descriptors if requested
	 */
	if (!noclose) {
		closeall(0);
		open("/dev/null", O_RDWR);
		dup2(0, STDOUT_FILENO);
		dup2(0, STDERR_FILENO);
	} else {
		/* 
		 * Otherwise, dup stdin, stdout, and stderr onto /dev/null
		 */
		int devnull = open("/dev/null", O_RDWR);
		if (devnull < 0)
			error("Unable to open /dev/null: %m");
		if (dup2(devnull, STDIN_FILENO) < 0)
			error("Unable to dup /dev/null onto stdin: %m");
		if (dup2(devnull, STDOUT_FILENO) < 0)
			error("Unable to dup /dev/null onto stdout: %m");
		if (dup2(devnull, STDERR_FILENO) < 0)
			error("Unable to dup /dev/null onto stderr: %m");
		if (close(devnull) < 0)
			error("Unable to close /dev/null: %m");
	}

	return 0;

}

/*
 * Read and return pid stored in pidfile.
 * Returns 0 if file doesn't exist or pid cannot be read.
 * If pidfd != NULL, the file will be kept open and the fd
 * returned.
 */
pid_t
read_pidfile(const char *pidfile, int *pidfd)
{
	int fd;
	FILE *fp = NULL;
	unsigned long pid;
	pid_t         lpid;

	if ((fd = open(pidfile, O_RDONLY)) < 0) 
		return ((pid_t) 0);

	if (!(fp = fdopen(fd, "r")) && (errno != ENOENT)) 
		error ("Unable to access old pidfile at `%s': %m", pidfile);

	if (fscanf(fp, "%lu", &pid) < 1) {
		error ("Possible corrupt pidfile `%s'", pidfile);
		return ((pid_t) 0);
	}

	if ((lpid = fd_is_read_lock_blocked(fd)) == (pid_t) 0) {
		verbose ("pidfile not locked, assuming no running daemon");
		return (lpid);
	}

	if (lpid != (pid_t) pid) 
		fatal ("pidfile locked by %lu but contains pid=%lu",
		       (unsigned long) lpid, (unsigned long) pid);

	if (pidfd != NULL)
		*pidfd = fd;
	else 
		(void) close(fd); /* Ignore errors */

	return (lpid);
}



int
create_pidfile(const char *pidfile)
{
	FILE *fp;

	xassert(pidfile != NULL);
	xassert(pidfile[0] == '/');

	if (!(fp = fopen(pidfile, "w"))) {
		error("Unable to open pidfile `%s': %m", pidfile);
		return -1;
	}

	if (fd_get_write_lock(fileno(fp)) < 0) {
		error ("Unable to lock pidfile `%s': %m", pidfile);
		goto error;
	}

	if (fprintf(fp, "%lu\n", (unsigned long) getpid()) == EOF) {
		error("Unable to write to pidfile `%s': %m", pidfile);
		goto error;
	}

	fflush(fp);
	
	/*
	 * if (fclose(fp) == EOF) {
         *	error("Unable to close pidfile `%s': %m", pidfile);
         *	goto error;
         *}
	 */
	return (fileno(fp));

  error:
	if (unlink(pidfile) < 0)
		error("Unable to remove pidfile `%s': %m", pidfile);
	return -1;
}

