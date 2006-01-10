/*****************************************************************************\
 *  kill_tree.c - Kill process tree based upon process IDs
 *	Used primarily for MPICH-GM
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Takao Hatazaki <takao.hatazaki@hp.com>
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
#include <dirent.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <strings.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>

#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/log.h"

typedef struct xpid_s {
	pid_t pid;
	int is_usercmd;
	char *cmd;
	struct xpid_s *next;
} xpid_t;

typedef struct xppid_s {
	pid_t ppid;
	xpid_t *list;
	struct xppid_s *next;
} xppid_t;

#define HASH_LEN 64

#define GET_HASH_IDX(ppid) ((ppid)%HASH_LEN)

static xpid_t *_alloc_pid(pid_t pid, int is_usercmd, char *cmd, xpid_t *next)
{
	xpid_t *new;

	new = (xpid_t *)xmalloc(sizeof(*new));
	new->pid = pid;
	new->is_usercmd = is_usercmd;
	new->cmd = xstrdup(cmd);
	new->next = next;
	return new;
}

static xppid_t *_alloc_ppid(pid_t ppid, pid_t pid, int is_usercmd, char *cmd,
			    xppid_t *next)
{
	xppid_t *new;

	new = xmalloc(sizeof(*new));
	new->ppid = ppid;
	new->list = _alloc_pid(pid, is_usercmd, cmd, NULL);
	new->next = next;
	return new;
}

static void _push_to_hashtbl(pid_t ppid, pid_t pid,
			     int is_usercmd, char *cmd, xppid_t **hashtbl)
{
	int idx;
	xppid_t *ppids, *newppid;
	xpid_t *newpid;

	idx = GET_HASH_IDX(ppid);
	ppids = hashtbl[idx];
	while (ppids) {
		if (ppids->ppid == ppid) {
			newpid = _alloc_pid(pid, is_usercmd, cmd, ppids->list);
			ppids->list = newpid;
			return;
		}
		ppids = ppids->next;
	}
	newppid = _alloc_ppid(ppid, pid, is_usercmd, cmd, hashtbl[idx]);
	hashtbl[idx] = newppid;    
}

static int get_myname(char *s)
{
	char path[PATH_MAX], rbuf[1024];
	int fd;

	sprintf(path, "/proc/%ld/stat", (long)getpid());
	if ((fd = open(path, O_RDONLY)) < 0) {
		error("Cannot open /proc/getpid()/stat");
		return -1;
	}
	if (read(fd, rbuf, 1024) <= 0) {
		error("Cannot read /proc/getpid()/stat");
		close(fd);
		return -1;
	}
	close(fd);
	if (sscanf(rbuf, "%*ld %s ", s) != 1) {
		error("Cannot get the command name from /proc/getpid()/stat");
		return -1;
	}
	return 0;
}

static xppid_t **_build_hashtbl()
{
	DIR *dir;
	struct dirent *de;
	char path[PATH_MAX], *endptr, *num, rbuf[1024];
	char myname[1024], cmd[1024];
	int fd;
	long pid, ppid;
	xppid_t **hashtbl;

	if ((dir = opendir("/proc")) == NULL) {
		error("opendir(/proc): %m");
		return NULL;
	}
	if (get_myname(myname) < 0) return NULL;
	debug3("Myname in build_hashtbl: %s", myname);

	hashtbl = (xppid_t **)xmalloc(HASH_LEN * sizeof(xppid_t *));

	while ((de = readdir(dir)) != NULL) {
		num = de->d_name;
		strtol(num, &endptr, 10);
		if (endptr == NULL || *endptr != 0)
			continue;
		sprintf(path, "/proc/%s/stat", num);
		if ((fd = open(path, O_RDONLY)) < 0) {
			continue;
		}
		if (read(fd, rbuf, 1024) <= 0) {
			close(fd);
			continue;
		}
		if (sscanf(rbuf, "%ld %s %*s %ld", &pid, cmd, &ppid) != 3) {
			close(fd);
			continue;
		}
		close(fd);

		/* Record cmd for debugging purpose */
		_push_to_hashtbl((pid_t)ppid, (pid_t)pid, 
				 strcmp(myname, cmd), cmd, hashtbl);
	}
	closedir(dir);
	return hashtbl;
}

static void _destroy_list(xpid_t *list)
{
	xpid_t *tmp;

	while (list) {
		tmp = list->next;
		xfree(list->cmd);
		xfree(list);
		list = tmp;
	}
}

static void _destroy_hashtbl(xppid_t **hashtbl)
{
	int i;
	xppid_t *ppid, *tmp;

	for (i=0; i<HASH_LEN; i++) {
		ppid = hashtbl[i];
		while (ppid) {
			_destroy_list(ppid->list);
			tmp = ppid->next;
			xfree(ppid);
			ppid = tmp;
		}
	}
	xfree(hashtbl);
}


static xpid_t *_get_list(int top, xpid_t *list, xppid_t **hashtbl)
{
	xppid_t *ppid;
	xpid_t *children;

	ppid = hashtbl[GET_HASH_IDX(top)];
	while (ppid) {
		if (ppid->ppid == top) {
			children = ppid->list;
			while (children) {
				list = _alloc_pid(children->pid,
						  children->is_usercmd,
						  children->cmd,
						  list);
				children = children->next;
			}
			children = ppid->list;
			while (children) {
				list = _get_list(children->pid, list, hashtbl);
				children = children->next;
			}
			break;
		}
		ppid = ppid->next;
	}
	return list;
}

static int _kill_proclist(xpid_t *list, int sig)
{
	int rc, rc0;

	rc = 0;
	while (list) {
		if (list->pid > 1) {
			if (! list->is_usercmd) {
				debug2("%ld %s is not a user command.  "
				       "Skipped sending signal %d",
				       (long)list->pid, list->cmd, sig);
			} else {
				verbose("Sending %d to %d %s",
					sig, list->pid, list->cmd);
				rc0 = kill(list->pid, sig);
				if (rc0) rc = errno; /* save the last error */
			}
		}
		list = list->next;
	}

	return rc;
}


/*
 * Some of processes may not be in the same process group
 * (e.g. GMPI processes).  So, find out the process tree,
 * then kill all that subtree.
 */
extern int kill_proc_tree(pid_t top, int sig)
{
	xpid_t *list;
	int rc = -1;
	xppid_t **hashtbl;

	if ((hashtbl = _build_hashtbl()) == NULL)
		return -1;
	
	list = _get_list(top, NULL, hashtbl);
	rc = _kill_proclist(list, sig);
	_destroy_hashtbl(hashtbl);
	_destroy_list(list);
	return rc;
}


/*
 * Return the pid of the process named "process_name" 
 * which is the ancestor of "process".
 */
extern pid_t find_ancestor(pid_t process, char *process_name)
{
	char path[PATH_MAX], rbuf[1024];
	int fd;
	long pid, ppid;

	pid = ppid = (long)process;
	do {
		if (ppid <= 1) {
			return 0;
		}

		sprintf(path, "/proc/%d/stat", ppid);
		if ((fd = open(path, O_RDONLY)) < 0) {
			return 0;
		}
		if (read(fd, rbuf, 1024) <= 0) {
			close(fd);
			return 0;
		}
		close(fd);
		if (sscanf(rbuf, "%ld %*s %*s %ld", &pid, &ppid) != 2) {
			return 0;
		}

		sprintf(path, "/proc/%d/cmdline", pid);
		if ((fd = open(path, O_RDONLY)) < 0) {
			continue;
		}
		if (read(fd, rbuf, 1024) <= 0) {
			close(fd);
			continue;
		}
		close(fd);
	} while (!strstr(rbuf, process_name));

	return pid;
}
