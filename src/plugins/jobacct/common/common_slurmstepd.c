/*****************************************************************************\
 *  jobacct_common.c - common functions for almost all jobacct plugins.
 *****************************************************************************
 *
 *  Copyright (C) 2005 Hewlett-Packard Development Company, L.P.
 *  Written by Danny Auble, <da@llnl.gov>
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
 *
 *  This file is patterned after jobcomp_linux.c, written by Morris Jette and
 *  Copyright (C) 2002 The Regents of the University of California.
\*****************************************************************************/

#include "jobacct_common.h"

bool jobacct_shutdown = false;
bool suspended = false;
List task_list = NULL;
pthread_mutex_t jobacct_lock = PTHREAD_MUTEX_INITIALIZER;

extern int common_endpoll()
{
	jobacct_shutdown = true;
       
	return SLURM_SUCCESS;
}
extern int common_add_task(pid_t pid, jobacct_id_t *jobacct_id)
{
	struct jobacctinfo *jobacct = common_alloc_jobacct(jobacct_id);
	
	slurm_mutex_lock(&jobacct_lock);
	if(pid <= 0) {
		error("invalid pid given (%d) for task acct", pid);
		goto error;
	} else if (!task_list) {
		error("no task list created!");
		goto error;
	}

	jobacct->pid = pid;
	jobacct->min_cpu = 0;
	debug2("adding task %u pid %d on node %uto jobacct", 
	       jobacct_id->taskid, pid, jobacct_id->nodeid);
	list_push(task_list, jobacct);
	slurm_mutex_unlock(&jobacct_lock);

	return SLURM_SUCCESS;
error:
	slurm_mutex_unlock(&jobacct_lock);
	common_free_jobacct(jobacct);
	return SLURM_ERROR;
}

extern struct jobacctinfo *common_stat_task(pid_t pid)
{
	struct jobacctinfo *jobacct = NULL;
	struct jobacctinfo *ret_jobacct = NULL;
	ListIterator itr = NULL;
	
	slurm_mutex_lock(&jobacct_lock);
	if (!task_list) {
		error("no task list created!");
		goto error;
	}

	itr = list_iterator_create(task_list);
	while((jobacct = list_next(itr))) { 
		if(jobacct->pid == pid)
			break;
	}
	list_iterator_destroy(itr);
	ret_jobacct = xmalloc(sizeof(struct jobacctinfo));
	memcpy(ret_jobacct, jobacct, sizeof(struct jobacctinfo));
error:
	slurm_mutex_unlock(&jobacct_lock);
	return ret_jobacct;
}

extern struct jobacctinfo *common_remove_task(pid_t pid)
{
	struct jobacctinfo *jobacct = NULL;
	struct jobacctinfo *ret_jobacct = NULL;
	ListIterator itr = NULL;

	slurm_mutex_lock(&jobacct_lock);
	if (!task_list) {
		error("no task list created!");
		goto error;
	}

	itr = list_iterator_create(task_list);
	while((jobacct = list_next(itr))) { 
		if(jobacct->pid == pid) {
			list_remove(itr);
			break;
		}
	}
	list_iterator_destroy(itr);
	if(jobacct) {
		debug2("removing task %u pid %d from jobacct", 
		       jobacct->max_vsize_id.taskid, jobacct->pid);
		ret_jobacct = xmalloc(sizeof(struct jobacctinfo));
		memcpy(ret_jobacct, jobacct, sizeof(struct jobacctinfo));
		common_free_jobacct(jobacct);
	} else {
		error("pid(%d) not being watched in jobacct!", pid);
	}
error:
	slurm_mutex_unlock(&jobacct_lock);
	return ret_jobacct;
}

extern void common_suspendpoll()
{
	if(suspended)
		suspended = false;
	else
		suspended = true;
}
