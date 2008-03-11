/*****************************************************************************\
 * job_scheduler.c - manage the scheduling of pending jobs in priority order
 *	Note there is a global job list (job_list)
 *****************************************************************************
 *  Copyright (C) 2002-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/common/node_select.h"
#include "src/common/slurm_accounting_storage.h"
#include "src/common/xassert.h"
#include "src/common/xstring.h"

#include "src/slurmctld/agent.h"
#include "src/slurmctld/job_scheduler.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/node_scheduler.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/srun_comm.h"
#include "src/slurmctld/assoc_mgr.h"

#define _DEBUG 0
#define MAX_RETRIES 10

static void _depend_list_del(void *dep_ptr);
static char **_xduparray(uint16_t size, char ** array);

/* 
 * build_job_queue - build (non-priority ordered) list of pending jobs
 * OUT job_queue - pointer to job queue
 * RET number of entries in job_queue
 * NOTE: the buffer at *job_queue must be xfreed by the caller
 */
extern int build_job_queue(struct job_queue **job_queue)
{
	ListIterator job_iterator;
	struct job_record *job_ptr = NULL;
	int job_buffer_size, job_queue_size;
	struct job_queue *my_job_queue;

	/* build list pending jobs */
	job_buffer_size = job_queue_size = 0;
	job_queue[0] = my_job_queue = NULL;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		xassert (job_ptr->magic == JOB_MAGIC);
		if ((job_ptr->job_state != JOB_PENDING)    ||
		    (job_ptr->job_state &  JOB_COMPLETING) ||
		    (job_ptr->priority == 0))	/* held */
			continue;
		if (!job_independent(job_ptr))	/* waiting for other job */
			continue;
		if (job_buffer_size <= job_queue_size) {
			job_buffer_size += 200;
			xrealloc(my_job_queue, job_buffer_size *
				 sizeof(struct job_queue));
		}
		my_job_queue[job_queue_size].job_ptr  = job_ptr;
		my_job_queue[job_queue_size].job_priority = job_ptr->priority;
		my_job_queue[job_queue_size].part_priority = 
						job_ptr->part_ptr->priority;
		job_queue_size++;
	}
	list_iterator_destroy(job_iterator);

	job_queue[0] = my_job_queue;
	return job_queue_size;
}

/*
 * job_is_completing - Determine if jobs are in the process of completing.
 * RET - True of any job is in the process of completing
 * NOTE: This function can reduce resource fragmentation, which is a 
 * critical issue on Elan interconnect based systems.
 */
extern bool job_is_completing(void)
{
	bool completing = false;
	ListIterator job_iterator;
	struct job_record *job_ptr = NULL;
	time_t recent = time(NULL) - (slurmctld_conf.kill_wait + 2);

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		if ((job_ptr->job_state & JOB_COMPLETING) &&
		    (job_ptr->end_time >= recent)) {
			completing = true;
			break;
		}
	}
	list_iterator_destroy(job_iterator);

	return completing;
}

/*
 * set_job_elig_time - set the eligible time for pending jobs once their
 *      dependencies are lifted (in job->details->begin_time)
 */
extern void set_job_elig_time(void)
{
	struct job_record *job_ptr = NULL;
	struct part_record *part_ptr = NULL;
	ListIterator job_iterator;
	slurmctld_lock_t job_write_lock =
		{ READ_LOCK, WRITE_LOCK, WRITE_LOCK, READ_LOCK };
	time_t now = time(NULL);

	lock_slurmctld(job_write_lock);
	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		part_ptr = job_ptr->part_ptr;
		if (job_ptr->job_state != JOB_PENDING)
			continue;
		if (part_ptr == NULL)
			continue;
		if ((job_ptr->details == NULL) || job_ptr->details->begin_time)
			continue;
		if (part_ptr->state_up == 0)
			continue;
		if ((job_ptr->time_limit != NO_VAL) &&
		    (job_ptr->time_limit > part_ptr->max_time))
			continue;
		if ((job_ptr->details->max_nodes != 0) &&
		    ((job_ptr->details->max_nodes < part_ptr->min_nodes) ||
		     (job_ptr->details->min_nodes > part_ptr->max_nodes)))
			continue;
		if (!job_independent(job_ptr))
			continue;
		job_ptr->details->begin_time = now;
	}
	list_iterator_destroy(job_iterator);
	unlock_slurmctld(job_write_lock);
}

/* Test of part_ptr can still run jobs or if its nodes have
 * already been reserved by higher priority jobs (those in
 * the failed_parts array) */
static bool _failed_partition(struct part_record *part_ptr,
			      struct part_record **failed_parts, 
			      int failed_part_cnt)
{
	int i;

	for (i = 0; i < failed_part_cnt; i++) {
		if (failed_parts[i] == part_ptr)
			return true;
	}
	return false;
}

#ifndef HAVE_BG
/* Add a partition to the failed_parts array, reserving its nodes
 * from use by lower priority jobs. Also flags all partitions with
 * nodes overlapping this partition. */
static void _add_failed_partition(struct part_record *failed_part_ptr, 
			      struct part_record **failed_parts, 
			      int *failed_part_cnt)
{
	int count = *failed_part_cnt;
	ListIterator part_iterator;
	struct part_record *part_ptr;

	failed_parts[count++] = failed_part_ptr;

	/* We also need to add partitions that have overlapping nodes */
	part_iterator = list_iterator_create(part_list);
	while ((part_ptr = (struct part_record *) list_next(part_iterator))) {
		if ((part_ptr == failed_part_ptr) ||
		    (_failed_partition(part_ptr, failed_parts, count)) ||
		    (!bit_overlap(part_ptr->node_bitmap, 
				  failed_part_ptr->node_bitmap)))
			continue;
		failed_parts[count++] = part_ptr;
	}
	list_iterator_destroy(part_iterator);

	*failed_part_cnt = count;
}
#endif

/* 
 * schedule - attempt to schedule all pending jobs
 *	pending jobs for each partition will be scheduled in priority  
 *	order until a request fails
 * RET count of jobs scheduled
 * Note: We re-build the queue every time. Jobs can not only be added 
 *	or removed from the queue, but have their priority or partition 
 *	changed with the update_job RPC. In general nodes will be in priority 
 *	order (by submit time), so the sorting should be pretty fast.
 */
extern int schedule(void)
{
	struct job_queue *job_queue;
	int i, error_code, failed_part_cnt = 0, job_queue_size, job_cnt = 0;
	struct job_record *job_ptr;
	struct part_record **failed_parts = NULL;
	/* Locks: Read config, write job, write node, read partition */
	slurmctld_lock_t job_write_lock =
	    { READ_LOCK, WRITE_LOCK, WRITE_LOCK, READ_LOCK };
#ifdef HAVE_BG
	char *ionodes = NULL;
	char tmp_char[256];
#endif
	static bool wiki_sched = false;
	static bool wiki_sched_test = false;
	DEF_TIMERS;

	START_TIMER;
	/* don't bother trying to avoid fragmentation with sched/wiki */
	if (!wiki_sched_test) {
		char *sched_type = slurm_get_sched_type();
		if ((strcmp(sched_type, "sched/wiki") == 0)
		||  (strcmp(sched_type, "sched/wiki2") == 0))
			wiki_sched = true;
		xfree(sched_type);
		wiki_sched_test = true;
	}

	lock_slurmctld(job_write_lock);
	/* Avoid resource fragmentation if important */
	if ((!wiki_sched) && switch_no_frag() && job_is_completing()) {
		unlock_slurmctld(job_write_lock);
		debug("schedule() returning, some job still completing");
		return SLURM_SUCCESS;
	}
	debug("Running job scheduler");
	job_queue_size = build_job_queue(&job_queue);
	if (job_queue_size == 0) {
		unlock_slurmctld(job_write_lock);
		return SLURM_SUCCESS;
	}
	sort_job_queue(job_queue, job_queue_size);

	failed_parts = xmalloc(sizeof(struct part_record *) * 
			       list_count(part_list));

	for (i = 0; i < job_queue_size; i++) {
		job_ptr = job_queue[i].job_ptr;
		if (job_ptr->priority == 0)	/* held */
			continue;
		if (_failed_partition(job_ptr->part_ptr, failed_parts, 
				      failed_part_cnt)) {
			continue;
		}
		if (validate_assoc_id(acct_db_conn, job_ptr->assoc_id)) {
			/* NOTE: This only happens if a user's account is 
			 * disabled between when the job was submitted and 
			 * the time we consider running it. It should be 
			 * very rare. */
			info("schedule: JobId=%u has invalid account",
				job_ptr->job_id);
			last_job_update = time(NULL);
			job_ptr->job_state = JOB_FAILED;
			job_ptr->exit_code = 1;
			job_ptr->state_reason = FAIL_BANK_ACCOUNT;
			job_ptr->start_time = job_ptr->end_time = time(NULL);
			job_completion_logger(job_ptr);
			delete_job_details(job_ptr);
			continue;
		}

		error_code = select_nodes(job_ptr, false, NULL);
		if (error_code == ESLURM_NODES_BUSY) {
#ifndef HAVE_BG 	/* keep trying to schedule jobs in partition */
			/* While we use static partitiioning on Blue Gene, 
			 * each job can be scheduled independently without 
			 * impacting other jobs with different characteristics
			 * (e.g. node-use [virtual or coprocessor] or conn-type
			 * [mesh, torus, or nav]). Because of this we sort and 
			 * then try to schedule every pending job. This does 
			 * increase the overhead of this job scheduling cycle, 
			 * but the only way to effectively avoid this is to 
			 * define each SLURM partition as containing a 
			 * single Blue Gene job partition type (e.g. 
			 * group all Blue Gene job partitions of type 
			 * 2x2x2 coprocessor mesh into a single SLURM
			 * partition, say "co-mesh-222") */
			_add_failed_partition(job_ptr->part_ptr, failed_parts,
				      &failed_part_cnt);
#endif
		} else if (error_code == SLURM_SUCCESS) {	
			/* job initiated */
			last_job_update = time(NULL);
#ifdef HAVE_BG
			select_g_get_jobinfo(job_ptr->select_jobinfo, 
					     SELECT_DATA_IONODES, 
					     &ionodes);
			if(ionodes) {
				sprintf(tmp_char,"%s[%s]",
						job_ptr->nodes,
						ionodes);
			} else {
				sprintf(tmp_char,"%s",job_ptr->nodes);
			}
			info("schedule: JobId=%u BPList=%s",
			     job_ptr->job_id, tmp_char);
			xfree(ionodes);
#else
			info("schedule: JobId=%u NodeList=%s",
			     job_ptr->job_id, job_ptr->nodes);
#endif
			if (job_ptr->batch_flag)
				launch_job(job_ptr);
			else
				srun_allocate(job_ptr->job_id);
			job_cnt++;
		} else if (error_code !=
		           ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE) {
			info("schedule: JobId=%u non-runnable: %s",
				job_ptr->job_id, 
				slurm_strerror(error_code));
			last_job_update = time(NULL);
			job_ptr->job_state = JOB_FAILED;
			job_ptr->exit_code = 1;
			job_ptr->state_reason = FAIL_BAD_CONSTRAINTS;
			job_ptr->start_time = job_ptr->end_time = time(NULL);
			job_completion_logger(job_ptr);
			delete_job_details(job_ptr);
		}
	}

	xfree(failed_parts);
	xfree(job_queue);
	unlock_slurmctld(job_write_lock);
	END_TIMER2("schedule");
	return job_cnt;
}


/* 
 * sort_job_queue - sort job_queue in decending priority order
 * IN job_queue_size - count of elements in the job queue
 * IN/OUT job_queue - pointer to sorted job queue
 */
extern void sort_job_queue(struct job_queue *job_queue, int job_queue_size)
{
	int i, j, top_prio_inx;
	struct job_record *tmp_job_ptr;
	uint32_t top_job_prio,  tmp_job_prio;
	uint16_t top_part_prio, tmp_part_prio;

	for (i = 0; i < job_queue_size; i++) {
		top_prio_inx  = i;
		top_job_prio  = job_queue[i].job_priority;
		top_part_prio = job_queue[i].part_priority;

		for (j = (i + 1); j < job_queue_size; j++) {
			if (top_part_prio > job_queue[j].part_priority)
				continue;
			if ((top_part_prio == job_queue[j].part_priority) &&
			    (top_job_prio  >= job_queue[j].job_priority))
				continue;

			top_prio_inx  = j;
			top_job_prio  = job_queue[j].job_priority;
			top_part_prio = job_queue[j].part_priority;
		}
		if (top_prio_inx == i)
			continue;	/* in correct order */

		/* swap records at top_prio_inx and i */
		tmp_job_ptr   = job_queue[i].job_ptr;
		tmp_job_prio  = job_queue[i].job_priority;
		tmp_part_prio = job_queue[i].part_priority;

		job_queue[i].job_ptr       = job_queue[top_prio_inx].job_ptr;
		job_queue[i].job_priority  = job_queue[top_prio_inx].job_priority;
		job_queue[i].part_priority = job_queue[top_prio_inx].part_priority;

		job_queue[top_prio_inx].job_ptr       = tmp_job_ptr;
		job_queue[top_prio_inx].job_priority  = tmp_job_prio;
		job_queue[top_prio_inx].part_priority = tmp_part_prio;

	}
}

/*
 * launch_job - send an RPC to a slurmd to initiate a batch job 
 * IN job_ptr - pointer to job that will be initiated
 */
extern void launch_job(struct job_record *job_ptr)
{
	batch_job_launch_msg_t *launch_msg_ptr;
	agent_arg_t *agent_arg_ptr;
	struct node_record *node_ptr;

	node_ptr = find_first_node_record(job_ptr->node_bitmap);
	if (node_ptr == NULL)
		return;

	/* Initialization of data structures */
	launch_msg_ptr = (batch_job_launch_msg_t *)
				xmalloc(sizeof(batch_job_launch_msg_t));
	launch_msg_ptr->job_id = job_ptr->job_id;
	launch_msg_ptr->step_id = NO_VAL;
	launch_msg_ptr->uid = job_ptr->user_id;
	launch_msg_ptr->gid = job_ptr->group_id;
	launch_msg_ptr->nprocs = job_ptr->details->num_tasks;
	launch_msg_ptr->nodes = xstrdup(job_ptr->nodes);
	launch_msg_ptr->overcommit = job_ptr->details->overcommit;
	launch_msg_ptr->open_mode  = job_ptr->details->open_mode;
	launch_msg_ptr->acctg_freq = job_ptr->details->acctg_freq;

	if (make_batch_job_cred(launch_msg_ptr)) {
		error("aborting batch job %u", job_ptr->job_id);
		/* FIXME: This is a kludge, but this event indicates a serious
		 * problem with OpenSSH and should never happen. We are
		 * too deep into the job launch to gracefully clean up. */
		job_ptr->end_time    = time(NULL);
		job_ptr->time_limit = 0;
		xfree(launch_msg_ptr->nodes);
		xfree(launch_msg_ptr);
		return;
	}

	launch_msg_ptr->err = xstrdup(job_ptr->details->err);
	launch_msg_ptr->in = xstrdup(job_ptr->details->in);
	launch_msg_ptr->out = xstrdup(job_ptr->details->out);
	launch_msg_ptr->work_dir = xstrdup(job_ptr->details->work_dir);
	launch_msg_ptr->argc = job_ptr->details->argc;
	launch_msg_ptr->argv = _xduparray(job_ptr->details->argc,
					job_ptr->details->argv);
	launch_msg_ptr->script = get_job_script(job_ptr);
	launch_msg_ptr->environment =
	    get_job_env(job_ptr, &launch_msg_ptr->envc);

	launch_msg_ptr->num_cpu_groups = job_ptr->num_cpu_groups;
	launch_msg_ptr->cpus_per_node  = xmalloc(sizeof(uint32_t) *
			job_ptr->num_cpu_groups);
	memcpy(launch_msg_ptr->cpus_per_node, job_ptr->cpus_per_node,
			(sizeof(uint32_t) * job_ptr->num_cpu_groups));
	launch_msg_ptr->cpu_count_reps  = xmalloc(sizeof(uint32_t) *
			job_ptr->num_cpu_groups);
	memcpy(launch_msg_ptr->cpu_count_reps, job_ptr->cpu_count_reps,
			(sizeof(uint32_t) * job_ptr->num_cpu_groups));

	launch_msg_ptr->select_jobinfo = select_g_copy_jobinfo(
			job_ptr->select_jobinfo);

	agent_arg_ptr = (agent_arg_t *) xmalloc(sizeof(agent_arg_t));
	agent_arg_ptr->node_count = 1;
	agent_arg_ptr->retry = 0;
	agent_arg_ptr->hostlist = hostlist_create(node_ptr->name);
	agent_arg_ptr->msg_type = REQUEST_BATCH_JOB_LAUNCH;
	agent_arg_ptr->msg_args = (void *) launch_msg_ptr;

	/* Launch the RPC via agent */
	agent_queue_request(agent_arg_ptr);
}

static char **
_xduparray(uint16_t size, char ** array)
{
	int i;
	char ** result;

	if (size == 0)
		return (char **)NULL;

	result = (char **) xmalloc(sizeof(char *) * size);
	for (i=0; i<size; i++)
		result[i] = xstrdup(array[i]);
	return result;
}

/*
 * make_batch_job_cred - add a job credential to the batch_job_launch_msg
 * IN/OUT launch_msg_ptr - batch_job_launch_msg in which job_id, step_id, 
 *                         uid and nodes have already been set 
 * RET 0 or error code
 */
extern int make_batch_job_cred(batch_job_launch_msg_t *launch_msg_ptr)
{
	slurm_cred_arg_t cred_arg;

	cred_arg.jobid     = launch_msg_ptr->job_id;
	cred_arg.stepid    = launch_msg_ptr->step_id;
	cred_arg.uid       = launch_msg_ptr->uid;
	cred_arg.hostlist  = launch_msg_ptr->nodes;
	cred_arg.alloc_lps_cnt = 0;
	cred_arg.alloc_lps = NULL;

	launch_msg_ptr->cred = slurm_cred_create(slurmctld_config.cred_ctx,
			 &cred_arg);

	if (launch_msg_ptr->cred)
		return SLURM_SUCCESS;
	error("slurm_cred_create failure for batch job %u", cred_arg.jobid);
	return SLURM_ERROR;
}

static void _depend_list_del(void *dep_ptr)
{
	xfree(dep_ptr);
}

/* Print a job's dependency information based upon job_ptr->depend_list */
extern void print_job_dependency(struct job_record *job_ptr)
{
	ListIterator depend_iter;
	struct depend_spec *dep_ptr;
	char *dep_str;

	info("Dependency information for job %u", job_ptr->job_id);
	if ((job_ptr->details == NULL) ||
	    (job_ptr->details->depend_list == NULL))
		return;

	depend_iter = list_iterator_create(job_ptr->details->depend_list);
	if (!depend_iter)
		fatal("list_iterator_create memory allocation failure");
	while ((dep_ptr = list_next(depend_iter))) {
		if      (dep_ptr->depend_type == SLURM_DEPEND_AFTER)
			dep_str = "after";
		else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER_ANY)
			dep_str = "afterany";
		else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER_NOT_OK)
			dep_str = "afternotok";
		else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER_OK)
			dep_str = "afterok";
		else
			dep_str = "unknown";
		info("  %s:%u", dep_str, dep_ptr->job_id);
	}
	list_iterator_destroy(depend_iter);
}

/*
 * Determine if a job's dependencies are met
 * RET: 0 = no dependencies
 *      1 = dependencies remain
 *      2 = failure (job completion code not per dependency), delete the job
 */
extern int test_job_dependency(struct job_record *job_ptr)
{
	ListIterator depend_iter;
	struct depend_spec *dep_ptr;
	bool failure = false;

	if ((job_ptr->details == NULL) ||
	    (job_ptr->details->depend_list == NULL))
		return 0;

	depend_iter = list_iterator_create(job_ptr->details->depend_list);
	if (!depend_iter)
		fatal("list_iterator_create memory allocation failure");
	while ((dep_ptr = list_next(depend_iter))) {
		if (dep_ptr->job_ptr->job_id != dep_ptr->job_id) {
			/* job is gone, dependency lifted */
			list_delete_item(depend_iter);
		} else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER) {
			if (!IS_JOB_PENDING(dep_ptr->job_ptr))
				list_delete_item(depend_iter);
			else
				break;
		} else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER_ANY) {
			if (IS_JOB_FINISHED(dep_ptr->job_ptr))
				list_delete_item(depend_iter);
			else
				break;
		} else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER_NOT_OK) {
			if (!IS_JOB_FINISHED(dep_ptr->job_ptr))
				break;
			if ((dep_ptr->job_ptr->job_state & (~JOB_COMPLETING))
			    != JOB_COMPLETE)
				list_delete_item(depend_iter);
			else {
				failure = true;
				break;
			}
		} else if (dep_ptr->depend_type == SLURM_DEPEND_AFTER_OK) {
			if (!IS_JOB_FINISHED(dep_ptr->job_ptr))
				break;
			if ((dep_ptr->job_ptr->job_state & (~JOB_COMPLETING))
			    == JOB_COMPLETE)
				list_delete_item(depend_iter);
			else {
				failure = true;
				break;
			}
		} else
			failure = true;
	}
	list_iterator_destroy(depend_iter);

	if (failure)
		return 2;
	if (dep_ptr)
		return 1;
	return 0;
}

/*
 * Parse a job dependency string and use it to establish a "depend_spec" 
 * list of dependencies. We accept both old format (a single job ID) and
 * new format (e.g. "afterok:123:124,after:128").
 * IN job_ptr - job record to have dependency and depend_list updated
 * IN new_depend - new dependency description
 * RET returns an error code from slurm_errno.h
 */
extern int update_job_dependency(struct job_record *job_ptr, char *new_depend)
{
	int rc = SLURM_SUCCESS;
	uint16_t depend_type = 0;
	uint32_t job_id = 0;
	char *tok = new_depend, *sep_ptr, *sep_ptr2;
	List new_depend_list = NULL;
	struct depend_spec *dep_ptr;
	struct job_record *dep_job_ptr;
	char dep_buf[32];

	if (job_ptr->details == NULL)
		return EINVAL;

	/* Clear dependencies on NULL or empty dependency input */
	if ((new_depend == NULL) || (new_depend[0] == '\0')) {
		xfree(job_ptr->details->dependency);
		if (job_ptr->details->depend_list)
			list_destroy(job_ptr->details->depend_list);
		return rc;

	}

	new_depend_list = list_create(_depend_list_del);
	/* validate new dependency string */
	while (rc == SLURM_SUCCESS) {
		sep_ptr = strchr(tok, ':');
		if ((sep_ptr == NULL) && (job_id == 0)) {
			job_id = strtol(tok, &sep_ptr, 10);
			if ((sep_ptr == NULL) || (sep_ptr[0] != '\0') ||
			    (job_id <= 0) || (job_id == job_ptr->job_id)) {
				rc = EINVAL;
				break;
			}
			/* old format, just a single job_id */
			dep_job_ptr = find_job_record(job_id);
			if (!dep_job_ptr)	/* assume already done */
				break;
			snprintf(dep_buf, sizeof(dep_buf), "afterany:%u", job_id);
			new_depend = dep_buf;
			dep_ptr = xmalloc(sizeof(struct depend_spec));
			dep_ptr->depend_type = SLURM_DEPEND_AFTER_ANY;
			dep_ptr->job_id = job_id;
			dep_ptr->job_ptr = dep_job_ptr;
			if (!list_append(new_depend_list, dep_ptr))
				fatal("list_append memory allocation failure");
			break;
		}

		if      (strncasecmp(tok, "afternotok", 10) == 0)
			depend_type = SLURM_DEPEND_AFTER_NOT_OK;
		else if (strncasecmp(tok, "afterany", 8) == 0)
			depend_type = SLURM_DEPEND_AFTER_ANY;
		else if (strncasecmp(tok, "afterok", 7) == 0)
			depend_type = SLURM_DEPEND_AFTER_OK;
		else if (strncasecmp(tok, "after", 5) == 0)
			depend_type = SLURM_DEPEND_AFTER;
		else {
			rc = EINVAL;
			break;
		}
		sep_ptr++;	/* skip over ":" */
		while (rc == SLURM_SUCCESS) {
			job_id = strtol(sep_ptr, &sep_ptr2, 10);
			if ((sep_ptr2 == NULL) || 
			    (job_id <= 0) || (job_id == job_ptr->job_id) ||
			    ((sep_ptr2[0] != '\0') && (sep_ptr2[0] != ',') && 
			     (sep_ptr2[0] != ':'))) {
				rc = EINVAL;
				break;
			}
			dep_job_ptr = find_job_record(job_id);
			if (dep_job_ptr) {	/* job still active */
				dep_ptr = xmalloc(sizeof(struct depend_spec));
				dep_ptr->depend_type = depend_type;
				dep_ptr->job_id = job_id;
				dep_ptr->job_ptr = dep_job_ptr;
				if (!list_append(new_depend_list, dep_ptr)) {
					fatal("list_append memory allocation "
						"failure");
				}
			}
			if (sep_ptr2[0] != ':')
				break;
			sep_ptr = sep_ptr2 + 1;	/* skip over ":" */
		}
		if (sep_ptr2[0] == ',')
			tok = sep_ptr2 + 1;
		else
			break;
	}

	if (rc == SLURM_SUCCESS) {
		xfree(job_ptr->details->dependency);
		job_ptr->details->dependency = xstrdup(new_depend);
		if (job_ptr->details->depend_list)
			list_destroy(job_ptr->details->depend_list);
		job_ptr->details->depend_list = new_depend_list;
#if _DEBUG
		print_job_dependency(job_ptr);
#endif
	} else {
		list_destroy(new_depend_list);
	}
	return rc;
}

/* Determine if a pending job will run using only the specified nodes
 * (in job_desc_msg->req_nodes), build response message and return 
 * SLURM_SUCCESS on success. Otherwise return an error code. Caller 
 * must free response message */
extern int job_start_data(job_desc_msg_t *job_desc_msg, 
			  will_run_response_msg_t **resp)
{
	struct job_record *job_ptr;
	struct part_record *part_ptr;
	bitstr_t *avail_bitmap = NULL;
	uint32_t min_nodes, max_nodes, req_nodes;
	int rc = SLURM_SUCCESS;

	job_ptr = find_job_record(job_desc_msg->job_id);
	if (job_ptr == NULL)
		return ESLURM_INVALID_JOB_ID;

	part_ptr = job_ptr->part_ptr;
	if (part_ptr == NULL)
		return ESLURM_INVALID_PARTITION_NAME;

	if ((job_ptr->details == NULL) ||
	    (job_ptr->job_state != JOB_PENDING))
		return ESLURM_DISABLED;

	if ((job_desc_msg->req_nodes == NULL) || 
	    (job_desc_msg->req_nodes == '\0')) {
		/* assume all nodes available to job for testing */
		avail_bitmap = bit_copy(avail_node_bitmap);
	} else if (node_name2bitmap(job_desc_msg->req_nodes, false, 
				    &avail_bitmap) != 0) {
		return ESLURM_INVALID_NODE_NAME;
	} else {
		/* Only consider nodes that are not DOWN or DRAINED */
		bit_and(avail_bitmap, avail_node_bitmap);
	}

	if (job_req_node_filter(job_ptr, avail_bitmap))
		rc = ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE;
	if (job_ptr->details->exc_node_bitmap) {
		bitstr_t *exc_node_mask = NULL;
		exc_node_mask = bit_copy(job_ptr->details->exc_node_bitmap);
		if (exc_node_mask == NULL)
			fatal("bit_copy malloc failure");
		bit_not(exc_node_mask);
		bit_and(avail_bitmap, exc_node_mask);
		FREE_NULL_BITMAP(exc_node_mask);
	}
	if (job_ptr->details->req_node_bitmap) {
		if (!bit_super_set(job_ptr->details->req_node_bitmap, 
				   avail_bitmap)) {
			rc = ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE;
		}
	}

	if (rc == SLURM_SUCCESS) {
		min_nodes = MAX(job_ptr->details->min_nodes, 
				part_ptr->min_nodes);
		if (job_ptr->details->max_nodes == 0)
			max_nodes = part_ptr->max_nodes;
		else
			max_nodes = MIN(job_ptr->details->max_nodes, 
					part_ptr->max_nodes);
		max_nodes = MIN(max_nodes, 500000);	/* prevent overflows */
		if (job_ptr->details->max_nodes)
			req_nodes = max_nodes;
		else
			req_nodes = min_nodes;

		rc = select_g_job_test(job_ptr, avail_bitmap,
				min_nodes, max_nodes, req_nodes, 
				SELECT_MODE_WILL_RUN);
	}

	if (rc == SLURM_SUCCESS) {
		will_run_response_msg_t *resp_data;
		resp_data = xmalloc(sizeof(will_run_response_msg_t));
		resp_data->job_id     = job_ptr->job_id;
		resp_data->start_time = job_ptr->start_time;
		resp_data->node_list  = bitmap2node_name(avail_bitmap);
		FREE_NULL_BITMAP(avail_bitmap);
		*resp = resp_data;
		return SLURM_SUCCESS;
	} else {
		FREE_NULL_BITMAP(avail_bitmap);
		return ESLURM_REQUESTED_NODE_CONFIG_UNAVAILABLE;
	}

}
