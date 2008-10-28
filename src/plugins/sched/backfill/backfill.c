/*****************************************************************************\
 *  backfill.c - simple backfill scheduler plugin. 
 *
 *  If a partition is does not have root only access and nodes are not shared
 *  then raise the priority of pending jobs if doing so does not adversely
 *  effect the expected initiation of any higher priority job. We do not alter
 *  a job's required or excluded node list, so this is a conservative 
 *  algorithm.
 *
 *  For example, consider a cluster "lx[01-08]" with one job executing on 
 *  nodes "lx[01-04]". The highest priority pending job requires five nodes 
 *  including "lx05". The next highest priority pending job requires any 
 *  three nodes. Without explicitly forcing the second job to use nodes 
 *  "lx[06-08]", we can't start it without possibly delaying the higher 
 *  priority job.
 *****************************************************************************
 *  Copyright (C) 2003-2007 The Regents of the University of California.
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
 *  LLNL-CODE-402394.
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "slurm/slurm.h"
#include "slurm/slurm_errno.h"

#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/common/node_select.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"

#include "src/slurmctld/acct_policy.h"
#include "src/slurmctld/job_scheduler.h"
#include "src/slurmctld/licenses.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/node_scheduler.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/srun_comm.h"
#include "backfill.h"

typedef struct node_space_map {
	time_t begin_time;
	time_t end_time;
	bitstr_t *avail_bitmap;
	int next;	/* next record, by time, zero termination */
} node_space_map_t;
int backfilled_jobs = 0;

/*********************** local variables *********************/
static bool new_work      = false;
static bool stop_backfill = false;
static pthread_mutex_t thread_flag_mutex = PTHREAD_MUTEX_INITIALIZER;

/* Backfill scheduling has considerable overhead, 
 *	so only attempt it every BACKFILL_INTERVAL seconds.
 * Much of the scheduling for BlueGene happens through backfill,
 *	so we run it more frequently. */
#ifndef BACKFILL_INTERVAL
#  ifdef HAVE_BG
#    define BACKFILL_INTERVAL	5
#  else
#    define BACKFILL_INTERVAL	10
#  endif
#endif

/* Set __DEBUG to get detailed logging for this thread without 
 * detailed logging for the entire slurmctld daemon */
#define __DEBUG			0

/* Do not attempt to build job/resource/time record for
 * more than MAX_BACKFILL_JOB_CNT records */
#define MAX_BACKFILL_JOB_CNT	100

/* Do not build job/resource/time record for more than this 
 * far in the future, in seconds, currently one day */
#define BACKFILL_WINDOW		(24 * 60 * 60)

/*********************** local functions *********************/
static void _add_reservation(uint32_t start_time, uint32_t end_reserve, 
			     bitstr_t *res_bitmap, 
			     node_space_map_t *node_space, 
			     int *node_space_recs);
static void _attempt_backfill(void);
static void _diff_tv_str(struct timeval *tv1,struct timeval *tv2,
		char *tv_str, int len_tv_str);
static bool _more_work(void);
static int _start_job(struct job_record *job_ptr, bitstr_t *avail_bitmap);

#if __DEBUG
/* Log resource allocate table */
static void _dump_node_space_table(node_space_map_t *node_space_ptr)
{
	int i = 0;
	char begin_buf[32], end_buf[32], *node_list;

	info("=========================================");
	while (1) {
		slurm_make_time_str(&node_space_ptr[i].begin_time,
				    begin_buf, sizeof(begin_buf));
		slurm_make_time_str(&node_space_ptr[i].end_time,
				    end_buf, sizeof(end_buf));
		node_list = bitmap2node_name(node_space_ptr[i].avail_bitmap);
		info("Begin:%s End:%s Nodes:%s", begin_buf, end_buf, node_list);
		xfree(node_list);
		if ((i = node_space_ptr[i].next) == 0)
			break;
	}
	info("=========================================");
}
#endif

/*
 * _diff_tv_str - build a string showing the time difference between two times
 * IN tv1 - start of event
 * IN tv2 - end of event
 * OUT tv_str - place to put delta time in format "usec=%ld"
 * IN len_tv_str - size of tv_str in bytes
 */
static void _diff_tv_str(struct timeval *tv1,struct timeval *tv2,
		char *tv_str, int len_tv_str)
{
	long delta_t;
	delta_t  = (tv2->tv_sec  - tv1->tv_sec) * 1000000;
	delta_t +=  tv2->tv_usec - tv1->tv_usec;
	snprintf(tv_str, len_tv_str, "usec=%ld", delta_t);
}

/* Terminate backfill_agent */
extern void stop_backfill_agent(void)
{
	stop_backfill = true;
}


/* backfill_agent - detached thread periodically attempts to backfill jobs */
extern void *backfill_agent(void *args)
{
	struct timeval tv1, tv2;
	char tv_str[20];
	time_t now;
	int i, iter;
	static time_t last_backfill_time = 0;
	/* Read config, and partitions; Write jobs and nodes */
	slurmctld_lock_t all_locks = {
		READ_LOCK, WRITE_LOCK, WRITE_LOCK, READ_LOCK };

	while (!stop_backfill) {
		iter = (BACKFILL_CHECK_SEC * 1000000) /
		       STOP_CHECK_USEC;
		for (i=0; ((i<iter) && (!stop_backfill)); i++) {
			/* test stop_backfill every 0.1 sec for
			 * 2.0 secs to avoid running continuously */
			usleep(STOP_CHECK_USEC);
		}
		
		now = time(NULL);
		/* Avoid resource fragmentation if important */
		if (job_is_completing())
			continue;
		if ((difftime(now, last_backfill_time) < BACKFILL_INTERVAL) ||
		    stop_backfill || (!_more_work()))
			continue;
		last_backfill_time = now;

		gettimeofday(&tv1, NULL);
		lock_slurmctld(all_locks);
		_attempt_backfill();
		unlock_slurmctld(all_locks);
		gettimeofday(&tv2, NULL);
		_diff_tv_str(&tv1, &tv2, tv_str, 20);
#if __DEBUG
		info("backfill: completed, %s", tv_str);
#endif
	}
	return NULL;
}

static void _attempt_backfill(void)
{
	bool filter_root = false;
	struct job_queue *job_queue = NULL;
	int i, j,job_queue_size, node_space_recs;
	struct job_record *job_ptr;
	struct part_record *part_ptr;
	uint32_t end_time, end_reserve, time_limit;
	uint32_t min_nodes, max_nodes, req_nodes;
	uint16_t orig_shared;
	bitstr_t *avail_bitmap = NULL, *tmp_bitmap;
	time_t now = time(NULL);
	node_space_map_t node_space[MAX_BACKFILL_JOB_CNT + 2];

	if (slurm_get_root_filter())
		filter_root = true;

	job_queue_size = build_job_queue(&job_queue);
	if (job_queue_size == 0)
		return;

	sort_job_queue(job_queue, job_queue_size);

	node_space[0].begin_time = now;
	node_space[0].end_time = now + BACKFILL_WINDOW;
	node_space[0].avail_bitmap = bit_copy(avail_node_bitmap);
	node_space[0].next = 0;
	node_space_recs = 1;
#if __DEBUG
	_dump_node_space_table(node_space);
#endif

	for (i = 0; i < job_queue_size; i++) {
		job_ptr = job_queue[i].job_ptr;
		part_ptr = job_ptr->part_ptr;
#if __DEBUG
		info("backfill test for job %u", job_ptr->job_id);
#endif

		if (part_ptr == NULL) {
			part_ptr = find_part_record(job_ptr->partition);
			xassert(part_ptr);
			job_ptr->part_ptr = part_ptr;
			error("partition pointer reset for job %u, part %s",
			      job_ptr->job_id, job_ptr->partition);
		}
		if ((part_ptr->state_up == 0) ||
		    (part_ptr->node_bitmap == NULL))
		 	continue;
		if ((part_ptr->root_only) && filter_root)
			continue;

		if (!acct_policy_job_runnable(job_ptr))
			continue;
		if (license_job_test(job_ptr) != SLURM_SUCCESS)
			continue;
		if (job_independent(job_ptr) != true)
			continue;

		/* Determine minimum and maximum node counts */
		min_nodes = MAX(job_ptr->details->min_nodes,
				part_ptr->min_nodes);
		if (job_ptr->details->max_nodes == 0)
			max_nodes = part_ptr->max_nodes;
		else
			max_nodes = MIN(job_ptr->details->max_nodes,
					part_ptr->max_nodes);
		max_nodes = MIN(max_nodes, 500000);     /* prevent overflows */
		if (job_ptr->details->max_nodes)
			req_nodes = max_nodes;
		else
			req_nodes = min_nodes;
		if (min_nodes > max_nodes) {
			/* job's min_nodes exceeds partition's max_nodes */
			continue;
		}

		/* Determine job's expected completion time */
		if (job_ptr->time_limit == NO_VAL) {
			if (part_ptr->max_time == INFINITE)
				time_limit = 365 * 24 * 60; /* one year */
			else
				time_limit = part_ptr->max_time;
		} else {
			if (part_ptr->max_time == INFINITE)
				time_limit = job_ptr->time_limit;
			else
				time_limit = MIN(job_ptr->time_limit,
						 part_ptr->max_time);
		}
		end_time = (time_limit * 60) + now;

		/* Identify usable nodes for this job */
		FREE_NULL_BITMAP(avail_bitmap);
		avail_bitmap = bit_copy(part_ptr->node_bitmap);
		bit_and(avail_bitmap, up_node_bitmap);
		for (j=0; ; ) {
			if (node_space[j].begin_time <= end_time) {
				bit_and(avail_bitmap, 
					node_space[j].avail_bitmap);
			} else
				break;
			if ((j = node_space[j].next) == 0)
				break;
		}
		if (job_req_node_filter(job_ptr, avail_bitmap))
			continue;	/* problem with features */
		if (job_ptr->details->exc_node_bitmap) {
			bit_not(job_ptr->details->exc_node_bitmap);
			bit_and(avail_bitmap, 
				job_ptr->details->exc_node_bitmap);
			bit_not(job_ptr->details->exc_node_bitmap);
		}
		if ((job_ptr->details->req_node_bitmap) &&
		    (!bit_super_set(job_ptr->details->req_node_bitmap,
				    avail_bitmap)))
			continue;	/* required nodes missing */
		if (bit_set_count(avail_bitmap) < min_nodes)
			continue;	/* insufficient nodes remain */

		/* Try to schedule the job. First on dedicated nodes
		 * then on shared nodes (if so configured). */
		orig_shared = job_ptr->details->shared;
		job_ptr->details->shared = 0;
		tmp_bitmap = bit_copy(avail_bitmap);
		j = select_g_job_test(job_ptr, avail_bitmap, min_nodes,
				      max_nodes, req_nodes,
				      SELECT_MODE_WILL_RUN);
		job_ptr->details->shared = orig_shared;
		if ((j != SLURM_SUCCESS) && (orig_shared != 0)) {
			FREE_NULL_BITMAP(avail_bitmap);
			avail_bitmap= tmp_bitmap;
			j = select_g_job_test(job_ptr, avail_bitmap, min_nodes,
					      max_nodes, req_nodes,
					      SELECT_MODE_WILL_RUN);
		} else
			FREE_NULL_BITMAP(tmp_bitmap);
		if (j != SLURM_SUCCESS)
			continue;	/* not runable */

		if ((job_ptr->start_time <= now) &&
		    (_start_job(job_ptr, avail_bitmap) != SLURM_SUCCESS)) {
			/* Planned to start job, but something bad happended */
			break;
		}
		if (job_ptr->start_time > (now + BACKFILL_WINDOW)) {
			/* Starts too far in the future to worry about */
			continue;
		}

		if (node_space_recs == MAX_BACKFILL_JOB_CNT) {
			/* Already have too many jobs to deal with */
			break;
		}

		/*
		 * Add reservation to scheduling table
		 */
		end_reserve = job_ptr->start_time + (time_limit * 60);
		bit_not(avail_bitmap);
		_add_reservation(job_ptr->start_time, end_reserve, 
				 avail_bitmap, node_space, &node_space_recs);
#if __DEBUG
		_dump_node_space_table(node_space);
#endif
	}
	FREE_NULL_BITMAP(avail_bitmap);

	for (i=0; ; ) {
		bit_free(node_space[i].avail_bitmap);
		if ((i = node_space[i].next) == 0)
			break;
	}
	xfree(job_queue);
}

static int _start_job(struct job_record *job_ptr, bitstr_t *avail_bitmap)
{
	int rc;
	bitstr_t *orig_exc_nodes = NULL;
	static uint32_t fail_jobid = 0;

	if (job_ptr->details->exc_node_bitmap)
		orig_exc_nodes = job_ptr->details->exc_node_bitmap;
	job_ptr->details->exc_node_bitmap = bit_copy(avail_bitmap);
	bit_not(job_ptr->details->exc_node_bitmap);

	rc = select_nodes(job_ptr, false, NULL);
	bit_free(job_ptr->details->exc_node_bitmap);
	job_ptr->details->exc_node_bitmap = orig_exc_nodes;
	if (rc == SLURM_SUCCESS) {	
		/* job initiated */
		last_job_update = time(NULL);
		info("backfill: Started JobId=%u on %s",
		     job_ptr->job_id, job_ptr->nodes);
		if (job_ptr->batch_flag)
			launch_job(job_ptr);
		else
			srun_allocate(job_ptr->job_id);
		backfilled_jobs++;
#if __DEBUG
		info("backfill: Jobs backfilled: %d", backfilled_jobs);
#endif
	} else if (job_ptr->job_id != fail_jobid) {
		char *node_list = bitmap2node_name(avail_bitmap);
		/* This happens when a job has sharing disabled and
		 * a selected node is still completing some job, 
		 * which should be a temporary situation. */
		verbose("backfill: Failed to start JobId=%u on %s: %s",
			job_ptr->job_id, node_list, slurm_strerror(rc));
		xfree(node_list);
		fail_jobid = job_ptr->job_id;
	} else {
		debug3("backfill: Failed to start JobId=%u", job_ptr->job_id);
	}

	return rc;
}

/* trigger the attempt of a backfill */
extern void run_backfill (void)
{
	pthread_mutex_lock( &thread_flag_mutex );
	new_work = true;
	pthread_mutex_unlock( &thread_flag_mutex );
}

/* Report if any changes occurred to job, node or partition information */
static bool _more_work (void)
{
	bool rc;
	static time_t backfill_job_time  = (time_t) 0;
	static time_t backfill_node_time = (time_t) 0;
	static time_t backfill_part_time = (time_t) 0;

	pthread_mutex_lock( &thread_flag_mutex );
	if ( (backfill_job_time  == last_job_update ) &&
	     (backfill_node_time == last_node_update) &&
	     (backfill_part_time == last_part_update) &&
	     (new_work == false) ) {
		rc = false;
	} else {
		backfill_job_time  = last_job_update;
		backfill_node_time = last_node_update;
		backfill_part_time = last_part_update;
		new_work = false;
		rc = true;
	}
	pthread_mutex_unlock( &thread_flag_mutex );
	return rc;
}

/* Create a reservation for a job in the future */
static void _add_reservation(uint32_t start_time, uint32_t end_reserve, 
			     bitstr_t *res_bitmap, 
			     node_space_map_t *node_space, 
			     int *node_space_recs)
{
	int i, j;

	for (j=0; ; ) {
		if (node_space[j].end_time > start_time) {
			/* insert start entry record */
			i = *node_space_recs;
			node_space[i].begin_time = start_time;
			node_space[i].end_time = node_space[j].end_time;
			node_space[j].end_time = start_time;
			node_space[i].avail_bitmap = 
				bit_copy(node_space[j].avail_bitmap);
			node_space[i].next = node_space[j].next;
			node_space[j].next = i;
			(*node_space_recs)++;
			break;
		}
		if (node_space[j].end_time == start_time) {
			/* no need to insert start entry record */
			break;
		}
		if ((j = node_space[j].next) == 0)
			break;
	}

	for (j=0; ; ) {
		if (node_space[j].begin_time >= start_time)
			bit_and(node_space[j].avail_bitmap, res_bitmap);
		if ((j = node_space[j].next) == 0)
			break;
	}
}
