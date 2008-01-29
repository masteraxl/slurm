/*****************************************************************************\
 *  select_linear.c - node selection plugin for simple one-dimensional 
 *  address space. Selects nodes for a job so as to minimize the number 
 *  of sets of consecutive nodes using a best-fit algorithm.
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004-2006 The Regents of the University of California.
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
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  endif
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <slurm/slurm.h>
#include <slurm/slurm_errno.h>

#include "src/common/list.h"
#include "src/common/log.h"
#include "src/common/node_select.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_resource_info.h"
#include "src/common/xassert.h"
#include "src/common/xmalloc.h"

#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/proc_req.h"
#include "src/plugins/select/linear/select_linear.h"

#define SELECT_DEBUG	0
#define NO_SHARE_LIMIT	0xfffe

static int _add_job_to_nodes(struct node_cr_record *node_cr_ptr,
			     struct job_record *job_ptr, char *pre_err, 
			     int suspended);
static void _cr_job_list_del(void *x);
static int  _cr_job_list_sort(void *x, void *y);
static void _dump_node_cr(struct node_cr_record *node_cr_ptr);
static struct node_cr_record *_dup_node_cr(struct node_cr_record *node_cr_ptr);
static int  _find_job_mate(struct job_record *job_ptr, bitstr_t *bitmap,
			   uint32_t min_nodes, uint32_t max_nodes,
			   uint32_t req_nodes);
static void _free_node_cr(struct node_cr_record *node_cr_ptr);
static void _init_node_cr(void);
static int _job_count_bitmap(struct node_cr_record *node_cr_ptr,
			     struct job_record *job_ptr, 
			     bitstr_t * bitmap, bitstr_t * jobmap, int job_cnt);
static int _job_test(struct job_record *job_ptr, bitstr_t *bitmap,
		     uint32_t min_nodes, uint32_t max_nodes, 
		     uint32_t req_nodes);
static int _rm_job_from_nodes(struct node_cr_record *node_cr_ptr,
			      struct job_record *job_ptr, char *pre_err, 
			      int remove_all);
static int _will_run_test(struct job_record *job_ptr, bitstr_t *bitmap,
			  uint32_t min_nodes, uint32_t max_nodes, 
			  int max_share, uint32_t req_nodes);

/*
 * These variables are required by the generic plugin interface.  If they
 * are not found in the plugin, the plugin loader will ignore it.
 *
 * plugin_name - a string giving a human-readable description of the
 * plugin.  There is no maximum length, but the symbol must refer to
 * a valid string.
 *
 * plugin_type - a string suggesting the type of the plugin or its
 * applicability to a particular form of data or method of data handling.
 * If the low-level plugin API is used, the contents of this string are
 * unimportant and may be anything.  SLURM uses the higher-level plugin
 * interface which requires this string to be of the form
 *
 *	<application>/<method>
 *
 * where <application> is a description of the intended application of
 * the plugin (e.g., "select" for SLURM node selection) and <method>
 * is a description of how this plugin satisfies that application.  SLURM will
 * only load select plugins if the plugin_type string has a 
 * prefix of "select/".
 *
 * plugin_version - an unsigned 32-bit integer giving the version number
 * of the plugin.  If major and minor revisions are desired, the major
 * version number may be multiplied by a suitable magnitude constant such
 * as 100 or 1000.  Various SLURM versions will likely require a certain
 * minimum versions for their plugins as the node selection API matures.
 */
const char plugin_name[]       	= "Linear node selection plugin";
const char plugin_type[]       	= "select/linear";
const uint32_t plugin_version	= 90;

static struct node_record *select_node_ptr = NULL;
static int select_node_cnt = 0;
static uint16_t select_fast_schedule;
static uint16_t cr_type;

static struct node_cr_record *node_cr_ptr = NULL;
static pthread_mutex_t cr_mutex = PTHREAD_MUTEX_INITIALIZER;

#ifdef HAVE_XCPU
#define XCPU_POLL_TIME 120
static pthread_t xcpu_thread = 0;
static pthread_mutex_t thread_flag_mutex = PTHREAD_MUTEX_INITIALIZER;
static int agent_fini = 0;

static void *xcpu_agent(void *args)
{
	int i;
	static time_t last_xcpu_test;
	char reason[12], clone_path[128], down_node_list[512];
	struct stat buf;
	time_t now;

	last_xcpu_test = time(NULL) + XCPU_POLL_TIME;
	while (!agent_fini) {
		now = time(NULL);

		if (difftime(now, last_xcpu_test) >= XCPU_POLL_TIME) {
			debug3("Running XCPU node state test");
			down_node_list[0] = '\0';

			for (i=0; i<select_node_cnt; i++) {
				snprintf(clone_path, sizeof(clone_path), 
					"%s/%s/xcpu/clone", XCPU_DIR, 
					select_node_ptr[i].name);
				if (stat(clone_path, &buf) == 0)
					continue;
				error("stat %s: %m", clone_path);
				if ((strlen(select_node_ptr[i].name) +
				     strlen(down_node_list) + 2) <
				    sizeof(down_node_list)) {
					if (down_node_list[0] != '\0')
						strcat(down_node_list,",");
					strcat(down_node_list, 
						select_node_ptr[i].name);
				} else
					error("down_node_list overflow");
			}
			if (down_node_list[0]) {
				char time_str[32];
				slurm_make_time_str(&now, time_str, 	
					sizeof(time_str));
				snprintf(reason, sizeof(reason),
					"select_linear: Can not stat XCPU "
					"[SLURM@%s]", time_str);
				slurm_drain_nodes(down_node_list, reason);
			}
			last_xcpu_test = now;
		}

		sleep(1);
	}
	return NULL;
}

static int _init_status_pthread(void)
{
	pthread_attr_t attr;

	slurm_mutex_lock( &thread_flag_mutex );
	if ( xcpu_thread ) {
		debug2("XCPU thread already running, not starting "
			"another");
		slurm_mutex_unlock( &thread_flag_mutex );
		return SLURM_ERROR;
	}

	slurm_attr_init( &attr );
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_DETACHED );
	pthread_create( &xcpu_thread, &attr, xcpu_agent, NULL);
	slurm_mutex_unlock( &thread_flag_mutex );
	slurm_attr_destroy( &attr );

	return SLURM_SUCCESS;
}

static int _fini_status_pthread(void)
{
	int i, rc = SLURM_SUCCESS;

	slurm_mutex_lock( &thread_flag_mutex );
	if ( xcpu_thread ) {
		agent_fini = 1;
		for (i=0; i<4; i++) {
			if (pthread_kill(xcpu_thread, 0)) {
				xcpu_thread = 0;
				break;
			}
			sleep(1);
		}
		if ( xcpu_thread ) {
			error("could not kill XCPU agent thread");
			rc = SLURM_ERROR;
		}
	}
	slurm_mutex_unlock( &thread_flag_mutex );
	return rc;
}
#endif

static bool _enough_nodes(int avail_nodes, int rem_nodes, 
		uint32_t min_nodes, uint32_t req_nodes)
{
	int needed_nodes;

	if (req_nodes > min_nodes)
		needed_nodes = rem_nodes + min_nodes - req_nodes;
	else
		needed_nodes = rem_nodes;

	return(avail_nodes >= needed_nodes);
}

/*
 * init() is called when the plugin is loaded, before any other functions
 * are called.  Put global initialization here.
 */
extern int init ( void )
{
	int rc = SLURM_SUCCESS;
#ifdef HAVE_XCPU
	rc = _init_status_pthread();
#endif
#ifdef HAVE_BG
	error("%s is incompatable with BlueGene", plugin_name);
	fatal("Use SelectType=select/bluegene");
#endif
	cr_type = (select_type_plugin_info_t)
			slurmctld_conf.select_type_param;
	return rc;
}

extern int fini ( void )
{
	int rc = SLURM_SUCCESS;
#ifdef HAVE_XCPU
	rc = _fini_status_pthread();
#endif
	slurm_mutex_lock(&cr_mutex);
	_free_node_cr(node_cr_ptr);
	node_cr_ptr = NULL;
	slurm_mutex_unlock(&cr_mutex);
	return rc;
}

/*
 * The remainder of this file implements the standard SLURM 
 * node selection API.
 */

extern int select_p_state_save(char *dir_name)
{
	return SLURM_SUCCESS;
}

extern int select_p_state_restore(char *dir_name)
{
	return SLURM_SUCCESS;
}

extern int select_p_job_init(List job_list)
{
	return SLURM_SUCCESS;
}

extern int select_p_node_init(struct node_record *node_ptr, int node_cnt)
{
	if (node_ptr == NULL) {
		error("select_p_node_init: node_ptr == NULL");
		return SLURM_ERROR;
	}

	if (node_cnt < 0) {
		error("select_p_node_init: node_cnt < 0");
		return SLURM_ERROR;
	}

	/* NOTE: We free the consumable resources info here, but
	 * can't rebuild it since the partition and node structures
	 * have not yet had node bitmaps reset. */
	slurm_mutex_lock(&cr_mutex);
	_free_node_cr(node_cr_ptr);
	node_cr_ptr = NULL;
	slurm_mutex_unlock(&cr_mutex);

	select_node_ptr = node_ptr;
	select_node_cnt = node_cnt;
	select_fast_schedule = slurm_get_fast_schedule();

	return SLURM_SUCCESS;
}

extern int select_p_block_init(List part_list)
{
	return SLURM_SUCCESS;
}

/*
 * _get_avail_cpus - Get the number of "available" cpus on a node
 *	given this number given the number of cpus_per_task and
 *	maximum sockets, cores, threads.  Note that the value of
 *	cpus is the lowest-level logical processor (LLLP).
 * IN job_ptr - pointer to job being scheduled
 * IN index - index of node's configuration information in select_node_ptr
 */
static uint16_t _get_avail_cpus(struct job_record *job_ptr, int index)
{
	struct node_record *node_ptr;
	uint16_t avail_cpus;
	uint16_t cpus, sockets, cores, threads;
	uint16_t cpus_per_task = 1;
	uint16_t ntasks_per_node = 0, ntasks_per_socket = 0, ntasks_per_core = 0;
	uint16_t max_sockets = 0xffff, max_cores = 0xffff, max_threads = 0xffff;
	multi_core_data_t *mc_ptr = NULL;
	int min_sockets = 0, min_cores = 0;

	if (job_ptr->details == NULL)
		return (uint16_t) 0;

	if (job_ptr->details->cpus_per_task)
		cpus_per_task = job_ptr->details->cpus_per_task;
	if (job_ptr->details->ntasks_per_node)
		ntasks_per_node = job_ptr->details->ntasks_per_node;
	if ((mc_ptr = job_ptr->details->mc_ptr)) {
		max_sockets       = mc_ptr->max_sockets;
		max_cores         = mc_ptr->max_cores;
		max_threads       = mc_ptr->max_threads;
		ntasks_per_socket = mc_ptr->ntasks_per_socket;
		ntasks_per_core   = mc_ptr->ntasks_per_core;
	}

	node_ptr = &(select_node_ptr[index]);
	if (select_fast_schedule) { /* don't bother checking each node */
		cpus    = node_ptr->config_ptr->cpus;
		sockets = node_ptr->config_ptr->sockets;
		cores   = node_ptr->config_ptr->cores;
		threads = node_ptr->config_ptr->threads;
	} else {
		cpus    = node_ptr->cpus;
		sockets = node_ptr->sockets;
		cores   = node_ptr->cores;
		threads = node_ptr->threads;
	}

#if 0
	info("host %s User_ sockets %u cores %u threads %u ", 
	     node_ptr->name, max_sockets, max_cores, max_threads);

	info("host %s HW_ cpus %u sockets %u cores %u threads %u ", 
	     node_ptr->name, cpus, sockets, cores, threads);
#endif

	avail_cpus = slurm_get_avail_procs(
			max_sockets, max_cores, max_threads, 
			min_sockets, min_cores, cpus_per_task,
			ntasks_per_node, ntasks_per_socket, ntasks_per_core,
	    		&cpus, &sockets, &cores, &threads, NULL, 
			SELECT_TYPE_INFO_NONE,
			job_ptr->job_id, node_ptr->name);

#if 0
	debug3("avail_cpus index %d = %d (out of %d %d %d %d)",
				index, avail_cpus,
				cpus, sockets, cores, threads);
#endif
	return(avail_cpus);
}

/*
 * select_p_job_test - Given a specification of scheduling requirements, 
 *	identify the nodes which "best" satisfy the request.
 * 	"best" is defined as either single set of consecutive nodes satisfying 
 *	the request and leaving the minimum number of unused nodes OR 
 *	the fewest number of consecutive node sets
 * IN/OUT job_ptr - pointer to job being considered for initiation,
 *                  set's start_time when job expected to start
 * IN/OUT bitmap - usable nodes are set on input, nodes not required to 
 *	satisfy the request are cleared, other left set
 * IN min_nodes - minimum count of nodes
 * IN req_nodes - requested (or desired) count of nodes
 * IN max_nodes - maximum count of nodes (0==don't care)
 * IN mode - SELECT_MODE_RUN_NOW: try to schedule job now
 *           SELECT_MODE_TEST_ONLY: test if job can ever run
 *           SELECT_MODE_WILL_RUN: determine when and where job can run
 * RET zero on success, EINVAL otherwise
 * globals (passed via select_p_node_init): 
 *	node_record_count - count of nodes configured
 *	node_record_table_ptr - pointer to global node table
 * NOTE: the job information that is considered for scheduling includes:
 *	req_node_bitmap: bitmap of specific nodes required by the job
 *	contiguous: allocated nodes must be sequentially located
 *	num_procs: minimum number of processors required by the job
 * NOTE: bitmap must be a superset of the job's required at the time that 
 *	select_p_job_test is called
 */
extern int select_p_job_test(struct job_record *job_ptr, bitstr_t *bitmap,
			uint32_t min_nodes, uint32_t max_nodes, 
			uint32_t req_nodes, int mode)
{
	bitstr_t *orig_map;
	int i, j, rc = EINVAL, prev_cnt = -1;
	int min_share = 0, max_share = 0;
	uint32_t save_mem = 0;

	xassert(bitmap);
	if (job_ptr->details == NULL)
		return EINVAL;

	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL) {
		_init_node_cr();
		if (node_cr_ptr == NULL) {
			slurm_mutex_unlock(&cr_mutex);
			error("select_p_job_test: node_cr_ptr not initialized");
			return SLURM_ERROR;
		}
	}

	if (bit_set_count(bitmap) < min_nodes) {
		slurm_mutex_unlock(&cr_mutex);
		return EINVAL;
	}

	if (mode != SELECT_MODE_TEST_ONLY) {
		if (job_ptr->details->shared == 1) {
			max_share = job_ptr->part_ptr->max_share & 
					~SHARED_FORCE;
		} else	/* ((shared == 0) || (shared == (uint16_t) NO_VAL)) */
			max_share = 0;
	}

	if (mode == SELECT_MODE_WILL_RUN) {
		rc = _will_run_test(job_ptr, bitmap, min_nodes, max_nodes,
				    max_share, req_nodes);
		slurm_mutex_unlock(&cr_mutex);
		return rc;
	} else if (mode == SELECT_MODE_TEST_ONLY) {
		min_share = max_share = NO_SHARE_LIMIT;
		save_mem = job_ptr->details->job_min_memory;
		job_ptr->details->job_min_memory = 0;
	}

	orig_map = bit_copy(bitmap);
	for (i=min_share; i<=max_share; i++) {
		j = _job_count_bitmap(node_cr_ptr, job_ptr, orig_map, bitmap, i);
		if ((j == prev_cnt) || (j < min_nodes))
			continue;
		prev_cnt = j;
		if ((mode == SELECT_MODE_RUN_NOW) && (i > 0)) {
			/* We need to share. 
			 * Try to find suitable job to share nodes with. */
			rc = _find_job_mate(job_ptr, bitmap, 
					    min_nodes, max_nodes, req_nodes);
			if (rc == SLURM_SUCCESS)
				break;
		}
		rc = _job_test(job_ptr, bitmap, min_nodes, max_nodes, 
			       req_nodes);
		if (rc == SLURM_SUCCESS)
			break;
		continue;
	}
	bit_free(orig_map);
	slurm_mutex_unlock(&cr_mutex);
	if (save_mem)
		job_ptr->details->job_min_memory = save_mem;
	return rc;
}

/*
 * Set the bits in 'jobmap' that correspond to bits in the 'bitmap'
 * that are running 'job_cnt' jobs or less, and clear the rest.
 */
static int _job_count_bitmap(struct node_cr_record *node_cr_ptr,
			     struct job_record *job_ptr, 
			     bitstr_t * bitmap, bitstr_t * jobmap, int job_cnt)
{
	int i, count = 0, total_jobs;
	struct part_cr_record *part_cr_ptr;
	uint32_t job_memory = 0;

	xassert(node_cr_ptr);
	if (job_ptr->details->job_min_memory  && (cr_type == CR_MEMORY))
		job_memory = job_ptr->details->job_min_memory;

	for (i = 0; i < node_record_count; i++) {
		if (!bit_test(bitmap, i)) {
			bit_clear(jobmap, i);
			continue;
		}

		if (select_fast_schedule) {
			if ((node_cr_ptr[i].alloc_memory + job_memory) >
			     node_record_table_ptr[i].config_ptr->real_memory) {
				bit_clear(jobmap, i);
				continue;
			}
		} else {
			if ((node_cr_ptr[i].alloc_memory + job_memory) >
			     node_record_table_ptr[i].real_memory) {
				bit_clear(jobmap, i);
				continue;
			}
		}

		if ((job_cnt != NO_SHARE_LIMIT) &&
		    (node_cr_ptr[i].exclusive_jobid != 0)) {
			/* already reserved by some exclusive job */
			bit_clear(jobmap, i);
			continue;
		}

		total_jobs = 0;
		part_cr_ptr = node_cr_ptr[i].parts;
		while (part_cr_ptr) {
			if (job_cnt == 0)
				total_jobs += part_cr_ptr->run_job_cnt;
			else if (part_cr_ptr->part_ptr == job_ptr->part_ptr) {
				total_jobs += part_cr_ptr->run_job_cnt;
				break;
			}
			part_cr_ptr = part_cr_ptr->next;
		}
		if ((job_cnt != 0) && (part_cr_ptr == NULL)) {
			error("_job_count_bitmap: could not find "
				"partition %s for node %s",
				job_ptr->part_ptr->name,
				node_record_table_ptr[i].name);
		}
		if (total_jobs <= job_cnt) {
			bit_set(jobmap, i);
			count++;
		} else {
			bit_clear(jobmap, i);
		}

	}
	return count;
}

/* _find_job_mate - does most of the real work for select_p_job_test(), 
 *	in trying to find a suitable job to mate this one with. This is 
 *	a pretty simple algorithm now, but could try to match the job 
 *	with multiple jobs that add up to the proper size or a single 
 *	job plus a few idle nodes. */
static int _find_job_mate(struct job_record *job_ptr, bitstr_t *bitmap,
			  uint32_t min_nodes, uint32_t max_nodes,
			  uint32_t req_nodes)
{
	ListIterator job_iterator;
	struct job_record *job_scan_ptr;

	job_iterator = list_iterator_create(job_list);
	while ((job_scan_ptr = (struct job_record *) list_next(job_iterator))) {
		if ((job_scan_ptr->part_ptr == job_ptr->part_ptr) &&
		    (job_scan_ptr->job_state == JOB_RUNNING) &&
		    (job_scan_ptr->node_cnt == req_nodes) &&
		    (job_scan_ptr->total_procs >= job_ptr->num_procs) &&
		    bit_super_set(job_scan_ptr->node_bitmap, bitmap)) {
			bit_and(bitmap, job_scan_ptr->node_bitmap);
			return SLURM_SUCCESS;
		}
	}
	list_iterator_destroy(job_iterator);
	return EINVAL;
}

/* _job_test - does most of the real work for select_p_job_test(), which 
 *	pretty much just handles load-leveling and max_share logic */
static int _job_test(struct job_record *job_ptr, bitstr_t *bitmap,
			uint32_t min_nodes, uint32_t max_nodes, 
			uint32_t req_nodes)
{
	int i, index, error_code = EINVAL, sufficient;
	int *consec_nodes;	/* how many nodes we can add from this 
				 * consecutive set of nodes */
	int *consec_cpus;	/* how many nodes we can add from this 
				 * consecutive set of nodes */
	int *consec_start;	/* where this consecutive set starts (index) */
	int *consec_end;	/* where this consecutive set ends (index) */
	int *consec_req;	/* are nodes from this set required 
				 * (in req_bitmap) */
	int consec_index, consec_size;
	int rem_cpus, rem_nodes;	/* remaining resources desired */
	int best_fit_nodes, best_fit_cpus, best_fit_req;
	int best_fit_location = 0, best_fit_sufficient;
	int avail_cpus;

	if ((job_ptr->details->req_node_bitmap) &&
	    (!bit_super_set(job_ptr->details->req_node_bitmap, bitmap)))
		return error_code;

	consec_index = 0;
	consec_size  = 50;	/* start allocation for 50 sets of 
				 * consecutive nodes */
	consec_cpus  = xmalloc(sizeof(int) * consec_size);
	consec_nodes = xmalloc(sizeof(int) * consec_size);
	consec_start = xmalloc(sizeof(int) * consec_size);
	consec_end   = xmalloc(sizeof(int) * consec_size);
	consec_req   = xmalloc(sizeof(int) * consec_size);


	/* Build table with information about sets of consecutive nodes */
	consec_cpus[consec_index] = consec_nodes[consec_index] = 0;
	consec_req[consec_index] = -1;	/* no required nodes here by default */
	rem_cpus = job_ptr->num_procs;
	if (req_nodes > min_nodes)
		rem_nodes = req_nodes;
	else
		rem_nodes = min_nodes;

	for (index = 0; index < select_node_cnt; index++) {
		if (bit_test(bitmap, index)) {
			if (consec_nodes[consec_index] == 0)
				consec_start[consec_index] = index;

			avail_cpus = _get_avail_cpus(job_ptr, index);

			if (job_ptr->details->req_node_bitmap
			&&  bit_test(job_ptr->details->req_node_bitmap, index)
			&&  (max_nodes > 0)) {
				if (consec_req[consec_index] == -1) {
					/* first required node in set */
					consec_req[consec_index] = index;
				}
				rem_cpus -= avail_cpus;
				rem_nodes--;
				max_nodes--;
			} else {	 /* node not required (yet) */
				bit_clear(bitmap, index); 
				consec_cpus[consec_index] += avail_cpus;
				consec_nodes[consec_index]++;
			}
		} else if (consec_nodes[consec_index] == 0) {
			consec_req[consec_index] = -1;
			/* already picked up any required nodes */
			/* re-use this record */
		} else {
			consec_end[consec_index] = index - 1;
			if (++consec_index >= consec_size) {
				consec_size *= 2;
				xrealloc(consec_cpus,
					 sizeof(int) * consec_size);
				xrealloc(consec_nodes,
					 sizeof(int) * consec_size);
				xrealloc(consec_start,
					 sizeof(int) * consec_size);
				xrealloc(consec_end,
					 sizeof(int) * consec_size);
				xrealloc(consec_req,
					 sizeof(int) * consec_size);
			}
			consec_cpus[consec_index] = 0;
			consec_nodes[consec_index] = 0;
			consec_req[consec_index] = -1;
		}
	}
	if (consec_nodes[consec_index] != 0)
		consec_end[consec_index++] = index - 1;

#if SELECT_DEBUG
	/* don't compile this, slows things down too much */
	debug3("rem_cpus=%d, rem_nodes=%d", rem_cpus, rem_nodes);
	for (i = 0; i < consec_index; i++) {
		if (consec_req[i] != -1)
			debug3
			    ("start=%s, end=%s, nodes=%d, cpus=%d, req=%s",
			     select_node_ptr[consec_start[i]].name,
			     select_node_ptr[consec_end[i]].name,
			     consec_nodes[i], consec_cpus[i],
			     select_node_ptr[consec_req[i]].name);
		else
			debug3("start=%s, end=%s, nodes=%d, cpus=%d",
			       select_node_ptr[consec_start[i]].name,
			       select_node_ptr[consec_end[i]].name,
			       consec_nodes[i], consec_cpus[i]);
	}
#endif

	/* accumulate nodes from these sets of consecutive nodes until */
	/*   sufficient resources have been accumulated */
	while (consec_index && (max_nodes > 0)) {
		best_fit_cpus = best_fit_nodes = best_fit_sufficient = 0;
		best_fit_req = -1;	/* first required node, -1 if none */
		for (i = 0; i < consec_index; i++) {
			if (consec_nodes[i] == 0)
				continue;
			sufficient = (consec_cpus[i] >= rem_cpus)
			&& _enough_nodes(consec_nodes[i], rem_nodes,
					 min_nodes, req_nodes);

			/* if first possibility OR */
			/* contains required nodes OR */
			/* first set large enough for request OR */
			/* tightest fit (less resource waste) OR */
			/* nothing yet large enough, but this is biggest */
			if ((best_fit_nodes == 0) ||	
			    ((best_fit_req == -1) && (consec_req[i] != -1)) ||
			    (sufficient && (best_fit_sufficient == 0)) ||
			    (sufficient && (consec_cpus[i] < best_fit_cpus)) ||	
			    ((sufficient == 0) && 
			     (consec_cpus[i] > best_fit_cpus))) {
				best_fit_cpus = consec_cpus[i];
				best_fit_nodes = consec_nodes[i];
				best_fit_location = i;
				best_fit_req = consec_req[i];
				best_fit_sufficient = sufficient;
			}
		}
		if (best_fit_nodes == 0)
			break;
		if (job_ptr->details->contiguous && 
		    ((best_fit_cpus < rem_cpus) ||
		     (!_enough_nodes(best_fit_nodes, rem_nodes, 
				     min_nodes, req_nodes))))
			break;	/* no hole large enough */
		if (best_fit_req != -1) {
			/* This collection of nodes includes required ones
			 * select nodes from this set, first working up
			 * then down from the required nodes */
			for (i = best_fit_req;
			     i <= consec_end[best_fit_location]; i++) {
				if ((max_nodes <= 0)
				||  ((rem_nodes <= 0) && (rem_cpus <= 0)))
					break;
				if (bit_test(bitmap, i))
					continue;
				bit_set(bitmap, i);
				rem_nodes--;
				max_nodes--;
				avail_cpus = _get_avail_cpus(job_ptr, i);
				rem_cpus -= avail_cpus;
			}
			for (i = (best_fit_req - 1);
			     i >= consec_start[best_fit_location]; i--) {
				if ((max_nodes <= 0)
				||  ((rem_nodes <= 0) && (rem_cpus <= 0)))
					break;
				if (bit_test(bitmap, i)) 
					continue;
				bit_set(bitmap, i);
				rem_nodes--;
				max_nodes--;
				avail_cpus = _get_avail_cpus(job_ptr, i);
				rem_cpus -= avail_cpus;
			}
		} else {
			for (i = consec_start[best_fit_location];
			     i <= consec_end[best_fit_location]; i++) {
				if ((max_nodes <= 0)
				||  ((rem_nodes <= 0) && (rem_cpus <= 0)))
					break;
				if (bit_test(bitmap, i))
					continue;
				bit_set(bitmap, i);
				rem_nodes--;
				max_nodes--;
				avail_cpus = _get_avail_cpus(job_ptr, i);
				rem_cpus -= avail_cpus;
			}
		}
		if (job_ptr->details->contiguous || 
		    ((rem_nodes <= 0) && (rem_cpus <= 0))) {
			error_code = SLURM_SUCCESS;
			break;
		}
		consec_cpus[best_fit_location] = 0;
		consec_nodes[best_fit_location] = 0;
	}

	if (error_code && (rem_cpus <= 0)
	&&  _enough_nodes(0, rem_nodes, min_nodes, req_nodes)) {
		error_code = SLURM_SUCCESS;
	}

	xfree(consec_cpus);
	xfree(consec_nodes);
	xfree(consec_start);
	xfree(consec_end);
	xfree(consec_req);
	return error_code;
}

extern int select_p_job_begin(struct job_record *job_ptr)
{
	int rc = SLURM_SUCCESS;
#ifdef HAVE_XCPU
	int i;
	char clone_path[128];

	xassert(job_ptr);
	xassert(job_ptr->node_bitmap);

	for (i=0; i<select_node_cnt; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		snprintf(clone_path, sizeof(clone_path), 
			"%s/%s/xcpu/clone", XCPU_DIR, 
			select_node_ptr[i].name);
		if (chown(clone_path, (uid_t)job_ptr->user_id, 
				(gid_t)job_ptr->group_id)) {
			error("chown %s: %m", clone_path);
			rc = SLURM_ERROR;
		} else {
			debug("chown %s to %u", clone_path, 
				job_ptr->user_id);
		}
	}
#endif
	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL)
		_init_node_cr();
	_add_job_to_nodes(node_cr_ptr, job_ptr, "select_p_job_begin", 1);
	slurm_mutex_unlock(&cr_mutex);
	return rc;
}

extern int select_p_job_fini(struct job_record *job_ptr)
{
	int rc = SLURM_SUCCESS;
#ifdef HAVE_XCPU
	int i;
	char clone_path[128];

	for (i=0; i<select_node_cnt; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		snprintf(clone_path, sizeof(clone_path), 
			"%s/%s/xcpu/clone", XCPU_DIR, 
			select_node_ptr[i].name);
		if (chown(clone_path, (uid_t)0, (gid_t)0)) {
			error("chown %s: %m", clone_path);
			rc = SLURM_ERROR;
		} else {
			debug("chown %s to 0", clone_path);
		}
	}
#endif
	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL)
		_init_node_cr();
	_rm_job_from_nodes(node_cr_ptr, job_ptr, "select_p_job_fini", 1);
	slurm_mutex_unlock(&cr_mutex);
	return rc;
}

extern int select_p_job_suspend(struct job_record *job_ptr)
{
	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL)
		_init_node_cr();
	_rm_job_from_nodes(node_cr_ptr, job_ptr, "select_p_job_suspend", 0);
	slurm_mutex_unlock(&cr_mutex);
	return SLURM_SUCCESS;
}

extern int select_p_job_resume(struct job_record *job_ptr)
{
	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL)
		_init_node_cr();
	_add_job_to_nodes(node_cr_ptr, job_ptr, "select_p_job_resume", 0);
	slurm_mutex_unlock(&cr_mutex);
	return SLURM_SUCCESS;
}

extern int select_p_job_ready(struct job_record *job_ptr)
{
	if (job_ptr->job_state != JOB_RUNNING)
		return 0;

	return 1;
}

extern int select_p_pack_node_info(time_t last_query_time, Buf *buffer_ptr)
{
	/* This function is always invalid on normal Linux clusters */
	return SLURM_ERROR;
}

extern int select_p_get_select_nodeinfo (struct node_record *node_ptr, 
                                         enum select_data_info info,
                                         void *data)
{
       return SLURM_SUCCESS;
}

extern int select_p_update_nodeinfo (struct job_record *job_ptr)
{
       return SLURM_SUCCESS;
}

extern int select_p_update_block (update_part_msg_t *part_desc_ptr)
{
	return SLURM_SUCCESS;
}

extern int select_p_update_sub_node (update_part_msg_t *part_desc_ptr)
{
	return SLURM_SUCCESS;
}
extern int select_p_get_extra_jobinfo (struct node_record *node_ptr, 
                                       struct job_record *job_ptr, 
                                       enum select_data_info info,
                                       void *data)
{
	int rc = SLURM_SUCCESS;
	uint16_t *tmp_16;

	xassert(job_ptr);
	xassert(job_ptr->magic == JOB_MAGIC);

	switch (info) {
	case SELECT_AVAIL_CPUS:
		tmp_16 = (uint16_t *) data;

		if ((job_ptr->details->cpus_per_task > 1)
		||  (job_ptr->details->mc_ptr)) {
			int index = (node_ptr - node_record_table_ptr);
			*tmp_16 = _get_avail_cpus(job_ptr, index);
		} else {
			if (slurmctld_conf.fast_schedule) {
				*tmp_16 = node_ptr->config_ptr->cpus;
			} else {
				*tmp_16 = node_ptr->cpus;
			}
		}
		break;
	default:
		error("select_g_get_extra_jobinfo info %d invalid", info);
		rc = SLURM_ERROR;
		break;
	}
	
	return rc;
}

extern int select_p_get_info_from_plugin (enum select_data_info info,
					  void *data)
{
	return SLURM_SUCCESS;
}

extern int select_p_update_node_state (int index, uint16_t state)
{
	return SLURM_SUCCESS;
}

extern int select_p_alter_node_cnt(enum select_node_cnt type, void *data)
{
	return SLURM_SUCCESS;
}

extern int select_p_reconfigure(void)
{
	slurm_mutex_lock(&cr_mutex);
	_free_node_cr(node_cr_ptr);
	node_cr_ptr = NULL;
	_init_node_cr();
	slurm_mutex_unlock(&cr_mutex);

	return SLURM_SUCCESS;
}

/*
 * deallocate resources that were assigned to this job 
 *
 * if remove_all = 0: the job has been suspended, so just deallocate CPUs
 * if remove_all = 1: deallocate all resources
 */
static int _rm_job_from_nodes(struct node_cr_record *node_cr_ptr,
			      struct job_record *job_ptr, char *pre_err, 
			      int remove_all)
{
	int i, rc = SLURM_SUCCESS;
	struct part_cr_record *part_cr_ptr;
	uint32_t job_memory = 0;

	if (node_cr_ptr == NULL) {
		error("%s: node_cr_ptr not initialized", pre_err);
		return SLURM_ERROR;
	}

	if (remove_all && job_ptr->details && 
	    job_ptr->details->job_min_memory && (cr_type == CR_MEMORY))
		job_memory = job_ptr->details->job_min_memory;

	for (i = 0; i < select_node_cnt; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		if (node_cr_ptr[i].alloc_memory >= job_memory)
			node_cr_ptr[i].alloc_memory -= job_memory;
		else {
			node_cr_ptr[i].alloc_memory = 0;
			error("%s: memory underflow for node %s",
				pre_err, node_record_table_ptr[i].name);
		}
		if (node_cr_ptr[i].exclusive_jobid == job_ptr->job_id)
			node_cr_ptr[i].exclusive_jobid = 0;
		part_cr_ptr = node_cr_ptr[i].parts;
		while (part_cr_ptr) {
			if (part_cr_ptr->part_ptr != job_ptr->part_ptr) {
				part_cr_ptr = part_cr_ptr->next;
				continue;
			}
			if (part_cr_ptr->run_job_cnt > 0)
				part_cr_ptr->run_job_cnt--;
			else {
				error("%s: run_job_cnt underflow for node %s",
					pre_err, node_record_table_ptr[i].name);
			}
			if (remove_all) {
				if (part_cr_ptr->tot_job_cnt > 0)
					part_cr_ptr->tot_job_cnt--;
				else {
					error("%s: tot_job_cnt underflow "
						"for node %s",
						node_record_table_ptr[i].name);
				}
				if ((part_cr_ptr->tot_job_cnt == 0) &&
				    (part_cr_ptr->run_job_cnt)) {
					part_cr_ptr->run_job_cnt = 0;
					error("%s: run_job_count out of sync "
						"for node %s",
						node_record_table_ptr[i].name);
				}
			}
			break;
		}
		if (part_cr_ptr == NULL) {
			error("%s: could not find partition %s for node %s",
				pre_err, job_ptr->part_ptr->name,
				node_record_table_ptr[i].name);
			rc = SLURM_ERROR;
		}
	}

	return rc;
}

/*
 * allocate resources to the given job
 *
 * if alloc_all = 0: the job has been suspended, so just re-allocate CPUs
 * if alloc_all = 1: allocate all resources (CPUs and memory)
 */
static int _add_job_to_nodes(struct node_cr_record *node_cr_ptr,
			     struct job_record *job_ptr, char *pre_err, 
			     int alloc_all)
{
	int i, rc = SLURM_SUCCESS, exclusive = 0;
	struct part_cr_record *part_cr_ptr;
	uint32_t job_memory = 0;

	if (node_cr_ptr == NULL) {
		error("%s: node_cr_ptr not initialized", pre_err);
		return SLURM_ERROR;
	}

	if (alloc_all && job_ptr->details && 
	    job_ptr->details->job_min_memory && (cr_type == CR_MEMORY))
		job_memory = job_ptr->details->job_min_memory;
	if (job_ptr->details->shared == 0)
		exclusive = 1;

	for (i = 0; i < select_node_cnt; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		node_cr_ptr[i].alloc_memory += job_memory;
		if (exclusive) {
			if (node_cr_ptr[i].exclusive_jobid) {
				error("select/linear: conflicting exclusive "
				      "jobs %u and %u on %s",
				      job_ptr->job_id, 
				      node_cr_ptr[i].exclusive_jobid,
				      node_record_table_ptr[i].name);
			}
			node_cr_ptr[i].exclusive_jobid = job_ptr->job_id;
		}

		part_cr_ptr = node_cr_ptr[i].parts;
		while (part_cr_ptr) {
			if (part_cr_ptr->part_ptr != job_ptr->part_ptr) {
				part_cr_ptr = part_cr_ptr->next;
				continue;
			}
			if (alloc_all)
				part_cr_ptr->tot_job_cnt++;
			part_cr_ptr->run_job_cnt++;
			break;
		}
		if (part_cr_ptr == NULL) {
			error("%s: could not find partition %s for node %s",
				pre_err, job_ptr->part_ptr->name,
				node_record_table_ptr[i].name);
			rc = SLURM_ERROR;
		}
	}

	return rc;
}

static void _free_node_cr(struct node_cr_record *node_cr_ptr)
{
	int i;
	struct part_cr_record *part_cr_ptr1, *part_cr_ptr2;

	if (node_cr_ptr == NULL)
		return;

	for (i = 0; i < select_node_cnt; i++) {
		part_cr_ptr1 = node_cr_ptr[i].parts;
		while (part_cr_ptr1) {
			part_cr_ptr2 = part_cr_ptr1->next;
			xfree(part_cr_ptr1);
			part_cr_ptr1 = part_cr_ptr2;
		}
	}
	xfree(node_cr_ptr);
}

static inline void _dump_node_cr(struct node_cr_record *node_cr_ptr)
{
#if SELECT_DEBUG
	int i;
	struct part_cr_record *part_cr_ptr;

	if (node_cr_ptr == NULL)
		return;

	for (i = 0; i < select_node_cnt; i++) {
		info("Node:%s exclusive:%u alloc_mem:%u", 
			node_record_table_ptr[i].name,
			node_cr_ptr[i].exclusive_jobid,
			node_cr_ptr[i].alloc_memory);

		part_cr_ptr = node_cr_ptr[i].parts;
		while (part_cr_ptr) {
			info("  Part:%s run:%u tot:%u", 
				part_cr_ptr->part_ptr->name,
				part_cr_ptr->run_job_cnt,
				part_cr_ptr->tot_job_cnt);
			part_cr_ptr = part_cr_ptr->next;
		}
	}
#endif
}

static struct node_cr_record *_dup_node_cr(struct node_cr_record *node_cr_ptr)
{
	int i;
	struct node_cr_record *new_node_cr_ptr;
	struct part_cr_record *part_cr_ptr, *new_part_cr_ptr;

	if (node_cr_ptr == NULL)
		return NULL;

	new_node_cr_ptr = xmalloc(select_node_cnt * 
				  sizeof(struct node_cr_record));

	for (i = 0; i < select_node_cnt; i++) {
		new_node_cr_ptr[i].alloc_memory = node_cr_ptr[i].alloc_memory;
		new_node_cr_ptr[i].exclusive_jobid = 
				node_cr_ptr[i].exclusive_jobid;
		part_cr_ptr = node_cr_ptr[i].parts;
		while (part_cr_ptr) {
			new_part_cr_ptr = xmalloc(sizeof(struct part_cr_record));
			new_part_cr_ptr->part_ptr    = part_cr_ptr->part_ptr;
			new_part_cr_ptr->run_job_cnt = part_cr_ptr->run_job_cnt;
			new_part_cr_ptr->tot_job_cnt = part_cr_ptr->tot_job_cnt;
			new_part_cr_ptr->next 	     = new_node_cr_ptr[i].parts;
			new_node_cr_ptr[i].parts     = new_part_cr_ptr;
			part_cr_ptr = part_cr_ptr->next;
		}
	}
	return new_node_cr_ptr;
}

static void _init_node_cr(void)
{
	struct part_record *part_ptr;
	struct part_cr_record *part_cr_ptr;
	ListIterator part_iterator;
	struct job_record *job_ptr;
	ListIterator job_iterator;
	uint32_t job_memory;
	int exclusive, i;

	if (node_cr_ptr)
		return;

	node_cr_ptr = xmalloc(select_node_cnt * sizeof(struct node_cr_record));

	/* build partition records */
	part_iterator = list_iterator_create(part_list);
	while ((part_ptr = (struct part_record *) list_next(part_iterator))) {
		for (i = 0; i < select_node_cnt; i++) {
			if (part_ptr->node_bitmap == NULL)
				break;
			if (!bit_test(part_ptr->node_bitmap, i))
				continue;
			part_cr_ptr = xmalloc(sizeof(struct part_cr_record));
			part_cr_ptr->next = node_cr_ptr[i].parts;
			part_cr_ptr->part_ptr = part_ptr;
			node_cr_ptr[i].parts = part_cr_ptr;
		}
		
	}
	list_iterator_destroy(part_iterator);

	/* record running and suspended jobs in node_cr_records */
	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		if ((job_ptr->job_state != JOB_RUNNING) &&
		    (job_ptr->job_state != JOB_SUSPENDED))
			continue;

		if (job_ptr->details && 
		    job_ptr->details->job_min_memory && (cr_type == CR_MEMORY))
			job_memory = job_ptr->details->job_min_memory;
		else
			job_memory = 0;
		if (job_ptr->details->shared == 0)
			exclusive = 1;
		else
			exclusive = 0;

		for (i = 0; i < select_node_cnt; i++) {
			if (job_ptr->node_bitmap == NULL)
				break;
			if (!bit_test(job_ptr->node_bitmap, i))
				continue;
			if (exclusive) {
				if (node_cr_ptr[i].exclusive_jobid) {
					error("select/linear: conflicting "
				 	      "exclusive jobs %u and %u on %s",
				 	      job_ptr->job_id, 
				 	      node_cr_ptr[i].exclusive_jobid,
				 	      node_record_table_ptr[i].name);
				}
				node_cr_ptr[i].exclusive_jobid = job_ptr->job_id;
			}
			node_cr_ptr[i].alloc_memory += job_memory;
			part_cr_ptr = node_cr_ptr[i].parts;
			while (part_cr_ptr) {
				if (part_cr_ptr->part_ptr != job_ptr->part_ptr) {
					part_cr_ptr = part_cr_ptr->next;
					continue;
				}
				part_cr_ptr->tot_job_cnt++;
				if (job_ptr->job_state == JOB_RUNNING)
					part_cr_ptr->run_job_cnt++;
				break;
			}
			if (part_cr_ptr == NULL) {
				error("_init_node_cr: could not find "
					"partition %s for node %s",
					job_ptr->part_ptr->name,
					node_record_table_ptr[i].name);
			}
		}
	}
	list_iterator_destroy(job_iterator);
	_dump_node_cr(node_cr_ptr);
}

/* Determine where and when the job at job_ptr can begin execution by updating 
 * a scratch node_cr_record structure to reflect each job terminating at the 
 * end of its time limit and use this to show where and when the job at job_ptr 
 * will begin execution. Used by Moab for backfill scheduling. */
static int _will_run_test(struct job_record *job_ptr, bitstr_t *bitmap,
			  uint32_t min_nodes, uint32_t max_nodes, 
			  int max_share, uint32_t req_nodes)
{
	struct node_cr_record *exp_node_cr;
	struct job_record *tmp_job_ptr, **tmp_job_pptr;
	List cr_job_list;
	ListIterator job_iterator;
	bitstr_t *orig_map;
	int i, rc = SLURM_ERROR;

	orig_map = bit_copy(bitmap);

	/* Try to run with currently available nodes */
	i = _job_count_bitmap(node_cr_ptr, job_ptr, orig_map, bitmap, max_share);
	if (i >= min_nodes) {
		rc = _job_test(job_ptr, bitmap, min_nodes, max_nodes, 
			       req_nodes);
		if (rc == SLURM_SUCCESS) {
			bit_free(orig_map);
			job_ptr->start_time = time(NULL);
			return SLURM_SUCCESS;
		}
	}

	/* Job is still pending. Simulate termination of jobs one at a time 
	 * to determine when and where the job can start. */
	exp_node_cr = _dup_node_cr(node_cr_ptr);
	if (exp_node_cr == NULL) {
		bit_free(orig_map);
		return SLURM_ERROR;
	}

	/* Build list of running jobs */
	cr_job_list = list_create(_cr_job_list_del);
	job_iterator = list_iterator_create(job_list);
	while ((tmp_job_ptr = (struct job_record *) list_next(job_iterator))) {
		if (tmp_job_ptr->job_state != JOB_RUNNING)
			continue;
		if (tmp_job_ptr->end_time == 0) {
			error("Job %u has zero end_time", tmp_job_ptr->job_id);
			continue;
		}
		tmp_job_pptr = xmalloc(sizeof(struct job_record *));
		*tmp_job_pptr = tmp_job_ptr;
		list_append(cr_job_list, tmp_job_pptr);
	}
	list_iterator_destroy(job_iterator);
	list_sort(cr_job_list, _cr_job_list_sort);

	/* Remove the running jobs one at a time from exp_node_cr and try
	 * scheduling the pending job after each one */
	job_iterator = list_iterator_create(cr_job_list);
	while ((tmp_job_pptr = (struct job_record **) list_next(job_iterator))) {
		tmp_job_ptr = *tmp_job_pptr;
		_rm_job_from_nodes(exp_node_cr, tmp_job_ptr,
				   "_will_run_test", 1);
		i = _job_count_bitmap(exp_node_cr, job_ptr, orig_map, bitmap, 
				      max_share);
		if (i < min_nodes)
			continue;
		rc = _job_test(job_ptr, bitmap, min_nodes, max_nodes, 
			       req_nodes);
		if (rc != SLURM_SUCCESS)
			continue;
		job_ptr->start_time = tmp_job_ptr->end_time;
		break;
	}
	list_iterator_destroy(job_iterator);
	list_destroy(cr_job_list);
	_free_node_cr(exp_node_cr);
	bit_free(orig_map);
	return rc;
}

static void _cr_job_list_del(void *x)
{
	xfree(x);
}
static int  _cr_job_list_sort(void *x, void *y)
{
	struct job_record **job1_pptr = (struct job_record **) x;
	struct job_record **job2_pptr = (struct job_record **) y;
	return (int) difftime(job1_pptr[0]->end_time, job2_pptr[0]->end_time);
}

extern int select_p_step_begin(struct step_record *step_ptr)
{
	slurm_step_layout_t *step_layout = step_ptr->step_layout;
	int i, node_inx = -1;
	uint32_t avail_mem, step_mem;

info("step_begin: mem:%u", step_ptr->mem_per_task);
	xassert(step_ptr->job_ptr);
	xassert(step_ptr->job_ptr->details);
	xassert(step_layout);
	xassert(step_ptr->step_node_bitmap);
	xassert(step_layout->node_cnt == 
		bit_set_count(step_ptr->step_node_bitmap));

	/* Don't track step memory use if job has reserved memory OR
	 * job has whole node OR we don't track memory usage */
	if (step_ptr->job_ptr->details->job_min_memory || 
	    (step_ptr->job_ptr->details->shared == 0) ||
	    (cr_type != CR_MEMORY))
		return SLURM_SUCCESS;

	/* test if there is sufficient memory */
	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL)
		_init_node_cr();
	for (i = 0; i < select_node_cnt; i++) {
		if (bit_test(step_ptr->step_node_bitmap, i) == 0)
			continue;
		node_inx++;
		step_mem = step_layout->tasks[node_inx] * step_ptr->mem_per_task;
		if (select_fast_schedule)
			avail_mem = node_record_table_ptr[i].
				    config_ptr->real_memory;
		else
			avail_mem = node_record_table_ptr[i].real_memory;
info("alloc %u need %u avail %u", node_cr_ptr[i].alloc_memory, step_mem, avail_mem);
		if ((node_cr_ptr[i].alloc_memory + step_mem) > avail_mem) {
			slurm_mutex_unlock(&cr_mutex);
			return SLURM_ERROR;	/* no room */
		}
	}

	/* reserve the memory */
	node_inx = -1;
	for (i = 0; i < select_node_cnt; i++) {
		if (bit_test(step_ptr->step_node_bitmap, i) == 0)
			continue;
		node_inx++;
		step_mem = step_layout->tasks[node_inx] * step_ptr->mem_per_task;
		node_cr_ptr[i].alloc_memory += step_mem;
	}
	slurm_mutex_unlock(&cr_mutex);
	return SLURM_SUCCESS;
}

extern int select_p_step_fini(struct step_record *step_ptr)
{
	slurm_step_layout_t *step_layout = step_ptr->step_layout;
	int i, node_inx = -1;
	uint32_t step_mem;

info("step_fini: mem:%u", step_ptr->mem_per_task);
	xassert(step_ptr->job_ptr);
	xassert(step_ptr->job_ptr->details);
	xassert(step_layout);
	xassert(step_ptr->step_node_bitmap);
	xassert(step_layout->node_cnt == 
		bit_set_count(step_ptr->step_node_bitmap));

	/* Don't track step memory use if job has reserved memory OR
	 * job has whole node OR we don't track memory usage */
	if (step_ptr->job_ptr->details->job_min_memory || 
	    (step_ptr->job_ptr->details->shared == 0) ||
	    (cr_type != CR_MEMORY))
		return SLURM_SUCCESS;

	/* release the memory */
	slurm_mutex_lock(&cr_mutex);
	if (node_cr_ptr == NULL)
		_init_node_cr();
	for (i = 0; i < select_node_cnt; i++) {
		if (bit_test(step_ptr->step_node_bitmap, i) == 0)
			continue;
		node_inx++;
		step_mem = step_layout->tasks[node_inx] * step_ptr->mem_per_task;
		if (node_cr_ptr[i].alloc_memory >= step_mem)
			node_cr_ptr[i].alloc_memory -= step_mem;
		else {
			node_cr_ptr[i].alloc_memory = 0;
			error("select/linear: alloc_memory underflow on %s",
				node_record_table_ptr[i].name);
		}
	}
	slurm_mutex_unlock(&cr_mutex);
	return SLURM_SUCCESS;
}
