/*****************************************************************************\
 *  node_scheduler.h - definitions of functions in node_scheduler.c
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette@llnl.gov> et. al.
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

#ifndef _HAVE_NODE_SCHEDULER_H
#define _HAVE_NODE_SCHEDULER_H

/*
 * allocate_nodes - change state of specified nodes to NODE_STATE_ALLOCATED
 * IN job_ptr - job being allocated resources
 * globals: node_record_count - number of nodes in the system
 *	node_record_table_ptr - pointer to global node table
 *	last_node_update - last update time of node table
 */
extern void allocate_nodes(struct job_record *job_ptr);

/*
 * build_node_details - set cpu counts and addresses for allocated nodes:
 *	cpu_count_reps, cpus_per_node, node_addr, node_cnt, num_cpu_groups
 * IN job_ptr - pointer to a job record
 */
extern void build_node_details(struct job_record *job_ptr);

/*
 * count_cpus - report how many cpus are associated with the identified nodes 
 * IN bitmap - map of nodes to tally
 * RET cpu count
 * globals: node_record_count - number of nodes configured
 *	node_record_table_ptr - pointer to global node table
 */
extern int count_cpus(unsigned *bitmap);

/*
 * deallocate_nodes - for a given job, deallocate its nodes and make
 *	their state NODE_STATE_COMPLETING
 * IN job_ptr - pointer to terminating job (already in some COMPLETING state)
 * IN timeout - true if job exhausted time limit, send REQUEST_KILL_TIMELIMIT
 *	RPC instead of REQUEST_TERMINATE_JOB
 * IN suspended - true if job was already suspended (node's job_run_cnt
 *	already decremented);
 * globals: node_record_count - number of nodes in the system
 *	node_record_table_ptr - pointer to global node table
 */
extern void deallocate_nodes(struct job_record *job_ptr, bool timeout,
		bool suspended);

/*
 * re_kill_job - for a given job, deallocate its nodes for a second time, 
 *	basically a cleanup for failed deallocate() calls
 * IN job_ptr - pointer to terminating job (already in some COMPLETING state)
 * globals: node_record_count - number of nodes in the system
 *	node_record_table_ptr - pointer to global node table
 */
extern void re_kill_job(struct job_record *job_ptr);

/*
 * select_nodes - select and allocate nodes to a specific job
 * IN job_ptr - pointer to the job record
 * IN test_only - if set do not allocate nodes, just confirm they  
 *	could be allocated now
 * RET 0 on success, ESLURM code from slurm_errno.h otherwise
 * globals: list_part - global list of partition info
 *	default_part_loc - pointer to default partition 
 *	config_list - global list of node configuration info
 * Notes: The algorithm is
 *	1) Build a table (node_set_ptr) of nodes with the requisite 
 *	   configuration. Each table entry includes their weight, 
 *	   node_list, features, etc.
 *	2) Call _pick_best_nodes() to select those nodes best satisfying 
 *	   the request, (e.g. best-fit or other criterion)
 *	3) Call allocate_nodes() to perform the actual allocation
 */
extern int select_nodes(struct job_record *job_ptr, bool test_only);

#endif /* !_HAVE_NODE_SCHEDULER_H*/
