/*****************************************************************************\
 *  bg_job_place.c - blue gene job placement (e.g. base block selection)
 *  functions.
 *
 *  $Id$ 
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Dan Phung <phung4@llnl.gov> and Morris Jette <jette1@llnl.gov>
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

#include "src/common/node_select.h"
#include "bluegene.h"

#define _DEBUG 0

#define SWAP(a,b,t)	\
_STMT_START {		\
	(t) = (a);	\
	(a) = (b);	\
	(b) = (t);	\
} _STMT_END

static int  _find_best_block_match(struct job_record* job_ptr,
				   bitstr_t* slurm_block_bitmap,
				   int min_nodes, int max_nodes,
				   int spec, bg_record_t** found_bg_record,
				   bool test_only);
static void _rotate_geo(uint16_t *req_geometry, int rot_cnt);

/* Rotate a 3-D geometry array through its six permutations */
static void _rotate_geo(uint16_t *req_geometry, int rot_cnt)
{
	uint16_t tmp;

	switch (rot_cnt) {
		case 0:		/* ABC -> ACB */
		case 2:		/* CAB -> CBA */
		case 4:		/* BCA -> BAC */
			SWAP(req_geometry[Y], req_geometry[Z], tmp);
			break;
		case 1:		/* ACB -> CAB */
		case 3:		/* CBA -> BCA */
		case 5:		/* BAC -> ABC */
			SWAP(req_geometry[X], req_geometry[Y], tmp);
			break;
	}
}

pthread_mutex_t create_dynamic_mutex = PTHREAD_MUTEX_INITIALIZER;

/*
 * finds the best match for a given job request 
 * 
 * IN - int spec right now holds the place for some type of
 * specification as to the importance of certain job params, for
 * instance, geometry, type, size, etc.
 * 
 * OUT - block_id of matched block, NULL otherwise
 * returns 1 for error (no match)
 * 
 */
static int _find_best_block_match(struct job_record* job_ptr, 
				  bitstr_t* slurm_block_bitmap, 
				  int min_nodes, int max_nodes,
				  int spec, bg_record_t** found_bg_record, 
				  bool test_only)
{
	ListIterator itr;
	ListIterator itr2;
	bg_record_t *record = NULL;
	bg_record_t *found_record = NULL;
	uint16_t req_geometry[BA_SYSTEM_DIMENSIONS];
	uint16_t conn_type, rotate, target_size = 0;
	uint32_t req_procs = job_ptr->num_procs;
	uint32_t proc_cnt;
	ba_request_t request; 
	int i;
	int rot_cnt = 0;
	int created = 0;
	int found = 0;
	int max_procs = (uint16_t) NO_VAL;
	List lists_of_lists = NULL;
	List temp_list = NULL;
	char tmp_char[256];
	bitstr_t* tmp_bitmap = NULL;

	if(!bg_list) {
		error("_find_best_block_match: There is no bg_list");
		return SLURM_ERROR;
	}
	
	select_g_get_jobinfo(job_ptr->select_jobinfo,
			     SELECT_DATA_CONN_TYPE, &conn_type);
	select_g_get_jobinfo(job_ptr->select_jobinfo,
			     SELECT_DATA_GEOMETRY, &req_geometry);
	select_g_get_jobinfo(job_ptr->select_jobinfo,
			     SELECT_DATA_ROTATE, &rotate);
	select_g_get_jobinfo(job_ptr->select_jobinfo,
			     SELECT_DATA_MAX_PROCS, &max_procs);

	if(req_geometry[0] != (uint16_t)NO_VAL)
		for (i=0; i<BA_SYSTEM_DIMENSIONS; i++)
			target_size *= (uint16_t)req_geometry[i];
	if (target_size == 0) {	/* no geometry specified */
		target_size = min_nodes;
		req_geometry[X] = (uint16_t)NO_VAL;
	}
	/* this is where we should have the control flow depending on
	 * the spec arguement */

	*found_bg_record = NULL;
try_again:	
	debug("got here");
	slurm_mutex_lock(&block_state_mutex);
	debug("number of blocks to check: %d state %d", 
	      list_count(bg_list),
	      test_only);
     	itr = list_iterator_create(bg_list);
	while ((record = (bg_record_t*) list_next(itr))) {		
		/* If test_only we want to fall through to tell the 
		   scheduler that it is runnable just not right now. 
		*/
		debug3("%s job_running = %d", 
		       record->bg_block_id, record->job_running);
		/*partition is being destroyed (-2), 
		  or is messed up some how (-3) ignore it*/
		if(record->job_running < -1)
			continue;
		else if((record->job_running != -1) 
		   && !test_only) {
			debug("block %s in use by %s job %d", 
			      record->bg_block_id,
			      record->user_name,
			      record->job_running);
			found = 1;
			continue;
		}
		
		/* Check processor count */
		proc_cnt = record->bp_count * record->cpus_per_bp;
		debug3("asking for %d-%d looking at %d", 
		      req_procs, max_procs, proc_cnt);
		if ((proc_cnt < req_procs)
		    || (max_procs != (uint16_t) NO_VAL 
			&& proc_cnt > max_procs)) {
			/* We use the proccessor count per partition here
			   mostly to see if we can run on a smaller partition. 
			 */
			convert_to_kilo(proc_cnt, tmp_char);
			debug("block %s CPU count (%s) not suitable",
			      record->bg_block_id, 
			      tmp_char);
			continue;
		}

		/*
		 * check that the number of nodes is suitable
		 */
 		debug3("asking for %d-%d bps looking at %d", 
		      min_nodes, max_nodes, record->bp_count);
		if ((record->bp_count < min_nodes)
		    ||  (max_nodes != 0 && record->bp_count > max_nodes)
		    ||  (record->bp_count < target_size)) {
			convert_to_kilo(record->node_cnt, tmp_char);
			debug("block %s node count (%s) not suitable",
			      record->bg_block_id,
			      tmp_char);
			continue;
		}
		
		/*
		 * Next we check that this block's bitmap is within 
		 * the set of nodes which the job can use. 
		 * Nodes not available for the job could be down,
		 * drained, allocated to some other job, or in some 
		 * SLURM block not available to this job.
		 */
		if (!bit_super_set(record->bitmap, slurm_block_bitmap)) {
			debug("bg block %s has nodes not usable by this job",
			      record->bg_block_id);
			continue;
		}

		/*
		 * Insure that any required nodes are in this BG block
		 */
		if (job_ptr->details->req_node_bitmap
		    && (!bit_super_set(job_ptr->details->req_node_bitmap,
				       record->bitmap))) {
			debug("bg block %s lacks required nodes",
				record->bg_block_id);
			continue;
		}
				
		/* Make sure no other partitions are under this partition 
		   are booted and running jobs
		*/
		itr2 = list_iterator_create(bg_list);
		while ((found_record = (bg_record_t*)
			list_next(itr2)) != NULL) {
			if ((!found_record->bg_block_id)
			    || (record->job_running == -2)
			    || (!strcmp(record->bg_block_id, 
					found_record->bg_block_id)))
				continue;
			if(blocks_overlap(record, found_record)) {
				if(!test_only
				   && bluegene_layout_mode == LAYOUT_OVERLAP) {
					if(!created && record->state 
						!= RM_PARTITION_READY)
						break;
					else if(created == 1 
						&& found_record->state 
						!= RM_PARTITION_FREE) {
						break;
					} 
				}
				if(!test_only
				   && ((found_record->job_running > -1)
				   || (found_record->job_running == -3))) {
					if(found_record->job_running > -1)
						debug("can't use %s, there is "
						      "a job (%d) running on "
						      "an overlapping "
						      "block %s", 
						      record->bg_block_id,
						      found_record->
						      job_running,
						      found_record->
						      bg_block_id);
					else
						error("can't use %s, "
						      "overlapping block %s "
						      "is in an error state.",
						      record->bg_block_id,
						      found_record->
						      bg_block_id);
					if(bluegene_layout_mode == 
					   LAYOUT_DYNAMIC) {
						temp_list = list_create(NULL);
						list_push(temp_list, record);
						free_block_list(temp_list);
						num_block_to_free++;
						list_destroy(temp_list);
					} 
					break;
				} 
			} 
		}
		list_iterator_destroy(itr2);
		if(found_record) {
			found = 1;
			continue;
		} 
				
			
		/***********************************************/
		/* check the connection type specified matches */
		/***********************************************/
		if ((conn_type != record->conn_type)
		    && (conn_type != SELECT_NAV)) {
			debug("bg block %s conn-type not usable asking for %s "
			      "record is %s", 
			      record->bg_block_id,
			      convert_conn_type(conn_type),
			      convert_conn_type(record->conn_type));
			continue;
		} 

		/*****************************************/
		/* match up geometry as "best" possible  */
		/*****************************************/
		if (req_geometry[X] == (uint16_t)NO_VAL)
			;	/* Geometry not specified */
		else {	/* match requested geometry */
			bool match = false;
			rot_cnt = 0;	/* attempt six rotations  */

			for (rot_cnt=0; rot_cnt<6; rot_cnt++) {		
				if ((record->geo[X] >= req_geometry[X])
				&&  (record->geo[Y] >= req_geometry[Y])
				&&  (record->geo[Z] >= req_geometry[Z])) {
					match = true;
					break;
				}
				if (!rotate)
					break;
				_rotate_geo(req_geometry, rot_cnt);
			}

			if (!match) 
				continue;	/* Not usable */
		}
		*found_bg_record = record;
		break;
	}
	list_iterator_destroy(itr);

	if(bluegene_layout_mode == LAYOUT_OVERLAP 
	   &&!test_only && created<2 && !*found_bg_record) {
		created++;
		slurm_mutex_unlock(&block_state_mutex);
		goto try_again;
	}
		
	if(!found && test_only && bluegene_layout_mode == LAYOUT_DYNAMIC) {
		slurm_mutex_unlock(&block_state_mutex);
		
		for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
			request.start[i] = 0;
			
		for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
			request.geometry[i] = req_geometry[i];
			
		request.save_name = NULL;
		request.elongate_geos = NULL;
		request.size = target_size;
		request.procs = req_procs;
		request.conn_type = conn_type;
		request.rotate = rotate;
		request.elongate = true;
		request.start_req=0;
		debug("trying with all free blocks");
		if(create_dynamic_block(&request, NULL) == SLURM_ERROR) {
			error("this job will never run on "
			      "this system");
			xfree(request.save_name);
			return SLURM_ERROR;
		} else {
			if(!request.save_name) {
				error("no name returned from "
				      "create_dynamic_block");
				return SLURM_ERROR;
			} 
			slurm_conf_lock();
			sprintf(tmp_char, "%s%s\0", 
				slurmctld_conf.node_prefix, request.save_name);
			slurm_conf_unlock();
			if (node_name2bitmap(tmp_char, 
					     false, 
					     &tmp_bitmap)) {
				fatal("Unable to convert nodes %s to bitmap", 
				      request.save_name);
			}
			
			bit_and(slurm_block_bitmap, tmp_bitmap);
			bit_free(tmp_bitmap);
			xfree(request.save_name);
			return SLURM_SUCCESS;
		}
	} else if(!*found_bg_record 
		  && !created 
		  && bluegene_layout_mode == LAYOUT_DYNAMIC) {
		debug2("going to create %d", target_size);
		slurm_mutex_unlock(&block_state_mutex);
		lists_of_lists = list_create(NULL);
		list_append(lists_of_lists, bg_list);
		list_append(lists_of_lists, bg_booted_block_list);
		list_append(lists_of_lists, bg_job_block_list);
		itr = list_iterator_create(lists_of_lists);
		while ((temp_list = (List)list_next(itr)) != NULL) {
			created++;

			for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
				request.start[i] = 0;
			
			for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
				request.geometry[i] = req_geometry[i];
			
			request.save_name = NULL;
			request.elongate_geos = NULL;
			request.size = target_size;
			request.procs = req_procs;
			request.conn_type = conn_type;
			request.rotate = rotate;
			request.elongate = true;
			request.start_req=0;
			/* 1- try empty space
			   2- we see if we can create one in the 
			   unused bps
			   3- see if we can create one in the non 
			   job running bps
			*/
			debug("trying with %d", created);
			if(create_dynamic_block(&request, temp_list) 
			   == SLURM_SUCCESS) {
				list_iterator_destroy(itr);
				list_destroy(lists_of_lists);
				lists_of_lists = NULL;
				goto try_again;
			}
		}
		list_iterator_destroy(itr);
		if(lists_of_lists)
			list_destroy(lists_of_lists);
		slurm_mutex_lock(&block_state_mutex);		
	}
	
	/* set the bitmap and do other allocation activities */
	if (*found_bg_record) {
		format_node_name(*found_bg_record, tmp_char);
	
		debug("_find_best_block_match %s <%s>", 
			(*found_bg_record)->bg_block_id, 
			tmp_char);
		bit_and(slurm_block_bitmap, (*found_bg_record)->bitmap);
		slurm_mutex_unlock(&block_state_mutex);
		return SLURM_SUCCESS;
	}
		
	slurm_mutex_unlock(&block_state_mutex);
	debug("_find_best_block_match none found");
	return SLURM_ERROR;
}

/*
 * Try to find resources for a given job request
 * IN job_ptr - pointer to job record in slurmctld
 * IN/OUT bitmap - nodes availble for assignment to job, clear those not to
 *	be used
 * IN min_nodes, max_nodes  - minimum and maximum number of nodes to allocate
 *	to this job (considers slurm block limits)
 * IN test_only - if true, only test if ever could run, not necessarily now
 * RET - SLURM_SUCCESS if job runnable now, error code otherwise
 */
extern int submit_job(struct job_record *job_ptr, bitstr_t *slurm_block_bitmap,
		      int min_nodes, int max_nodes, bool test_only)
{
	int spec = 1; /* this will be like, keep TYPE a priority, etc,  */
	bg_record_t* record = NULL;
	char buf[100];
	int i, rc = SLURM_SUCCESS;
	uint16_t geo[BA_SYSTEM_DIMENSIONS];
	uint16_t tmp16 = (uint16_t)NO_VAL;
	
	
	select_g_sprint_jobinfo(job_ptr->select_jobinfo, buf, sizeof(buf), 
				SELECT_PRINT_MIXED);
	debug("bluegene:submit_job: %s nodes=%d-%d", 
	      buf, 
	      min_nodes, 
	      max_nodes);
	if(bluegene_layout_mode == LAYOUT_DYNAMIC)
		slurm_mutex_lock(&create_dynamic_mutex);
	
	rc = _find_best_block_match(job_ptr, slurm_block_bitmap, min_nodes, 
				    max_nodes, spec, &record, test_only);
	
	if (rc == SLURM_SUCCESS) {
		if(!record) {
			debug2("can run, but block not made");
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_BLOCK_ID,
					     "unassigned");
			if(job_ptr->num_procs < bluegene_bp_node_cnt) {
				i = procs_per_node/job_ptr->num_procs;
				info("divide by %d",i);
			} else 
				i = 1;
			min_nodes *= bluegene_bp_node_cnt/i;
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_NODE_CNT,
					     &min_nodes);
			
			for (i=0; i<BA_SYSTEM_DIMENSIONS; i++)
				geo[i] = 0;
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_GEOMETRY, 
					     &geo);
			
		} else {
			/* set the block id and info about block */
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_BLOCK_ID, 
					     record->bg_block_id);
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_QUARTER, 
					     &record->quarter);
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_NODECARD, 
					     &record->nodecard);
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_NODE_CNT, 
					     &record->node_cnt);
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_GEOMETRY, 
					     &record->geo);
			tmp16 = record->conn_type;
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_CONN_TYPE, 
					     &tmp16);
		}
		if(test_only) {
			select_g_set_jobinfo(job_ptr->select_jobinfo,
					     SELECT_DATA_BLOCK_ID,
					     "unassigned");
		} 
	}
	if(bluegene_layout_mode == LAYOUT_DYNAMIC)
		slurm_mutex_unlock(&create_dynamic_mutex);
	
	return rc;
}
