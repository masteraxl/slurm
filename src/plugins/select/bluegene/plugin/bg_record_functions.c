/*****************************************************************************\
 *  bg_record_functions.c - header for creating blocks in a static environment.
 *
 *  $Id: bg_record_functions.c 12954 2008-01-04 20:37:49Z da $
 *****************************************************************************
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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

#include "bluegene.h"
#include "dynamic_block.h"

#include "src/common/uid.h"
#include "src/slurmctld/trigger_mgr.h"

/* some local functions */
#ifdef HAVE_BG
static int  _addto_node_list(bg_record_t *bg_record, int *start, int *end);
static int  _ba_node_cmpf_inc(ba_node_t *node_a, ba_node_t *node_b);
#endif

extern void print_bg_record(bg_record_t* bg_record)
{
	if (!bg_record) {
		error("print_bg_record, record given is null");
		return;
	}
#if _DEBUG
	info(" bg_record: ");
	if (bg_record->bg_block_id)
		info("\tbg_block_id: %s", bg_record->bg_block_id);
	info("\tnodes: %s", bg_record->nodes);
	info("\tsize: %d BPs %u Nodes %d cpus", 
	     bg_record->bp_count,
	     bg_record->node_cnt,
	     bg_record->cpus_per_bp * bg_record->bp_count);
	info("\tgeo: %ux%ux%u", bg_record->geo[X], bg_record->geo[Y], 
	     bg_record->geo[Z]);
	info("\tconn_type: %s", convert_conn_type(bg_record->conn_type));
	info("\tnode_use: %s", convert_node_use(bg_record->node_use));
	if (bg_record->bitmap) {
		char bitstring[BITSIZE];
		bit_fmt(bitstring, BITSIZE, bg_record->bitmap);
		info("\tbitmap: %s", bitstring);
	}
#else
{
	char tmp_char[256];
	format_node_name(bg_record, tmp_char, sizeof(tmp_char));
	info("Record: BlockID:%s Nodes:%s Conn:%s",
	     bg_record->bg_block_id, tmp_char,
	     convert_conn_type(bg_record->conn_type));
}
#endif
}

extern void destroy_bg_record(void *object)
{
	bg_record_t* bg_record = (bg_record_t*) object;

	if (bg_record) {
		xfree(bg_record->bg_block_id);
		xfree(bg_record->nodes);
		xfree(bg_record->ionodes);
		xfree(bg_record->user_name);
		xfree(bg_record->target_name);
		if(bg_record->bg_block_list) {
			list_destroy(bg_record->bg_block_list);
			bg_record->bg_block_list = NULL;
		}
		FREE_NULL_BITMAP(bg_record->bitmap);
		FREE_NULL_BITMAP(bg_record->ionode_bitmap);

		xfree(bg_record->blrtsimage);
		xfree(bg_record->linuximage);
		xfree(bg_record->mloaderimage);
		xfree(bg_record->ramdiskimage);

		xfree(bg_record);
	}
}

extern int block_exist_in_list(List my_list, bg_record_t *bg_record)
{
	ListIterator itr = list_iterator_create(my_list);
	bg_record_t *found_record = NULL;
	int rc = 0;

	while ((found_record = (bg_record_t *) list_next(itr)) != NULL) {
		/* check for full node bitmap compare */
		if(bit_equal(bg_record->bitmap, found_record->bitmap)
		   && bit_equal(bg_record->ionode_bitmap,
				found_record->ionode_bitmap)) {
			if(bg_record->ionodes)
				debug3("This block %s[%s] "
				       "is already in the list %s",
				       bg_record->nodes,
				       bg_record->ionodes,
				       found_record->bg_block_id);
			else
				debug3("This block %s "
				       "is already in the list %s",
				       bg_record->nodes,
				       found_record->bg_block_id);
				
			rc = 1;
			break;
		}
	}
	list_iterator_destroy(itr);
	return rc;
}

extern void process_nodes(bg_record_t *bg_record, bool startup)
{
#ifdef HAVE_BG
	int j=0, number;
	int diff=0;
	int largest_diff=-1;
	int best_start[BA_SYSTEM_DIMENSIONS];
	int start[BA_SYSTEM_DIMENSIONS];
	int end[BA_SYSTEM_DIMENSIONS];
	ListIterator itr;
	ba_node_t* ba_node = NULL;
	
	if(!bg_record->bg_block_list 
	   || !list_count(bg_record->bg_block_list)) {
		if(!bg_record->bg_block_list) {
			bg_record->bg_block_list =
				list_create(destroy_ba_node);
		}
		memset(&best_start, 0, sizeof(best_start));
		bg_record->bp_count = 0;
		if((bg_record->conn_type == SELECT_SMALL) && (!startup))
			error("We shouldn't be here there could be some "
			      "badness if we use this logic %s",
			      bg_record->nodes);
		while (bg_record->nodes[j] != '\0') {
			if ((bg_record->nodes[j] == '['
			     || bg_record->nodes[j] == ',')
			    && (bg_record->nodes[j+8] == ']' 
				|| bg_record->nodes[j+8] == ',')
			    && (bg_record->nodes[j+4] == 'x'
				|| bg_record->nodes[j+4] == '-')) {
				j++;
				number = xstrntol(bg_record->nodes + j,
						  NULL, BA_SYSTEM_DIMENSIONS,
						  HOSTLIST_BASE);
				start[X] = number / 
					(HOSTLIST_BASE * HOSTLIST_BASE);
				start[Y] = (number % 
					    (HOSTLIST_BASE * HOSTLIST_BASE))
					/ HOSTLIST_BASE;
				start[Z] = (number % HOSTLIST_BASE);
				j += 4;
				number = xstrntol(bg_record->nodes + j,
						NULL, 3, HOSTLIST_BASE);
				end[X] = number /
					(HOSTLIST_BASE * HOSTLIST_BASE);
				end[Y] = (number 
					  % (HOSTLIST_BASE * HOSTLIST_BASE))
					/ HOSTLIST_BASE;
				end[Z] = (number % HOSTLIST_BASE);
				j += 3;
				diff = end[X]-start[X];
				if(diff > largest_diff) {
					best_start[X] = start[X];
					best_start[Y] = start[Y];
					best_start[Z] = start[Z];
					debug3("start is now %dx%dx%d",
					       best_start[X],
					       best_start[Y],
					       best_start[Z]);
					largest_diff = diff;
				}
				bg_record->bp_count += _addto_node_list(
					bg_record, 
					start, 
					end);
				if(bg_record->nodes[j] != ',')
					break;
				j--;
			} else if((bg_record->nodes[j] >= '0'
				   && bg_record->nodes[j] <= '9')
				  || (bg_record->nodes[j] >= 'A'
				      && bg_record->nodes[j] <= 'Z')) {
				
				number = xstrntol(bg_record->nodes + j,
						  NULL, BA_SYSTEM_DIMENSIONS,
						  HOSTLIST_BASE);
				start[X] = number / 
					(HOSTLIST_BASE * HOSTLIST_BASE);
				start[Y] = (number % 
					    (HOSTLIST_BASE * HOSTLIST_BASE))
					/ HOSTLIST_BASE;
				start[Z] = (number % HOSTLIST_BASE);
				j+=3;
				diff = 0;
				if(diff > largest_diff) {
					best_start[X] = start[X];
					best_start[Y] = start[Y];
					best_start[Z] = start[Z];
					debug3("start is now %dx%dx%d",
					       best_start[X],
					       best_start[Y],
					       best_start[Z]);
					largest_diff = diff;
				}
				bg_record->bp_count += _addto_node_list(
					bg_record, 
					start, 
					start);
				if(bg_record->nodes[j] != ',')
					break;
				j--;
			}
			j++;
		}
		if(largest_diff == -1) 
			fatal("No hostnames given here");

		bg_record->start[X] = best_start[X];
		bg_record->start[Y] = best_start[Y];
		bg_record->start[Z] = best_start[Z];
		debug2("start is %dx%dx%d",
		       bg_record->start[X],
		       bg_record->start[Y],
		       bg_record->start[Z]);
	}
	
	bg_record->geo[X] = 0;
	bg_record->geo[Y] = 0;
	bg_record->geo[Z] = 0;
	end[X] = -1;
	end[Y] = -1;
	end[Z] = -1;

	list_sort(bg_record->bg_block_list, (ListCmpF) _ba_node_cmpf_inc);

	itr = list_iterator_create(bg_record->bg_block_list);
	while ((ba_node = list_next(itr)) != NULL) {
		if(!ba_node->used)
			continue;
		debug4("%c%c%c is included in this block",
		       alpha_num[ba_node->coord[X]],
		       alpha_num[ba_node->coord[Y]],
		       alpha_num[ba_node->coord[Z]]);
		       
		if(ba_node->coord[X]>end[X]) {
			bg_record->geo[X]++;
			end[X] = ba_node->coord[X];
		}
		if(ba_node->coord[Y]>end[Y]) {
			bg_record->geo[Y]++;
			end[Y] = ba_node->coord[Y];
		}
		if(ba_node->coord[Z]>end[Z]) {
			bg_record->geo[Z]++;
			end[Z] = ba_node->coord[Z];
		}
	}
	list_iterator_destroy(itr);
	debug3("geo = %c%c%c bp count is %d\n",
	       alpha_num[bg_record->geo[X]],
	       alpha_num[bg_record->geo[Y]],
	       alpha_num[bg_record->geo[Z]],
	       bg_record->bp_count);

	if ((bg_record->geo[X] == DIM_SIZE[X])
	    && (bg_record->geo[Y] == DIM_SIZE[Y])
	    && (bg_record->geo[Z] == DIM_SIZE[Z])) {
		bg_record->full_block = 1;	
	}	
	
/* #ifndef HAVE_BG_FILES */
/* 	max_dim[X] = MAX(max_dim[X], end[X]); */
/* 	max_dim[Y] = MAX(max_dim[Y], end[Y]); */
/* 	max_dim[Z] = MAX(max_dim[Z], end[Z]); */
/* #endif */
   
	if (node_name2bitmap(bg_record->nodes, 
			     false, 
			     &bg_record->bitmap)) {
		fatal("1 Unable to convert nodes %s to bitmap", 
		      bg_record->nodes);
	}
#endif
	return;
}

/* 
 * NOTE: This function does not do a mutex lock so if you are copying the
 * main bg_list you need to lock 'block_state_mutex' before calling
 */
extern List copy_bg_list(List in_list) 
{
	bg_record_t *bg_record = NULL;
	bg_record_t *new_record = NULL;
	List out_list = list_create(destroy_bg_record);
	ListIterator itr = list_iterator_create(in_list);

	while ((bg_record = (bg_record_t *) list_next(itr))) { 
		new_record = xmalloc(sizeof(bg_record_t));
		new_record->original = bg_record;
		copy_bg_record(bg_record, new_record);
		list_append(out_list, new_record);
	}

	list_iterator_destroy(itr);
	
	return out_list;	
}

extern void copy_bg_record(bg_record_t *fir_record, bg_record_t *sec_record)
{
	int i;
	ListIterator itr = NULL;
	ba_node_t *ba_node = NULL, *new_ba_node = NULL;
	
	if(!fir_record || !sec_record) {
		error("copy_bg_record: "
		      "given a null for either first record or second record");
		return;
	}

	xfree(sec_record->bg_block_id);
	sec_record->bg_block_id = xstrdup(fir_record->bg_block_id);
	xfree(sec_record->nodes);
	sec_record->nodes = xstrdup(fir_record->nodes);
	xfree(sec_record->ionodes);
	sec_record->ionodes = xstrdup(fir_record->ionodes);
	xfree(sec_record->user_name);
	sec_record->user_name = xstrdup(fir_record->user_name);
	xfree(sec_record->target_name);
	sec_record->target_name = xstrdup(fir_record->target_name);

	xfree(sec_record->blrtsimage);
	sec_record->blrtsimage = xstrdup(fir_record->blrtsimage);
	xfree(sec_record->linuximage);
	sec_record->linuximage = xstrdup(fir_record->linuximage);
	xfree(sec_record->mloaderimage);
	sec_record->mloaderimage = xstrdup(fir_record->mloaderimage);
	xfree(sec_record->ramdiskimage);
	sec_record->ramdiskimage = xstrdup(fir_record->ramdiskimage);

	sec_record->user_uid = fir_record->user_uid;
	sec_record->state = fir_record->state;
	sec_record->conn_type = fir_record->conn_type;
	sec_record->node_use = fir_record->node_use;
	sec_record->bp_count = fir_record->bp_count;
	sec_record->switch_count = fir_record->switch_count;
	sec_record->boot_state = fir_record->boot_state;
	sec_record->boot_count = fir_record->boot_count;
	sec_record->full_block = fir_record->full_block;

	for(i=0;i<BA_SYSTEM_DIMENSIONS;i++) {
		sec_record->geo[i] = fir_record->geo[i];
		sec_record->start[i] = fir_record->start[i];
	}

	FREE_NULL_BITMAP(sec_record->bitmap);
	if(fir_record->bitmap 
	   && (sec_record->bitmap = bit_copy(fir_record->bitmap)) == NULL) {
		error("Unable to copy bitmap for %s", fir_record->nodes);
		sec_record->bitmap = NULL;
	}
	FREE_NULL_BITMAP(sec_record->ionode_bitmap);
	if(fir_record->ionode_bitmap 
	   && (sec_record->ionode_bitmap
	       = bit_copy(fir_record->ionode_bitmap)) == NULL) {
		error("Unable to copy ionode_bitmap for %s",
		      fir_record->nodes);
		sec_record->ionode_bitmap = NULL;
	}
	if(sec_record->bg_block_list)
		list_destroy(sec_record->bg_block_list);
	sec_record->bg_block_list = list_create(destroy_ba_node);
	if(fir_record->bg_block_list) {
		itr = list_iterator_create(fir_record->bg_block_list);
		while((ba_node = list_next(itr))) {
			new_ba_node = ba_copy_node(ba_node);
			list_push(sec_record->bg_block_list, new_ba_node);
		}
		list_iterator_destroy(itr);
	}
	sec_record->job_running = fir_record->job_running;
	sec_record->job_ptr = fir_record->job_ptr;
	sec_record->cpus_per_bp = fir_record->cpus_per_bp;
	sec_record->node_cnt = fir_record->node_cnt;
	sec_record->quarter = fir_record->quarter;
	sec_record->nodecard = fir_record->nodecard;
}

extern bg_record_t *find_bg_record_in_list(List my_list, char *bg_block_id)
{
	ListIterator itr;
	bg_record_t *bg_record = NULL;
		
	if(!bg_block_id)
		return NULL;
			
	if(my_list) {
		slurm_mutex_lock(&block_state_mutex);
		itr = list_iterator_create(my_list);
		while ((bg_record = (bg_record_t *) list_next(itr)) != NULL) {
			if(bg_record->bg_block_id)
				if (!strcmp(bg_record->bg_block_id, 
					    bg_block_id))
					break;
		}
		list_iterator_destroy(itr);
		slurm_mutex_unlock(&block_state_mutex);
		if(bg_record)
			return bg_record;
		else
			return NULL;
	} else {
		error("find_bg_record_in_list: no list");
		return NULL;
	}
	
}

/* All changes to the bg_list target_name must 
   be done before this function is called. 
   also slurm_conf_lock() must be called before calling this
   function along with slurm_conf_unlock() afterwards.		
*/
extern int update_block_user(bg_record_t *bg_record, int set) 
{
	if(!bg_record->target_name) {
		error("Must set target_name to run update_block_user.");
		return -1;
	}
	if(!bg_record->user_name) {
		error("No user_name");
		bg_record->user_name = xstrdup(slurmctld_conf.slurm_user_name);
	}
#ifdef HAVE_BG_FILES
	int rc=0;	
	if(set) {
		if((rc = remove_all_users(bg_record->bg_block_id, 
					  bg_record->target_name))
		   == REMOVE_USER_ERR) {
			error("1 Something happened removing "
			      "users from block %s", 
			      bg_record->bg_block_id);
			return -1;
		} else if (rc == REMOVE_USER_NONE) {
			if (strcmp(bg_record->target_name, 
				   slurmctld_conf.slurm_user_name)) {
				info("Adding user %s to Block %s",
				     bg_record->target_name, 
				     bg_record->bg_block_id);
				
				if ((rc = bridge_add_block_user(
					     bg_record->bg_block_id, 
					     bg_record->target_name)) 
				    != STATUS_OK) {
					error("bridge_add_block_user"
					      "(%s,%s): %s", 
					      bg_record->bg_block_id, 
					      bg_record->target_name,
					      bg_err_str(rc));
					return -1;
				} 
			}
		}
	}
#endif
	
	if(strcmp(bg_record->target_name, bg_record->user_name)) {
		uid_t pw_uid;
		xfree(bg_record->user_name);
		bg_record->user_name = xstrdup(bg_record->target_name);
		pw_uid = uid_from_string(bg_record->user_name);
		if(pw_uid == (uid_t) -1) {
			error("No such user: %s", bg_record->user_name);
			return -1;
		} else {
			bg_record->user_uid = pw_uid; 
		}		
		return 1;
	}
	
	return 0;
}

/* If any nodes in node_list are drained, draining, or down, 
 *   then just return
 *   else drain all of the nodes
 * This function lets us drain an entire bgblock only if 
 * we have not already identified a specific node as bad. */
extern void drain_as_needed(bg_record_t *bg_record, char *reason)
{
	bool needed = true;
	hostlist_t hl;
	char *host = NULL;
	char bg_down_node[128];

	if(bg_record->job_running > NO_JOB_RUNNING)
		slurm_fail_job(bg_record->job_running);			

	/* small blocks */
	if(bg_record->cpus_per_bp != procs_per_node) {
		debug2("small block");
		goto end_it;
	}
	
	/* at least one base partition */
	hl = hostlist_create(bg_record->nodes);
	if (!hl) {
		slurm_drain_nodes(bg_record->nodes, reason);
		return;
	}
	while ((host = hostlist_shift(hl))) {
		if (node_already_down(bg_down_node)) {
			needed = false;
			free(host);
			break;
		}
		free(host);
	}
	hostlist_destroy(hl);
	
	if (needed) {
		slurm_drain_nodes(bg_record->nodes, reason);
	}
end_it:
	while(bg_record->job_running > NO_JOB_RUNNING) {
		debug2("block %s is still running job %d",
		       bg_record->bg_block_id, bg_record->job_running);
		sleep(1);
	}
	
	slurm_mutex_lock(&block_state_mutex);
	error("Setting Block %s to ERROR state.", bg_record->bg_block_id);
	bg_record->job_running = BLOCK_ERROR_STATE;
	bg_record->state = RM_PARTITION_ERROR;
	slurm_mutex_unlock(&block_state_mutex);
	trigger_block_error();
	return;
}

extern int set_ionodes(bg_record_t *bg_record)
{
	int i = 0;
	int start_bit = 0;
	int size = 0;
	char bitstring[BITSIZE];
	
	if(!bg_record)
		return SLURM_ERROR;
	/* set the bitmap blank here if it is a full node we don't
	   want anything set we also don't want the bg_record->ionodes set.
	*/
	bg_record->ionode_bitmap = bit_alloc(bluegene_numpsets);
	if(bg_record->quarter == (uint16_t)NO_VAL) {
		return SLURM_SUCCESS;
	}

	start_bit = bluegene_quarter_ionode_cnt*bg_record->quarter;
	
	if(bg_record->nodecard != (uint16_t)NO_VAL
	   && bluegene_nodecard_ionode_cnt) {
		start_bit += bluegene_nodecard_ionode_cnt*bg_record->nodecard;
		size = bluegene_nodecard_ionode_cnt;
	} else
		size = bluegene_quarter_ionode_cnt;
	size += start_bit;

	if(size == start_bit) {
		error("start bit is the same as the end bit %d", size);
		return SLURM_ERROR;
	}
	for(i=start_bit; i<size; i++)
		bit_set(bg_record->ionode_bitmap, i);
	
	bit_fmt(bitstring, BITSIZE, bg_record->ionode_bitmap);
	bg_record->ionodes = xstrdup(bitstring);

	return SLURM_SUCCESS;
}

extern int add_bg_record(List records, List used_nodes, blockreq_t *blockreq)
{
	bg_record_t *bg_record = NULL;
	bg_record_t *found_record = NULL;
	ba_node_t *ba_node = NULL;
	ListIterator itr;
	uid_t pw_uid;
	int i, len;
	int small_size = 0;
	int small_count = 0;
	uint16_t quarter = 0;
	uint16_t nodecard = 0;
	int node_cnt = 0;
	
	if(!records) {
		fatal("add_bg_record: no records list given");
	}
	bg_record = (bg_record_t*) xmalloc(sizeof(bg_record_t));
	
	slurm_conf_lock();
	bg_record->user_name = 
		xstrdup(slurmctld_conf.slurm_user_name);
	bg_record->target_name = 
		xstrdup(slurmctld_conf.slurm_user_name);
	slurm_conf_unlock();
	pw_uid = uid_from_string(bg_record->user_name);
	if(pw_uid == (uid_t) -1) {
		error("No such user: %s", bg_record->user_name);
	} else {
		bg_record->user_uid = pw_uid;
	}

	bg_record->bg_block_list = list_create(destroy_ba_node);
	if(used_nodes) {
		if(copy_node_path(used_nodes, &bg_record->bg_block_list)
		   == SLURM_ERROR)
			error("couldn't copy the path for the allocation");
		bg_record->bp_count = list_count(used_nodes);
	}
	bg_record->quarter = (uint16_t)NO_VAL;
	bg_record->nodecard = (uint16_t)NO_VAL;
	if(set_ionodes(bg_record) == SLURM_ERROR) {
		fatal("add_bg_record: problem creating ionodes");
	}
	/* bg_record->boot_state = 0; 	Implicit */
	/* bg_record->state = 0;	Implicit */
	debug2("asking for %s %d %d %s", 
	       blockreq->block, blockreq->quarters, blockreq->nodecards,
	       convert_conn_type(blockreq->conn_type));
	len = strlen(blockreq->block);
	i=0;
	while(i<len 
	      && blockreq->block[i] != '[' 
	      && (blockreq->block[i] < '0' || blockreq->block[i] > 'Z'
		  || (blockreq->block[i] > '9' && blockreq->block[i] < 'A')))
		i++;
	
	if(i<len) {
		len -= i;
		slurm_conf_lock();
		len += strlen(slurmctld_conf.node_prefix)+1;
		bg_record->nodes = xmalloc(len);
		snprintf(bg_record->nodes, len, "%s%s", 
			slurmctld_conf.node_prefix, blockreq->block+i);
		slurm_conf_unlock();
			
	} else 
		fatal("BPs=%s is in a weird format", blockreq->block); 
	
	process_nodes(bg_record, false);
	
	bg_record->node_use = SELECT_COPROCESSOR_MODE;
	bg_record->conn_type = blockreq->conn_type;
	bg_record->cpus_per_bp = procs_per_node;
	bg_record->node_cnt = bluegene_bp_node_cnt * bg_record->bp_count;
	bg_record->job_running = NO_JOB_RUNNING;

	if(blockreq->blrtsimage)
		bg_record->blrtsimage = xstrdup(blockreq->blrtsimage);
	else
		bg_record->blrtsimage = xstrdup(default_blrtsimage);

	if(blockreq->linuximage)
		bg_record->linuximage = xstrdup(blockreq->linuximage);
	else
		bg_record->linuximage = xstrdup(default_linuximage);

	if(blockreq->mloaderimage)
		bg_record->mloaderimage = xstrdup(blockreq->mloaderimage);
	else
		bg_record->mloaderimage = xstrdup(default_mloaderimage);

	if(blockreq->ramdiskimage)
		bg_record->ramdiskimage = xstrdup(blockreq->ramdiskimage);
	else
		bg_record->ramdiskimage = xstrdup(default_ramdiskimage);
		
	if(bg_record->conn_type != SELECT_SMALL) {
		/* this needs to be an append so we keep things in the
		   order we got them, they will be sorted later */
		list_append(records, bg_record);
		/* this isn't a correct list so we need to set it later for
		   now we just used it to be the bp number */
		if(!used_nodes) {
			debug4("we didn't get a request list so we are "
			       "destroying this bp list");
			list_destroy(bg_record->bg_block_list);
			bg_record->bg_block_list = NULL;
		}
	} else {
		debug("adding a small block");
		/* if the ionode cnt for nodecards is 0 then don't
		   allow a nodecard allocation 
		*/
		if(!bluegene_nodecard_ionode_cnt) {
			if(blockreq->nodecards) 
				fatal("There is an error in your "
				      "bluegene.conf file.\n"
				      "Can't create a 32 node block with "
				      "Numpsets=%u. (Try setting it to 64)",
				      bluegene_numpsets);
		}

		if(blockreq->nodecards==0 && blockreq->quarters==0) {
			info("No specs given for this small block, "
			     "I am spliting this block into 4 quarters");
			blockreq->quarters=4;
		}		

		i = (blockreq->nodecards*bluegene_nodecard_node_cnt) + 
			(blockreq->quarters*bluegene_quarter_node_cnt);
		if(i != bluegene_bp_node_cnt)
			fatal("There is an error in your bluegene.conf file.\n"
			      "I am unable to request %d nodes consisting of "
			      "%u nodecards and\n%u quarters in one "
			      "base partition with %u nodes.", 
			      i, bluegene_bp_node_cnt, 
			      blockreq->nodecards, blockreq->quarters);
		small_count = blockreq->nodecards+blockreq->quarters; 
		
		/* Automatically create 4-way split if 
		 * conn_type == SELECT_SMALL in bluegene.conf
		 * Here we go through each node listed and do the same thing
		 * for each node.
		 */
		itr = list_iterator_create(bg_record->bg_block_list);
		while ((ba_node = list_next(itr)) != NULL) {
			/* break base partition up into 16 parts */
			small_size = 16;
			node_cnt = 0;
			quarter = 0;
			nodecard = 0;
			for(i=0; i<small_count; i++) {
				if(i == blockreq->nodecards) {
					/* break base partition 
					   up into 4 parts */
					small_size = 4;
				}
									
				if(small_size == 4)
					nodecard = (uint16_t)NO_VAL;
				else
					nodecard = i%4; 
				found_record = create_small_record(bg_record,
								   quarter,
								   nodecard);
								 
				/* this needs to be an append so we
				   keep things in the order we got
				   them, they will be sorted later */
				list_append(records, found_record);
				node_cnt += bluegene_bp_node_cnt/small_size;
				if(node_cnt == 128) {
					node_cnt = 0;
					quarter++;
				}
			}
		}
		list_iterator_destroy(itr);
		destroy_bg_record(bg_record);
	} 
	return SLURM_SUCCESS;
}

extern int format_node_name(bg_record_t *bg_record, char *buf, int buf_size)
{
	if(bg_record->ionodes) {
		snprintf(buf, buf_size, "%s[%s]",
			bg_record->nodes,
			bg_record->ionodes);
	} else {
		snprintf(buf, buf_size, "%s", bg_record->nodes);
	}
	return SLURM_SUCCESS;
}

/************************* local functions ***************************/

#ifdef HAVE_BG
static int _addto_node_list(bg_record_t *bg_record, int *start, int *end)
{
	int node_count=0;
	int x,y,z;
	char node_name_tmp[255];
	ba_node_t *ba_node = NULL;

	if ((start[X] < 0) || (start[Y] < 0) || (start[Z] < 0)) {
		fatal("bluegene.conf starting coordinate is invalid: %d%d%d",
		      start[X], start[Y], start[Z]);
	}
	if ((end[X] >= DIM_SIZE[X]) || (end[Y] >= DIM_SIZE[Y])
	||  (end[Z] >= DIM_SIZE[Z])) {
		fatal("bluegene.conf matrix size exceeds space defined in " 
		      "slurm.conf %c%c%cx%d%d%d => %c%c%c",
		      alpha_num[start[X]], alpha_num[start[Y]],
		      alpha_num[start[Z]], 
		      end[X], end[Y], end[Z], 
		      alpha_num[DIM_SIZE[X]], alpha_num[DIM_SIZE[Y]], 
		      alpha_num[DIM_SIZE[Z]]);
	}
	debug3("adding bps: %c%c%cx%c%c%c",
	       alpha_num[start[X]], alpha_num[start[Y]], alpha_num[start[Z]],
	       alpha_num[end[X]], alpha_num[end[Y]], alpha_num[end[Z]]);
	debug3("slurm.conf:    %c%c%c",
	       alpha_num[DIM_SIZE[X]], alpha_num[DIM_SIZE[Y]],
	       alpha_num[DIM_SIZE[Z]]); 
	
	for (x = start[X]; x <= end[X]; x++) {
		for (y = start[Y]; y <= end[Y]; y++) {
			for (z = start[Z]; z <= end[Z]; z++) {
				slurm_conf_lock();
				snprintf(node_name_tmp, sizeof(node_name_tmp),
					 "%s%c%c%c", 
					 slurmctld_conf.node_prefix,
					 alpha_num[x], alpha_num[y],
					 alpha_num[z]);		
				slurm_conf_unlock();
				ba_node = ba_copy_node(
					&ba_system_ptr->grid[x][y][z]);
				ba_node->used = 1;
				list_append(bg_record->bg_block_list, ba_node);
				node_count++;
			}
		}
	}
	return node_count;
}

static int _ba_node_cmpf_inc(ba_node_t *node_a, ba_node_t *node_b)
{
	if (node_a->coord[X] < node_b->coord[X])
		return -1;
	else if (node_a->coord[X] > node_b->coord[X])
		return 1;
	
	if (node_a->coord[Y] < node_b->coord[Y])
		return -1;
	else if (node_a->coord[Y] > node_b->coord[Y])
		return 1;

	if (node_a->coord[Z] < node_b->coord[Z])
		return -1;
	else if (node_a->coord[Z] > node_b->coord[Z])
		return 1;

	error("You have the node %c%c%c in the list twice",
	      alpha_num[node_a->coord[X]],
	      alpha_num[node_a->coord[Y]],
	      alpha_num[node_a->coord[Z]]); 
	return 0;
}
#endif //HAVE_BG


