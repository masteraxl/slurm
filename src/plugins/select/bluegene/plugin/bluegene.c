/*****************************************************************************\
 *  bluegene.c - blue gene node configuration processing module. 
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <auble1@llnl.gov> et. al.
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
#include "defined_block.h"
#include <stdio.h>

#define MMCS_POLL_TIME 120	/* poll MMCS for down switches and nodes 
				 * every 120 secs */
#define BG_POLL_TIME 0	        /* poll bg blocks every 3 secs */

#define _DEBUG 0

char* bg_conf = NULL;

/* Global variables */
rm_BGL_t *bg = NULL;

List bg_list = NULL;			/* total list of bg_record entries */
List bg_curr_block_list = NULL;  	/* current bg blocks in bluegene.conf*/
List bg_job_block_list = NULL;  	/* jobs running in these blocks */
List bg_booted_block_list = NULL;  	/* blocks that are booted */
List bg_freeing_list = NULL;  	        /* blocks that being freed */

List bg_blrtsimage_list = NULL;
List bg_linuximage_list = NULL;
List bg_mloaderimage_list = NULL;
List bg_ramdiskimage_list = NULL;
char *default_blrtsimage = NULL, *default_linuximage = NULL;
char *default_mloaderimage = NULL, *default_ramdiskimage = NULL;
char *bridge_api_file = NULL; 
bg_layout_t bluegene_layout_mode = NO_VAL;
uint16_t bluegene_numpsets = 0;
uint16_t bluegene_bp_node_cnt = 0;
uint16_t bluegene_quarter_node_cnt = 0;
uint16_t bluegene_quarter_ionode_cnt = 0;
uint16_t bluegene_nodecard_node_cnt = 0;
uint16_t bluegene_nodecard_ionode_cnt = 0;
uint16_t bridge_api_verb = 0;
bool agent_fini = false;
time_t last_bg_update;
pthread_mutex_t block_state_mutex = PTHREAD_MUTEX_INITIALIZER;
int num_block_to_free = 0;
int num_block_freed = 0;
int blocks_are_created = 0;
int num_unused_cpus = 0;

pthread_mutex_t freed_cnt_mutex = PTHREAD_MUTEX_INITIALIZER;
List bg_free_block_list = NULL;  	/* blocks to be deleted */
List bg_destroy_block_list = NULL;       /* blocks to be destroyed */
int free_cnt = 0;
int destroy_cnt = 0;

#ifndef HAVE_BG_FILES
# if BA_SYSTEM_DIMENSIONS==3
int max_dim[BA_SYSTEM_DIMENSIONS] = { 0, 0, 0 };
# else
int max_dim[BA_SYSTEM_DIMENSIONS] = { 0 };
# endif
#endif


static void _set_bg_lists();
static int  _validate_config_nodes(List *bg_found_block_list);
static int  _bg_record_cmpf_inc(bg_record_t *rec_a, bg_record_t *rec_b);
static int _delete_old_blocks(List bg_found_block_list);
static char *_get_bg_conf(void);
static int  _reopen_bridge_log(void);

/* Initialize all plugin variables */
extern int init_bg(void)
{
#ifdef HAVE_BG_FILES
	int rc;
	rm_size3D_t bp_size;
	
	info("Attempting to contact MMCS");
	if ((rc = bridge_get_bg(&bg)) != STATUS_OK) {
		fatal("init_bg: rm_get_BGL(): %s", bg_err_str(rc));
		return SLURM_ERROR;
	}
	
	if ((rc = bridge_get_data(bg, RM_Msize, &bp_size)) != STATUS_OK) {
		fatal("init_bg: rm_get_data(): %s", bg_err_str(rc));
		return SLURM_ERROR;
	}
	verbose("BlueGene configured with %d x %d x %d base blocks",
		bp_size.X, bp_size.Y, bp_size.Z);
	DIM_SIZE[X]=bp_size.X;
	DIM_SIZE[Y]=bp_size.Y;
	DIM_SIZE[Z]=bp_size.Z;
#endif
	ba_init(NULL);
	info("BlueGene plugin loaded successfully");

	return SLURM_SUCCESS;
}

/* Purge all plugin variables */
extern void fini_bg(void)
{
#ifdef HAVE_BG_FILES
	int rc;
#endif
	_set_bg_lists();
	
	if (bg_list) {
		list_destroy(bg_list);
		bg_list = NULL;
	}	
	if (bg_curr_block_list) {
		list_destroy(bg_curr_block_list);
		bg_curr_block_list = NULL;
	}	
	if (bg_job_block_list) {
		list_destroy(bg_job_block_list);
		bg_job_block_list = NULL;
		num_unused_cpus = 0;
	}
	if (bg_booted_block_list) {
		list_destroy(bg_booted_block_list);
		bg_booted_block_list = NULL;
	}

	/* wait for the free threads to finish up don't destroy the
	 * bg_free_block_list here */
	while(free_cnt > 0)
		usleep(1000);
	/* wait for the destroy threads to finish up don't destroy the
	 * bg_destroy_block_list here */
	while(destroy_cnt > 0)
		usleep(1000);
	
	if(bg_blrtsimage_list) {
		list_destroy(bg_blrtsimage_list);
		bg_blrtsimage_list = NULL;
	}
	
	if(bg_linuximage_list) {
		list_destroy(bg_linuximage_list);
		bg_linuximage_list = NULL;
	}
	
	if(bg_mloaderimage_list) {
		list_destroy(bg_mloaderimage_list);
		bg_mloaderimage_list = NULL;
	}

	if(bg_ramdiskimage_list) {
		list_destroy(bg_ramdiskimage_list);
		bg_ramdiskimage_list = NULL;
	}

	xfree(default_blrtsimage);
	xfree(default_linuximage);
	xfree(default_mloaderimage);
	xfree(default_ramdiskimage);
	xfree(bridge_api_file);
	xfree(bg_conf);
	
#ifdef HAVE_BG_FILES
	if(bg)
		if ((rc = bridge_free_bg(bg)) != STATUS_OK)
			error("bridge_free_BGL(): %s", bg_err_str(rc));
#endif	
	ba_fini();
}

/* 
 * block_state_mutex should be locked before calling this function
 */
extern bool blocks_overlap(bg_record_t *rec_a, bg_record_t *rec_b)
{
	bitstr_t *my_bitmap = NULL;

	if((rec_a->bp_count > 1) && (rec_b->bp_count > 1)) {
		/* Test for conflicting passthroughs */
		reset_ba_system(false);
		check_and_set_node_list(rec_a->bg_block_list);
		if(check_and_set_node_list(rec_b->bg_block_list)
		   == SLURM_ERROR) 
			return true;
	}
	
	
	my_bitmap = bit_copy(rec_a->bitmap);
	bit_and(my_bitmap, rec_b->bitmap);
	if (bit_ffs(my_bitmap) == -1) {
		FREE_NULL_BITMAP(my_bitmap);
		return false;
	}
	FREE_NULL_BITMAP(my_bitmap);
		
	if(rec_a->quarter != (uint16_t) NO_VAL) {
		if(rec_b->quarter == (uint16_t) NO_VAL)
			return true;
		else if(rec_a->quarter != rec_b->quarter)
			return false;
		if(rec_a->nodecard != (uint16_t) NO_VAL) {
			if(rec_b->nodecard == (uint16_t) NO_VAL)
				return true;
			else if(rec_a->nodecard
				!= rec_b->nodecard)
				return false;
		}
	}
	
	return true;
}

extern int remove_all_users(char *bg_block_id, char *user_name) 
{
	int returnc = REMOVE_USER_NONE;
#ifdef HAVE_BG_FILES
	char *user;
	rm_partition_t *block_ptr = NULL;
	int rc, i, user_count;

	if ((rc = bridge_get_block(bg_block_id,  &block_ptr)) != STATUS_OK) {
		if(rc == INCONSISTENT_DATA
		   && bluegene_layout_mode == LAYOUT_DYNAMIC)
			return REMOVE_USER_FOUND;
			
		error("bridge_get_block(%s): %s", 
		      bg_block_id, 
		      bg_err_str(rc));
		return REMOVE_USER_ERR;
	}	
	
	if((rc = bridge_get_data(block_ptr, RM_PartitionUsersNum, 
				 &user_count)) 
	   != STATUS_OK) {
		error("bridge_get_data(RM_PartitionUsersNum): %s", 
		      bg_err_str(rc));
		returnc = REMOVE_USER_ERR;
		user_count = 0;
	} else
		debug2("got %d users for %s",user_count, bg_block_id);
	for(i=0; i<user_count; i++) {
		if(i) {
			if ((rc = bridge_get_data(block_ptr, 
						  RM_PartitionNextUser, 
						  &user)) 
			    != STATUS_OK) {
				error("bridge_get_data"
				      "(RM_PartitionNextUser): %s", 
				      bg_err_str(rc));
				returnc = REMOVE_USER_ERR;
				break;
			}
		} else {
			if ((rc = bridge_get_data(block_ptr, 
						  RM_PartitionFirstUser, 
						  &user)) 
			    != STATUS_OK) {
				error("bridge_get_data"
				      "(RM_PartitionFirstUser): %s",
				      bg_err_str(rc));
				returnc = REMOVE_USER_ERR;
				break;
			}
		}
		if(!user) {
			error("No user was returned from database");
			continue;
		}
		if(!strcmp(user, slurmctld_conf.slurm_user_name)) {
			free(user);
			continue;
		}

		if(user_name) {
			if(!strcmp(user, user_name)) {
				returnc = REMOVE_USER_FOUND;
				free(user);
				continue;
			}
		}
		
		info("Removing user %s from Block %s", user, bg_block_id);
		if ((rc = bridge_remove_block_user(bg_block_id, user)) 
		    != STATUS_OK) {
			debug("user %s isn't on block %s",
			      user, 
			      bg_block_id);
		}
		free(user);
	}
	if ((rc = bridge_free_block(block_ptr)) != STATUS_OK) {
		error("bridge_free_block(): %s", bg_err_str(rc));
	}
#endif
	return returnc;
}

/* if SLURM_ERROR you will need to fail the job with
   slurm_fail_job(bg_record->job_running);
*/

extern int set_block_user(bg_record_t *bg_record) 
{
	int rc = 0;
	debug("resetting the boot state flag and "
	      "counter for block %s.",
	      bg_record->bg_block_id);
	bg_record->boot_state = 0;
	bg_record->boot_count = 0;
	slurm_conf_lock();
	if((rc = update_block_user(bg_record, 1)) == 1) {
		last_bg_update = time(NULL);
		rc = SLURM_SUCCESS;
	} else if (rc == -1) {
		error("Unable to add user name to block %s. "
		      "Cancelling job.",
		      bg_record->bg_block_id);
		rc = SLURM_ERROR;
	}	
	xfree(bg_record->target_name);
	bg_record->target_name = 
		xstrdup(slurmctld_conf.slurm_user_name);
	slurm_conf_unlock();	
	return rc;
}

extern char* convert_conn_type(rm_connection_type_t conn_type)
{
	switch (conn_type) {
	case (SELECT_MESH): 
		return "MESH"; 
	case (SELECT_TORUS): 
		return "TORUS"; 
	case (SELECT_SMALL): 
		return "SMALL"; 
	case (SELECT_NAV):
		return "NAV";
	default:
		break;
	}
	return "";
}

extern char* convert_node_use(rm_partition_mode_t pt)
{
	switch (pt) {
	case (SELECT_COPROCESSOR_MODE): 
		return "COPROCESSOR"; 
	case (SELECT_VIRTUAL_NODE_MODE): 
		return "VIRTUAL"; 
	default:
		break;
	}
	return "";
}

/** 
 * sort the partitions by increasing size
 */
extern void sort_bg_record_inc_size(List records){
	if (records == NULL)
		return;
	list_sort(records, (ListCmpF) _bg_record_cmpf_inc);
	last_bg_update = time(NULL);
}

/*
 * bluegene_agent - detached thread periodically updates status of
 * bluegene nodes. 
 * 
 * NOTE: I don't grab any locks here because slurm_drain_nodes grabs
 * the necessary locks.
 */
extern void *bluegene_agent(void *args)
{
	static time_t last_mmcs_test;
	static time_t last_bg_test;
	int rc;

	last_mmcs_test = time(NULL) + MMCS_POLL_TIME;
	last_bg_test = time(NULL) + BG_POLL_TIME;
	while (!agent_fini) {
		time_t now = time(NULL);

		if (difftime(now, last_bg_test) >= BG_POLL_TIME) {
			if (agent_fini)		/* don't bother */
				return NULL;	/* quit now */
			if(blocks_are_created) {
				last_bg_test = now;
				if((rc = update_block_list()) == 1) {
					slurm_mutex_lock(&block_state_mutex);
					last_bg_update = now;
					slurm_mutex_unlock(&block_state_mutex);
				} else if(rc == -1)
					error("Error with update_block_list");
				if(bluegene_layout_mode == LAYOUT_DYNAMIC) {
					if((rc = update_freeing_block_list())
					   == -1)
						error("Error with "
						      "update_block_list 2");
				}
			}
		}

		if (difftime(now, last_mmcs_test) >= MMCS_POLL_TIME) {
			if (agent_fini)		/* don't bother */
				return NULL;	/* quit now */
			last_mmcs_test = now;
			test_mmcs_failures();	/* can run for a while */
		}	
				
		sleep(1);
	}
	return NULL;
}

/* must set the protecting mutex if any before this function is called */

extern int remove_from_bg_list(List my_bg_list, bg_record_t *bg_record)
{
	bg_record_t *found_record = NULL;
	ListIterator itr;
	int rc = SLURM_ERROR;

	if(!bg_record)
		return rc;

	//slurm_mutex_lock(&block_state_mutex);	
	itr = list_iterator_create(my_bg_list);
	while ((found_record = (bg_record_t *) list_next(itr)) != NULL) {
		if(found_record)
			if(bg_record == found_record) {
				list_remove(itr);
				rc = SLURM_SUCCESS;
				break;
			}
	}
	list_iterator_destroy(itr);
	//slurm_mutex_unlock(&block_state_mutex);

	return rc;
}

/* This is here to remove from the orignal list when dealing with
 * copies like above all locks need to be set.  This function does not
 * free anything you must free it when you are done */
extern bg_record_t *find_and_remove_org_from_bg_list(List my_list, 
						     bg_record_t *bg_record)
{
	ListIterator itr = list_iterator_create(my_list);
	bg_record_t *found_record = NULL;

	while ((found_record = (bg_record_t *) list_next(itr)) != NULL) {
		/* check for full node bitmap compare */
		if(bit_equal(bg_record->bitmap, found_record->bitmap)
		   && bit_equal(bg_record->ionode_bitmap,
				found_record->ionode_bitmap)) {
			if(!strcmp(bg_record->bg_block_id,
				   found_record->bg_block_id)) {
				list_remove(itr);
				debug2("got the block");
				break;
			}
		}
	}
	list_iterator_destroy(itr);
	return found_record;
}

/* This is here to remove from the orignal list when dealing with
 * copies like above all locks need to be set */
extern bg_record_t *find_org_in_bg_list(List my_list, bg_record_t *bg_record)
{
	ListIterator itr = list_iterator_create(my_list);
	bg_record_t *found_record = NULL;

	while ((found_record = (bg_record_t *) list_next(itr)) != NULL) {
		/* check for full node bitmap compare */
		if(bit_equal(bg_record->bitmap, found_record->bitmap)
		   && bit_equal(bg_record->ionode_bitmap,
				found_record->ionode_bitmap)) {
			
			if(!strcmp(bg_record->bg_block_id,
				   found_record->bg_block_id)) {
				debug2("got the block");
				break;
			}
		}
	}
	list_iterator_destroy(itr);
	return found_record;
}

extern int bg_free_block(bg_record_t *bg_record)
{
#ifdef HAVE_BG_FILES
	int rc;
#endif
	if(!bg_record) {
		error("bg_free_block: there was no bg_record");
		return SLURM_ERROR;
	}
	
	while (1) {
		if(!bg_record) {
			error("bg_free_block: there was no bg_record");
			return SLURM_ERROR;
		}
		
		slurm_mutex_lock(&block_state_mutex);			
		if (bg_record->state != (uint16_t)NO_VAL
		    && bg_record->state != RM_PARTITION_FREE 
		    && bg_record->state != RM_PARTITION_DEALLOCATING) {
#ifdef HAVE_BG_FILES
			debug2("bridge_destroy %s",bg_record->bg_block_id);
			
			rc = bridge_destroy_block(bg_record->bg_block_id);
			if (rc != STATUS_OK) {
				if(rc == PARTITION_NOT_FOUND) {
					debug("block %s is not found",
					      bg_record->bg_block_id);
					break;
				} else if(rc == INCOMPATIBLE_STATE) {
					debug2("bridge_destroy_partition"
					       "(%s): %s State = %d",
					       bg_record->bg_block_id, 
					       bg_err_str(rc), 
					       bg_record->state);
				} else {
					error("bridge_destroy_partition"
					      "(%s): %s State = %d",
					      bg_record->bg_block_id, 
					      bg_err_str(rc), 
					      bg_record->state);
				}
			}
#else
			bg_record->state = RM_PARTITION_FREE;	
#endif
		}
		
		if ((bg_record->state == RM_PARTITION_FREE)
		    ||  (bg_record->state == RM_PARTITION_ERROR)) {
			break;
		}
		slurm_mutex_unlock(&block_state_mutex);			
		sleep(3);
	}
	remove_from_bg_list(bg_booted_block_list, bg_record);
	slurm_mutex_unlock(&block_state_mutex);			
		
	return SLURM_SUCCESS;
}

/* Free multiple blocks in parallel */
extern void *mult_free_block(void *args)
{
	bg_record_t *bg_record = NULL;
	
	slurm_mutex_lock(&freed_cnt_mutex);
	if ((bg_freeing_list == NULL) 
	    && ((bg_freeing_list = list_create(destroy_bg_record)) == NULL))
		fatal("malloc failure in bg_freeing_list");
	slurm_mutex_unlock(&freed_cnt_mutex);
	
	/*
	 * Don't just exit when there is no work left. Creating 
	 * pthreads from within a dynamically linked object (plugin)
	 * causes large memory leaks on some systems that seem 
	 * unavoidable even from detached pthreads.
	 */
	while (!agent_fini) {
		slurm_mutex_lock(&freed_cnt_mutex);
		bg_record = list_dequeue(bg_free_block_list);
		slurm_mutex_unlock(&freed_cnt_mutex);
		if (!bg_record) {
			usleep(100000);
			continue;
		}
		if(bg_record->job_ptr) {
			info("We are freeing a block (%s) that "
			     "has job %u(%u), This should never happen.\n",
			     bg_record->bg_block_id,
			     bg_record->job_ptr->job_id, 
			     bg_record->job_running);
			term_jobs_on_block(bg_record->bg_block_id);
		}
		debug("freeing the block %s.", bg_record->bg_block_id);
		bg_free_block(bg_record);	
		debug("done\n");
		slurm_mutex_lock(&freed_cnt_mutex);
		num_block_freed++;
		slurm_mutex_unlock(&freed_cnt_mutex);
	}
	slurm_mutex_lock(&freed_cnt_mutex);
	free_cnt--;
	if(bg_freeing_list) {
		list_destroy(bg_freeing_list);
		bg_freeing_list = NULL;
	}
	if(free_cnt == 0) {
		list_destroy(bg_free_block_list);
		bg_free_block_list = NULL;
	}
	slurm_mutex_unlock(&freed_cnt_mutex);
	return NULL;
}

/* destroy multiple blocks in parallel */
extern void *mult_destroy_block(void *args)
{
	bg_record_t *bg_record = NULL;

#ifdef HAVE_BG_FILES
	int rc;
#endif
	slurm_mutex_lock(&freed_cnt_mutex);
	if ((bg_freeing_list == NULL) 
	    && ((bg_freeing_list = list_create(destroy_bg_record)) == NULL))
		fatal("malloc failure in bg_freeing_list");
	slurm_mutex_unlock(&freed_cnt_mutex);
	
	/*
	 * Don't just exit when there is no work left. Creating 
	 * pthreads from within a dynamically linked object (plugin)
	 * causes large memory leaks on some systems that seem 
	 * unavoidable even from detached pthreads.
	 */
	while (!agent_fini) {
		slurm_mutex_lock(&freed_cnt_mutex);
		bg_record = list_dequeue(bg_destroy_block_list);
		slurm_mutex_unlock(&freed_cnt_mutex);
		if (!bg_record) {
			usleep(100000);
			continue;
		}
		slurm_mutex_lock(&block_state_mutex);
		remove_from_bg_list(bg_list, bg_record);
		list_push(bg_freeing_list, bg_record);
		
		/* 
		 * we only are sorting this so when we send it to a
		 * tool such as smap it will be in a nice order
		 */
		sort_bg_record_inc_size(bg_freeing_list);
		if(remove_from_bg_list(bg_job_block_list, bg_record) 
		   == SLURM_SUCCESS) {
			num_unused_cpus += 
				bg_record->bp_count*bg_record->cpus_per_bp;
		}
		slurm_mutex_unlock(&block_state_mutex);
		debug3("removing the jobs on block %s\n",
		       bg_record->bg_block_id);
		term_jobs_on_block(bg_record->bg_block_id);
		
		debug2("destroying %s", (char *)bg_record->bg_block_id);
		if(bg_free_block(bg_record) == SLURM_ERROR) {
			debug("there was an error");
			goto already_here;
		}
		debug2("done destroying");
		slurm_mutex_lock(&block_state_mutex);
		remove_from_bg_list(bg_freeing_list, bg_record);
		slurm_mutex_unlock(&block_state_mutex);
								
#ifdef HAVE_BG_FILES
		debug2("removing from database %s", 
		       (char *)bg_record->bg_block_id);
		
		rc = bridge_remove_block(bg_record->bg_block_id);
		if (rc != STATUS_OK) {
			if(rc == PARTITION_NOT_FOUND) {
				debug("block %s is not found",
				      bg_record->bg_block_id);
			} else {
				error("1 rm_remove_partition(%s): %s",
				      bg_record->bg_block_id,
				      bg_err_str(rc));
			}
		} else
			debug2("done %s", 
			       (char *)bg_record->bg_block_id);
#endif
		destroy_bg_record(bg_record);
		debug2("destroyed");
		
	already_here:
		slurm_mutex_lock(&freed_cnt_mutex);
		num_block_freed++;
		slurm_mutex_unlock(&freed_cnt_mutex);
				
	}
	slurm_mutex_lock(&freed_cnt_mutex);
	destroy_cnt--;
	if(bg_freeing_list) {
		list_destroy(bg_freeing_list);
		bg_freeing_list = NULL;
	}
	if(destroy_cnt == 0) {
		list_destroy(bg_destroy_block_list);
		bg_destroy_block_list = NULL;
	}
	slurm_mutex_unlock(&freed_cnt_mutex);

	return NULL;
}

extern int free_block_list(List delete_list)
{
	bg_record_t *found_record = NULL;
	int retries;
	List *block_list = NULL;
	int *count = NULL;
	pthread_attr_t attr_agent;
	pthread_t thread_agent;

	if(!delete_list || !list_count(delete_list))
		return SLURM_SUCCESS;

	/* set up which list to push onto */
	if(bluegene_layout_mode == LAYOUT_DYNAMIC) {
		block_list = &bg_destroy_block_list;
		count = &destroy_cnt;
	} else {
		block_list = &bg_free_block_list;
		count = &free_cnt;
	}
	
	slurm_mutex_lock(&freed_cnt_mutex);
	
	if ((*block_list == NULL) 
	    && ((*block_list = list_create(NULL)) == NULL))
		fatal("malloc failure in free_block_list");
	
	while ((found_record = (bg_record_t*)list_pop(delete_list)) != NULL) {
		/* push job onto queue in a FIFO */
		debug3("adding %s to be freed", found_record->bg_block_id);
		if (list_push(*block_list, found_record) == NULL)
			fatal("malloc failure in _block_op/list_push");
		
		/* already running MAX_AGENTS we don't really need more 
		   since they don't end until we shut down the controller */
		if (*count > MAX_AGENT_COUNT) 
			continue;
		
		(*count)++;
		
		slurm_attr_init(&attr_agent);
		if (pthread_attr_setdetachstate(
			    &attr_agent, 
			    PTHREAD_CREATE_DETACHED))
			error("pthread_attr_setdetachstate error %m");
		retries = 0;
		if(bluegene_layout_mode == LAYOUT_DYNAMIC) {
			while (pthread_create(&thread_agent, 
					      &attr_agent, 
					      mult_destroy_block,
					      NULL)) {
				error("pthread_create "
				      "error %m");
				if (++retries > MAX_PTHREAD_RETRIES)
					fatal("Can't create "
					      "pthread");
				/* sleep and retry */
				usleep(1000);	
			}
		} else {
			while (pthread_create(&thread_agent, 
					      &attr_agent, 
					      mult_free_block, 
					      NULL)) {
				error("pthread_create "
				      "error %m");
				if (++retries > MAX_PTHREAD_RETRIES)
					fatal("Can't create "
					      "pthread");
				/* sleep and retry */
				usleep(1000);	
			}
		}
		slurm_attr_destroy(&attr_agent);
	}
	slurm_mutex_unlock(&freed_cnt_mutex);
	return SLURM_SUCCESS;
}

/*
 * Read and process the bluegene.conf configuration file so to interpret what
 * blocks are static/dynamic, torus/mesh, etc.
 */
extern int read_bg_conf(void)
{
	int i;
	int count = 0;
	s_p_hashtbl_t *tbl = NULL;
	char *layout = NULL;
	blockreq_t **blockreq_array = NULL;
	image_t **image_array = NULL;
	image_t *image = NULL;
	static time_t last_config_update = (time_t) 0;
	struct stat config_stat;
	ListIterator itr = NULL;
	/* found bg blocks already on system */
	List bg_found_block_list = list_create(NULL);
	
	debug("Reading the bluegene.conf file");

	/* check if config file has changed */
	if (!bg_conf)
		bg_conf = _get_bg_conf();
	if (stat(bg_conf, &config_stat) < 0)
		fatal("can't stat bluegene.conf file %s: %m", bg_conf);
	if (last_config_update) {
		_reopen_bridge_log();
		if(last_config_update == config_stat.st_mtime)
			debug("%s unchanged", bg_conf);
		else {
			info("Restart slurmctld for %s changes to take effect", 
			     bg_conf);
		}
		last_config_update = config_stat.st_mtime; 
		return SLURM_SUCCESS;
	}
	last_config_update = config_stat.st_mtime; 

	/* initialization */
	/* bg_conf defined in bg_node_alloc.h */
	tbl = s_p_hashtbl_create(bg_conf_file_options);
	
	if(s_p_parse_file(tbl, bg_conf) == SLURM_ERROR)
		fatal("something wrong with opening/reading bluegene "
		      "conf file");
	
	_set_bg_lists();	
	if (s_p_get_array((void ***)&image_array, 
			  &count, "AltBlrtsImage", tbl)) {
		for (i = 0; i < count; i++) {
			list_append(bg_blrtsimage_list, image_array[i]);
			image_array[i] = NULL;
		}
	}
	if (!s_p_get_string(&default_blrtsimage, "BlrtsImage", tbl)) {
		if(!list_count(bg_blrtsimage_list))
			fatal("BlrtsImage not configured "
			      "in bluegene.conf");
		itr = list_iterator_create(bg_blrtsimage_list);
		image = list_next(itr);
		image->def = true;
		list_iterator_destroy(itr);
		default_blrtsimage = xstrdup(image->name);
		info("Warning: using %s as the default BlrtsImage.  "
		     "If this isn't correct please set BlrtsImage",
		     default_blrtsimage); 
	} else {
		debug3("default BlrtsImage %s", default_blrtsimage);
		image = xmalloc(sizeof(image_t));
		image->name = xstrdup(default_blrtsimage);
		image->def = true;
		image->groups = NULL;
		/* we want it to be first */
		list_push(bg_blrtsimage_list, image);
	}
		
	if (s_p_get_array((void ***)&image_array, 
			  &count, "AltLinuxImage", tbl)) {
		for (i = 0; i < count; i++) {
			list_append(bg_linuximage_list, image_array[i]);
			image_array[i] = NULL;
		}
	}
	if (!s_p_get_string(&default_linuximage, "LinuxImage", tbl)) {
		if(!list_count(bg_linuximage_list))
			fatal("LinuxImage not configured "
			      "in bluegene.conf");
		itr = list_iterator_create(bg_linuximage_list);
		image = list_next(itr);
		image->def = true;
		list_iterator_destroy(itr);
		default_linuximage = xstrdup(image->name);
		info("Warning: using %s as the default LinuxImage.  "
		     "If this isn't correct please set LinuxImage",
		     default_linuximage); 
	} else {
		debug3("default LinuxImage %s", default_linuximage);
		image = xmalloc(sizeof(image_t));
		image->name = xstrdup(default_linuximage);
		image->def = true;
		image->groups = NULL;
		/* we want it to be first */
		list_push(bg_linuximage_list, image);		
	}

	if (s_p_get_array((void ***)&image_array, 
			  &count, "AltMloaderImage", tbl)) {
		for (i = 0; i < count; i++) {
			list_append(bg_mloaderimage_list, image_array[i]);
			image_array[i] = NULL;
		}
	}
	if (!s_p_get_string(&default_mloaderimage,
			    "MloaderImage", tbl)) {
		if(!list_count(bg_mloaderimage_list))
			fatal("MloaderImage not configured "
			      "in bluegene.conf");
		itr = list_iterator_create(bg_mloaderimage_list);
		image = list_next(itr);
		image->def = true;
		list_iterator_destroy(itr);
		default_mloaderimage = xstrdup(image->name);
		info("Warning: using %s as the default MloaderImage.  "
		     "If this isn't correct please set MloaderImage",
		     default_mloaderimage); 
	} else {
		debug3("default MloaderImage %s", default_mloaderimage);
		image = xmalloc(sizeof(image_t));
		image->name = xstrdup(default_mloaderimage);
		image->def = true;
		image->groups = NULL;
		/* we want it to be first */
		list_push(bg_mloaderimage_list, image);		
	}

	if (s_p_get_array((void ***)&image_array, 
			  &count, "AltRamDiskImage", tbl)) {
		for (i = 0; i < count; i++) {
			list_append(bg_ramdiskimage_list, image_array[i]);
			image_array[i] = NULL;
		}
	}
	if (!s_p_get_string(&default_ramdiskimage,
			    "RamDiskImage", tbl)) {
		if(!list_count(bg_ramdiskimage_list))
			fatal("RamDiskImage not configured "
			      "in bluegene.conf");
		itr = list_iterator_create(bg_ramdiskimage_list);
		image = list_next(itr);
		image->def = true;
		list_iterator_destroy(itr);
		default_ramdiskimage = xstrdup(image->name);
		info("Warning: using %s as the default RamDiskImage.  "
		     "If this isn't correct please set RamDiskImage",
		     default_ramdiskimage); 
	} else {
		debug3("default RamDiskImage %s", default_ramdiskimage);
		image = xmalloc(sizeof(image_t));
		image->name = xstrdup(default_ramdiskimage);
		image->def = true;
		image->groups = NULL;
		/* we want it to be first */
		list_push(bg_ramdiskimage_list, image);		
	}

	if (!s_p_get_uint16(&bluegene_numpsets, "Numpsets", tbl))
		fatal("Warning: Numpsets not configured in bluegene.conf");
	if (!s_p_get_uint16(&bridge_api_verb, "BridgeAPIVerbose", tbl))
		info("Warning: BridgeAPIVerbose not configured "
		     "in bluegene.conf");
	if (!s_p_get_string(&bridge_api_file, "BridgeAPILogFile", tbl)) 
		info("BridgeAPILogFile not configured in bluegene.conf");
	else
		_reopen_bridge_log();
	if (!s_p_get_string(&layout, "LayoutMode", tbl)) {
		info("Warning: LayoutMode was not specified in bluegene.conf "
		     "defaulting to STATIC partitioning");
		bluegene_layout_mode = LAYOUT_STATIC;
	} else {
		if(!strcasecmp(layout,"STATIC")) 
			bluegene_layout_mode = LAYOUT_STATIC;
		else if(!strcasecmp(layout,"OVERLAP")) 
			bluegene_layout_mode = LAYOUT_OVERLAP;
		else if(!strcasecmp(layout,"DYNAMIC")) 
			bluegene_layout_mode = LAYOUT_DYNAMIC;
		else {
			fatal("I don't understand this LayoutMode = %s", 
			      layout);
		}
		xfree(layout);
	}
	if (!s_p_get_uint16(
		    &bluegene_bp_node_cnt, "BasePartitionNodeCnt", tbl)) {
		error("BasePartitionNodeCnt not configured in bluegene.conf "
		      "defaulting to 512 as BasePartitionNodeCnt");
		bluegene_bp_node_cnt = 512;
		bluegene_quarter_node_cnt = 128;
	} else {
		if(bluegene_bp_node_cnt<=0)
			fatal("You should have more than 0 nodes "
			      "per base partition");

		bluegene_quarter_node_cnt = bluegene_bp_node_cnt/4;
	}

	if (!s_p_get_uint16(
		    &bluegene_nodecard_node_cnt, "NodeCardNodeCnt", tbl)) {
		error("NodeCardNodeCnt not configured in bluegene.conf "
		      "defaulting to 32 as NodeCardNodeCnt");
		bluegene_nodecard_node_cnt = 32;
	}
	
	if(bluegene_nodecard_node_cnt<=0)
		fatal("You should have more than 0 nodes per nodecard");

	if(bluegene_numpsets) {
		bluegene_quarter_ionode_cnt = bluegene_numpsets/4;
		bluegene_nodecard_ionode_cnt = bluegene_quarter_ionode_cnt/4;
		if((int)bluegene_nodecard_ionode_cnt < 1) {
			bluegene_nodecard_ionode_cnt = 0;
		}
	} else {
		fatal("your numpsets is 0");
	}

	/* add blocks defined in file */
	if(bluegene_layout_mode != LAYOUT_DYNAMIC) {
		if (!s_p_get_array((void ***)&blockreq_array, 
				   &count, "BPs", tbl)) {
			info("WARNING: no blocks defined in bluegene.conf, "
			     "only making full system block");
			create_full_system_block(NULL);
		}
		
		for (i = 0; i < count; i++) {
			add_bg_record(bg_list, NULL, blockreq_array[i]);
		}
	}
//#if 0	
	/* Check to see if the configs we have are correct */
	if (_validate_config_nodes(&bg_found_block_list) == SLURM_ERROR) { 
		_delete_old_blocks(bg_found_block_list);
	}
//#endif
	/* looking for blocks only I created */
	if(bluegene_layout_mode == LAYOUT_DYNAMIC) {
		init_wires();
		info("No blocks created until jobs are submitted");
	} else {
		if (create_defined_blocks(bluegene_layout_mode,
					  bg_found_block_list) 
		    == SLURM_ERROR) {
			/* error in creating the static blocks, so
			 * blocks referenced by submitted jobs won't
			 * correspond to actual slurm blocks.
			 */
			fatal("Error, could not create the static blocks");
			return SLURM_ERROR;
		}
	} 
	
	slurm_mutex_lock(&block_state_mutex);
	list_destroy(bg_curr_block_list);
	bg_curr_block_list = NULL;
	if(bg_found_block_list) {
		list_destroy(bg_found_block_list);
		bg_found_block_list = NULL;
	}
	last_bg_update = time(NULL);
	blocks_are_created = 1;
	sort_bg_record_inc_size(bg_list);
	slurm_mutex_unlock(&block_state_mutex);
	debug("Blocks have finished being created.");
	s_p_hashtbl_destroy(tbl);

	return SLURM_SUCCESS;
}


static void _set_bg_lists()
{
	slurm_mutex_lock(&block_state_mutex);
	if(bg_booted_block_list) 
		list_destroy(bg_booted_block_list);
	bg_booted_block_list = list_create(NULL);
	if(bg_job_block_list) 
		list_destroy(bg_job_block_list);
	bg_job_block_list = list_create(NULL);	
	num_unused_cpus = 
		DIM_SIZE[X] * DIM_SIZE[Y] * DIM_SIZE[Z] * procs_per_node;
	if(bg_curr_block_list)
		list_destroy(bg_curr_block_list);	
	bg_curr_block_list = list_create(destroy_bg_record);
	
	
	if(bg_list) 
		list_destroy(bg_list);
	bg_list = list_create(destroy_bg_record);

	slurm_mutex_unlock(&block_state_mutex);	
	
	if(bg_blrtsimage_list)
		list_destroy(bg_blrtsimage_list);
	bg_blrtsimage_list = list_create(destroy_image);
	if(bg_linuximage_list)
		list_destroy(bg_linuximage_list);
	bg_linuximage_list = list_create(destroy_image);
	if(bg_mloaderimage_list)
		list_destroy(bg_mloaderimage_list);
	bg_mloaderimage_list = list_create(destroy_image);
	if(bg_ramdiskimage_list)
		list_destroy(bg_ramdiskimage_list);
	bg_ramdiskimage_list = list_create(destroy_image);
	
}

/*
 * _validate_config_nodes - Match slurm configuration information with
 *                          current BG block configuration.
 * IN/OUT bg_found_block_list - if NULL is created and then any blocks
 *                              found on the system are then pushed on.
 * RET - SLURM_SUCCESS if they match, else an error 
 * code. Writes bg_block_id into bg_list records.
 */

static int _validate_config_nodes(List *bg_found_block_list)
{
	int rc = SLURM_ERROR;
#ifdef HAVE_BG_FILES
	bg_record_t* bg_record = NULL;	
	bg_record_t* init_bg_record = NULL;
	bg_record_t* full_system_bg_record = NULL;	
	int full_created = 0;
	ListIterator itr_conf;
	ListIterator itr_curr;
	rm_partition_mode_t node_use;
	char tmp_char[256];

	/* read current bg block info into bg_curr_block_list */
	if (read_bg_blocks() == SLURM_ERROR)
		return SLURM_ERROR;
	
	if(!bg_recover) 
		return SLURM_ERROR;

	if(!bg_curr_block_list)
		return SLURM_ERROR;

	itr_curr = list_iterator_create(bg_curr_block_list);
	while ((init_bg_record = list_next(itr_curr))) 
		if(init_bg_record->full_block) 
			full_system_bg_record = init_bg_record;	
	
	if(!*bg_found_block_list)
		(*bg_found_block_list) = list_create(NULL);

	
	itr_conf = list_iterator_create(bg_list);
	while ((bg_record = (bg_record_t*) list_next(itr_conf))) {
		/* translate hostlist to ranged 
		   string for consistent format
		   search here 
		*/
		node_use = SELECT_COPROCESSOR_MODE; 
		list_iterator_reset(itr_curr);
		while ((init_bg_record = list_next(itr_curr))) {
			if (strcasecmp(bg_record->nodes, 
				       init_bg_record->nodes))
				continue; /* wrong nodes */
			if (bg_record->conn_type 
			    != init_bg_record->conn_type)
				continue; /* wrong conn_type */
			if(bg_record->quarter !=
			   init_bg_record->quarter)
				continue; /* wrong quart */
			if(bg_record->nodecard !=
			   init_bg_record->nodecard)
				continue; /* wrong nodecard */
			if(bg_record->blrtsimage &&
			   strcasecmp(bg_record->blrtsimage,
				      init_bg_record->blrtsimage)) 
				continue;
			if(bg_record->linuximage &&
			   strcasecmp(bg_record->linuximage,
				      init_bg_record->linuximage))
				continue;
			if(bg_record->mloaderimage &&
			   strcasecmp(bg_record->mloaderimage,
				      init_bg_record->mloaderimage))
				continue;
			if(bg_record->ramdiskimage &&
			   strcasecmp(bg_record->ramdiskimage,
				      init_bg_record->ramdiskimage))
				continue;
		       			
			copy_bg_record(init_bg_record, 
				       bg_record);
			break;
		}
			
		if (!bg_record->bg_block_id) {
			format_node_name(bg_record, tmp_char,
					 sizeof(tmp_char));
			info("Block found in bluegene.conf to be "
			     "created: Nodes:%s", 
			     tmp_char);
			rc = SLURM_ERROR;
		} else {
			if(bg_record->full_block)
				full_created = 1;

			list_push(*bg_found_block_list, bg_record);
			format_node_name(bg_record, tmp_char,
					 sizeof(tmp_char));
			info("Existing: BlockID:%s Nodes:%s Conn:%s",
			     bg_record->bg_block_id, 
			     tmp_char,
			     convert_conn_type(bg_record->conn_type));
			if(((bg_record->state == RM_PARTITION_READY)
			    || (bg_record->state == RM_PARTITION_CONFIGURING))
			   && !block_exist_in_list(bg_booted_block_list, 
						   bg_record))
				list_push(bg_booted_block_list, bg_record);
		}
	}		
	list_iterator_destroy(itr_conf);
	list_iterator_destroy(itr_curr);
	if(bluegene_layout_mode == LAYOUT_DYNAMIC)
		goto finished;

	if(!full_created && full_system_bg_record) {
		bg_record = xmalloc(sizeof(bg_record_t));
		copy_bg_record(full_system_bg_record, bg_record);
		list_append(bg_list, bg_record);
		list_push(*bg_found_block_list, bg_record);
		format_node_name(bg_record, tmp_char, sizeof(tmp_char));
		info("Existing: BlockID:%s Nodes:%s Conn:%s",
		     bg_record->bg_block_id, 
		     tmp_char,
		     convert_conn_type(bg_record->conn_type));
		if(((bg_record->state == RM_PARTITION_READY)
		    || (bg_record->state == RM_PARTITION_CONFIGURING))
		   && !block_exist_in_list(bg_booted_block_list, 
					   bg_record))
			list_push(bg_booted_block_list, bg_record);
	}
		
finished:
	if(list_count(bg_list) == list_count(bg_curr_block_list))
		rc = SLURM_SUCCESS;
	
#endif

	return rc;
}

/* 
 * Comparator used for sorting blocks smallest to largest
 * 
 * returns: -1: rec_a >rec_b   0: rec_a == rec_b   1: rec_a < rec_b
 * 
 */
static int _bg_record_cmpf_inc(bg_record_t* rec_a, bg_record_t* rec_b)
{
	int size_a = rec_a->node_cnt;
	int size_b = rec_b->node_cnt;
	if (size_a < size_b)
		return -1;
	else if (size_a > size_b)
		return 1;
	if(rec_a->nodes && rec_b->nodes) {
		size_a = strcmp(rec_a->nodes, rec_b->nodes);
		if (size_a < 0)
			return -1;
		else if (size_a > 0)
			return 1;
	}
	if (rec_a->quarter < rec_b->quarter)
		return -1;
	else if (rec_a->quarter > rec_b->quarter)
		return 1;

	if(rec_a->nodecard < rec_b->nodecard)
		return -1;
	else if(rec_a->nodecard > rec_b->nodecard)
		return 1;

	return 0;
}

static int _delete_old_blocks(List bg_found_block_list)
{
#ifdef HAVE_BG_FILES
	ListIterator itr_curr, itr_found;
	bg_record_t *found_record = NULL, *init_record = NULL;
	pthread_attr_t attr_agent;
	pthread_t thread_agent;
	int retries;
	List bg_destroy_list = list_create(NULL);

	info("removing unspecified blocks");
	if(!bg_recover) {
		if(bg_curr_block_list) {
			itr_curr = list_iterator_create(bg_curr_block_list);
			while ((init_record = 
				(bg_record_t*)list_next(itr_curr))) {
				list_remove(itr_curr);
				list_push(bg_destroy_list, init_record);
			}
			list_iterator_destroy(itr_curr);
		} else {
			error("_delete_old_blocks: "
			      "no bg_curr_block_list 1");
			list_destroy(bg_destroy_list);
			return SLURM_ERROR;
		}
	} else {
		if(bg_curr_block_list) {
			itr_curr = list_iterator_create(bg_curr_block_list);
			while ((init_record = list_next(itr_curr))) {
				if(bg_found_block_list) {
					itr_found = list_iterator_create(
						bg_found_block_list);
					while ((found_record 
						= list_next(itr_found)) 
					       != NULL) {
						if (!strcmp(init_record->
							    bg_block_id, 
							    found_record->
							    bg_block_id)) {
							/* don't delete 
							   this one 
							*/
							break;	
						}
					}
					list_iterator_destroy(itr_found);
				} else {
					error("_delete_old_blocks: "
					      "no bg_found_block_list");
					list_iterator_destroy(itr_curr);
					list_destroy(bg_destroy_list);
					return SLURM_ERROR;
				}
				if(found_record == NULL) {
					list_remove(itr_curr);
					list_push(bg_destroy_list, 
						  init_record);
				}
			}		
			list_iterator_destroy(itr_curr);
		} else {
			error("_delete_old_blocks: "
			      "no bg_curr_block_list 2");
			list_destroy(bg_destroy_list);
			return SLURM_ERROR;
		}
	}

	slurm_mutex_lock(&freed_cnt_mutex);
	if ((bg_destroy_block_list == NULL) 
	    && ((bg_destroy_block_list = list_create(NULL)) == NULL))
		fatal("malloc failure in block_list");
		
	itr_curr = list_iterator_create(bg_destroy_list);
	while ((init_record = (bg_record_t*) list_next(itr_curr))) {
		list_push(bg_destroy_block_list, init_record);
		num_block_to_free++;
		if (destroy_cnt > MAX_AGENT_COUNT) 
			continue;
		
		destroy_cnt++;

		slurm_attr_init(&attr_agent);
		if (pthread_attr_setdetachstate(&attr_agent, 
						PTHREAD_CREATE_DETACHED))
			error("pthread_attr_setdetachstate error %m");
		
		retries = 0;
		while (pthread_create(&thread_agent, 
				      &attr_agent, 
				      mult_destroy_block, 
				      NULL)) {
			error("pthread_create "
			      "error %m");
			if (++retries > MAX_PTHREAD_RETRIES)
				fatal("Can't create "
				      "pthread");
			/* sleep and retry */
			usleep(1000);	
		}
		slurm_attr_destroy(&attr_agent);
	}
	list_iterator_destroy(itr_curr);
	slurm_mutex_unlock(&freed_cnt_mutex);
	list_destroy(bg_destroy_list);
		
	retries=30;
	while(num_block_to_free > num_block_freed) {
		update_freeing_block_list();
		if(retries==30) {
			info("Waiting for old blocks to be "
			     "freed.  Have %d of %d",
			     num_block_freed, 
			     num_block_to_free);
			retries=0;
		}
		retries++;
		sleep(1);
	}
	
	info("I am done deleting");
#endif	
	return SLURM_SUCCESS;
}

static char *_get_bg_conf(void)
{
	char *val = getenv("SLURM_CONF");
	char *rc;
	int i;

	if (!val)
		return xstrdup(BLUEGENE_CONFIG_FILE);

	/* Replace file name on end of path */
	i = strlen(val) - strlen("slurm.conf") + strlen("bluegene.conf") + 1;
	rc = xmalloc(i);
	strcpy(rc, val);
	val = strrchr(rc, (int)'/');
	if (val)	/* absolute path */
		val++;
	else		/* not absolute path */
		val = rc;
	strcpy(val, "bluegene.conf");
	return rc;
}

static int _reopen_bridge_log(void)
{
	int rc = SLURM_SUCCESS;

	if (bridge_api_file == NULL)
		return rc;
	
#ifdef HAVE_BG_FILES
	rc = bridge_set_log_params(bridge_api_file, bridge_api_verb);
#endif
	debug3("Bridge api file set to %s, verbose level %d\n", 
	       bridge_api_file, bridge_api_verb);
	
	return rc;
}

