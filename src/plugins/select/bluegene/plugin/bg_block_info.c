/*****************************************************************************\
 *  bg_block_info.c - bluegene block information from the db2 database.
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004-2006 The Regents of the University of California.
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
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  endif
#  if WITH_PTHREADS
#    include <pthread.h>
#  endif
#endif

#include <signal.h>
#include <unistd.h>

#include <slurm/slurm_errno.h>

#include <pwd.h>
#include <sys/types.h>
#include "src/common/hostlist.h"
#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/common/node_select.h"
#include "src/common/uid.h"
#include "src/common/xstring.h"
#include "src/slurmctld/proc_req.h"
#include "src/api/job_info.h"
#include "bluegene.h"

#define _DEBUG 0
#define RETRY_BOOT_COUNT 3

#ifdef HAVE_BG_FILES
static int  _block_is_deallocating(bg_record_t *bg_record);

static int _block_is_deallocating(bg_record_t *bg_record)
{
	int jobid = bg_record->job_running;
	char *user_name = NULL;
	slurm_conf_lock();
	user_name = xstrdup(slurmctld_conf.slurm_user_name);
	if(remove_all_users(bg_record->bg_block_id, NULL) 
	   == REMOVE_USER_ERR) {
		error("Something happened removing "
		      "users from block %s", 
		      bg_record->bg_block_id);
	} 
	slurm_conf_unlock();
	
	if(bg_record->target_name 
	   && bg_record->user_name) {
		if(!strcmp(bg_record->target_name, user_name)) {
			if(strcmp(bg_record->target_name, 
				  bg_record->user_name)) {
				error("Block %s was in a ready state "
				      "for user %s but is being freed. "
				      "Job %d was lost.",
				      bg_record->bg_block_id,
				      bg_record->user_name,
				      jobid);
				slurm_mutex_unlock(&block_state_mutex);
				if(jobid > -1)
					slurm_fail_job(jobid);
				if(remove_from_bg_list(bg_job_block_list, 
						       bg_record) 
				   == SLURM_SUCCESS) {
					slurm_mutex_lock(&block_state_mutex);
					num_unused_cpus += bg_record->bp_count
						*bg_record->cpus_per_bp;
				} else {
					slurm_mutex_lock(&block_state_mutex);
				}
			} else {
				debug("Block %s was in a ready state "
				      "but is being freed. No job running.",
				      bg_record->bg_block_id);
			}
		} else {
			error("State went to free on a boot "
			      "for block %s.",
			      bg_record->bg_block_id);
		}
		slurm_mutex_unlock(&block_state_mutex);
		remove_from_bg_list(bg_booted_block_list, bg_record);
		slurm_mutex_lock(&block_state_mutex);
	} else if(bg_record->user_name) {
		error("Target Name was not set "
		      "not set for block %s.",
		      bg_record->bg_block_id);
		bg_record->target_name = 
			xstrdup(bg_record->user_name);
	} else {
		error("Target Name and User Name are "
		      "not set for block %s.",
		      bg_record->bg_block_id);
		bg_record->user_name = xstrdup(user_name);
		bg_record->target_name = 
			xstrdup(bg_record->user_name);
	}

	xfree(user_name);
			
	return SLURM_SUCCESS;
}
#endif


/*
 * check to see if block is ready to execute.  Meaning
 * User is added to the list of users able to run, and no one 
 * else is running on the block.
 *
 * NOTE: This happens in parallel with srun and slurmd spawning
 * the job. A prolog script is expected to defer initiation of
 * the job script until the BG block is available for use.
 */
extern int block_ready(struct job_record *job_ptr)
{
	int rc = 1;
	char *block_id = NULL;
	bg_record_t *bg_record = NULL;
	
	rc = select_g_get_jobinfo(job_ptr->select_jobinfo,
				  SELECT_DATA_BLOCK_ID, &block_id);
	if (rc == SLURM_SUCCESS) {
		bg_record = find_bg_record_in_list(bg_list, block_id);
		slurm_mutex_lock(&block_state_mutex);
		
		if(bg_record) {
			if(bg_record->job_running != job_ptr->job_id) {
				rc = 0;
			} else if ((bg_record->user_uid == job_ptr->user_id)
				   && (bg_record->state 
				       == RM_PARTITION_READY)) {
				rc = 1;
			} else if (bg_record->user_uid != job_ptr->user_id)
				rc = 0;
			else
				rc = READY_JOB_ERROR;	/* try again */
		} else {
			error("block_ready: block %s not in bg_list.",
			      block_id);
			rc = READY_JOB_FATAL;	/* fatal error */
		}
		slurm_mutex_unlock(&block_state_mutex);
		xfree(block_id);
	} else
		rc = READY_JOB_ERROR;
	return rc;
}				

/* Pack all relevent information about a block */
extern void pack_block(bg_record_t *bg_record, Buf buffer)
{
	packstr(bg_record->nodes, buffer);
	packstr(bg_record->user_name, buffer);
	packstr(bg_record->bg_block_id, buffer);
	pack16((uint16_t)bg_record->state, buffer);
	pack16((uint16_t)bg_record->conn_type, buffer);
	pack16((uint16_t)bg_record->node_use, buffer);	
	pack16((uint16_t)bg_record->quarter, buffer);	
	pack16((uint16_t)bg_record->nodecard, buffer);	
	pack32((uint32_t)bg_record->node_cnt, buffer);	
}

extern int update_block_list()
{
	int updated = 0;
#ifdef HAVE_BG_FILES
	int rc;
	rm_partition_t *block_ptr = NULL;
	rm_partition_mode_t node_use;
	rm_partition_state_t state;
	char *name = NULL;
	bg_record_t *bg_record = NULL;
	time_t now;
	int skipped_dealloc = 0;
	ListIterator itr = NULL;
	
	if(!bg_list) 
		return updated;
	
	slurm_mutex_lock(&block_state_mutex);
	itr = list_iterator_create(bg_list);
	while ((bg_record = (bg_record_t *) list_next(itr)) != NULL) {
		if(!bg_record->bg_block_id)
			continue;
		name = bg_record->bg_block_id;
		if ((rc = bridge_get_block_info(name, &block_ptr)) 
		    != STATUS_OK) {
			if(rc == INCONSISTENT_DATA
			   && bluegene_layout_mode == LAYOUT_DYNAMIC)
				continue;
			
			error("bridge_get_block_info(%s): %s", 
			      name, 
			      bg_err_str(rc));
			continue;
		}
				
		if ((rc = bridge_get_data(block_ptr, RM_PartitionMode,
					  &node_use))
		    != STATUS_OK) {
			error("bridge_get_data(RM_PartitionMode): %s",
			      bg_err_str(rc));
			updated = -1;
			goto next_block;
		} else if(bg_record->node_use != node_use) {
			debug("node_use of Block %s was %d "
			      "and now is %d",
			      bg_record->bg_block_id, 
			      bg_record->node_use, 
			      node_use);
			bg_record->node_use = node_use;
			updated = 1;
		}
		
		if ((rc = bridge_get_data(block_ptr, RM_PartitionState,
					  &state))
		    != STATUS_OK) {
			error("bridge_get_data(RM_PartitionState): %s",
			      bg_err_str(rc));
			updated = -1;
			goto next_block;
		} else if(bg_record->job_running != -3 
			  //plugin set error
			  && bg_record->state != state) {
			debug("state of Block %s was %d and now is %d",
			      bg_record->bg_block_id, 
			      bg_record->state, 
			      state);
			/* 
			   check to make sure block went 
			   through freeing correctly 
			*/
			if(bg_record->state 
			   != RM_PARTITION_DEALLOCATING
			   && state == RM_PARTITION_FREE)
				skipped_dealloc = 1;

			bg_record->state = state;

			if(bg_record->state 
			   == RM_PARTITION_DEALLOCATING) {
				_block_is_deallocating(bg_record);
			} else if(skipped_dealloc) {
				_block_is_deallocating(bg_record);
				skipped_dealloc = 0;
			} else if(bg_record->state 
				  == RM_PARTITION_CONFIGURING)
				bg_record->boot_state = 1;
			updated = 1;
		}

		/* check the boot state */
		debug3("boot state for block %s is %d",
		       bg_record->bg_block_id,
		       bg_record->boot_state);
		if(bg_record->boot_state == 1) {
			switch(bg_record->state) {
			case RM_PARTITION_CONFIGURING:
				debug3("checking to make sure user %s "
				       "is the user.",
				       bg_record->target_name);
				slurm_conf_lock();
				if(update_block_user(bg_record, 0) 
				   == 1)
					last_bg_update = time(NULL);
				slurm_conf_unlock();
				break;
			case RM_PARTITION_ERROR:
				error("block in an error state");
			case RM_PARTITION_FREE:
				if(bg_record->boot_count 
				   < RETRY_BOOT_COUNT) {
					slurm_mutex_unlock(
						&block_state_mutex);
					if((rc = boot_block(bg_record))
					   != SLURM_SUCCESS) {
						updated = -1;
					}
					slurm_mutex_lock(
						&block_state_mutex);
					debug("boot count for block "
					      "%s is %d",
					      bg_record->bg_block_id,
					      bg_record->boot_count);
					bg_record->boot_count++;
				} else {
					char reason[128], time_str[32];

					error("Couldn't boot Block %s "
					      "for user %s",
					      bg_record->bg_block_id, 
					      bg_record->target_name);
					slurm_mutex_unlock(
						&block_state_mutex);
					
					now = time(NULL);
					slurm_make_time_str(&now, time_str,
							    sizeof(time_str));
					snprintf(reason, 
						 sizeof(reason),
						 "update_block_list: "
						 "Boot fails "
						 "[SLURM@%s]", 
						 time_str);
					drain_as_needed(bg_record, reason);
					slurm_mutex_lock(&block_state_mutex);
					bg_record->boot_state = 0;
					bg_record->boot_count = 0;
				}
				break;
			case RM_PARTITION_READY:
				debug("block %s is ready.",
				      bg_record->bg_block_id);
				set_block_user(bg_record); 	
				break;
			default:
				debug("Hey the state of the "
				      "Block is %d doing nothing.",
				      bg_record->state);
				break;
			}
		}
	next_block:
		if ((rc = bridge_free_block(block_ptr)) 
		    != STATUS_OK) {
			error("bridge_free_block(): %s", 
			      bg_err_str(rc));
		}				
	}
	list_iterator_destroy(itr);
	slurm_mutex_unlock(&block_state_mutex);
		
#endif
	return updated;
}
