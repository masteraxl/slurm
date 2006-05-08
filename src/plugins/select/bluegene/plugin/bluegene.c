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
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#include "bluegene.h"
#include <stdio.h>

#define BUFSIZE 4096
#define BITSIZE 128
#define MMCS_POLL_TIME 120	/* poll MMCS for down switches and nodes 
				 * every 120 secs */
#define BG_POLL_TIME 0	        /* poll bg blocks every 3 secs */

#define _DEBUG 0

char* bg_conf = NULL;

/* Global variables */
rm_BGL_t *bg;

List bg_list = NULL;			/* total list of bg_record entries */
List bg_curr_block_list = NULL;  	/* current bg blocks in bluegene.conf*/
List bg_found_block_list = NULL;  	/* found bg blocks already on system */
List bg_job_block_list = NULL;  	/* jobs running in these blocks */
List bg_booted_block_list = NULL;  	/* blocks that are booted */
List bg_freeing_list = NULL;  	        /* blocks that being freed */

char *bluegene_blrts = NULL, *bluegene_linux = NULL, *bluegene_mloader = NULL;
char *bluegene_ramdisk = NULL, *bridge_api_file = NULL; 
bg_layout_t bluegene_layout_mode = NO_VAL;
uint16_t bluegene_numpsets = 0;
uint16_t bluegene_bp_node_cnt = 0;
uint16_t bluegene_quarter_node_cnt = 0;
uint16_t bluegene_nodecard_node_cnt = 0;
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

#ifdef HAVE_BG_FILES
  static int _update_bg_record_state(List bg_destroy_list);
#else
# if BA_SYSTEM_DIMENSIONS==3
    int max_dim[BA_SYSTEM_DIMENSIONS] = { 0, 0, 0 };
# else
    int max_dim[BA_SYSTEM_DIMENSIONS] = { 0 };
# endif
#endif

/* some local functions */
#ifdef HAVE_BG
static int  _addto_node_list(bg_record_t *bg_record, int *start, int *end);
#endif

static void _set_bg_lists();
static int  _validate_config_nodes(void);
static int  _bg_record_cmpf_inc(bg_record_t *rec_a, bg_record_t *rec_b);
static int _delete_old_blocks(void);
static char *_get_bg_conf(void);
static int _add_block_db(bg_record_t *bg_record, int *block_inx);
static int _split_block(bg_record_t *bg_record, int procs, int *block_inx);
static int _breakup_blocks(ba_request_t *request, List my_block_list, 
			   int *block_inx);
static bg_record_t *_create_small_record(bg_record_t *bg_record, 
					 uint16_t quarter, uint16_t nodecard);
static int _add_bg_record(List records, blockreq_t *blockreq);
static int  _reopen_bridge_log(void);

/* Initialize all plugin variables */
extern int init_bg(void)
{
#ifdef HAVE_BG_FILES
	int rc;
	rm_size3D_t bp_size;
	
	info("Attempting to contact MMCS");
	slurm_mutex_lock(&api_file_mutex);
	if ((rc = rm_set_serial(BG_SERIAL)) != STATUS_OK) {
		slurm_mutex_unlock(&api_file_mutex);
		fatal("init_bg: rm_set_serial(): %s", bg_err_str(rc));
		return SLURM_ERROR;
	}
	
	if ((rc = rm_get_BGL(&bg)) != STATUS_OK) {
		slurm_mutex_unlock(&api_file_mutex);
		fatal("init_bg: rm_get_BGL(): %s", bg_err_str(rc));
		return SLURM_ERROR;
	}
	slurm_mutex_unlock(&api_file_mutex);
	
	if ((rc = rm_get_data(bg, RM_Msize, &bp_size)) != STATUS_OK) {
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
	if (bg_found_block_list) {
		list_destroy(bg_found_block_list);
		bg_found_block_list = NULL;
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
	xfree(bluegene_blrts);
	xfree(bluegene_linux);
	xfree(bluegene_mloader);
	xfree(bluegene_ramdisk);
	xfree(bridge_api_file);
	xfree(bluegene_layout_mode);
	xfree(bg_conf);
	
#ifdef HAVE_BG_FILES
	if(bg)
		if ((rc = rm_free_BGL(bg)) != STATUS_OK)
			error("rm_free_BGL(): %s", bg_err_str(rc));
#endif	
	ba_fini();
}

extern void print_bg_record(bg_record_t* bg_record)
{
	char tmp_char[256];

	if (!bg_record) {
		error("print_bg_record, record given is null");
		return;
	}
#if _DEBUG
	info(" bg_record: ");
	if (bg_record->bg_block_id)
		info("\tbg_block_id: %s", bg_record->bg_block_id);
	info("\tnodes: %s", bg_record->nodes);
	info("\tsize: %d BPs %d Nodes %d cpus", 
	     bg_record->bp_count,
	     bg_record->node_cnt,
	     bg_record->cpus_per_bp * bg_record->bp_count);
	info("\tgeo: %dx%dx%d", bg_record->geo[X], bg_record->geo[Y], 
	     bg_record->geo[Z]);
	info("\tlifecycle: %s", convert_lifecycle(bg_record->block_lifecycle));
	info("\tconn_type: %s", convert_conn_type(bg_record->conn_type));
	info("\tnode_use: %s", convert_node_use(bg_record->node_use));
	if (bg_record->hostlist) {
		char buffer[BUFSIZE];
		hostlist_ranged_string(bg_record->hostlist, BUFSIZE, buffer);
		info("\thostlist %s", buffer);
	}
	if (bg_record->bitmap) {
		char bitstring[BITSIZE];
		bit_fmt(bitstring, BITSIZE, bg_record->bitmap);
		info("\tbitmap: %s", bitstring);
	}
#else
	format_node_name(bg_record, tmp_char);
	info("Record: BlockID:%s Nodes:%s Conn:%s",
	     bg_record->bg_block_id, tmp_char,
	     convert_conn_type(bg_record->conn_type));
#endif
}

extern void destroy_bg_record(void *object)
{
	bg_record_t* bg_record = (bg_record_t*) object;

	if (bg_record) {
		xfree(bg_record->bg_block_id);
		bg_record->bg_block_id = NULL;
		xfree(bg_record->nodes);
		xfree(bg_record->user_name);
		xfree(bg_record->target_name);
		if(bg_record->bg_block_list)
			list_destroy(bg_record->bg_block_list);
		if(bg_record->hostlist)
			hostlist_destroy(bg_record->hostlist);
		if(bg_record->bitmap)
			bit_free(bg_record->bitmap);
		
		xfree(bg_record);
		bg_record = NULL;
	}
}

extern int block_exist_in_list(List my_list, bg_record_t *bg_record)
{
	ListIterator itr = list_iterator_create(my_list);
	bg_record_t *found_record = NULL;
	int rc = 0;

	while ((found_record = (bg_record_t *) list_next(itr)) != NULL) {
		if(bit_equal(bg_record->bitmap, found_record->bitmap)
		   && (bg_record->quarter == found_record->quarter)
		   && (bg_record->nodecard == found_record->nodecard)){
			debug2("This partition %s %d %d "
			       "already exists here %s",
			       bg_record->nodes,
			       bg_record->quarter,
			       bg_record->nodecard,
			       found_record->bg_block_id);
			rc = 1;
			break;
		}
	}
	list_iterator_destroy(itr);
	return rc;
}

extern void process_nodes(bg_record_t *bg_record)
{
#ifdef HAVE_BG
	int j=0, number;
	int start[BA_SYSTEM_DIMENSIONS];
	int end[BA_SYSTEM_DIMENSIONS];
	ListIterator itr;
	ba_node_t* ba_node = NULL;
	
	bg_record->bp_count = 0;
				
	while (bg_record->nodes[j] != '\0') {
		if ((bg_record->nodes[j] == '['
		     || bg_record->nodes[j] == ',')
		    && (bg_record->nodes[j+8] == ']' 
			|| bg_record->nodes[j+8] == ',')
		    && (bg_record->nodes[j+4] == 'x'
			|| bg_record->nodes[j+4] == '-')) {
			j++;
			number = atoi(bg_record->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j += 4;
			number = atoi(bg_record->nodes + j);
			end[X] = number / 100;
			end[Y] = (number % 100) / 10;
			end[Z] = (number % 10);
			j += 3;
			if(!bg_record->bp_count) {
				bg_record->start[X] = start[X];
				bg_record->start[Y] = start[Y];
				bg_record->start[Z] = start[Z];
				debug2("start is %d%d%d",
				      bg_record->start[X],
				      bg_record->start[Y],
				      bg_record->start[Z]);
			}
			bg_record->bp_count += _addto_node_list(bg_record, 
								start, 
								end);
			if(bg_record->nodes[j] != ',')
				break;
			j--;
		} else if((bg_record->nodes[j] < 58 
			   && bg_record->nodes[j] > 47)) {
					
			number = atoi(bg_record->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j+=3;
			if(!bg_record->bp_count) {
				bg_record->start[X] = start[X];
				bg_record->start[Y] = start[Y];
				bg_record->start[Z] = start[Z];
				debug2("start is %d%d%d",
				      bg_record->start[X],
				      bg_record->start[Y],
				      bg_record->start[Z]);
			}
			bg_record->bp_count += _addto_node_list(bg_record, 
								 start, 
								 start);
			if(bg_record->nodes[j] != ',')
				break;
		}
		j++;
	}
	j=0;
	bg_record->geo[X] = 0;
	bg_record->geo[Y] = 0;
	bg_record->geo[Z] = 0;
	end[X] = -1;
	end[Y] = -1;
	end[Z] = -1;
	
	itr = list_iterator_create(bg_record->bg_block_list);
	while ((ba_node = list_next(itr)) != NULL) {
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
	debug3("geo = %d%d%d\n",
	       bg_record->geo[X],
	       bg_record->geo[Y],
	       bg_record->geo[Z]);
	
#ifndef HAVE_BG_FILES
	max_dim[X] = MAX(max_dim[X], end[X]);
	max_dim[Y] = MAX(max_dim[Y], end[Y]);
	max_dim[Z] = MAX(max_dim[Z], end[Z]);
#endif
   
	if (node_name2bitmap(bg_record->nodes, 
			     false, 
			     &bg_record->bitmap)) {
		fatal("Unable to convert nodes %s to bitmap", 
		      bg_record->nodes);
	}
#endif
	return;
}

extern void copy_bg_record(bg_record_t *fir_record, bg_record_t *sec_record)
{
	int i;

	xfree(sec_record->bg_block_id);
	sec_record->bg_block_id = xstrdup(fir_record->bg_block_id);
	xfree(sec_record->nodes);
	sec_record->nodes = xstrdup(fir_record->nodes);
	xfree(sec_record->user_name);
	sec_record->user_name = xstrdup(fir_record->user_name);
	xfree(sec_record->target_name);
	sec_record->target_name = xstrdup(fir_record->target_name);
	sec_record->user_uid = fir_record->user_uid;
	sec_record->block_lifecycle = fir_record->block_lifecycle;
	sec_record->state = fir_record->state;
	sec_record->conn_type = fir_record->conn_type;
	sec_record->node_use = fir_record->node_use;
	sec_record->bp_count = fir_record->bp_count;
	sec_record->switch_count = fir_record->switch_count;
	sec_record->boot_state = fir_record->boot_state;
	sec_record->boot_count = fir_record->boot_count;

	for(i=0;i<BA_SYSTEM_DIMENSIONS;i++) {
		sec_record->geo[i] = fir_record->geo[i];
		sec_record->start[i] = fir_record->start[i];
	}

	if(sec_record->bitmap)
		bit_free(sec_record->bitmap);
	if(fir_record->bitmap 
	   && (sec_record->bitmap = bit_copy(fir_record->bitmap)) == NULL) {
		error("Unable to copy bitmap for %s", fir_record->nodes);
		sec_record->bitmap = NULL;
	}
	sec_record->job_running = fir_record->job_running;
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
*/
extern int update_block_user(bg_record_t *bg_record, int set) 
{
	struct passwd *pw_ent = NULL;

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
				
				slurm_mutex_lock(&api_file_mutex);
				if ((rc = rm_add_part_user(
					     bg_record->bg_block_id, 
					     bg_record->target_name)) 
				    != STATUS_OK) {
					slurm_mutex_unlock(&api_file_mutex);
					error("rm_add_part_user(%s,%s): %s", 
					      bg_record->bg_block_id, 
					      bg_record->target_name,
					      bg_err_str(rc));
					return -1;
				} 
				slurm_mutex_unlock(&api_file_mutex);
			}
		}
	}
#endif
	
	if(strcmp(bg_record->target_name, bg_record->user_name)) {
		xfree(bg_record->user_name);
		bg_record->user_name = xstrdup(bg_record->target_name);
		if((pw_ent = getpwnam(bg_record->user_name)) == NULL) {
			error("getpwnam(%s): %m", bg_record->user_name);
			return -1;
		} else {
			bg_record->user_uid = pw_ent->pw_uid; 
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

	if(bg_record->job_running > -1)
		slurm_fail_job(bg_record->job_running);			

	/* small blocks */
	if(bg_record->cpus_per_bp != procs_per_node) {
		debug2("small block");
		while(bg_record->job_running > -1) 
			sleep(1);
		slurm_mutex_lock(&block_state_mutex);
		bg_record->job_running = -3;
		bg_record->state = RM_PARTITION_ERROR;
		slurm_mutex_unlock(&block_state_mutex);
		return;
	}
	
	/* at least one base partition */
	hl = hostlist_create(bg_record->nodes);
	if (!hl) {
		slurm_drain_nodes(bg_record->nodes, reason);
		return;
	}
	while ((host = hostlist_shift(hl))) {
		if (node_already_down(host)) {
			needed = false;
			free(host);
			break;
		}
		free(host);
	}
	hostlist_destroy(hl);
	
	if (needed)
		slurm_drain_nodes(bg_record->nodes, reason);
}

extern int format_node_name(bg_record_t *bg_record, char tmp_char[])
{
	if(bg_record->quarter != (uint16_t)NO_VAL) {
		if(bg_record->nodecard != (uint16_t)NO_VAL) {
			sprintf(tmp_char,"%s.%d.%d",
				bg_record->nodes,
				bg_record->quarter,
				bg_record->nodecard);
		} else {
			sprintf(tmp_char,"%s.%d",
				bg_record->nodes,
				bg_record->quarter);
		}
	} else {
		sprintf(tmp_char,"%s",bg_record->nodes);
	}
	return SLURM_SUCCESS;
}

extern bool blocks_overlap(bg_record_t *rec_a, bg_record_t *rec_b)
{
	bitstr_t *my_bitmap = NULL;
	int rc;
	if(rec_a->bp_count > 1 && rec_a->bp_count > 1) {
		reset_ba_system();
		load_block_wiring(rec_a->bg_block_id);
		rc = load_block_wiring(rec_b->bg_block_id);
		if(rc == SLURM_ERROR)
			return true;
	}
	my_bitmap = bit_copy(rec_a->bitmap);
	bit_and(my_bitmap, rec_b->bitmap);
	if (bit_ffs(my_bitmap) == -1) {
		bit_free(my_bitmap);
		return false;
	}
	bit_free(my_bitmap);
		
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

	slurm_mutex_lock(&api_file_mutex);
	if ((rc = rm_get_partition(bg_block_id,  &block_ptr)) != STATUS_OK) {
		slurm_mutex_unlock(&api_file_mutex);
		if(rc == INCONSISTENT_DATA
		   && bluegene_layout_mode == LAYOUT_DYNAMIC)
			return REMOVE_USER_FOUND;
			
		error("rm_get_partition(%s): %s", 
		      bg_block_id, 
		      bg_err_str(rc));
		return REMOVE_USER_ERR;
	}	
	slurm_mutex_unlock(&api_file_mutex);
	
	if((rc = rm_get_data(block_ptr, RM_PartitionUsersNum, &user_count)) 
	   != STATUS_OK) {
		error("rm_get_data(RM_PartitionUsersNum): %s", 
		      bg_err_str(rc));
		returnc = REMOVE_USER_ERR;
		user_count = 0;
	} else
		debug2("got %d users for %s",user_count, bg_block_id);
	for(i=0; i<user_count; i++) {
		if(i) {
			if ((rc = rm_get_data(block_ptr, 
					      RM_PartitionNextUser, 
					      &user)) 
			    != STATUS_OK) {
				error("rm_get_data(RM_PartitionNextUser): %s", 
				      bg_err_str(rc));
				returnc = REMOVE_USER_ERR;
				break;
			}
		} else {
			if ((rc = rm_get_data(block_ptr, 
					      RM_PartitionFirstUser, 
					      &user)) 
			    != STATUS_OK) {
				error("rm_get_data(RM_PartitionFirstUser): %s",
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
		slurm_mutex_lock(&api_file_mutex);
		if ((rc = rm_remove_part_user(bg_block_id, user)) 
		    != STATUS_OK) {
			debug("user %s isn't on block %s",
			      user, 
			      bg_block_id);
		}
		slurm_mutex_unlock(&api_file_mutex);
		free(user);
	}
	if ((rc = rm_free_partition(block_ptr)) != STATUS_OK) {
		error("rm_free_partition(): %s", bg_err_str(rc));
	}
#endif
	return returnc;
}

extern void set_block_user(bg_record_t *bg_record) 
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
	} else if (rc == -1) {
		error("Unable to add user name to block %s. "
		      "Cancelling job.",
		      bg_record->bg_block_id);
		(void) slurm_fail_job(bg_record->job_running);
	}	
	xfree(bg_record->target_name);
	bg_record->target_name = 
		xstrdup(slurmctld_conf.slurm_user_name);
	slurm_conf_unlock();			
}

extern char* convert_lifecycle(lifecycle_type_t lifecycle)
{
	if (lifecycle == DYNAMIC)
		return "DYNAMIC";
	else 
		return "STATIC";
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
	slurm_mutex_lock(&block_state_mutex);
	list_sort(records, (ListCmpF) _bg_record_cmpf_inc);
	last_bg_update = time(NULL);
	slurm_mutex_unlock(&block_state_mutex);
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

/*
 * create_defined_blocks - create the static blocks that will be used
 * for scheduling, all partitions must be able to be created and booted
 * at once.  
 * IN - int overlapped, 1 if partitions are to be overlapped, 0 if they are
 * static.
 * RET - success of fitting all configurations
 */
extern int create_defined_blocks(bg_layout_t overlapped)
{
	int rc = SLURM_SUCCESS;

	ListIterator itr;
	bg_record_t *bg_record = NULL;

#ifdef HAVE_BG_FILES
	bg_record_t *found_record = NULL;
	char *name = NULL;
	int geo[BA_SYSTEM_DIMENSIONS];
	int i;
	ListIterator itr_found;
	init_wires();
#endif
	slurm_mutex_lock(&block_state_mutex);
	reset_ba_system();
	
#ifdef HAVE_BG_FILES
	if(bg_list) {
		itr = list_iterator_create(bg_list);
		while ((bg_record = (bg_record_t *) list_next(itr)) 
		       != NULL) {
			if(bg_found_block_list) {
				itr_found = list_iterator_create(
					bg_found_block_list);
				while ((found_record = (bg_record_t*) 
					list_next(itr_found)) != NULL) {
/* 					info("%s.%d.%d ?= %s.%d.%d\n", */
/* 					     bg_record->nodes, */
/* 					     bg_record->quarter, */
/* 					     bg_record->nodecard, */
/* 					     found_record->nodes, */
/* 					     found_record->quarter, */
/* 					     found_record->nodecard); */
					
					if ((bit_equal(bg_record->bitmap, 
						       found_record->bitmap))
					    && (bg_record->quarter ==
						found_record->quarter)
					    && (bg_record->nodecard ==
						found_record->nodecard)) {
						/* don't reboot this one */
						break;	
					}
				}
				list_iterator_destroy(itr_found);
			} else {
				error("create_defined_blocks: "
				      "no bg_found_block_list 1");
			}
			if(bg_record->bp_count>0 
			   && !bg_record->full_block
			   && bg_record->cpus_per_bp == procs_per_node) {
				if(overlapped == LAYOUT_OVERLAP)
					reset_ba_system();
				for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
					geo[i] = bg_record->geo[i];
				debug2("adding %s %d%d%d %d%d%d",
				       bg_record->nodes,
				       bg_record->start[X],
				       bg_record->start[Y],
				       bg_record->start[Z],
				       geo[X],
				       geo[Y],
				       geo[Z]);	
				if(bg_record->bg_block_id) {
					rc = SLURM_ERROR;
					rc = load_block_wiring(
						bg_record->bg_block_id);
				}
				if(rc == -1) {
					name = set_bg_block(NULL,
							    bg_record->start, 
							    geo, 
							    bg_record->
							    conn_type);
					if(!name) {			
						debug("I was unable to make "
						      "the requested block.");
						slurm_mutex_unlock(
							&block_state_mutex);
						return SLURM_ERROR;
					}
					xfree(name);
				} else if (rc == SLURM_ERROR) {
					debug("something happened in the "
					      "load of %s", 
					      bg_record->bg_block_id);
					slurm_mutex_unlock(
						&block_state_mutex);
					return SLURM_ERROR;
				}
			}
			if(found_record == NULL) {
				if((rc = configure_block(bg_record)) 
				   == SLURM_ERROR) {
					list_iterator_destroy(itr);
					slurm_mutex_unlock(&block_state_mutex);
					return rc;
				}
				print_bg_record(bg_record);
			}
		}
		list_iterator_destroy(itr);
	} else {
		error("create_defined_blocks: no bg_list 2");
		slurm_mutex_unlock(&block_state_mutex);
		return SLURM_ERROR;
	}
#endif
	slurm_mutex_unlock(&block_state_mutex);
	create_full_system_block();

	sort_bg_record_inc_size(bg_list);

	
#ifndef HAVE_BG_FILES
	char tmp_char[256];
	static int block_inx = 0;
	if(bg_list) {
		slurm_mutex_lock(&block_state_mutex);
		itr = list_iterator_create(bg_list);
		while ((bg_record = (bg_record_t*) list_next(itr))) {
			if (bg_record->bg_block_id)
				continue;
			bg_record->bg_block_id = xmalloc(sizeof(char)*8);
			snprintf(bg_record->bg_block_id, 8, "RMP%d", 
				 block_inx++);
			format_node_name(bg_record, tmp_char);
			info("Record: BlockID:%s Nodes:%s Conn:%s",
			     bg_record->bg_block_id, tmp_char,
			     convert_conn_type(bg_record->conn_type));
		}
		list_iterator_destroy(itr);
		slurm_mutex_unlock(&block_state_mutex);
	} else {
		error("create_defined_blocks: no bg_list 4");
		return SLURM_ERROR;
	}
	
	
#endif	/* not have HAVE_BG_FILES */
	

#ifdef _PRINT_BLOCKS_AND_EXIT
	if(bg_list) {
		itr = list_iterator_create(bg_list);
		debug("\n\n");
		while ((found_record = (bg_record_t *) list_next(itr)) 
		       != NULL) {
			print_bg_record(found_record);
		}
		list_iterator_destroy(itr);
	} else {
		error("create_defined_blocks: no bg_list 5");
	}
 	exit(0);
#endif	/* _PRINT_BLOCKS_AND_EXIT */
	rc = SLURM_SUCCESS;
	//exit(0);
	return rc;
}



/*
 * create_dynamic_block - create a new block to be used for a new
 * job allocation.  This will be added to the booted and job bg_lists.
 * RET - success of fitting configuration in the running system.
 */
extern int create_dynamic_block(ba_request_t *request, List my_block_list)
{
	int rc = SLURM_SUCCESS;
	static int block_inx = 0;
	
	ListIterator itr;
	bg_record_t *bg_record = NULL;
	List results = list_create(NULL);
	uint16_t num_quarter=0, num_nodecard=0;
	char *name = NULL;
	bitstr_t *my_bitmap = NULL;
	int geo[BA_SYSTEM_DIMENSIONS];
	int i;
	blockreq_t blockreq;
			
	slurm_mutex_lock(&block_state_mutex);
	reset_ba_system();
		
	if(my_block_list) {
		itr = list_iterator_create(my_block_list);
		while ((bg_record = (bg_record_t *) list_next(itr)) != NULL) {
			if(!my_bitmap) {
				my_bitmap = 
					bit_alloc(bit_size(bg_record->bitmap));
			}
				
			if(bg_record->job_running != -2 
			   && !bit_super_set(bg_record->bitmap, my_bitmap)) {
				bit_or(my_bitmap, bg_record->bitmap);
	
				for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
					geo[i] = bg_record->geo[i];
				debug2("adding %s %d%d%d %d%d%d",
				       bg_record->nodes,
				       bg_record->start[X],
				       bg_record->start[Y],
				       bg_record->start[Z],
				       geo[X],
				       geo[Y],
				       geo[Z]);

				if(bg_record->bg_block_id) {
					rc = SLURM_ERROR;
					rc = load_block_wiring(
						bg_record->bg_block_id);
				}
				if(rc == -1) {
					name = set_bg_block(NULL,
							    bg_record->start, 
							    geo, 
							    bg_record->
							    conn_type);
					if(!name) {
						debug("I was unable to make "
						      "the requested block.");
						bit_free(my_bitmap);
						slurm_mutex_unlock(
							&block_state_mutex);
						return SLURM_ERROR;
					}
					xfree(name);
				} else if (rc == SLURM_ERROR) {
					debug("something happened in the "
					      "load of %s", 
					      bg_record->bg_block_id);
					slurm_mutex_unlock(
						&block_state_mutex);
					return SLURM_ERROR;
				}
			}
		}
		list_iterator_destroy(itr);
		if(my_bitmap)
			bit_free(my_bitmap);
	} else {
		debug("No list was given");
	}
	
	if(request->size==1 && request->procs < bluegene_bp_node_cnt) {
		request->conn_type = SELECT_SMALL;
		if(request->procs == (procs_per_node/16)) {
			num_nodecard=4;
			num_quarter=3;
		} else {
			num_quarter=4;
		}
		if(_breakup_blocks(request, my_block_list, &block_inx) 
		   == SLURM_SUCCESS)
			goto finished;
	}
	
	if(request->conn_type == SELECT_NAV)
		request->conn_type = SELECT_TORUS;
	
	if(!new_ba_request(request)) {
		error("Problems with request for size %d geo %dx%dx%d", 
		      request->size,
		      request->geometry[X], 
		      request->geometry[Y], 
		      request->geometry[Z]);
		rc = SLURM_ERROR;
		goto finished;
	} 
	
	if(!list_count(bg_list) || !my_block_list) {
		bg_record = NULL;
		goto no_list;
	}

	/*Try to put block starting in the smallest of the exisiting blocks*/
	itr = list_iterator_create(bg_list);
	while ((bg_record = (bg_record_t *) list_next(itr)) != NULL) {
		request->rotate_count = 0;
		request->elongate_count = 1;
		
		if(bg_record->job_running == -1 
		   && (bg_record->quarter == (uint16_t) NO_VAL
		       || (bg_record->quarter == 0 
			   && (bg_record->nodecard == (uint16_t) NO_VAL
			       || bg_record->nodecard == 0)))) {
			
			for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
				request->start[i] = bg_record->start[i];
			debug2("allocating %d%d%d %d",
			       bg_record->nodes,
			       request->start[X],
			       request->start[Y],
			       request->start[Z],
			       request->size);
			request->start_req = 1;
			rc = SLURM_SUCCESS;
			if (!allocate_block(request, results)){
				debug2("allocate failure for size %d "
				       "base partitions", 
				       request->size);
				rc = SLURM_ERROR;
			} else 
				break;
		}
	}
	list_iterator_destroy(itr);

no_list:
	if(!bg_record) {
		request->start_req = 0;
		rc = SLURM_SUCCESS;
		if (!allocate_block(request, results)){
			debug("allocate failure for size %d base partitions", 
			      request->size);
			rc = SLURM_ERROR;
		}
	}

	if(rc == SLURM_ERROR || !my_block_list) {
		goto finished;
	}
	/*set up bg_record(s) here */
	list_destroy(results);
	results = list_create(destroy_bg_record);
	
	blockreq.block = request->save_name;
	blockreq.conn_type = request->conn_type;
	blockreq.nodecards = num_nodecard;
	blockreq.quarters = num_quarter;

	_add_bg_record(results, &blockreq);

	while((bg_record = (bg_record_t *) list_pop(results)) != NULL) {
		if(block_exist_in_list(bg_list, bg_record))
			destroy_bg_record(bg_record);
		else {
			if(_add_block_db(bg_record, &block_inx) == SLURM_ERROR)
				goto finished;
			list_push(bg_list, bg_record);
			print_bg_record(bg_record);
		}
	}

finished:
	if(my_block_list)
		xfree(request->save_name);
	if(request->elongate_geos)
		list_destroy(request->elongate_geos);
	if(results)
		list_destroy(results);
	
	slurm_mutex_unlock(&block_state_mutex);
	sort_bg_record_inc_size(bg_list);
	
	return rc;
}

extern int create_full_system_block()
{
	int rc = SLURM_SUCCESS;
	ListIterator itr;
	bg_record_t *bg_record = NULL;
	char *name = NULL;
	List records = NULL;
	int geo[BA_SYSTEM_DIMENSIONS];
	int i;
	blockreq_t blockreq;

	/* Here we are adding a block that in for the entire machine 
	   just in case it isn't in the bluegene.conf file.
	*/
	slurm_mutex_lock(&block_state_mutex);
	
#ifdef HAVE_BG_FILES
	geo[X] = DIM_SIZE[X] - 1;
	geo[Y] = DIM_SIZE[Y] - 1;
	geo[Z] = DIM_SIZE[Z] - 1;
#else
	geo[X] = max_dim[X];
	geo[Y] = max_dim[Y];
	geo[Z] = max_dim[Z];
#endif
	slurm_conf_lock();
	name = xmalloc(sizeof(char)*(10+strlen(slurmctld_conf.node_prefix)));
	if((geo[X] == 0) && (geo[Y] == 0) && (geo[Z] == 0))
		sprintf(name, "%s000", slurmctld_conf.node_prefix);
	else
		sprintf(name, "%s[000x%d%d%d]",
			slurmctld_conf.node_prefix,
			geo[X], geo[Y], geo[Z]);
	slurm_conf_unlock();
			
	if(bg_found_block_list) {
		itr = list_iterator_create(bg_found_block_list);
		while ((bg_record = (bg_record_t *) list_next(itr)) != NULL) {
			if (!strcmp(name, bg_record->nodes)) {
				xfree(name);
				list_iterator_destroy(itr);
				/* don't create total already there */
				goto no_total;	
			}
		}
		list_iterator_destroy(itr);
	} else {
		error("create_full_system_block: no bg_found_block_list 2");
	}
	
	if(bg_list) {
		itr = list_iterator_create(bg_list);
		while ((bg_record = (bg_record_t *) list_next(itr)) 
		       != NULL) {
			if (!strcmp(name, bg_record->nodes)) {
				xfree(name);
				list_iterator_destroy(itr);
				/* don't create total already there */
				goto no_total;	
			}
		}
		list_iterator_destroy(itr);
	} else {
		xfree(name);
		error("create_overlapped_blocks: no bg_list 3");
		rc = SLURM_ERROR;
		goto no_total;
	}

	records = list_create(destroy_bg_record);
	blockreq.block = name;
	blockreq.conn_type = SELECT_TORUS;
	blockreq.nodecards = 0;
	blockreq.quarters = 0;
	_add_bg_record(records, &blockreq);
	xfree(name);
	
	bg_record = (bg_record_t *) list_pop(records);
	if(!bg_record) {
		error("Nothing was returned from full system create");
			rc = SLURM_ERROR;
			goto no_total;
	}
	reset_ba_system();
	for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) 
		geo[i] = bg_record->geo[i];
	debug2("adding %s %d%d%d %d%d%d",
	       bg_record->nodes,
	       bg_record->start[X],
	       bg_record->start[Y],
	       bg_record->start[Z],
	       geo[X],
	       geo[Y],
	       geo[Z]);
	
	name = set_bg_block(NULL,
			    bg_record->start, 
			    geo, 
			    bg_record->conn_type);		
	if(!name) {
		error("I was unable to make the "
		      "requested block.");
		rc = SLURM_ERROR;
		destroy_bg_record(bg_record);
		goto no_total;
	}
	xfree(name);
	
#ifdef HAVE_BG_FILES
	if((rc = configure_block(bg_record)) == SLURM_ERROR) {
		error("unable to configure block in api");
		destroy_bg_record(bg_record);
		goto no_total;
	}
#endif	/* HAVE_BG_FILES */
	list_push(bg_list, bg_record);

no_total:
	if(records)
		list_destroy(records);
	slurm_mutex_unlock(&block_state_mutex);
	return rc;
}

extern int remove_from_bg_list(List my_bg_list, bg_record_t *bg_record)
{
	bg_record_t *found_record = NULL;
	ListIterator itr;
	int rc = SLURM_ERROR;

	slurm_mutex_lock(&block_state_mutex);	
	itr = list_iterator_create(my_bg_list);
	while ((found_record = (bg_record_t *) list_next(itr)) != NULL) {
		if(found_record && bg_record)
			if(bg_record == found_record) {
				list_remove(itr);
				rc = SLURM_SUCCESS;
				break;
			}
	}
	list_iterator_destroy(itr);
	slurm_mutex_unlock(&block_state_mutex);

	return rc;
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
		if (bg_record->state != NO_VAL
		    && bg_record->state != RM_PARTITION_FREE 
		    && bg_record->state != RM_PARTITION_DEALLOCATING) {
#ifdef HAVE_BG_FILES
			debug2("pm_destroy %s",bg_record->bg_block_id);
			
			slurm_mutex_lock(&api_file_mutex);
			rc = pm_destroy_partition(bg_record->bg_block_id);
			slurm_mutex_unlock(&api_file_mutex);
			if (rc != STATUS_OK) {
				if(rc == PARTITION_NOT_FOUND) {
					debug("block %s is not found",
					      bg_record->bg_block_id);
					break;
				} else if(rc == INCOMPATIBLE_STATE) {
					debug2("pm_destroy_partition(%s): %s "
					      "State = %d",
					      bg_record->bg_block_id, 
					      bg_err_str(rc), 
					      bg_record->state);
				} else {
					error("pm_destroy_partition(%s): %s "
					      "State = %d",
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
	slurm_mutex_unlock(&block_state_mutex);			
	remove_from_bg_list(bg_booted_block_list, bg_record);
		
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
	slurm_mutex_unlock(&freed_cnt_mutex);	
	return NULL;
}

/* destroy multiple blocks in parallel */
extern void *mult_destroy_block(void *args)
{
	bg_record_t *bg_record = NULL;
	bg_record_t *found_record = NULL;
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
		if(bg_record->job_running == -2) {
			slurm_mutex_unlock(&block_state_mutex);
			goto already_here;
		}
		bg_record->job_running = -2;
		slurm_mutex_unlock(&block_state_mutex);
			
		slurm_mutex_lock(&freed_cnt_mutex);
		if(find_bg_record_in_list(bg_freeing_list, 
					  bg_record->bg_block_id)) {
			slurm_mutex_unlock(&freed_cnt_mutex);
			goto already_here;	
		}

		found_record = xmalloc(sizeof(bg_record_t));
		found_record->bg_block_id = xstrdup(bg_record->bg_block_id);
		list_push(bg_freeing_list, found_record);
		slurm_mutex_unlock(&freed_cnt_mutex);

		debug2("removing the jobs on block %s\n",
		      bg_record->bg_block_id);
		term_jobs_on_block(bg_record->bg_block_id);
		
		debug2("destroying %s", (char *)bg_record->bg_block_id);
		if(bg_free_block(bg_record) == SLURM_ERROR) {
			debug("there was an error");
			goto already_here;
		}
		debug2("done destroying");
		remove_from_bg_list(bg_list, bg_record);
		
#ifdef HAVE_BG_FILES
		debug2("removing from database %s", 
		       (char *)found_record->bg_block_id);

		slurm_mutex_lock(&api_file_mutex);
		rc = rm_remove_partition(found_record->bg_block_id);
		if (rc != STATUS_OK) {
			if(rc == PARTITION_NOT_FOUND) {
				debug("block %s is not found",
				      found_record->bg_block_id);
			} else {
				error("1 rm_remove_partition(%s): %s",
				      found_record->bg_block_id,
				      bg_err_str(rc));
			}
		} else
			debug2("done");
		slurm_mutex_unlock(&api_file_mutex);	
#endif
		slurm_mutex_lock(&block_state_mutex);
		if(blocks_are_created)
			destroy_bg_record(bg_record);
		destroy_bg_record(found_record);
		slurm_mutex_unlock(&block_state_mutex);
		slurm_mutex_lock(&freed_cnt_mutex);
		remove_from_bg_list(bg_freeing_list, found_record);
		slurm_mutex_unlock(&freed_cnt_mutex);
		
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
	/* already running MAX_AGENTS we don't really need more 
	   since they never end */
	
	while ((found_record = (bg_record_t*)list_pop(delete_list)) != NULL) {
		/* push job onto queue in a FIFO */
		if (list_push(*block_list, found_record) == NULL)
			fatal("malloc failure in _block_op/list_push");
		
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
	static time_t last_config_update = (time_t) 0;
	struct stat config_stat;
		
	debug("Reading the bluegene.conf file");

	/* check if config file has changed */
	if (!bg_conf)
		bg_conf = _get_bg_conf();
	if (stat(bg_conf, &config_stat) < 0)
		fatal("can't stat bluegene.conf file %s: %m", bg_conf);
	if (last_config_update) {
		_reopen_bridge_log();
		last_config_update = config_stat.st_mtime; 
		if(last_config_update == config_stat.st_mtime)
			debug("bluegene.conf unchanged");
		else
			debug("bluegene.conf changed, doing nothing");
		return SLURM_SUCCESS;
	}
	last_config_update = config_stat.st_mtime; 

	/* initialization */
	/* bg_conf defined in bg_node_alloc.h */
	tbl = s_p_hashtbl_create(bg_conf_file_options);
	
	if(s_p_parse_file(tbl, bg_conf) == SLURM_ERROR)
		fatal("something wrong with opening/reading bluegene "
		      "conf file");
	
	if (!s_p_get_string(&bluegene_blrts, "BlrtsImage", tbl)) 
		fatal("BlrtsImage not configured in bluegene.conf");
	if (!s_p_get_string(&bluegene_linux, "LinuxImage", tbl)) 
		fatal("LinuxImage not configured in bluegene.conf");
	if (!s_p_get_string(&bluegene_mloader, "MloaderImage", tbl)) 
		fatal("MloaderImage not configured in bluegene.conf");
	if (!s_p_get_string(&bluegene_ramdisk, "RamDiskImage", tbl)) 
		fatal("RamDiskImage not configured in bluegene.conf");
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
		bluegene_quarter_node_cnt = bluegene_bp_node_cnt/4;
	}
	if (!s_p_get_uint16(
		    &bluegene_nodecard_node_cnt, "NodeCardNodeCnt", tbl)) {
		error("NodeCardNodeCnt not configured in bluegene.conf "
		      "defaulting to 32 as NodeCardNodeCnt");
		bluegene_nodecard_node_cnt = 32;
	}
	    
	_set_bg_lists();	
	/* add blocks defined in file */
	if(bluegene_layout_mode != LAYOUT_DYNAMIC) {
		if (!s_p_get_array((void ***)&blockreq_array, 
				   &count, "BPs", tbl)) {
			info("WARNING: no blocks defined in bluegene.conf, "
			     "only making full system block");
			create_full_system_block();
		}
		
		for (i = 0; i < count; i++) {
			_add_bg_record(bg_list, blockreq_array[i]);
		}
	}
//#if 0	
	/* Check to see if the configs we have are correct */
	if (_validate_config_nodes() == SLURM_ERROR) { 
		_delete_old_blocks();
	}
//#endif
	/* looking for blocks only I created */
	if(bluegene_layout_mode == LAYOUT_DYNAMIC) {
		init_wires();
		info("No blocks created until jobs are submitted");
	} else {
		if (create_defined_blocks(bluegene_layout_mode) 
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
	list_destroy(bg_found_block_list);
	bg_found_block_list = NULL;
	last_bg_update = time(NULL);
	blocks_are_created = 1;
	slurm_mutex_unlock(&block_state_mutex);
	sort_bg_record_inc_size(bg_list);
	debug("Blocks have finished being created.");
	s_p_hashtbl_destroy(tbl);

	return SLURM_SUCCESS;
}

#ifdef HAVE_BG_FILES
static int _update_bg_record_state(List bg_destroy_list)
{
	rm_partition_state_flag_t block_state = PARTITION_ALL_FLAG;
	char *name = NULL;
	rm_partition_list_t *block_list = NULL;
	int j, rc, func_rc = SLURM_SUCCESS, num_blocks = 0;
	rm_partition_state_t state;
	rm_partition_t *block_ptr = NULL;
	ListIterator itr;
	bg_record_t* bg_record = NULL;	

	if(!bg_destroy_list) {
		return SLURM_SUCCESS;
	}
	
	slurm_mutex_lock(&api_file_mutex);
	if ((rc = rm_get_partitions_info(block_state, &block_list))
	    != STATUS_OK) {
		slurm_mutex_unlock(&api_file_mutex);
		error("1 rm_get_partitions_info(): %s", bg_err_str(rc));
		return SLURM_ERROR; 
	}
	
	if ((rc = rm_get_data(block_list, RM_PartListSize, &num_blocks))
	    != STATUS_OK) {
		error("rm_get_data(RM_PartListSize): %s", bg_err_str(rc));
		func_rc = SLURM_ERROR;
		num_blocks = 0;
	}
		
	for (j=0; j<num_blocks; j++) {
		if (j) {
			if ((rc = rm_get_data(block_list, 
					      RM_PartListNextPart, 
					      &block_ptr)) 
			    != STATUS_OK) {
				error("rm_get_data(RM_PartListNextPart): %s",
				      bg_err_str(rc));
				func_rc = SLURM_ERROR;
				break;
			}
		} else {
			if ((rc = rm_get_data(block_list, 
					      RM_PartListFirstPart, 
					      &block_ptr)) 
			    != STATUS_OK) {
				error("rm_get_data(RM_PartListFirstPart: %s",
				      bg_err_str(rc));
				func_rc = SLURM_ERROR;
				break;
			}
		}
		if ((rc = rm_get_data(block_ptr, 
				      RM_PartitionID, 
				      &name))
		    != STATUS_OK) {
			error("rm_get_data(RM_PartitionID): %s", 
			      bg_err_str(rc));
			func_rc = SLURM_ERROR;
			break;
		}
		if (!name) {
			error("RM_Partition is NULL");
			continue;
		}
		
		slurm_mutex_lock(&block_state_mutex);
		itr = list_iterator_create(bg_destroy_list);
		while ((bg_record = (bg_record_t*) list_next(itr))) {	
			if(!bg_record->bg_block_id) 
				continue;
			if(strcmp(bg_record->bg_block_id, name)) {
				continue;		
			}
		       
			if ((rc = rm_get_data(block_ptr, 
					      RM_PartitionState, 
					      &state))
			    != STATUS_OK) {
				error("rm_get_data(RM_PartitionState): %s",
				      bg_err_str(rc));
			} else if(bg_record->state != state) {
				debug("state of Block %s was %d "
				      "and now is %d",
				      name, bg_record->state, state);
				bg_record->state = state;
			}
			break;
		}
		list_iterator_destroy(itr);
		slurm_mutex_unlock(&block_state_mutex);
			
		free(name);
	}
	
	if ((rc = rm_free_partition_list(block_list)) != STATUS_OK) {
		error("rm_free_partition_list(): %s", bg_err_str(rc));
	}
	slurm_mutex_unlock(&api_file_mutex);
	
	return func_rc;
}
#endif /* HAVE_BG_FILES */

#ifdef HAVE_BG
static int _addto_node_list(bg_record_t *bg_record, int *start, int *end)
{
	int node_count=0;
	int x,y,z;
	char node_name_tmp[255];
	debug3("%d%d%dx%d%d%d",
	     start[X],
	     start[Y],
	     start[Z],
	     end[X],
	     end[Y],
	     end[Z]);
	debug3("%d%d%d",
	     DIM_SIZE[X],
	     DIM_SIZE[Y],
	     DIM_SIZE[Z]);
	     
	assert(end[X] < DIM_SIZE[X]);
	assert(start[X] >= 0);
	assert(end[Y] < DIM_SIZE[Y]);
	assert(start[Y] >= 0);
	assert(end[Z] < DIM_SIZE[Z]);
	assert(start[Z] >= 0);
	
	for (x = start[X]; x <= end[X]; x++) {
		for (y = start[Y]; y <= end[Y]; y++) {
			for (z = start[Z]; z <= end[Z]; z++) {
				slurm_conf_lock();
				sprintf(node_name_tmp, "%s%d%d%d", 
					slurmctld_conf.node_prefix,
					x, y, z);		
				slurm_conf_unlock();
				list_append(bg_record->bg_block_list, 
					    &ba_system_ptr->grid[x][y][z]);
				node_count++;
			}
		}
	}
	return node_count;
}
#endif //HAVE_BG

static void _set_bg_lists()
{
	slurm_mutex_lock(&block_state_mutex);
	if (bg_found_block_list)
		list_destroy(bg_found_block_list);
	bg_found_block_list = list_create(NULL);
	if (bg_booted_block_list) 
		list_destroy(bg_booted_block_list);
	bg_booted_block_list = list_create(NULL);
	if (bg_job_block_list) 
		list_destroy(bg_job_block_list);
	bg_job_block_list = list_create(NULL);	
	num_unused_cpus = 
		DIM_SIZE[X] * DIM_SIZE[Y] * DIM_SIZE[Z] * procs_per_node;
	if (bg_curr_block_list)
		list_destroy(bg_curr_block_list);	
	bg_curr_block_list = list_create(destroy_bg_record);
	
	
	if (bg_list) 
		list_destroy(bg_list);
	bg_list = list_create(destroy_bg_record);
	slurm_mutex_unlock(&block_state_mutex);		
}

/*
 * Match slurm configuration information with current BG block 
 * configuration. Return SLURM_SUCCESS if they match, else an error 
 * code. Writes bg_block_id into bg_list records.
 */

static int _validate_config_nodes(void)
{
	int rc = SLURM_ERROR;
#ifdef HAVE_BG_FILES
	bg_record_t* bg_record = NULL;	
	bg_record_t* init_bg_record = NULL;
	ListIterator itr_conf;
	ListIterator itr_curr;
	rm_partition_mode_t node_use;
	char tmp_char[256];
	/* read current bg block info into bg_curr_block_list */
	if (read_bg_blocks() == SLURM_ERROR)
		return SLURM_ERROR;
	
	if(!bg_recover) 
		return SLURM_ERROR;
	
	itr_conf = list_iterator_create(bg_list);
	while ((bg_record = (bg_record_t*) list_next(itr_conf))) {
		/* translate hostlist to ranged 
		   string for consistent format
		   search here 
		*/
		node_use = SELECT_COPROCESSOR_MODE; 
		itr_curr = list_iterator_create(bg_curr_block_list);
		while ((init_bg_record = (bg_record_t*) 
			list_next(itr_curr)) 
		       != NULL) {
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
			copy_bg_record(init_bg_record, 
				       bg_record);
			break;
		}
		list_iterator_destroy(itr_curr);
			
		if (!bg_record->bg_block_id) {
			format_node_name(bg_record, tmp_char);	
			info("Block found in bluegene.conf to be "
			     "created: Nodes:%s", 
			     tmp_char);
			rc = SLURM_ERROR;
		} else {
			list_push(bg_found_block_list, bg_record);
			format_node_name(bg_record, tmp_char);
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
	if(bluegene_layout_mode == LAYOUT_DYNAMIC)
		goto finished;
		
	itr_curr = list_iterator_create(bg_curr_block_list);
	while ((init_bg_record = (bg_record_t*) list_next(itr_curr)) 
	       != NULL) {
		debug3("%s %d %d%d%d %d%d%d",
		       init_bg_record->bg_block_id, 
		       init_bg_record->bp_count, 
		       init_bg_record->geo[X],
		       init_bg_record->geo[Y],
		       init_bg_record->geo[Z],
		       DIM_SIZE[X],
		       DIM_SIZE[Y],
		       DIM_SIZE[Z]);
		if ((init_bg_record->geo[X] == DIM_SIZE[X])
		    && (init_bg_record->geo[Y] == DIM_SIZE[Y])
		    && (init_bg_record->geo[Z] == DIM_SIZE[Z]))
		{
			bg_record = (bg_record_t*) 
				xmalloc(sizeof(bg_record_t));
			list_push(bg_list, bg_record);
			list_push(bg_found_block_list, bg_record);
			copy_bg_record(init_bg_record, bg_record);
			bg_record->full_block = 1;
			debug("full system %s",
			      bg_record->bg_block_id);
			format_node_name(bg_record, tmp_char);
			info("Existing: BlockID:%s Nodes:%s Conn:%s",
			     bg_record->bg_block_id, 
			     tmp_char,
			     convert_conn_type(bg_record->conn_type));
			if(((bg_record->state == RM_PARTITION_READY)
			    || (bg_record->state == RM_PARTITION_CONFIGURING))
			   && !block_exist_in_list(bg_booted_block_list, 
						    bg_record))
				list_push(bg_booted_block_list, bg_record);
			break;
		}
	}
	list_iterator_destroy(itr_curr);
		
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
	size_a = strcmp(rec_a->nodes, rec_b->nodes);
	if (size_a < 0)
		return -1;
	else if (size_a > 0)
		return 1;
	
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

static int _delete_old_blocks(void)
{
#ifdef HAVE_BG_FILES
	ListIterator itr_curr, itr_found;
	bg_record_t *found_record = NULL, *init_record = NULL;
	pthread_attr_t attr_agent;
	pthread_t thread_agent;
	int retries;
	List bg_destroy_list = list_create(NULL);

	num_block_to_free = 0;
	num_block_freed = 0;

	if(!bg_recover) {
		if(bg_curr_block_list) {
			itr_curr = list_iterator_create(bg_curr_block_list);
			while ((init_record = 
				(bg_record_t*)list_next(itr_curr))) {
				list_push(bg_destroy_list, init_record);
			}
			list_iterator_destroy(itr_curr);
		} else {
			error("_delete_old_blocks: "
			      "no bg_curr_block_list 1");
			return SLURM_ERROR;
		}
	} else {
		if(bg_curr_block_list) {
			itr_curr = list_iterator_create(bg_curr_block_list);
			while ((init_record = (bg_record_t*) 
				list_next(itr_curr))) {
				if(bg_found_block_list) {
					itr_found = list_iterator_create(
						bg_found_block_list);
					while ((found_record = (bg_record_t*) 
						list_next(itr_found)) 
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
					return SLURM_ERROR;
				}
				if(found_record == NULL) {
					list_push(bg_destroy_list, 
						  init_record);
				}
			}		
			list_iterator_destroy(itr_curr);
		} else {
			error("_delete_old_blocks: "
			      "no bg_curr_block_list 2");
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
		
	retries=30;
	while(num_block_to_free != num_block_freed) {
		_update_bg_record_state(bg_destroy_list);
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
	
	list_destroy(bg_destroy_list);
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

static int _add_block_db(bg_record_t *bg_record, int *block_inx)
{
#ifdef HAVE_BG_FILES
	if(configure_block(bg_record) == SLURM_ERROR) {
		xfree(bg_record);
		error("unable to configure block in api");
		return SLURM_ERROR;
	}
#else
	bg_record->bg_block_id = xmalloc(sizeof(char)*8);
	snprintf(bg_record->bg_block_id, 8, "RMP%d", 
		 (*block_inx)++);
#endif
	return SLURM_SUCCESS;
}
static int _split_block(bg_record_t *bg_record, int procs, int *block_inx) 
{
	bg_record_t *found_record = NULL;
	bool full_bp = false; 
	int small_count = 0;
	int small_size = 0;
	uint16_t num_nodecard = 0, num_quarter = 0;
	int i;
	int node_cnt = 0;
	uint16_t quarter = 0;
	uint16_t nodecard = 0;

	if(bg_record->quarter == (uint16_t) NO_VAL)
		full_bp = true;
	
	if(procs == (procs_per_node/16)) {
		num_nodecard=4;
		if(full_bp)
			num_quarter=3;
	} else if(full_bp) {
		num_quarter = 4;
	} else {
		error("you asked for something that was already this size");
		return SLURM_ERROR;
	}
	debug2("asking for %d 32s from a %d block",
	     num_nodecard, bg_record->node_cnt);
	small_count = num_nodecard+num_quarter; 

	/* break base partition up into 16 parts */
	small_size = bluegene_bp_node_cnt/bluegene_nodecard_node_cnt;
	node_cnt = 0;
	if(!full_bp)
		quarter = bg_record->quarter;
	else
		quarter = 0;
	nodecard = 0;
	for(i=0; i<small_count; i++) {
		if(i == num_nodecard) {
			/* break base partition up into 4 parts */
			small_size = 4;
		}
		
		if(small_size == 4)
			nodecard = (uint16_t)NO_VAL;
		else
			nodecard = i%4; 
		found_record = _create_small_record(bg_record,
						    quarter,
						    nodecard);
		if(block_exist_in_list(bg_list, found_record)) {
			destroy_bg_record(found_record);
		} else {
			if(_add_block_db(found_record, block_inx) 
			   == SLURM_ERROR)
				return SLURM_ERROR;
			list_push(bg_list, found_record);
			print_bg_record(found_record);
		}
		node_cnt += bluegene_bp_node_cnt/small_size;
		if(node_cnt == 128) {
			node_cnt = 0;
			quarter++;
		}
	}
		
	return SLURM_SUCCESS;
}

static int _breakup_blocks(ba_request_t *request, List my_block_list, 
			   int *block_inx)
{
	int rc = SLURM_ERROR;
	bg_record_t *bg_record = NULL;
	ListIterator itr;
	int proc_cnt=0;
	int total_proc_cnt=0;
	uint16_t last_quarter = (uint16_t) NO_VAL;
	char tmp_char[256];
	
	debug2("proc count = %d size = %d",
	      request->procs, request->size);
	
	itr = list_iterator_create(bg_list);			
	while ((bg_record = (bg_record_t *) list_next(itr)) != NULL) {
		if(bg_record->job_running != -1)
			continue;
		if(bg_record->state != RM_PARTITION_FREE)
			continue;
		proc_cnt = bg_record->bp_count * 
			bg_record->cpus_per_bp;
		if(proc_cnt == request->procs) {
			debug2("found it here %s, %s",
			       bg_record->bg_block_id,
			       bg_record->nodes);
			request->save_name = xmalloc(sizeof(char) * 4);
			sprintf(request->save_name, "%d%d%d",
				bg_record->start[X],
				bg_record->start[Y],
				bg_record->start[Z]);
			rc = SLURM_SUCCESS;
			goto finished;
		}
		if(bg_record->node_cnt > bluegene_bp_node_cnt)
			continue;
		if(proc_cnt < request->procs) {
			if(last_quarter != bg_record->quarter){
				last_quarter = bg_record->quarter;
				total_proc_cnt = proc_cnt;
			} else {
				total_proc_cnt += proc_cnt;
			}
			debug2("1 got %d on quarter %d",
			       total_proc_cnt, last_quarter);
			if(total_proc_cnt == request->procs) {
				request->save_name = xmalloc(sizeof(char) * 4);
				sprintf(request->save_name, "%d%d%d",
					bg_record->start[X],
					bg_record->start[Y],
					bg_record->start[Z]);
				if(!my_block_list) {
					rc = SLURM_SUCCESS;
					goto finished;	
				}
						
				bg_record = _create_small_record(
					bg_record,
					last_quarter,
					(uint16_t) NO_VAL);
				if(block_exist_in_list(bg_list, bg_record))
					destroy_bg_record(bg_record);
				else {
					if(_add_block_db(bg_record, block_inx)
					   == SLURM_ERROR)
						return SLURM_ERROR;
					list_push(bg_list, bg_record);
					print_bg_record(bg_record);
				}
				rc = SLURM_SUCCESS;
				goto finished;	
			}
			continue;
		}
		break;
	}
	if(bg_record) {
		debug2("got one on the first pass");
		goto found_one;
	}
	list_iterator_reset(itr);
	last_quarter = (uint16_t) NO_VAL;
	while ((bg_record = (bg_record_t *) list_next(itr)) 
	       != NULL) {
		if(bg_record->job_running != -1)
			continue;
		proc_cnt = bg_record->bp_count * bg_record->cpus_per_bp;
		if(proc_cnt == request->procs) {
			debug2("found it here %s, %s",
			       bg_record->bg_block_id,
			       bg_record->nodes);
			request->save_name = xmalloc(sizeof(char) * 4);
			sprintf(request->save_name, "%d%d%d",
				bg_record->start[X],
				bg_record->start[Y],
				bg_record->start[Z]);
			rc = SLURM_SUCCESS;
			goto finished;
		} 

		if(bg_record->node_cnt > bluegene_bp_node_cnt)
			continue;
		if(proc_cnt < request->procs) {
			if(last_quarter != bg_record->quarter){
				last_quarter = bg_record->quarter;
				total_proc_cnt = proc_cnt;
			} else {
				total_proc_cnt += proc_cnt;
			}
			debug2("got %d on quarter %d",
			       total_proc_cnt, last_quarter);
			if(total_proc_cnt == request->procs) {
				request->save_name = xmalloc(sizeof(char) * 4);
				sprintf(request->save_name, "%d%d%d",
					bg_record->start[X],
					bg_record->start[Y],
					bg_record->start[Z]);
				if(!my_block_list) {
					rc = SLURM_SUCCESS;
					goto finished;	
				}
				bg_record = _create_small_record(
					bg_record,
					last_quarter,
					(uint16_t) NO_VAL);
				if(block_exist_in_list(bg_list, bg_record))
					destroy_bg_record(bg_record);
				else {
					if(_add_block_db(bg_record, block_inx)
					   == SLURM_ERROR)
						return SLURM_ERROR;
					list_push(bg_list, bg_record);
					print_bg_record(bg_record);
				}
				rc = SLURM_SUCCESS;
				goto finished;	
			}
			continue;
		}				
		break;
	}
found_one:
	if(bg_record) {
		format_node_name(bg_record, tmp_char);
			
		debug2("going to split %s, %s",
		      bg_record->bg_block_id,
		      tmp_char);
		request->save_name = xmalloc(sizeof(char) * 4);
		sprintf(request->save_name, "%d%d%d",
			bg_record->start[X],
			bg_record->start[Y],
			bg_record->start[Z]);
		if(!my_block_list) {
			rc = SLURM_SUCCESS;
			goto finished;	
		}
		_split_block(bg_record, request->procs,	block_inx);
		rc = SLURM_SUCCESS;
		goto finished;
	}
	
finished:
	list_iterator_destroy(itr);
		
	return rc;
}

static bg_record_t *_create_small_record(bg_record_t *bg_record, 
					 uint16_t quarter, uint16_t nodecard)
{
	bg_record_t *found_record = NULL;
	int small_size = 4;
	
	found_record = (bg_record_t*) xmalloc(sizeof(bg_record_t));
				
	found_record->job_running = -1;
	found_record->user_name = xstrdup(bg_record->user_name);
	found_record->user_uid = bg_record->user_uid;
	found_record->bg_block_list = list_create(NULL);
	found_record->hostlist = hostlist_create(NULL);
	found_record->nodes = xstrdup(bg_record->nodes);

	process_nodes(found_record);
				
	found_record->conn_type = SELECT_SMALL;
				
	found_record->node_use = SELECT_COPROCESSOR_MODE;
	if(nodecard != (uint16_t) NO_VAL)
		small_size = 16;
	found_record->cpus_per_bp = procs_per_node/small_size;
	found_record->node_cnt = bluegene_bp_node_cnt/small_size;
	found_record->quarter = quarter; 
	found_record->nodecard = nodecard;
			
	return found_record;
}

static int _add_bg_record(List records, blockreq_t *blockreq)
{
	bg_record_t *bg_record = NULL;
	bg_record_t *found_record = NULL;
	ba_node_t *ba_node = NULL;
	ListIterator itr;
	struct passwd *pw_ent = NULL;
	int i, len;
	int small_size = 0;
	int small_count = 0;
	uint16_t quarter = 0;
	uint16_t nodecard = 0;
	int node_cnt = 0;
	
	bg_record = (bg_record_t*) xmalloc(sizeof(bg_record_t));
	
	slurm_conf_lock();
	bg_record->user_name = 
		xstrdup(slurmctld_conf.slurm_user_name);
	bg_record->target_name = 
		xstrdup(slurmctld_conf.slurm_user_name);
	slurm_conf_unlock();
	if((pw_ent = getpwnam(bg_record->user_name)) == NULL) {
		error("getpwnam(%s): %m", bg_record->user_name);
	} else {
		bg_record->user_uid = pw_ent->pw_uid;
	}
	bg_record->bg_block_list = list_create(NULL);		
	bg_record->hostlist = hostlist_create(NULL);
	bg_record->quarter = (uint16_t)NO_VAL;
	bg_record->nodecard = (uint16_t)NO_VAL;
	/* bg_record->boot_state = 0; 	Implicit */
	/* bg_record->state = 0;	Implicit */
	debug2("asking for %s %d %d %s", 
	       blockreq->block, blockreq->quarters, blockreq->nodecards,
	       convert_conn_type(blockreq->conn_type));
	len = strlen(blockreq->block);
	i=0;
	while((blockreq->block[i] != '[' 
	       && (blockreq->block[i] > 57 || blockreq->block[i] < 48)) 
	      && (i<len)) 		
		i++;
	
	if(i<len) {
		len -= i;
		slurm_conf_lock();
		bg_record->nodes = xmalloc(sizeof(char)*
					   (len
					    +strlen(slurmctld_conf.node_prefix)
					    +1));
		
		sprintf(bg_record->nodes, "%s%s", 
			slurmctld_conf.node_prefix, blockreq->block+i);
		slurm_conf_unlock();
			
	} else 
		fatal("BPs=%s is in a weird format", blockreq->block); 
	
	process_nodes(bg_record);
	
	bg_record->node_use = SELECT_COPROCESSOR_MODE;
	bg_record->conn_type = blockreq->conn_type;
	bg_record->cpus_per_bp = procs_per_node;
	bg_record->node_cnt = bluegene_bp_node_cnt * bg_record->bp_count;
	bg_record->job_running = -1;
	
	if(bg_record->conn_type != SELECT_SMALL)
		list_push(records, bg_record);
	else {
		debug("adding a small block");
		if(blockreq->nodecards==0 && blockreq->quarters==0) {
			info("No specs given for this small block, "
			     "I am spliting this block into 4 quarters");
			blockreq->quarters=4;
		}
		i = (blockreq->nodecards*bluegene_nodecard_node_cnt) + 
			(blockreq->quarters*bluegene_quarter_node_cnt);
		if(i != bluegene_bp_node_cnt)
			fatal("There is an error in your bluegene.conf file.\n"
			      "I am unable to request %d nodes in one "
			      "base partition with %d nodes.", 
			      i, bluegene_bp_node_cnt);
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
				found_record = _create_small_record(bg_record,
								    quarter,
								    nodecard);
								 
				list_push(records, found_record);
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

static int _reopen_bridge_log(void)
{
	static FILE *fp = NULL;

	if (bridge_api_file == NULL)
		return SLURM_SUCCESS;

	slurm_mutex_lock(&api_file_mutex);
	if(fp)
		fclose(fp);
	fp = fopen(bridge_api_file, "a");
	
	if (fp == NULL) { 
		error("can't open file for bridgeapi.log at %s: %m", 
		      bridge_api_file);
		slurm_mutex_unlock(&api_file_mutex);
		return SLURM_ERROR;
	}

#ifdef HAVE_BG_FILES
	setSayMessageParams(fp, bridge_api_verb);
#else
	if (fprintf(fp, "bridgeapi.log to write here at level %d\n", 
		    bridge_api_verb) < 20) {
		error("can't write to bridgeapi.log: %m");
		slurm_mutex_unlock(&api_file_mutex);	
		return SLURM_ERROR;
	}
#endif
	slurm_mutex_unlock(&api_file_mutex);	
	
	return SLURM_SUCCESS;
}

