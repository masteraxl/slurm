/*****************************************************************************\
 *  bluegene.c - blue gene node configuration processing module. 
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

#define _DEBUG 0

char* bgl_conf = BLUEGENE_CONFIG_FILE;

/* Global variables */
rm_BGL_t *bgl;
List bgl_list = NULL;			/* list of bgl_record entries */
List bgl_curr_part_list = NULL;  	/* current bgl partitions */
List bgl_found_part_list = NULL;  	/* found bgl partitions */
char *bluegene_blrts = NULL, *bluegene_linux = NULL, *bluegene_mloader = NULL;
char *bluegene_ramdisk = NULL, *bridge_api_file = NULL;
char *change_numpsets = NULL;
int numpsets;
bool agent_fini = false;
int bridge_api_verb = 0;
		
/* some local functions */
#ifdef HAVE_BGL
static int  _addto_node_list(bgl_record_t *bgl_record, int *start, int *end);
#endif
static void _set_bgl_lists();
static int  _validate_config_nodes(void);
static int  _bgl_record_cmpf_inc(bgl_record_t* rec_a, bgl_record_t* rec_b);
static int  _parse_bgl_spec(char *in_line);
static void _process_nodes(bgl_record_t *bgl_record);
static int  _reopen_bridge_log(void);
static void _strip_13_10(char *line);

/* Initialize all plugin variables */
extern int init_bgl(void)
{
#ifdef HAVE_BGL_FILES
	int rc;
	rm_size3D_t bp_size;
	
	info("Attempting to contact MMCS");
	if ((rc = rm_set_serial(BGL_SERIAL)) != STATUS_OK) {
		fatal("init_bgl: rm_set_serial(): %s", bgl_err_str(rc));
		return SLURM_ERROR;
	}
	
	if ((rc = rm_get_BGL(&bgl)) != STATUS_OK) {
		fatal("init_bgl: rm_get_BGL(): %s", bgl_err_str(rc));
		return SLURM_ERROR;
	}

	if ((rc = rm_get_data(bgl, RM_Msize, &bp_size)) != STATUS_OK) {
		fatal("init_bgl: rm_get_data(): %s", bgl_err_str(rc));
		return SLURM_ERROR;
	}
	verbose("BlueGene configured with %d x %d x %d base partitions",
		bp_size.X, bp_size.Y, bp_size.Z);
	DIM_SIZE[X]=bp_size.X;
	DIM_SIZE[Y]=bp_size.Y;
	DIM_SIZE[Z]=bp_size.Z;
#endif
	pa_init(NULL);

	info("BlueGene plugin loaded successfully");

	return SLURM_SUCCESS;
}

/* Purge all plugin variables */
extern void fini_bgl(void)
{
	_set_bgl_lists();
	
	if (bgl_list) {
		list_destroy(bgl_list);
		bgl_list = NULL;
	}
	
	if (bgl_curr_part_list) {
		list_destroy(bgl_curr_part_list);
		bgl_curr_part_list = NULL;
	}
	
	if (bgl_found_part_list) {
		list_destroy(bgl_found_part_list);
		bgl_found_part_list = NULL;
	}

	xfree(bluegene_blrts);
	xfree(bluegene_linux);
	xfree(bluegene_mloader);
	xfree(bluegene_ramdisk);
	xfree(bridge_api_file);

#ifdef HAVE_BGL_FILES
	if(bgl)
		rm_free_BGL(bgl);
#endif	
	pa_fini();
}

extern void print_bgl_record(bgl_record_t* bgl_record)
{
	if (!bgl_record) {
		error("print_bgl_record, record given is null");
		return;
	}
#if _DEBUG
	info(" bgl_record: ");
	if (bgl_record->bgl_part_id)
		info("\tbgl_part_id: %s", bgl_record->bgl_part_id);
	info("\tnodes: %s", bgl_record->nodes);
	info("\tsize: %d", bgl_record->bp_count);
	info("\tgeo: %dx%dx%d", bgl_record->geo[X], bgl_record->geo[Y], 
	     bgl_record->geo[Z]);
	info("\tlifecycle: %s", convert_lifecycle(bgl_record->part_lifecycle));
	info("\tconn_type: %s", convert_conn_type(bgl_record->conn_type));
	info("\tnode_use: %s", convert_node_use(bgl_record->node_use));
	if (bgl_record->hostlist) {
		char buffer[BUFSIZE];
		hostlist_ranged_string(bgl_record->hostlist, BUFSIZE, buffer);
		info("\thostlist %s", buffer);
	}
	if (bgl_record->bitmap) {
		char bitstring[BITSIZE];
		bit_fmt(bitstring, BITSIZE, bgl_record->bitmap);
		info("\tbitmap: %s", bitstring);
	}
#else
	info("bgl_part_id=%s nodes=%s", bgl_record->bgl_part_id, 
	     bgl_record->nodes);
#endif
}

extern void destroy_bgl_record(void* object)
{
	bgl_record_t* bgl_record = (bgl_record_t*) object;

	if (bgl_record) {
		if(bgl_record->nodes) 
			xfree(bgl_record->nodes);
		if(bgl_record->owner_name)
			xfree(bgl_record->owner_name);
		if(bgl_record->bgl_part_list)
			list_destroy(bgl_record->bgl_part_list);
		if(bgl_record->hostlist)
			hostlist_destroy(bgl_record->hostlist);
		if(bgl_record->bitmap)
			bit_free(bgl_record->bitmap);
		if(bgl_record->bgl_part_id)
			xfree(bgl_record->bgl_part_id);
		
		xfree(bgl_record);
	}
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
extern void sort_bgl_record_inc_size(List records){
	if (records == NULL)
		return;
	list_sort(records, (ListCmpF) _bgl_record_cmpf_inc);
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

	last_mmcs_test = time(NULL) + MMCS_POLL_TIME;
	while (!agent_fini) {
		time_t now = time(NULL);

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
 * Convert a BGL API error code to a string
 * IN inx - error code from any of the BGL Bridge APIs
 * RET - string describing the error condition
 */
extern char *bgl_err_str(status_t inx)
{
#ifdef HAVE_BGL_FILES
	switch (inx) {
	case STATUS_OK:
		return "Status OK";
	case PARTITION_NOT_FOUND:
		return "Partition not found";
	case JOB_NOT_FOUND:
		return "Job not found";
	case BP_NOT_FOUND:
		return "Base partition not found";
	case SWITCH_NOT_FOUND:
		return "Switch not found";
	case JOB_ALREADY_DEFINED:
		return "Job already defined";
	case CONNECTION_ERROR:
		return "Connection error";
	case INTERNAL_ERROR:
		return "Internal error";
	case INVALID_INPUT:
		return "Invalid input";
	case INCOMPATIBLE_STATE:
		return "Incompatible state";
	case INCONSISTENT_DATA:
		return "Inconsistent data";
	}
#endif

	return "?";
}

/*
 * create_static_partitions - create the static partitions that will be used
 *   for scheduling.  
 * IN/OUT part_list - (global, from slurmctld): SLURM's partition 
 *   configurations. Fill in bgl_part_id 
 * RET - success of fitting all configurations
 */
extern int create_static_partitions(List part_list)
{
	int rc = SLURM_SUCCESS;

	ListIterator itr, itr_found;
	bgl_record_t *bgl_record, *found_record;
	
	reset_pa_system();
	
	itr = list_iterator_create(bgl_list);
	while ((bgl_record = (bgl_record_t *) list_next(itr)) != NULL) {
			
		if(bgl_record->bp_count>0 && bgl_record->node_use==SELECT_COPROCESSOR_MODE)
			set_bgl_part(bgl_record->bgl_part_list, 
				     bgl_record->bp_count, 
				     bgl_record->conn_type);
			
	}
	list_iterator_destroy(itr);
	

	itr = list_iterator_create(bgl_list);
	while ((bgl_record = (bgl_record_t *) list_next(itr)) != NULL) {
		itr_found = list_iterator_create(bgl_found_part_list);
		while ((found_record = (bgl_record_t*) list_next(itr_found)) != NULL) {
			if (!strcmp(bgl_record->nodes, found_record->nodes)) {
				break;	/* don't reboot this one */
			}
		}
		list_iterator_destroy(itr_found);
		if(found_record == NULL) {
#ifdef HAVE_BGL_FILES
			//bgl_record->node_use = SELECT_VIRTUAL_NODE_MODE;
			print_bgl_record(bgl_record);
			if((rc = configure_partition(bgl_record)) == SLURM_ERROR)
				return rc;
			
			/* Here we are adding some partitions manually because of the way
			   We want to run the system.  This will need to be changed for 
			   the real system because this is not going to work in the real
			   deal.
			*/
			bgl_record = (bgl_record_t *) list_next(itr);
			if(bgl_record == NULL)
				break;
			print_bgl_record(bgl_record);
			if((rc = configure_partition(bgl_record)) == SLURM_ERROR)
				return rc;
#endif
		}
	}
	list_iterator_destroy(itr);

	/* Here we are adding some partitions manually because of the way
	   We want to run the system.  This will need to be changed for 
	   the real system because this is not going to work in the real
	   deal.
	*/
	
#ifdef HAVE_BGL_FILES
	reset_pa_system();

	bgl_record = (bgl_record_t*) xmalloc(sizeof(bgl_record_t));
	
	bgl_record->nodes = xmalloc(sizeof(char)*13);
	memset(bgl_record->nodes, 0, 13);
	if(DIM_SIZE[X]==1 && DIM_SIZE[Y]==1 && DIM_SIZE[Z]==1)
		sprintf(bgl_record->nodes, "bgl000");
       	else
		sprintf(bgl_record->nodes, "bgl[000x%d%d%d]", DIM_SIZE[X]-1,  
			DIM_SIZE[Y]-1, DIM_SIZE[Z]-1);
       	itr = list_iterator_create(bgl_found_part_list);
	while ((found_record = (bgl_record_t *) list_next(itr)) != NULL) {
		if (!strcmp(bgl_record->nodes, found_record->nodes)) {
			goto no_total;	/* don't reboot this one */
		}
	}
	list_iterator_destroy(itr);
	bgl_record->bgl_part_list = list_create(NULL);			
	bgl_record->hostlist = hostlist_create(NULL);
	_process_nodes(bgl_record);
	list_push(bgl_list, bgl_record);
	
	bgl_record->conn_type = SELECT_TORUS;
	
	set_bgl_part(bgl_record->bgl_part_list, 
		     bgl_record->bp_count, 
		     bgl_record->conn_type);
	bgl_record->node_use = SELECT_COPROCESSOR_MODE;
	print_bgl_record(bgl_record);
	if((rc = configure_partition(bgl_record)) == SLURM_ERROR)
		return rc;
	
	found_record = (bgl_record_t*) xmalloc(sizeof(bgl_record_t));
	list_push(bgl_list, found_record);
			
	found_record->bgl_part_list = bgl_record->bgl_part_list;			
	found_record->hostlist = bgl_record->hostlist;
	found_record->nodes = xstrdup(bgl_record->nodes);
	
	found_record->bp_count = bgl_record->bp_count;
	found_record->switch_count = bgl_record->switch_count;
	found_record->geo[X] = bgl_record->geo[X];
	found_record->geo[Y] = bgl_record->geo[Y];
	found_record->geo[Z] = bgl_record->geo[Z];
			
	found_record->conn_type = bgl_record->conn_type;
	found_record->bitmap = bgl_record->bitmap;
	found_record->node_use = SELECT_VIRTUAL_NODE_MODE;
	print_bgl_record(found_record);
	if((rc = configure_partition(bgl_record)) == SLURM_ERROR)
		return rc;
	
no_total:
	
	rc = SLURM_SUCCESS;
#ifdef _PRINT_PARTS_AND_EXIT
 	itr = list_iterator_create(bgl_list);
	debug("\n\n");
 	while ((found_record = (bgl_record_t *) list_next(itr)) != NULL) {
 		print_bgl_record(found_record);
 	}
 	list_iterator_destroy(itr);
 	exit(0);
#endif	/* _PRINT_PARTS_AND_EXIT */
#endif	/* HAVE_BGL_FILES */

	return rc;
}

#ifdef HAVE_BGL_FILES
static rm_partition_state_t _get_state_partition(pm_partition_id_t part_id)
{
	rm_partition_state_t state = RM_PARTITION_NAV;
	rm_partition_t *part_ptr;
	int rc, j, num_parts;
	rm_partition_list_t *part_list;
	rm_partition_state_flag_t part_state = PARTITION_ALL_FLAG;
	char *name;
		
	if ((rc = rm_get_partitions_info(part_state, &part_list))
	    != STATUS_OK) {
		error("rm_get_partitions(): %s\n", bgl_err_str(rc));
		return state;	
	}

	if ((rc = rm_get_data(part_list, RM_PartListSize, &num_parts))
	    != STATUS_OK) {
		error("rm_get_data(RM_PartListSize): %s\n", bgl_err_str(rc));
		num_parts = 0;
	}

	for (j=0; j<num_parts; j++) {
		if (j) {
			if ((rc = rm_get_data(part_list,
					      RM_PartListNextPart, &part_ptr))
			    != STATUS_OK) {
				error("rm_get_data(RM_PartListNextPart): %s",
				      bgl_err_str(rc));
				break;
			}
		} else {
			if ((rc = rm_get_data(part_list,
					      RM_PartListFirstPart, &part_ptr))
			    != STATUS_OK) {
				error("rm_get_data(RM_PartListFirstPart): %s",
				      bgl_err_str(rc));
				break;
			}
		}
		
		if ((rc = rm_get_data(part_ptr, RM_PartitionID, &name)) 
		    != STATUS_OK) {
			error("rm_get_data(RM_PartitionID): %s",
			      bgl_err_str(rc));
			continue;
		}
		
		if (strcmp(part_id, name))
			continue;
		
		if ((rc = rm_get_data(part_ptr, RM_PartitionState, &state))
		    != STATUS_OK) {
			error("rm_get_data(RM_PartitionState): %s",
			      bgl_err_str(rc));
			state = RM_PARTITION_NAV;
		}
		break;
	}
	if ((rc = rm_free_partition_list(part_list)) != STATUS_OK)
		error("rm_free_partition_list(): %s", bgl_err_str(rc));
	
	return state;
}
#endif

extern int bgl_free_partition(pm_partition_id_t part_id)
{
#ifdef HAVE_BGL_FILES
	rm_partition_state_t state;
	int rc;

	while (1) {
		state = _get_state_partition(part_id);

		if (state != RM_PARTITION_FREE) {
			if ((rc = pm_destroy_partition(part_id)) != STATUS_OK) {
				if(rc == PARTITION_NOT_FOUND) {
					debug("partition %s is not found");
					break;
				}
				error("pm_destroy_partition(%s): %s",
				      part_id, bgl_err_str(rc));
			}
		}

		if ((state == RM_PARTITION_FREE)
		    ||  (state == RM_PARTITION_ERROR))
			break;
		sleep(3);
	}
#endif
	return SLURM_SUCCESS;
}

#ifdef HAVE_BGL
static int _addto_node_list(bgl_record_t *bgl_record, int *start, int *end)
{
	int node_count=0;
	int x,y,z;
	char node_name_tmp[7];
	
	assert(end[X] < DIM_SIZE[X]);
	assert(start[X] >= 0);
	assert(end[Y] < DIM_SIZE[Y]);
	assert(start[Y] >= 0);
	assert(end[Z] < DIM_SIZE[Z]);
	assert(start[Z] >= 0);
	
	for (x = start[X]; x <= end[X]; x++) {
		for (y = start[Y]; y <= end[Y]; y++) {
			for (z = start[Z]; z <= end[Z]; z++) {
				sprintf(node_name_tmp, "bgl%d%d%d", 
					x, y, z);		
				hostlist_push(bgl_record->hostlist, 
					      node_name_tmp);
				list_append(bgl_record->bgl_part_list, 
					    &pa_system_ptr->grid[x][y][z]);
				node_count++;
			}
		}
	}
	return node_count;
}
#endif

static void _set_bgl_lists()
{
	bgl_record_t *bgl_record;
	
	if (bgl_found_part_list) {
		while ((bgl_record = list_pop(bgl_found_part_list)) != NULL) {
		}
	} else
		bgl_found_part_list = list_create(NULL);
	
	if (bgl_curr_part_list){
		while ((bgl_record = list_pop(bgl_curr_part_list)) != NULL){
			destroy_bgl_record(bgl_record);
		}
	} else
		bgl_curr_part_list = list_create(destroy_bgl_record);
	
/* empty the old list before reading new data */
	if (bgl_list) {
		while ((bgl_record = list_pop(bgl_list)) != NULL) {
			destroy_bgl_record(bgl_record);		
		}
	} else
		bgl_list = list_create(destroy_bgl_record);
	
}

/*
 * Match slurm configuration information with current BGL partition 
 * configuration. Return SLURM_SUCCESS if they match, else an error 
 * code. Writes bgl_partition_id into bgl_list records.
 */

static int _validate_config_nodes(void)
{
	int rc = SLURM_ERROR;
#ifdef HAVE_BGL_FILES
	bgl_record_t* record;	/* records from configuration files */
	bgl_record_t* init_record;	/* records from actual BGL config */
	ListIterator itr_conf, itr_curr;
	rm_partition_mode_t node_use;
	
	/* read current bgl partition info into bgl_curr_part_list */
	if (read_bgl_partitions() == SLURM_ERROR)
		return SLURM_ERROR;
	
	if(!bgl_recover) 
		return SLURM_ERROR;
	
	itr_conf = list_iterator_create(bgl_list);
	while ((record = (bgl_record_t*) list_next(itr_conf))) {
		/* translate hostlist to ranged string for consistent format */
        	/* search here */
		node_use = SELECT_COPROCESSOR_MODE; 
		
		itr_curr = list_iterator_create(bgl_curr_part_list);
		while ((init_record = (bgl_record_t*) list_next(itr_curr)) 
		       != NULL) {
				
			if (strcasecmp(record->nodes, init_record->nodes))
				continue;	/* wrong nodes */
			if (record->conn_type != init_record->conn_type)
				continue;      /* wrong conn_type */
			if(record->node_use != init_record->node_use)
				continue;     /* wrong mode */
			record->bgl_part_id = xstrdup(init_record->bgl_part_id);
			break;
		}
		list_iterator_destroy(itr_curr);
		if (!record->bgl_part_id) {
			info("BGL PartitionID:NONE Nodes:%s", record->nodes);
			rc = SLURM_ERROR;
		} else {
			list_push(bgl_found_part_list, record);
			info("BGL PartitionID:%s Nodes:%s Conn:%s Mode:%s",
			     init_record->bgl_part_id, record->nodes,
			     convert_conn_type(record->conn_type),
			     convert_node_use(record->node_use));
		}
	}
	
	list_iterator_destroy(itr_conf);
	if(list_count(bgl_list) == list_count(bgl_curr_part_list))
		rc = SLURM_SUCCESS;
#endif

	return rc;
}

/* 
 * Comparator used for sorting partitions smallest to largest
 * 
 * returns: -1: rec_a >rec_b   0: rec_a == rec_b   1: rec_a < rec_b
 * 
 */
static int _bgl_record_cmpf_inc(bgl_record_t* rec_a, bgl_record_t* rec_b)
{
	if (rec_a->bp_count < rec_b->bp_count)
		return -1;
	else if (rec_a->bp_count > rec_b->bp_count)
		return 1;
	else
		return 0;
}


static int _delete_old_partitions(void)
{
#ifdef HAVE_BGL_FILES
	int rc;
	ListIterator itr_curr, itr_found;
	bgl_record_t *found_record, *init_record;
        
	if(!bgl_recover) {
		itr_curr = list_iterator_create(bgl_curr_part_list);
		while ((init_record = (bgl_record_t*) list_next(itr_curr))) {
			debug("removing the jobs on partition %s\n",
			      init_record->bgl_part_id);
			term_jobs_on_part(init_record->bgl_part_id);
			
			debug("destroying %s\n",(char *)init_record->bgl_part_id);
			bgl_free_partition(init_record->bgl_part_id);
			
			rc = rm_remove_partition(init_record->bgl_part_id);
			if (rc != STATUS_OK) {
				error("rm_remove_partition(%s): %s",
					init_record->bgl_part_id,
					bgl_err_str(rc));
			} else
				debug("done\n");
		}
		list_iterator_destroy(itr_curr);
	} else {
		itr_curr = list_iterator_create(bgl_curr_part_list);
		while ((init_record = (bgl_record_t*) list_next(itr_curr))) {
			itr_found = list_iterator_create(bgl_found_part_list);
			while ((found_record = (bgl_record_t*) list_next(itr_found)) 
			       != NULL) {
				if (!strcmp(init_record->bgl_part_id, 
					    found_record->bgl_part_id)) {
					break;	/* delete this one */
				}
			}
			list_iterator_destroy(itr_found);
			if(found_record == NULL) {
				debug("removing the jobs on partition %s\n",
				      init_record->bgl_part_id);
				term_jobs_on_part(init_record->bgl_part_id);
			
				debug("destroying %s\n",(char *)init_record->bgl_part_id);
				bgl_free_partition(init_record->bgl_part_id);
			
				rc = rm_remove_partition(init_record->
						bgl_part_id);
				if (rc != STATUS_OK) {
					error("rm_remove_partition(%s): %s",
						init_record->bgl_part_id,
						bgl_err_str(rc));
				} else
					debug("done\n");
			}
		}		
		list_iterator_destroy(itr_curr);
	}
#endif	
	return SLURM_SUCCESS;
}

/*
 * Read and process the bluegene.conf configuration file so to interpret what
 * partitions are static/dynamic, torus/mesh, etc.
 */
extern int read_bgl_conf(void)
{
	FILE *bgl_spec_file;	/* pointer to input data file */
	int line_num;		/* line number in input file */
	char in_line[BUFSIZE];	/* input line */
	int i, j, error_code = SLURM_SUCCESS;
	static time_t last_config_update = (time_t) 0;
	struct stat config_stat;

	debug("Reading the bluegene.conf file");

	/* check if config file has changed */
	if (!bgl_conf)
		fatal("bluegene.conf file not defined");
	if (stat(bgl_conf, &config_stat) < 0)
		fatal("can't stat bluegene.conf file %s: %m", bgl_conf);
	if (last_config_update
	    &&  (last_config_update == config_stat.st_mtime)) {
		debug("bluegene.conf unchanged");
		_reopen_bridge_log();
		return SLURM_SUCCESS;
	}
	last_config_update = config_stat.st_mtime; 

	/* initialization */
	/* bgl_conf defined in bgl_node_alloc.h */
	bgl_spec_file = fopen(bgl_conf, "r");
	if (bgl_spec_file == NULL)
		fatal("_read_bgl_conf error opening file %s, %m",
		      bgl_conf);
	
	_set_bgl_lists();	
	
	/* process the data file */
	line_num = 0;
	while (fgets(in_line, BUFSIZE, bgl_spec_file) != NULL) {
		line_num++;
		_strip_13_10(in_line);
		if (strlen(in_line) >= (BUFSIZE - 1)) {
			error("_read_bgl_config line %d, of input file %s "
			      "too long", line_num, bgl_conf);
			fclose(bgl_spec_file);
			return E2BIG;
		}

		/* everything after a non-escaped "#" is a comment */
		/* replace comment flag "#" with an end of string (NULL) */
		/* escape sequence "\#" translated to "#" */
		for (i = 0; i < BUFSIZE; i++) {
			if (in_line[i] == (char) NULL)
				break;
			if (in_line[i] != '#')
				continue;
			if ((i > 0) && (in_line[i - 1] == '\\')) {
				for (j = i; j < BUFSIZE; j++) {
					in_line[j - 1] = in_line[j];
				}
				continue;
			}
			in_line[i] = (char) NULL;
			break;
		}
		
		/* parse what is left, non-comments */
		/* partition configuration parameters */
		error_code = _parse_bgl_spec(in_line);
		
		/* report any leftover strings on input line */
		report_leftover(in_line, line_num);
	}
	fclose(bgl_spec_file);
		
	if (!bluegene_blrts)
		fatal("BlrtsImage not configured in bluegene.conf");
	if (!bluegene_linux)
		fatal("LinuxImage not configured in bluegene.conf");
	if (!bluegene_mloader)
		fatal("MloaderImage not configured in bluegene.conf");
	if (!bluegene_ramdisk)
		fatal("RamDiskImage not configured in bluegene.conf");
	if (!bridge_api_file)
		info("BridgeAPILogFile not configured in bluegene.conf");
	else
		_reopen_bridge_log();	
	if (!numpsets)
		info("Warning: Numpsets not configured in bluegene.conf");
	
	/* Check to see if the configs we have are correct */
	if (_validate_config_nodes() == SLURM_ERROR) { 
		_delete_old_partitions();
	}
	
	/* looking for partitions only I created */
	if (create_static_partitions(NULL) == SLURM_ERROR) {
		/* error in creating the static partitions, so
		 * partitions referenced by submitted jobs won't
		 * correspond to actual slurm partitions/bgl
		 * partitions.
		 */
		fatal("Error, could not create the static partitions");
		return SLURM_ERROR;
	}
	return error_code;
}

/* Explicitly strip out  new-line and carriage-return */
static void _strip_13_10(char *line)
{
	int len = strlen(line);
	int i;

	for(i=0;i<len;i++) {
		if(line[i]==13 || line[i]==10) {
			line[i] = '\0';
			return;
		}
	}
}

/*
 *
 * _parse_bgl_spec - parse the partition specification, build table and 
 *	set values
 * IN/OUT in_line - line from the configuration file, parsed keywords 
 *	and values replaced by blanks
 * RET 0 if no error, error code otherwise
 * Note: Operates on common variables
 * global: part_list - global partition list pointer
 *	default_part - default parameters for a partition
 */
static int _parse_bgl_spec(char *in_line)
{
	int error_code = SLURM_SUCCESS;
	char *nodes = NULL, *conn_type = NULL, *node_use = NULL;
	char *blrts_image = NULL,   *linux_image = NULL;
	char *mloader_image = NULL, *ramdisk_image = NULL;
	char *api_file = NULL;
	int pset_num=-1, api_verb=-1;
	bgl_record_t *bgl_record, *found_record;
	
	error_code = slurm_parser(in_line,
				  "BlrtsImage=", 's', &blrts_image,
				  "LinuxImage=", 's', &linux_image,
				  "MloaderImage=", 's', &mloader_image,
				  "Numpsets=", 'd', &pset_num,
				  "BridgeAPIVerbose=", 'd', &api_verb,
				  "BridgeAPILogFile=", 's', &api_file,
				  "Nodes=", 's', &nodes,
				  "RamDiskImage=", 's', &ramdisk_image,
				  "Type=", 's', &conn_type,
				  "Use=", 's', &node_use,
				  "END");

	if (error_code)
		return SLURM_ERROR;

	/* Process system-wide info */
	if (blrts_image) {
		xfree(bluegene_blrts);
		bluegene_blrts = blrts_image;
		blrts_image = NULL;	/* nothing left to xfree */
	}
	if (linux_image) {
		xfree(bluegene_linux);
		bluegene_linux = linux_image;
		linux_image = NULL;	/* nothing left to xfree */
	}
	if (mloader_image) {
		xfree(bluegene_mloader);
		bluegene_mloader = mloader_image;
		mloader_image = NULL;	/* nothing left to xfree */
	}
	if (ramdisk_image) {
		xfree(bluegene_ramdisk);
		bluegene_ramdisk = ramdisk_image;
		ramdisk_image = NULL;	/* nothing left to xfree */
	}
	if (api_file) {
		xfree(bridge_api_file);
		bridge_api_file = api_file;
		api_file = NULL;	/* nothing left to xfree */
	}
	
	if (pset_num > 0) {
		numpsets = pset_num;
	}
	if (api_verb >= 0) {
		bridge_api_verb = api_verb;
	}

	/* Process node information */
	if (!nodes)
		return SLURM_SUCCESS;	/* not partition line. */
	
	bgl_record = (bgl_record_t*) xmalloc(sizeof(bgl_record_t));
	list_push(bgl_list, bgl_record);
	
	bgl_record->bgl_part_list = list_create(NULL);			
	bgl_record->hostlist = hostlist_create(NULL);
	
	bgl_record->nodes = xstrdup(nodes);
	xfree(nodes);	/* pointer moved, nothing left to xfree */
	
	_process_nodes(bgl_record);
	
	if (!conn_type || !strcasecmp(conn_type,"TORUS"))
		bgl_record->conn_type = SELECT_TORUS;
	else
		bgl_record->conn_type = SELECT_MESH;
	
	if (conn_type)
		xfree(conn_type);	/* pointer moved, nothing left to xfree */

	if (node_use) {
		/* First we check to see if we only want one type of mode */
		if(!strcasecmp(conn_type,"COPROCESSOR"))
			bgl_record->node_use = SELECT_COPROCESSOR_MODE;
		else
			bgl_record->node_use = SELECT_VIRTUAL_NODE_MODE;
	} else {
		/* If not then we will make both. */

		/* this is here to make a co_proc and virtual partition 
		 * just like each other */

		bgl_record->node_use = SELECT_VIRTUAL_NODE_MODE;
			
		found_record = (bgl_record_t*) xmalloc(sizeof(bgl_record_t));
		list_push(bgl_list, found_record);
	
		found_record->bgl_part_list = bgl_record->bgl_part_list;			
		found_record->hostlist = bgl_record->hostlist;
		found_record->nodes = xstrdup(bgl_record->nodes);
	
		found_record->bp_count = bgl_record->bp_count;
		found_record->switch_count = bgl_record->switch_count;
		found_record->geo[X] = bgl_record->geo[X];
		found_record->geo[Y] = bgl_record->geo[Y];
		found_record->geo[Z] = bgl_record->geo[Z];
	
		found_record->conn_type = bgl_record->conn_type;
		found_record->bitmap = bgl_record->bitmap;
		found_record->node_use = SELECT_COPROCESSOR_MODE;
	
	}
#if _DEBUG
	debug("_parse_bgl_spec: added nodes=%s type=%s use=%s", 
	      bgl_record->nodes, 
	      convert_conn_type(bgl_record->conn_type), 
	      convert_node_use(bgl_record->node_use));
#endif
	return SLURM_SUCCESS;
}

static void _process_nodes(bgl_record_t *bgl_record)
{
#ifdef HAVE_BGL
	int j=0, number;
	int start[PA_SYSTEM_DIMENSIONS];
	int end[PA_SYSTEM_DIMENSIONS];
	char buffer[BUFSIZE];
	int funky=0;

	bgl_record->bp_count = 0;

	while (bgl_record->nodes[j] != '\0') {
		if ((bgl_record->nodes[j]   == '[')
		    && (bgl_record->nodes[j+8] == ']')
		    && ((bgl_record->nodes[j+4] == 'x')
			|| (bgl_record->nodes[j+4] == '-'))) {
			j++;
			number = atoi(bgl_record->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j += 4;
			number = atoi(bgl_record->nodes + j);
			end[X] = number / 100;
			end[Y] = (number % 100) / 10;
			end[Z] = (number % 10);
			j += 5;
			bgl_record->bp_count += _addto_node_list(bgl_record, 
								 start, 
								 end);
			if(bgl_record->nodes[j] != ',')
				break;
		} else if((bgl_record->nodes[j] < 58 && bgl_record->nodes[j] > 47) 
			  && bgl_record->nodes[j-1] != '[') {
					
			number = atoi(bgl_record->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j+=3;
			bgl_record->bp_count += _addto_node_list(bgl_record, 
								 start, 
								 start);
			if(bgl_record->nodes[j] != ',')
				break;	
		}
		j++;
	}
	hostlist_ranged_string(bgl_record->hostlist, BUFSIZE, buffer);
	if(strcmp(buffer,bgl_record->nodes)) {
		xfree(bgl_record->nodes);
		bgl_record->nodes = xstrdup(buffer);
	}
	j=0;
	while (bgl_record->nodes[j] != '\0') {
		if ((bgl_record->nodes[j]   == '[')
		    && (bgl_record->nodes[j+8] == ']')
		    && ((bgl_record->nodes[j+4] == 'x')
			|| (bgl_record->nodes[j+4] == '-'))) {
			j++;
			number = atoi(bgl_record->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j += 4;
			number = atoi(bgl_record->nodes + j);
			end[X] = number / 100;
			end[Y] = (number % 100) / 10;
			end[Z] = (number % 10);
			j += 5;			
			if(bgl_record->nodes[j] != ',') 
				break;
			funky=1;
		} else if((bgl_record->nodes[j] < 58 && bgl_record->nodes[j] > 47) 
			  && bgl_record->nodes[j-1] != '[') {
			number = atoi(bgl_record->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j+=3;
			if(bgl_record->nodes[j] != ',')
				break;	
			funky=1;
		}
		
		j++;
	}

	if(!funky) {
		bgl_record->geo[X] = (end[X] - start[X])+1;
		bgl_record->geo[Y] = (end[Y] - start[Y])+1;
		bgl_record->geo[Z] = (end[Z] - start[Z])+1;
	}
		
	if (node_name2bitmap(bgl_record->nodes, 
			     false, 
			     &bgl_record->bitmap)) {
		error("Unable to convert nodes %s to bitmap", 
		      bgl_record->nodes);
	}
#endif
	return;
}

static int _reopen_bridge_log(void)
{
	static FILE *fp = NULL;

	if (bridge_api_file == NULL)
		return SLURM_SUCCESS;

	if(fp)
		fclose(fp);
	fp = fopen(bridge_api_file,"a");
	if (fp == NULL) { 
		error("can't open file for bridgeapi.log at %s: %m", 
		      bridge_api_file);
		return SLURM_ERROR;
	}

#ifdef HAVE_BGL_FILES
	setSayMessageParams(fp, bridge_api_verb);
#else
	if (fprintf(fp, "bridgeapi.log to write here at level %d\n", 
			bridge_api_verb) < 20) {
		error("can't write to bridgeapi.log: %m");
		return SLURM_ERROR;
	}
#endif
		
	return SLURM_SUCCESS;
}
