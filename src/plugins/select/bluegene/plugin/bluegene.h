/*****************************************************************************\
 *  bluegene.h - header for blue gene configuration processing module. 
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Dan Phung <phung4@llnl.gov> and Danny Auble <da@llnl.gov>
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

#ifndef _BLUEGENE_H_
#define _BLUEGENE_H_

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdlib.h>
#include <sys/stat.h>
#include <pwd.h>

#include "src/common/bitstring.h"
#include "src/common/hostlist.h"
#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/slurmctld/slurmctld.h"
#include "../block_allocator/block_allocator.h"
#include "../wrap_rm_api.h"

typedef int lifecycle_type_t;

enum block_lifecycle {DYNAMIC, STATIC};

typedef enum bg_layout_type {
	LAYOUT_STATIC,  /* no overlaps, except for full system block
			   blocks never change */
	LAYOUT_OVERLAP, /* overlaps permitted, must be defined in 
			   bluegene.conf file */
	LAYOUT_DYNAMIC	/* slurm will make all blocks */
}bg_layout_t;

typedef struct bg_record {
	pm_partition_id_t bg_block_id;	/* ID returned from MMCS	*/
	char *nodes;			/* String of nodes in block */
	char *user_name;		/* user using the block */
	char *target_name;		/* when a block is freed this 
					   is the name of the user we 
					   want on the block */
	int full_block;                 /* wether or not block is the full
					   block */
	uid_t user_uid;   		/* Owner of block uid	*/
	lifecycle_type_t block_lifecycle;/* either STATIC or DYNAMIC	*/
	rm_partition_state_t state;   	/* the allocated block   */
	int start[BA_SYSTEM_DIMENSIONS];/* start node */
	int geo[BA_SYSTEM_DIMENSIONS];  /* geometry */
	rm_connection_type_t conn_type;	/* Mesh or Torus or NAV */
	rm_partition_mode_t node_use;	/* either COPROCESSOR or VIRTUAL */
	rm_partition_t *bg_block;       /* structure to hold info from db2 */
	List bg_block_list;             /* node list of blocks in block */
	hostlist_t hostlist;		/* expanded form of hosts */
	int bp_count;                   /* size */
	int switch_count;               /* number of switches used. */
	int boot_state;                 /* check to see if boot failed. 
					   -1 = fail, 
					   0 = not booting, 
					   1 = booting */
	int boot_count;                 /* number of attemts boot attempts */
	bitstr_t *bitmap;               /* bitmap to check the name 
					   of block */
	int job_running;                /* job id if there is a job running 
					   on the block */
	int cpus_per_bp;                /* count of cpus per base part */
	int node_cnt;                   /* count of nodes per block */
	int quarter;                    /* used for small blocks 
					   determine quarter of BP */
	int segment;                    /* used for small blocks 
					   determine segment of quarter */
} bg_record_t;

typedef struct {
	int source;
	int target;
} bg_conn_t;

typedef struct {
	int dim;
	List conn_list;
} bg_switch_t;

typedef struct {
	int *coord;
	int used;
	List switch_list;
} bg_bp_t;


/* Global variables */
extern rm_BGL_t *bg;
extern char *bluegene_blrts;
extern char *bluegene_linux;
extern char *bluegene_mloader;
extern char *bluegene_ramdisk;
extern char *bridge_api_file;
extern bg_layout_t bluegene_layout_mode;
extern int bluegene_numpsets;
extern int bluegene_mp_node_cnt;
extern int bluegene_nc_node_cnt;
extern ba_system_t *ba_system_ptr;
extern time_t last_bg_update;

extern List bg_curr_block_list; 	/* Initial bg block state */
extern List bg_list;			/* List of configured BG blocks */
extern List bg_job_block_list;  	/* jobs running in these blocks */
extern List bg_booted_block_list;  	/* blocks that are booted */

extern bool agent_fini;
extern pthread_mutex_t block_state_mutex;
extern int num_block_to_free;
extern int num_block_freed;
extern int blocks_are_created;
extern int procs_per_node;

#define MAX_PTHREAD_RETRIES  1
#define MAX_AGENT_COUNT      30

#include "bg_block_info.h"
#include "bg_job_place.h"
#include "bg_job_run.h"
#include "state_test.h"

/* bluegene.c */
/**********************************************/

/* Initialize all plugin variables */
extern int init_bg(void);

/* Purge all plugin variables */
extern void fini_bg(void);

/* Log a bg_record's contents */
extern void print_bg_record(bg_record_t *record);
extern void destroy_bg_record(void *object);
extern void copy_bg_record(bg_record_t *fir_record, bg_record_t *sec_record);

/* return bg_record from bg_list */
extern bg_record_t *find_bg_record(char *bg_block_id);

/* change username of a block bg_record_t target_name needs to be 
   updated before call of function. 
*/
extern int update_block_user(bg_record_t *bg_block_id); 
extern int format_node_name(bg_record_t *bg_record, char tmp_char[]);
extern bool blocks_overlap(bg_record_t *rec_a, bg_record_t *rec_b);


/* remove all users from a block but what is in user_name */
/* Note return codes */
#define REMOVE_USER_ERR  -1
#define REMOVE_USER_NONE  0
#define REMOVE_USER_FOUND 2
extern int remove_all_users(char *bg_block_id, char *user_name);
extern void set_block_user(bg_record_t *bg_record);

/* Return strings representing blue gene data types */
extern char *convert_lifecycle(lifecycle_type_t lifecycle);
extern char *convert_conn_type(rm_connection_type_t conn_type);
extern char *convert_node_use(rm_partition_mode_t pt);

/* sort a list of bg_records by size (node count) */
extern void sort_bg_record_inc_size(List records);

/* bluegene_agent - detached thread periodically tests status of bluegene 
 * nodes and switches */
extern void *bluegene_agent(void *args);

/*
 * Convert a BG API error code to a string
 * IN inx - error code from any of the BG Bridge APIs
 * RET - string describing the error condition
 */
extern char *bg_err_str(status_t inx);

/*
 * create_*_block(s) - functions for creating blocks that will be used
 *   for scheduling.
 * RET - success of fitting all configurations
 */
extern int create_defined_blocks(bg_layout_t overlapped);
extern int create_dynamic_block(ba_request_t *request, List my_block_list);
extern int create_full_system_block();

extern int bg_free_block(bg_record_t *bg_record);
extern int remove_from_bg_list(List my_bg_list, bg_record_t *bg_record);
extern void *mult_free_block(void *args);
extern void *mult_destroy_block(void *args);
extern int free_block_list(List delete_list);
extern int read_bg_conf(void);

/* block_sys.c */
/*****************************************************/
extern int configure_block(bg_record_t * bg_conf_record);
extern int read_bg_blocks();

/* bg_switch_connections.c */
/*****************************************************/
extern int configure_small_block(bg_record_t *bg_record);
extern int configure_block_switches(bg_record_t * bg_conf_record);


#endif /* _BLUEGENE_H_ */
 
