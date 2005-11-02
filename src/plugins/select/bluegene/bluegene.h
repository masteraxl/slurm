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
#include "src/partition_allocator/partition_allocator.h"
#include "src/plugins/select/bluegene/wrap_rm_api.h"

typedef int lifecycle_type_t;

enum part_lifecycle {DYNAMIC, STATIC};

typedef struct bgl_record {
	char *nodes;			/* String of nodes in partition */
	char *user_name;		/* user using the partition */
	char *target_name;		/* when a partition is freed this 
					   is the name of the user we 
					   want on the partition */
	uid_t user_uid;   		/* Owner of partition uid	*/
	pm_partition_id_t bgl_part_id;	/* ID returned from MMCS	*/
	lifecycle_type_t part_lifecycle;/* either STATIC or DYNAMIC	*/
	rm_partition_state_t state;   	/* the allocated partition   */
	int start[PA_SYSTEM_DIMENSIONS];/* start node */
	int geo[PA_SYSTEM_DIMENSIONS];  /* geometry */
	rm_connection_type_t conn_type;	/* Mesh or Torus or NAV */
	rm_partition_mode_t node_use;	/* either COPROCESSOR or VIRTUAL */
	rm_partition_t *bgl_part;       /* structure to hold info from db2 */
	List bgl_part_list;             /* node list of blocks in partition */
	hostlist_t hostlist;		/* expanded form of hosts */
	int bp_count;                   /* size */
	int switch_count;               /* number of switches used. */
	int boot_state;                 /* check to see if boot failed. 
					   -1 = fail, 
					   0 = not booting, 
					   1 = booting */
	int boot_count;                 /* number of attemts boot attempts */
	bitstr_t *bitmap;               /* bitmap to check the name 
					   of partition */
	int full_partition;             /* wether or not partition is the full
					   partition */
	int job_running;                /* signal if there is a job running 
					   on the partition */
	int cnodes_per_bp;              /* count of cnodes per Base part */
	int quarter;                    /* used for small partitions 
					   determine quarter of BP */
} bgl_record_t;

typedef struct {
	int source;
	int target;
} bgl_conn_t;

typedef struct {
	int dim;
	List conn_list;
} bgl_switch_t;

typedef struct {
	int *coord;
	int used;
	List switch_list;
} bgl_bp_t;


/* Global variables */
extern rm_BGL_t *bgl;
extern char *bluegene_blrts;
extern char *bluegene_linux;
extern char *bluegene_mloader;
extern char *bluegene_ramdisk;
extern char *bridge_api_file;
extern int numpsets;
extern pa_system_t *pa_system_ptr;
extern time_t last_bgl_update;
extern List bgl_curr_part_list; 	/* Initial bgl partition state */
extern List bgl_list;			/* List of configured BGL blocks */
extern bool agent_fini;
extern pthread_mutex_t part_state_mutex;
extern int num_part_to_free;
extern int num_part_freed;
extern int partitions_are_created;
extern int procs_per_node;
extern bgl_record_t *full_system_partition;


#define MAX_PTHREAD_RETRIES  1

#include "bgl_part_info.h"
#include "bgl_job_place.h"
#include "bgl_job_run.h"
#include "state_test.h"
/*
 * bgl_conf_record is used to store the elements read from the bluegene.conf
 * file and is loaded by init().
 */
/* typedef struct bgl_conf_record { */
/* 	char* nodes; */
/* 	rm_connection_type_t conn_type;/\* Mesh or Torus or NAV *\/ */
/* 	rm_partition_mode_t node_use; */
/* 	rm_partition_t *bgl_part; */
/* } bgl_conf_record_t; */



/* bluegene.c */
/**********************************************/

/* Initialize all plugin variables */
extern int init_bgl(void);

/* Purge all plugin variables */
extern void fini_bgl(void);

/* Log a bgl_record's contents */
extern void print_bgl_record(bgl_record_t* record);
extern void destroy_bgl_record(void* object);

/* return bgl_record from bgl_list */
extern bgl_record_t *find_bgl_record(char *bgl_part_id);

/* change username of a partition bgl_record_t target_name needs to be 
   updated before call of function. 
*/
extern int update_partition_user(bgl_record_t *bgl_part_id); 

/* remove all users from a partition but what is in user_name */
/* Note return codes */
#define REMOVE_USER_ERR  -1
#define REMOVE_USER_NONE  0
#define REMOVE_USER_FOUND 2
extern int remove_all_users(char *bgl_part_id, char *user_name);
extern void set_part_user(bgl_record_t *bgl_record);

/* Return strings representing blue gene data types */
extern char *convert_lifecycle(lifecycle_type_t lifecycle);
extern char *convert_conn_type(rm_connection_type_t conn_type);
extern char *convert_node_use(rm_partition_mode_t pt);

/* sort a list of bgl_records by size (node count) */
extern void sort_bgl_record_inc_size(List records);

/* bluegene_agent - detached thread periodically tests status of bluegene 
 * nodes and switches */
extern void *bluegene_agent(void *args);

/*
 * Convert a BGL API error code to a string
 * IN inx - error code from any of the BGL Bridge APIs
 * RET - string describing the error condition
 */
extern char *bgl_err_str(status_t inx);

/*
 * create_static_partitions - create the static partitions that will be used
 *   for scheduling.
 * IN/OUT part_list - (global, from slurmctld): SLURM's partition 
 *   configurations. Fill in bgl_part_id                 
 * RET - success of fitting all configurations
 */
extern int create_static_partitions(List part_list);

extern int bgl_free_partition(bgl_record_t *bgl_record);
extern void *mult_free_part(void *args);
extern void *mult_destroy_part(void *args);
extern int read_bgl_conf(void);

/* partition_sys.c */
/*****************************************************/
extern int configure_partition(bgl_record_t * bgl_conf_record);
extern int read_bgl_partitions();

/* bgl_switch_connections.c */
/*****************************************************/
extern int configure_small_partition(bgl_record_t *bgl_record);
extern int configure_partition_switches(bgl_record_t * bgl_conf_record);


#endif /* _BLUEGENE_H_ */
 
