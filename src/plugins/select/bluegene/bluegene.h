/*****************************************************************************\
 *  bluegene.h - header for blue gene configuration processing module. 
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

#include "src/common/bitstring.h"
#include "src/common/hostlist.h"
#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/slurmctld/slurmctld.h"
#include "src/partition_allocator/partition_allocator.h"

//#include "bgl_job_place.h"
//#include "bgl_job_run.h"

#ifdef HAVE_BGL_FILES
# include "rm_api.h"

/*
 * There is presently a huge amount of untested code to use the APIs.
 * Surround the code with "#ifdef USE_BGL_FILES". When it is confirmed 
 * to work, use "#ifdef HAVE_BGL_FILES" around the code using the APIs.
 */
/* #define USE_BGL_FILES 1 */

#else
  typedef char *   pm_partition_id_t;
  typedef int      rm_connection_type_t;
  typedef int      rm_partition_mode_t;
  typedef uint16_t rm_partition_t;
  typedef char *   rm_BGL_t;
  typedef char *   rm_component_id_t;
  typedef rm_component_id_t rm_bp_id_t;
  typedef int      rm_BP_state_t;
  typedef int      status_t;
  typedef int      rm_partition_state_t;
#endif

#define USER_NAME "da"

/* Global variables */
extern rm_BGL_t *bgl;			/* DB2 pointer */
extern char *bluegene_blrts;
extern char *bluegene_linux;
extern char *bluegene_mloader;
extern char *bluegene_ramdisk;
extern pa_system_t *pa_system_ptr;
extern int DIM_SIZE[PA_SYSTEM_DIMENSIONS];

extern List bgl_curr_part_list; 	/* Initial bgl partition state */
extern List bgl_list;			/* List of configured BGL blocks */
extern bool agent_fini;

typedef int lifecycle_type_t;
enum part_lifecycle {DYNAMIC, STATIC};

typedef struct {
	char *nodes;			/* String of nodes in partition */
	char *owner_name;		/* Owner of partition		*/
	pm_partition_id_t bgl_part_id;	/* ID returned from CMCS	*/
	lifecycle_type_t part_lifecycle;/* either STATIC or DYNAMIC	*/
	rm_partition_state_t state;   	/* the allocated partition   */
	int coord[SYSTEM_DIMENSIONS];   /* bottom left coordinates */
	rm_connection_type_t conn_type;	/* Mesh or Torus or NAV */
	rm_partition_mode_t node_use;	/* either COPROCESSOR or VIRTUAL */
	rm_partition_t *bgl_part;
	List bgl_part_list;
	hostlist_t hostlist;		/* expanded form of hosts */
	int bp_count;
	int switch_count;
	bitstr_t *bitmap;
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
void print_bgl_record(bgl_record_t* record);
void destroy_bgl_record(void* object);

/* Return strings representing blue gene data types */
char *convert_lifecycle(lifecycle_type_t lifecycle);
char *convert_conn_type(rm_connection_type_t conn_type);
char *convert_node_use(rm_partition_mode_t pt);

/* sort a list of bgl_records by size (node count) */
void sort_bgl_record_inc_size(List records);
/* void sort_bgl_record_dec_size(List records); */

/* bluegene_agent - detached thread periodically tests status of bluegene 
 * nodes and switches */
void *bluegene_agent(void *args);

/*
 * Convert a BGL API error code to a string
 * IN inx - error code from any of the BGL Bridge APIs
 * RET - string describing the error condition
 */
char *bgl_err_str(status_t inx);

/*
 * create_static_partitions - create the static partitions that will be used
 *   for scheduling.
 * IN/OUT part_list - (global, from slurmctld): SLURM's partition 
 *   configurations. Fill in bgl_part_id                 
 * RET - success of fitting all configurations
 */
int create_static_partitions(List part_list);

int read_bgl_conf();

/* partition_sys.c */
/*****************************************************/
int configure_partition(bgl_record_t * bgl_conf_record);
int read_bgl_partitions();

/* state_test.c */
/*****************************************************/

/* Test for nodes that are DOWN in BlueGene database, 
 * if so DRAIN them in SLURM */
void test_down_nodes(void);

/* Test for switches that are DOWN in BlueGene database,  
 * if so DRAIN them in SLURM and configure their base partition DOWN */
void test_down_switches(void);


/* bgl_job_place.c */
/*****************************************************/
int submit_job(struct job_record *job_ptr, bitstr_t *bitmap,
	       int min_nodes, int max_nodes);


/* bgl_job_run.c */
/*****************************************************/
/*
 * Perform any setup required to initiate a job
 * job_ptr IN - pointer to the job being initiated
 * RET - SLURM_SUCCESS or an error code 
 *
 * NOTE: This happens in parallel with srun and slurmd spawning 
 * the job. A prolog script is expected to defer initiation of 
 * the job script until the BGL block is available for use.
 */
int start_job(struct job_record *job_ptr);

/*
 * Perform any work required to terminate a jobs on a partition
 * bgl_part_id IN - partition name
 * RET - SLURM_SUCCESS or an error code
 *
 * NOTE: This happens when new partitions are created and we 
 * need to clean up jobs on them.
 */
int term_jobs_on_part(pm_partition_id_t bgl_part_id);

/* 
 * Perform any work required to terminate a job
 * job_ptr IN - pointer to the job being terminated
 * RET - SLURM_SUCCESS or an error code
 *
 * NOTE: This happens in parallel with srun and slurmd terminating
 * the job. Insure that this function, mpirun and the epilog can 
 * all deal with termination race conditions.
 */
int term_job(struct job_record *job_ptr);

/*
 * Synchronize BGL block state to that of currently active jobs.
 * This can recover from slurmctld crashes when partition ownership
 * changes were queued
 */
int sync_jobs(List job_list);

/* bgl_switch_connections.c */
/*****************************************************/
int configure_partition_switches(bgl_record_t * bgl_conf_record);

#endif /* _BLUEGENE_H_ */
