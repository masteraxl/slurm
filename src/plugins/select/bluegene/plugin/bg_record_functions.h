/*****************************************************************************\
 *  bg_record_functions.h - header for creating blocks in a static environment.
 *
 *  $Id: bg_record_functions.h 12954 2008-01-04 20:37:49Z da $
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

#ifndef _BLUEGENE_BG_RECORD_FUNCTIONS_H_
#define _BLUEGENE_BG_RECORD_FUNCTIONS_H_

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
#include "src/common/node_select.h"
#include "src/common/parse_time.h"
#include "src/slurmctld/slurmctld.h"

#include "../block_allocator/block_allocator.h"

typedef struct bg_record {
	pm_partition_id_t bg_block_id;	/* ID returned from MMCS	*/
	char *nodes;			/* String of nodes in block */
	char *ionodes; 		        /* String of ionodes in block
					 * NULL if not a small block*/
	char *user_name;		/* user using the block */
	char *target_name;		/* when a block is freed this 
					   is the name of the user we 
					   want on the block */
	int full_block;                 /* whether or not block is the full
					   block */
	int modifying;                  /* flag to say the block is
					   being modified or not at
					   job launch usually */
	uid_t user_uid;   		/* Owner of block uid	*/
	rm_partition_state_t state;     /* Current state of the block */
	int start[BA_SYSTEM_DIMENSIONS];/* start node */
	uint16_t geo[BA_SYSTEM_DIMENSIONS];  /* geometry */
	rm_connection_type_t conn_type;  /* MESH or Torus or NAV */
	rm_partition_mode_t node_use;	 /* either COPROCESSOR or VIRTUAL */
	rm_partition_t *bg_block;       /* structure to hold info from db2 */
	List bg_block_list;             /* node list of blocks in block */
	int bp_count;                   /* size */
	int switch_count;               /* number of switches used. */
	int boot_state;                 /* check to see if boot failed. 
					   -1 = fail, 
					   0 = not booting, 
					   1 = booting */
	int boot_count;                 /* number of attemts boot attempts */
	bitstr_t *bitmap;               /* bitmap to check the name 
					   of block */
	bitstr_t *ionode_bitmap;        /* for small blocks bitmap to
					   keep track which ionodes we
					   are on.  NULL if not a small block*/
	struct job_record *job_ptr;	/* pointer to job running on
					 * block or NULL if no job */
	int job_running;                /* job id of job running of if
					 * block is in an error state
					 * BLOCK_ERROR_STATE */
	int cpus_per_bp;                /* count of cpus per base part */
	uint32_t node_cnt;              /* count of nodes per block */
	uint16_t quarter;               /* used for small blocks 
					   determine quarter of BP */
	uint16_t nodecard;              /* used for small blocks 
					   determine nodecard of quarter */
	char *blrtsimage;              /* BlrtsImage for this block */
	char *linuximage;              /* LinuxImage for this block */
	char *mloaderimage;            /* mloaderImage for this block */
	char *ramdiskimage;            /* RamDiskImage for this block */
	struct bg_record *original;    /* if this is a copy this is a
					  pointer to the original */
} bg_record_t;

/* Log a bg_record's contents */
extern void print_bg_record(bg_record_t *record);
extern void destroy_bg_record(void *object);
extern int block_exist_in_list(List my_list, bg_record_t *bg_record);
extern void process_nodes(bg_record_t *bg_reord, bool startup);
extern List copy_bg_list(List in_list);
extern void copy_bg_record(bg_record_t *fir_record, bg_record_t *sec_record);

/* return bg_record from a bg_list */
extern bg_record_t *find_bg_record_in_list(List my_list, char *bg_block_id);

/* change username of a block bg_record_t target_name needs to be 
   updated before call of function. 
*/
extern int update_block_user(bg_record_t *bg_block_id, int set); 
extern void drain_as_needed(bg_record_t *bg_record, char *reason);

extern int set_ionodes(bg_record_t *bg_record);

extern int add_bg_record(List records, List used_nodes, blockreq_t *blockreq);

extern int format_node_name(bg_record_t *bg_record, char *buf, int buf_size);

#endif /* _BLUEGENE_BG_RECORD_FUNCTIONS_H_ */
