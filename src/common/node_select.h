/*****************************************************************************\
 *  node_select.h - Define node selection plugin functions.
 *
 * $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
 *  UCRL-CODE-217948.
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

#ifndef _NODE_SELECT_H 
#define _NODE_SELECT_H

#include "src/api/node_select_info.h"
#include "src/common/list.h"
#include "src/slurmctld/slurmctld.h"

/*****************************************\
 * GLOBAL SELECT STATE MANGEMENT FUNCIONS *
\*****************************************/

/*
 * Initialize context for node selection plugin
 */
extern int slurm_select_init(void);

/*
 * Terminate plugin and free all associated memory
 */
extern int slurm_select_fini(void);

/*
 * Save any global state information
 * IN dir_name - directory into which the data can be stored
 */
extern int select_g_state_save(char *dir_name);

/*
 * Initialize context for node selection plugin and
 * restore any global state information
 * IN dir_name - directory from which the data can be restored
 */
extern int select_g_state_restore(char *dir_name);

/*
 * Note re/initialization of node record data structure
 * IN node_ptr - current node data
 * IN node_count - number of node entries
 */
extern int select_g_node_init(struct node_record *node_ptr, int node_cnt);

/* 
 * Get select data from a specific node record
 * IN node_pts  - current node record
 * IN cr_info   - type of data to get from the node record
 * IN/OUT data  - the data to get from node record
 */
extern int select_g_get_select_nodeinfo (struct node_record *node_ptr, 
                                         enum select_data_info cr_info, 
					 void *data);

/* 
 * Update select data for a specific node record for a specific job 
 * IN cr_info   - type of data to update for a given job record
 * IN job_ptr - current job record
 */
extern int select_g_update_nodeinfo (struct job_record *job_ptr, 
				     enum select_data_info cr_info);

/* 
 * Get select data from a plugin
 * IN node_pts  - current node record
 * IN cr_info   - type of data to get from the node record (see enum select_data_info)
 * IN/OUT data  - the data to get from node record
 */
extern int select_g_get_info_from_plugin (enum select_data_info cr_info, 
					  void *data);

/* 
 * Alter the node count for a job given the type of system we are on
 * IN/OUT job_desc  - current job desc
 */
extern int select_g_alter_node_cnt (enum select_node_cnt type, void *data);

/*
 * Note re/initialization of partition record data structure
 * IN part_list - list of partition records
 */
extern int select_g_part_init(List part_list);

/* 
 * Note the initialization of job records, issued upon restart of 
 * slurmctld and used to synchronize any job state.
 */
extern int select_g_job_init(List job_list);

/******************************************************\
 * JOB-SPECIFIC SELECT CREDENTIAL MANAGEMENT FUNCIONS *
\******************************************************/

/*
 * Select the "best" nodes for given job from those available
 * IN job_ptr - pointer to job being considered for initiation
 * IN/OUT bitmap - map of nodes being considered for allocation on input,
 *                 map of nodes actually to be assigned on output
 * IN min_nodes - minimum number of nodes to allocate to job
 * IN max_nodes - maximum number of nodes to allocate to job 
 */
extern int select_g_job_test(struct job_record *job_ptr, bitstr_t *bitmap,
			     int min_nodes, int max_nodes, bool test_only);

/*
 * Note initiation of job is about to begin. Called immediately 
 * after select_g_job_test(). Executed from slurmctld.
 * IN job_ptr - pointer to job being initiated
 */
extern int select_g_job_begin(struct job_record *job_ptr);

/*
 * determine if job is ready to execute per the node select plugin
 * IN job_ptr - pointer to job being tested
 * RET -1 on error, 1 if ready to execute, 0 otherwise
 */
extern int select_g_job_ready(struct job_record *job_ptr);

/*
 * Note termination of job is starting. Executed from slurmctld.
 * IN job_ptr - pointer to job being terminated
 */
extern int select_g_job_fini(struct job_record *job_ptr);

/*
 * Suspend a job. Executed from slurmctld.
 * IN job_ptr - pointer to job being suspended
 * RET SLURM_SUCCESS or error code
 */
extern int select_g_job_suspend(struct job_record *job_ptr);

/*
 * Resume a job. Executed from slurmctld.
 * IN job_ptr - pointer to job being resumed
 * RET SLURM_SUCCESS or error code
 */
extern int select_g_job_resume(struct job_record *job_ptr);

/* allocate storage for a select job credential
 * OUT jobinfo - storage for a select job credential
 * RET         - slurm error code
 * NOTE: storage must be freed using select_g_free_jobinfo
 */
extern int select_g_alloc_jobinfo (select_jobinfo_t *jobinfo);

/* fill in a previously allocated select job credential
 * IN/OUT jobinfo  - updated select job credential
 * IN data_type - type of data to enter into job credential
 * IN data - the data to enter into job credential
 */
extern int select_g_set_jobinfo (select_jobinfo_t jobinfo,
				 enum select_data_type data_type, void *data);

/* get data from a select job credential
 * IN jobinfo  - updated select job credential
 * IN data_type - type of data to enter into job credential
 * OUT data - the data to get from job credential, caller must xfree
 *	data for data_tyep == SELECT_DATA_PART_ID
 */
extern int select_g_get_jobinfo (select_jobinfo_t jobinfo,
				 enum select_data_type data_type, void *data);

/* copy a select job credential
 * IN jobinfo - the select job credential to be copied
 * RET        - the copy or NULL on failure
 * NOTE: returned value must be freed using select_g_free_jobinfo
 */
extern select_jobinfo_t select_g_copy_jobinfo(select_jobinfo_t jobinfo);

/* free storage previously allocated for a select job credential
 * IN jobinfo  - the select job credential to be freed
 * RET         - slurm error code
 */
extern int select_g_free_jobinfo  (select_jobinfo_t *jobinfo);
 
/* 
 * Get selected data from a given node for a specific job. 
 * IN node_ptr  - current node record
 * IN job_ptr   - current job record
 * IN cr_info   - type of data to get from the node record
 * IN/OUT data  - the data to get from node record
 */
extern int select_g_get_extra_jobinfo (struct node_record *node_ptr, 
				       struct job_record *job_ptr, 
                                       enum select_data_info cr_info,
                                       void *data);

/* pack a select job credential into a buffer in machine independent form
 * IN jobinfo  - the select job credential to be saved
 * OUT buffer  - buffer with select credential appended
 * RET         - slurm error code
 */
extern int  select_g_pack_jobinfo  (select_jobinfo_t jobinfo, Buf buffer);

/* unpack a select job credential from a buffer
 * OUT jobinfo - the select job credential read
 * IN  buffer  - buffer with select credential read from current pointer loc
 * RET         - slurm error code
 * NOTE: returned value must be freed using select_g_free_jobinfo
 */
extern int  select_g_unpack_jobinfo(select_jobinfo_t jobinfo, Buf buffer);

/* write select job credential to a string
 * IN jobinfo - a select job credential
 * OUT buf    - location to write job credential contents
 * IN size    - byte size of buf
 * IN mode    - print mode, see enum select_print_mode
 * RET        - the string, same as buf
 */
extern char *select_g_sprint_jobinfo(select_jobinfo_t jobinfo,
				     char *buf, size_t size, int mode);

/******************************************************\
 * NODE-SELECT PLUGIN SPECIFIC INFORMATION FUNCTIONS  *
\******************************************************/

/* pack node-select plugin specific information into a buffer in 
 *	machine independent form
 * IN last_update_time - time of latest information consumer has
 * OUT buffer - location to hold the data, consumer must free
 * RET - slurm error code
 */
extern int select_g_pack_node_info(time_t last_query_time, Buf *buffer);
 
/* Unpack node select info from a buffer */
extern int select_g_unpack_node_info(node_select_info_msg_t **
				     node_select_info_msg_pptr, Buf buffer);

/* Free a node select information buffer */
extern int select_g_free_node_info(node_select_info_msg_t **
				   node_select_info_msg_pptr);

#endif /*__SELECT_PLUGIN_API_H__*/
