/****************************************************************************\
 *  slurmdbd_defs.h - definitions used for Slurm DBD RPCs
 *****************************************************************************
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
 *  UCRL-CODE-226842.
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

#ifndef _SLURMDBD_DEFS_H
#define _SLURMDBD_DEFS_H

#if HAVE_CONFIG_H 
#  include "config.h"
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  else
#    if HAVE_STDINT_H
#      include <stdint.h>
#    endif
#  endif			/* HAVE_INTTYPES_H */
#else				/* !HAVE_CONFIG_H */
#  include <inttypes.h>
#endif				/*  HAVE_CONFIG_H */

#include <slurm/slurm.h>

#include "src/common/pack.h"
#include "src/common/list.h"
#include "src/common/slurm_accounting_storage.h"

/* Increment SLURMDBD_VERSION if any of the RPCs change */
#define SLURMDBD_VERSION 01

/* SLURM DBD message types */
typedef enum {
	DBD_INIT = 1400,	/* Connection initialization		*/
	DBD_ADD_ACCOUNTS,       /* Add new account to the mix           */
	DBD_ADD_ACCOUNT_COORDS, /* Add new coordinatior to an account   */
	DBD_ADD_ASSOCS,         /* Add new association to the mix       */
	DBD_ADD_CLUSTERS,       /* Add new cluster to the mix           */
	DBD_ADD_USERS,          /* Add new user to the mix              */
	DBD_CLUSTER_PROCS,	/* Record total processors on cluster	*/
	DBD_GET_ACCOUNTS,	/* Get account information		*/
	DBD_GET_ASSOCS,         /* Get assocation information	*/
	DBD_GET_ASSOC_DAY,	/* Get assoc daily usage information	*/
	DBD_GET_ASSOC_HOUR,     /* Get assoc hourly usage information	*/
	DBD_GET_ASSOC_MONTH,	/* Get assoc monthly usage information	*/
	DBD_GET_CLUSTERS,	/* Get account information		*/
	DBD_GET_CLUSTER_HOUR,	/* Get cluster hourly usage information	*/
	DBD_GET_CLUSTER_DAY,	/* Get cluster hourly usage information	*/
	DBD_GET_CLUSTER_MONTH,	/* Get cluster hourly usage information	*/
	DBD_GET_JOBS,		/* Get job information			*/
	DBD_GET_USERS,  	/* Get account information		*/
	DBD_GOT_ACCOUNTS,	/* Response to DBD_GET_ACCOUNTS		*/
	DBD_GOT_ASSOCS, 	/* Response to DBD_GET_ASSOCS   	*/
	DBD_GOT_ASSOC_DAY,	/* Response to DBD_GET_ASSOC_DAY	*/
	DBD_GOT_ASSOC_HOUR,	/* Response to DBD_GET_ASSOC_HOUR	*/
	DBD_GOT_ASSOC_MONTH,	/* Response to DBD_GET_ASSOC_MONTH	*/
	DBD_GOT_CLUSTERS,	/* Response to DBD_GET_CLUSTERS		*/
	DBD_GOT_CLUSTER_DAY,	/* Response to DBD_GET_CLUSTER_DAY	*/
	DBD_GOT_CLUSTER_HOUR,	/* Response to DBD_GET_CLUSTER_HOUR	*/
	DBD_GOT_CLUSTER_MONTH,	/* Response to DBD_GET_CLUSTER_MONTH	*/
	DBD_GOT_JOBS,		/* Response to DBD_GET_JOBS		*/
	DBD_GOT_USERS,  	/* Response to DBD_GET_USERS		*/
	DBD_JOB_COMPLETE,	/* Record job completion 		*/
	DBD_JOB_START,		/* Record job starting			*/
	DBD_JOB_START_RC,	/* return db_index from job insertion 	*/
	DBD_JOB_SUSPEND,	/* Record job suspension		*/
	DBD_MODIFY_ACCOUNTS,    /* Modify existing account              */
	DBD_MODIFY_ASSOCS,      /* Modify existing association          */
	DBD_MODIFY_CLUSTERS,    /* Modify existing cluster              */
	DBD_MODIFY_USERS,       /* Modify existing user                 */
	DBD_MODIFY_USER_ADMIN_LEVEL,/* Modify existing user             */
	DBD_NODE_STATE,		/* Record node state transition		*/
	DBD_RC,			/* Return code from operation		*/
	DBD_REMOVE_ACCOUNTS,    /* Remove existing account              */
	DBD_REMOVE_ACCOUNT_COORDS,/* Remove existing coordinatior from
				   * an account */
	DBD_REMOVE_ASSOCS,      /* Remove existing association          */
	DBD_REMOVE_CLUSTERS,    /* Remove existing cluster              */
	DBD_REMOVE_USERS,        /* Remove existing user                 */
	DBD_STEP_COMPLETE,	/* Record step completion		*/
	DBD_STEP_START		/* Record step starting			*/
} slurmdbd_msg_type_t;

/*****************************************************************************\
 * Slurm DBD protocol data structures
\*****************************************************************************/

typedef struct slurmdbd_msg {
	uint16_t msg_type;	/* see slurmdbd_msg_type_t above */
	void * data;		/* pointer to a message type below */
} slurmdbd_msg_t;

typedef struct {
	char *acct;
	acct_user_cond_t *cond;
} dbd_acct_coord_msg_t;

typedef struct dbd_cluster_procs_msg {
	char *cluster_name;	/* name of cluster */
	uint32_t proc_count;	/* total processor count */
	time_t event_time;	/* time of transition */
} dbd_cluster_procs_msg_t;

typedef struct {
	void *cond; /* this could be anything based on the type types
		     * are defined in slurm_accounting_storage.h
		     * *_cond_t */
} dbd_cond_msg_t;

typedef struct {
	void *rec;
	time_t start;
	time_t end;
} dbd_usage_msg_t;

typedef struct dbd_get_jobs_msg {
	char *cluster_name; /* name of cluster to query */
	uint32_t gid;      /* group id */
	List selected_steps; /* List of jobacct_selected_step_t *'s */
	List selected_parts; /* List of char *'s */
	char *user;        /* user name */
} dbd_get_jobs_msg_t;

typedef struct dbd_init_msg {
	uint16_t version;	/* protocol version */
	uint32_t uid;		/* UID originating connection,
				 * filled by authtentication plugin*/
} dbd_init_msg_t;

typedef struct dbd_job_comp_msg {
	uint32_t assoc_id;	/* accounting association id needed to
				 * find job record in db */
	uint32_t db_index;	/* index into the db for this job */
	time_t   end_time;	/* job termintation time */
	uint32_t exit_code;	/* job exit code or signal */
	uint32_t job_id;	/* job ID */
	uint16_t job_state;	/* job state */
	char *   nodes;		/* hosts allocated to the job */
	time_t   start_time;	/* job start time */
	time_t   submit_time;	/* job submit time needed to find job
				 * record in db */
} dbd_job_comp_msg_t;

typedef struct dbd_job_start_msg {
	uint32_t alloc_cpus;	/* count of allocated processors */
	uint32_t assoc_id;	/* accounting association id */
	char *   block_id;      /* Bluegene block id */
	time_t   eligible_time;	/* time job becomes eligible to run */
	uint32_t gid;	        /* group ID */
	uint32_t job_id;	/* job ID */
	uint16_t job_state;	/* job state */
	char *   name;		/* job name */
	char *   nodes;		/* hosts allocated to the job */
	char *   partition;	/* partition job is running on */
	uint32_t priority;	/* job priority */
	uint32_t req_cpus;	/* count of req processors */
	time_t   start_time;	/* job start time */
	time_t   submit_time;	/* job submit time */
} dbd_job_start_msg_t;

typedef struct dbd_job_start_rc_msg {
	uint32_t db_index;	/* db_index */
	uint32_t return_code;
} dbd_job_start_rc_msg_t;

typedef struct dbd_job_suspend_msg {
	uint32_t assoc_id;	/* accounting association id needed
				 * to find job record in db */
	uint32_t db_index;	/* index into the db for this job */
	uint32_t job_id;	/* job ID needed to find job record
				 * in db */
	uint16_t job_state;	/* job state */
	time_t   submit_time;	/* job submit time needed to find job record
				 * in db */
	time_t   suspend_time;	/* job suspend or resume time */
} dbd_job_suspend_msg_t;

typedef struct {
	List my_list; /* this list could be of any type as long as it
			* is handled correctly on both ends */
} dbd_list_msg_t;

typedef struct {
	void *cond;
	void *rec;
} dbd_modify_msg_t;

#define DBD_NODE_STATE_DOWN  1
#define DBD_NODE_STATE_UP    2
typedef struct dbd_node_state_msg {
	char *cluster_name;	/* name of cluster */
	time_t event_time;	/* time of transition */
	char *hostlist;		/* name of hosts */
	uint16_t new_state;	/* new state of host, see DBD_NODE_STATE_* */
	char *reason;		/* explanation for the node's state */
} dbd_node_state_msg_t;

typedef struct dbd_rc_msg {
	char *   comment;	/* reason for failure */
	uint32_t return_code;
} dbd_rc_msg_t;

typedef struct dbd_step_comp_msg {
	uint32_t assoc_id;	/* accounting association id */
	uint32_t db_index;	/* index into the db for this job */
	time_t   end_time;	/* job termintation time */
	jobacctinfo_t *jobacct; /* status info */
	uint32_t job_id;	/* job ID */
	uint32_t req_uid;	/* requester user ID */
	time_t   start_time;	/* step start time */
	time_t   job_submit_time;/* job submit time needed to find job record
				  * in db */
	uint32_t step_id;	/* step ID */
	uint32_t total_procs;	/* count of allocated processors */
} dbd_step_comp_msg_t;

typedef struct dbd_step_start_msg {
	uint32_t assoc_id;	/* accounting association id */
	uint32_t db_index;	/* index into the db for this job */
	uint32_t job_id;	/* job ID */
	char *   name;		/* step name */
	char *   nodes;		/* hosts allocated to the step */
	time_t   start_time;	/* step start time */
	time_t   job_submit_time;/* job submit time needed to find job record
				  * in db */
	uint32_t step_id;	/* step ID */
	uint32_t total_procs;	/* count of allocated processors */
} dbd_step_start_msg_t;

/*****************************************************************************\
 * Slurm DBD message processing functions
\*****************************************************************************/

/* Open a socket connection to SlurmDbd using SlurmdbdAuthInfo specified */
extern int slurm_open_slurmdbd_conn(char *auth_info);

/* Close the SlurmDBD socket connection */
extern int slurm_close_slurmdbd_conn(void);

/* Send an RPC to the SlurmDBD. Do not wait for the reply. The RPC
 * will be queued and processed later if the SlurmDBD is not responding.
 * Returns SLURM_SUCCESS or an error code */
extern int slurm_send_slurmdbd_msg(slurmdbd_msg_t *req);

/* Send an RPC to the SlurmDBD and wait for an arbitrary reply message.
 * The RPC will not be queued if an error occurs.
 * The "resp" message must be freed by the caller.
 * Returns SLURM_SUCCESS or an error code */
extern int slurm_send_recv_slurmdbd_msg(slurmdbd_msg_t *req, 
					slurmdbd_msg_t *resp);

/* Send an RPC to the SlurmDBD and wait for the return code reply.
 * The RPC will not be queued if an error occurs.
 * Returns SLURM_SUCCESS or an error code */
extern int slurm_send_slurmdbd_recv_rc_msg(slurmdbd_msg_t *req, int *rc);

extern Buf pack_slurmdbd_msg(slurmdbd_msg_t *req);
extern int unpack_slurmdbd_msg(slurmdbd_msg_t *resp, Buf buffer);
/*****************************************************************************\
 * Free various SlurmDBD message structures
\*****************************************************************************/
void inline slurmdbd_free_acct_coord_msg(dbd_acct_coord_msg_t *msg);
void inline slurmdbd_free_cluster_procs_msg(dbd_cluster_procs_msg_t *msg);
void inline slurmdbd_free_cond_msg(slurmdbd_msg_type_t type,
				    dbd_cond_msg_t *msg);
void inline slurmdbd_free_get_jobs_msg(dbd_get_jobs_msg_t *msg);
void inline slurmdbd_free_init_msg(dbd_init_msg_t *msg);
void inline slurmdbd_free_job_complete_msg(dbd_job_comp_msg_t *msg);
void inline slurmdbd_free_job_start_msg(dbd_job_start_msg_t *msg);
void inline slurmdbd_free_job_start_rc_msg(dbd_job_start_rc_msg_t *msg);
void inline slurmdbd_free_job_suspend_msg(dbd_job_suspend_msg_t *msg);
void inline slurmdbd_free_list_msg(dbd_list_msg_t *msg);
void inline slurmdbd_free_modify_msg(slurmdbd_msg_type_t type,
				      dbd_modify_msg_t *msg);
void inline slurmdbd_free_node_state_msg(dbd_node_state_msg_t *msg);
void inline slurmdbd_free_rc_msg(dbd_rc_msg_t *msg);
void inline slurmdbd_free_step_complete_msg(dbd_step_comp_msg_t *msg);
void inline slurmdbd_free_step_start_msg(dbd_step_start_msg_t *msg);
void inline slurmdbd_free_usage_msg(slurmdbd_msg_type_t type,
				     dbd_usage_msg_t *msg);

/*****************************************************************************\
 * Pack various SlurmDBD message structures into a buffer
\*****************************************************************************/
void inline slurmdbd_pack_acct_coord_msg(dbd_acct_coord_msg_t *msg,
					  Buf buffer);
void inline slurmdbd_pack_cluster_procs_msg(dbd_cluster_procs_msg_t *msg,
					     Buf buffer);
void inline slurmdbd_pack_cond_msg(slurmdbd_msg_type_t type,
				    dbd_cond_msg_t *msg, Buf buffer);
void inline slurmdbd_pack_get_jobs_msg(dbd_get_jobs_msg_t *msg, Buf buffer);
void inline slurmdbd_pack_init_msg(dbd_init_msg_t *msg, Buf buffer,
				    char *auth_info);
void inline slurmdbd_pack_job_complete_msg(dbd_job_comp_msg_t *msg,
					    Buf buffer);
void inline slurmdbd_pack_job_start_msg(dbd_job_start_msg_t *msg,
					 Buf buffer);
void inline slurmdbd_pack_job_start_rc_msg(dbd_job_start_rc_msg_t *msg,
					    Buf buffer);
void inline slurmdbd_pack_job_suspend_msg(dbd_job_suspend_msg_t *msg,
					   Buf buffer);
void inline slurmdbd_pack_list_msg(slurmdbd_msg_type_t type,
				    dbd_list_msg_t *msg, Buf buffer);
void inline slurmdbd_pack_modify_msg(slurmdbd_msg_type_t type,
				      dbd_modify_msg_t *msg, Buf buffer);
void inline slurmdbd_pack_node_state_msg(dbd_node_state_msg_t *msg,
					  Buf buffer);
void inline slurmdbd_pack_rc_msg(dbd_rc_msg_t *msg, Buf buffer);
void inline slurmdbd_pack_step_complete_msg(dbd_step_comp_msg_t *msg,
					     Buf buffer);
void inline slurmdbd_pack_step_start_msg(dbd_step_start_msg_t *msg,
					  Buf buffer);
void inline slurmdbd_pack_usage_msg(slurmdbd_msg_type_t type,
				     dbd_usage_msg_t *msg, Buf buffer);

/*****************************************************************************\
 * Unpack various SlurmDBD message structures from a buffer
\*****************************************************************************/
int inline slurmdbd_unpack_acct_coord_msg(dbd_acct_coord_msg_t **msg,
					   Buf buffer);
int inline slurmdbd_unpack_cluster_procs_msg(dbd_cluster_procs_msg_t **msg,
					      Buf buffer);
int inline slurmdbd_unpack_cond_msg(slurmdbd_msg_type_t type,
				     dbd_cond_msg_t **msg, Buf buffer);
int inline slurmdbd_unpack_get_jobs_msg(dbd_get_jobs_msg_t **msg, Buf buffer);
int inline slurmdbd_unpack_init_msg(dbd_init_msg_t **msg, Buf buffer,
				     char *auth_info);
int inline slurmdbd_unpack_job_complete_msg(dbd_job_comp_msg_t **msg,
					     Buf buffer);
int inline slurmdbd_unpack_job_start_msg(dbd_job_start_msg_t **msg,
					  Buf buffer);
int inline slurmdbd_unpack_job_start_rc_msg(dbd_job_start_rc_msg_t **msg,
					     Buf buffer);
int inline slurmdbd_unpack_job_suspend_msg(dbd_job_suspend_msg_t **msg,
					    Buf buffer);
int inline slurmdbd_unpack_list_msg(slurmdbd_msg_type_t type,
				     dbd_list_msg_t **msg, Buf buffer);
int inline slurmdbd_unpack_modify_msg(slurmdbd_msg_type_t type,
				       dbd_modify_msg_t **msg, Buf buffer);
int inline slurmdbd_unpack_node_state_msg(dbd_node_state_msg_t **msg,
					   Buf buffer);
int inline slurmdbd_unpack_rc_msg(dbd_rc_msg_t **msg, Buf buffer);
int inline slurmdbd_unpack_step_complete_msg(dbd_step_comp_msg_t **msg,
					      Buf buffer);
int inline slurmdbd_unpack_step_start_msg(dbd_step_start_msg_t **msg,
					   Buf buffer);
int inline slurmdbd_unpack_usage_msg(slurmdbd_msg_type_t type,
				      dbd_usage_msg_t **msg,
				      Buf buffer);

#endif	/* !_SLURMDBD_DEFS_H */
