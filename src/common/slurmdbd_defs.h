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

#include "src/common/pack.h"

/* Increment SLURM_DBD_VERSION if any of the RPCs change */
#define SLURM_DBD_VERSION 01

/* SLURM DBD message types */
typedef enum {
	DBD_INIT = 1400,	/* Connection initialization		*/
	DBD_CLUSTER_PROCS,	/* Record tota processors on cluster	*/
	DBD_GET_JOBS,		/* Get job information			*/
	DBD_GOT_JOBS,		/* Response to DBD_GET_JOBS		*/
	DBD_JOB_COMPLETE,	/* Record job completion 		*/
	DBD_JOB_START,		/* Record job starting			*/
	DBD_JOB_SUSPEND,	/* Record job suspension		*/
	DBD_RC,			/* Return code from operation		*/
	DBD_NODE_STATE,		/* Record node state transition		*/
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

typedef struct dbd_cluster_procs_msg {
	char *cluster_name;	/* name of cluster */
	uint32_t proc_count;	/* total processor count */
	time_t event_time;	/* time of transition */
} dbd_cluster_procs_msg_t;

typedef struct dbd_get_jobs_msg {
	uint32_t job_count;	/* count of job ID filters */
	uint32_t *job_ids;	/* array of job ID filters */
} dbd_get_jobs_msg_t;

typedef struct dbd_got_jobs_msg {
	uint32_t job_count;	/* count of job IDs */
	uint32_t *job_ids;	/* array of job IDs */
} dbd_got_jobs_msg_t;

typedef struct dbd_init_msg {
	uint16_t version;	/* protocol version */
	uint32_t uid;		/* UID originating connection,
				 * filled by authtentication plugin*/
} dbd_init_msg_t;

typedef struct dbd_job_comp_msg {
	char *   account;
	time_t   end_time;
	uint32_t exit_code;
	uint32_t job_id;
	uint16_t job_state;
	char *   name;
	char *   nodes;
	uint32_t priority;
	time_t   start_time;
	uint32_t total_procs;
} dbd_job_comp_msg_t;

typedef struct dbd_job_start_msg {
	char *   account;
	uint32_t job_id;
	uint16_t job_state;
	char *   name;
	char *   nodes;
	uint32_t priority;
	time_t   start_time;
	uint32_t total_procs;
} dbd_job_start_msg_t;

typedef struct dbd_job_suspend_msg {
	uint32_t job_id;
} dbd_job_suspend_msg_t;

typedef struct dbd_rc_msg {
	uint32_t return_code;
} dbd_rc_msg_t;

#define DBD_NODE_STATE_DOWN  1
#define DBD_NODE_STATE_UP    2
typedef struct dbd_node_state_msg {
	char *hostlist;		/* name of hosts */
	char *reason;		/* explanation for the node's state */
	uint16_t new_state;	/* new state of host, see DBD_NODE_STATE_* */
	time_t event_time;	/* time of transition */
} dbd_node_state_msg_t;

typedef struct dbd_step_comp_msg {
	uint32_t job_id;
	uint32_t step_id;
} dbd_step_comp_msg_t;

typedef struct dbd_step_start_msg {
	uint32_t job_id;
	uint32_t step_id;
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

/*****************************************************************************\
 * Free various SlurmDBD message structures
\*****************************************************************************/
void inline slurm_dbd_free_cluster_procs_msg(dbd_cluster_procs_msg_t *msg);
void inline slurm_dbd_free_get_jobs_msg(dbd_get_jobs_msg_t *msg);
void inline slurm_dbd_free_got_jobs_msg(dbd_got_jobs_msg_t *msg);
void inline slurm_dbd_free_init_msg(dbd_init_msg_t *msg);
void inline slurm_dbd_free_job_complete_msg(dbd_job_comp_msg_t *msg);
void inline slurm_dbd_free_job_start_msg(dbd_job_start_msg_t *msg);
void inline slurm_dbd_free_job_suspend_msg(dbd_job_suspend_msg_t *msg);
void inline slurm_dbd_free_rc_msg(dbd_rc_msg_t *msg);
void inline slurm_dbd_free_node_state_msg(dbd_node_state_msg_t *msg);
void inline slurm_dbd_free_step_complete_msg(dbd_step_comp_msg_t *msg);
void inline slurm_dbd_free_step_start_msg(dbd_step_start_msg_t *msg);

/*****************************************************************************\
 * Pack various SlurmDBD message structures into a buffer
\*****************************************************************************/
void inline slurm_dbd_pack_cluster_procs_msg(dbd_cluster_procs_msg_t *msg,
								       Buf buffer);
void inline slurm_dbd_pack_get_jobs_msg(dbd_get_jobs_msg_t *msg,       Buf buffer);
void inline slurm_dbd_pack_got_jobs_msg(dbd_got_jobs_msg_t *msg,       Buf buffer);
void inline slurm_dbd_pack_init_msg(dbd_init_msg_t *msg,               Buf buffer,
				    char *auth_info);
void inline slurm_dbd_pack_job_complete_msg(dbd_job_comp_msg_t *msg,   Buf buffer);
void inline slurm_dbd_pack_job_start_msg(dbd_job_start_msg_t *msg,     Buf buffer);
void inline slurm_dbd_pack_job_suspend_msg(dbd_job_suspend_msg_t *msg, Buf buffer);
void inline slurm_dbd_pack_rc_msg(dbd_rc_msg_t *msg,                   Buf buffer);
void inline slurm_dbd_pack_node_state_msg(dbd_node_state_msg_t *msg,   Buf buffer);
void inline slurm_dbd_pack_step_complete_msg(dbd_step_comp_msg_t *msg, Buf buffer);
void inline slurm_dbd_pack_step_start_msg(dbd_step_start_msg_t *msg,   Buf buffer);

/*****************************************************************************\
 * Unpack various SlurmDBD message structures from a buffer
\*****************************************************************************/
int inline slurm_dbd_unpack_cluster_procs_msg(dbd_cluster_procs_msg_t **msg,
									 Buf buffer);
int inline slurm_dbd_unpack_get_jobs_msg(dbd_get_jobs_msg_t **msg,       Buf buffer);
int inline slurm_dbd_unpack_got_jobs_msg(dbd_got_jobs_msg_t **msg,       Buf buffer);
int inline slurm_dbd_unpack_init_msg(dbd_init_msg_t **msg,               Buf buffer,
				     char *auth_info);
int inline slurm_dbd_unpack_job_complete_msg(dbd_job_comp_msg_t **msg,   Buf buffer);
int inline slurm_dbd_unpack_job_start_msg(dbd_job_start_msg_t **msg,     Buf buffer);
int inline slurm_dbd_unpack_job_suspend_msg(dbd_job_suspend_msg_t **msg, Buf buffer);
int inline slurm_dbd_unpack_rc_msg(dbd_rc_msg_t **msg,                   Buf buffer);
int inline slurm_dbd_unpack_node_state_msg(dbd_node_state_msg_t **msg,   Buf buffer);
int inline slurm_dbd_unpack_step_complete_msg(dbd_step_comp_msg_t **msg, Buf buffer);
int inline slurm_dbd_unpack_step_start_msg(dbd_step_start_msg_t **msg,   Buf buffer);

#endif	/* !_SLURMDBD_DEFS_H */
