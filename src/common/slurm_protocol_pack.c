/****************************************************************************\
 *  slurm_protocol_pack.c - functions to pack and unpack structures for RPCs
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Kevin Tew <tew1@llnl.gov>, et. al.
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "src/api/slurm_pmi.h"
#include "src/common/bitstring.h"
#include "src/common/log.h"
#include "src/common/node_select.h"
#include "src/common/pack.h"
#include "src/common/slurm_auth.h"
#include "src/common/slurm_cred.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_protocol_defs.h"
#include "src/common/slurm_protocol_pack.h"
#include "src/common/switch.h"
#include "src/common/xmalloc.h"
#include "src/common/xassert.h"
#include "src/common/forward.h"

#define _pack_job_info_msg(msg,buf)		_pack_buffer_msg(msg,buf)
#define _pack_job_step_info_msg(msg,buf)	_pack_buffer_msg(msg,buf)
#define _pack_node_select_info_msg(msg,buf)	_pack_buffer_msg(msg,buf)
#define _pack_node_info_msg(msg,buf)		_pack_buffer_msg(msg,buf)

static void _pack_update_node_msg(update_node_msg_t * msg, Buf buffer);
static int _unpack_update_node_msg(update_node_msg_t ** msg, Buf buffer);

static void
 _pack_node_registration_status_msg(slurm_node_registration_status_msg_t *
				    msg, Buf buffer);
static int
 _unpack_node_registration_status_msg(slurm_node_registration_status_msg_t
				      ** msg, Buf buffer);

static void _pack_job_ready_msg(job_id_msg_t * msg, Buf buffer);
static int _unpack_job_ready_msg(job_id_msg_t ** msg_ptr, Buf buffer);

static void
 _pack_resource_allocation_response_msg(resource_allocation_response_msg_t *
					msg, Buf buffer);
static int
 _unpack_resource_allocation_response_msg(resource_allocation_response_msg_t
					  ** msg, Buf buffer);

static void
 _pack_resource_allocation_and_run_response_msg
    (resource_allocation_and_run_response_msg_t * msg, Buf buffer);
static int
 _unpack_resource_allocation_and_run_response_msg
    (resource_allocation_and_run_response_msg_t ** msg, Buf buffer);

static void _pack_submit_response_msg(submit_response_msg_t * msg,
				      Buf buffer);
static int _unpack_submit_response_msg(submit_response_msg_t ** msg,
				       Buf buffer);

static void _pack_node_info_request_msg(
			node_info_request_msg_t * msg, Buf buffer);

static int _unpack_node_info_request_msg(
			node_info_request_msg_t ** msg, Buf bufer);

static int _unpack_node_info_msg(node_info_msg_t ** msg, Buf buffer);
static int _unpack_node_info_members(node_info_t * node, Buf buffer);
static int _unpack_node_select_info_msg(node_select_info_msg_t ** msg, 
			Buf buffer);

static void _pack_update_partition_msg(update_part_msg_t * msg, Buf buffer);
static int _unpack_update_partition_msg(update_part_msg_t ** msg, Buf buffer);

static void _pack_delete_partition_msg(delete_part_msg_t * msg, Buf buffer);
static int _unpack_delete_partition_msg(delete_part_msg_t ** msg, Buf buffer);

static void _pack_job_step_create_request_msg(job_step_create_request_msg_t
					      * msg, Buf buffer);
static int _unpack_job_step_create_request_msg(
			job_step_create_request_msg_t ** msg, Buf buffer);

static void _pack_kill_job_msg(kill_job_msg_t * msg, Buf buffer);
static int _unpack_kill_job_msg(kill_job_msg_t ** msg, Buf buffer);

static void _pack_signal_job_msg(signal_job_msg_t * msg, Buf buffer);
static int _unpack_signal_job_msg(signal_job_msg_t ** msg, Buf buffer);

static void _pack_epilog_comp_msg(epilog_complete_msg_t * msg, Buf buffer);
static int  _unpack_epilog_comp_msg(epilog_complete_msg_t ** msg, Buf buffer);

static void _pack_update_job_time_msg(job_time_msg_t * msg, Buf buffer);
static int _unpack_update_job_time_msg(job_time_msg_t ** msg, Buf buffer);

static void  _pack_job_step_create_response_msg(
			job_step_create_response_msg_t * msg, Buf buffer);
static int _unpack_job_step_create_response_msg(
			job_step_create_response_msg_t ** msg, Buf buffer);

static void _pack_part_info_request_msg(part_info_request_msg_t * msg, 
			Buf buffer);

static int _unpack_part_info_request_msg(part_info_request_msg_t ** 
			msg, Buf buffer);

static void _pack_partition_info_msg(slurm_msg_t * msg, Buf buffer);
static int _unpack_partition_info_msg(partition_info_msg_t ** msg,
				      Buf buffer);
static int _unpack_partition_info_members(partition_info_t * part,
					  Buf buffer);

static void _pack_launch_tasks_request_msg(launch_tasks_request_msg_t *
					   msg, Buf buffer);
static int _unpack_launch_tasks_request_msg(launch_tasks_request_msg_t **
					    msg_ptr, Buf buffer);


static void _pack_spawn_task_request_msg(spawn_task_request_msg_t *
					   msg, Buf buffer);
static int _unpack_spawn_task_request_msg(spawn_task_request_msg_t **
					    msg_ptr, Buf buffer);

static void _pack_cancel_tasks_msg(kill_tasks_msg_t * msg, Buf buffer);
static int _unpack_cancel_tasks_msg(kill_tasks_msg_t ** msg_ptr, Buf buffer);

static void _pack_launch_tasks_response_msg(launch_tasks_response_msg_t *
					    msg, Buf buffer);
static int _unpack_launch_tasks_response_msg(launch_tasks_response_msg_t **
					     msg_ptr, Buf buffer);

static void _pack_shutdown_msg(shutdown_msg_t * msg, Buf buffer);
static int _unpack_shutdown_msg(shutdown_msg_t ** msg_ptr, Buf buffer);

static void _pack_reattach_tasks_request_msg(reattach_tasks_request_msg_t *,
					     Buf);
static int _unpack_reattach_tasks_request_msg(reattach_tasks_request_msg_t **,
					      Buf);

static void
 _pack_reattach_tasks_response_msg(reattach_tasks_response_msg_t *, Buf);
static int
 _unpack_reattach_tasks_response_msg(reattach_tasks_response_msg_t **, Buf);

static void _pack_task_exit_msg(task_exit_msg_t * msg, Buf buffer);
static int _unpack_task_exit_msg(task_exit_msg_t ** msg_ptr, Buf buffer);

static void _pack_old_job_desc_msg(old_job_alloc_msg_t * job_desc_ptr,
				   Buf buffer);
static int _unpack_old_job_desc_msg(old_job_alloc_msg_t **
				    job_desc_buffer_ptr, Buf buffer);

static void _pack_return_code_msg(return_code_msg_t * msg, Buf buffer);
static int _unpack_return_code_msg(return_code_msg_t ** msg, Buf buffer);

static void _pack_slurm_ctl_conf_msg(slurm_ctl_conf_info_msg_t * build_ptr,
				     Buf buffer);
static int _unpack_slurm_ctl_conf_msg(slurm_ctl_conf_info_msg_t **
				      build_buffer_ptr, Buf buffer);

static void _pack_job_info_request_msg(job_info_request_msg_t * 
		msg, Buf buffer);
static int _unpack_job_info_request_msg(job_info_request_msg_t** 
		msg, Buf buffer);

static void _pack_node_select_info_req_msg(node_info_select_request_msg_t *
		msg, Buf buffer);
static int _unpack_node_select_info_req_msg(node_info_select_request_msg_t **
		msg, Buf buffer);

static void _pack_job_step_info_req_msg(job_step_info_request_msg_t * msg,
					Buf buffer);
static int _unpack_job_step_info_req_msg(job_step_info_request_msg_t **
					 msg, Buf buffer);
static int _unpack_job_step_info_response_msg(job_step_info_response_msg_t
					      ** msg, Buf buffer);
static int _unpack_job_step_info_members(job_step_info_t * step, Buf buffer);

static void _pack_complete_job_allocation_msg(
	complete_job_allocation_msg_t * msg, Buf buffer);
static int _unpack_complete_job_allocation_msg(
	complete_job_allocation_msg_t ** msg_ptr, Buf buffer);
static void _pack_complete_batch_script_msg(
	complete_batch_script_msg_t * msg, Buf buffer);
static int _unpack_complete_batch_script_msg(
	complete_batch_script_msg_t ** msg_ptr, Buf buffer);
static void _pack_step_complete_msg(step_complete_msg_t * msg,
				    Buf buffer);
static int _unpack_step_complete_msg(step_complete_msg_t **
				     msg_ptr, Buf buffer);
static int _unpack_job_info_members(job_info_t * job, Buf buffer);

static void _pack_batch_job_launch_msg(batch_job_launch_msg_t * msg,
				       Buf buffer);
static int _unpack_batch_job_launch_msg(batch_job_launch_msg_t ** msg,
					Buf buffer);

static void _pack_job_desc_msg(job_desc_msg_t * job_desc_ptr, Buf buffer);
static int _unpack_job_desc_msg(job_desc_msg_t ** job_desc_buffer_ptr,
				Buf buffer);
static int _unpack_job_info_msg(job_info_msg_t ** msg, Buf buffer);

static void _pack_last_update_msg(last_update_msg_t * msg, Buf buffer);
static int _unpack_last_update_msg(last_update_msg_t ** msg, Buf buffer);

static void _pack_slurm_addr_array(slurm_addr * slurm_address,
				   uint16_t size_val, Buf buffer);
static int _unpack_slurm_addr_array(slurm_addr ** slurm_address,
				    uint16_t * size_val, Buf buffer);

static void _pack_ret_list(List ret_list, uint16_t size_val, Buf buffer);
static int _unpack_ret_list(List *ret_list, uint16_t size_val, Buf buffer);

static void _pack_job_id_request_msg(job_id_request_msg_t * msg, Buf buffer);
static int  _unpack_job_id_request_msg(job_id_request_msg_t ** msg, Buf buffer);

static void _pack_job_id_response_msg(job_id_response_msg_t * msg, Buf buffer);
static int  _unpack_job_id_response_msg(job_id_response_msg_t ** msg, 
					Buf buffer);

static void _pack_job_step_kill_msg(job_step_kill_msg_t * msg, Buf buffer);
static int  _unpack_job_step_kill_msg(job_step_kill_msg_t ** msg_ptr, 
				      Buf buffer);

static void _pack_srun_ping_msg(srun_ping_msg_t * msg, Buf buffer);
static int  _unpack_srun_ping_msg(srun_ping_msg_t ** msg_ptr, Buf buffer);

static void _pack_srun_node_fail_msg(srun_node_fail_msg_t * msg, Buf buffer);
static int  _unpack_srun_node_fail_msg(srun_node_fail_msg_t ** msg_ptr, 
					Buf buffer);

static void _pack_srun_timeout_msg(srun_timeout_msg_t * msg, Buf buffer);
static int  _unpack_srun_timeout_msg(srun_timeout_msg_t ** msg_ptr, Buf buffer);

static void _pack_checkpoint_msg(checkpoint_msg_t *msg, Buf buffer);
static int  _unpack_checkpoint_msg(checkpoint_msg_t **msg_ptr, Buf buffer);

static void _pack_checkpoint_resp_msg(checkpoint_resp_msg_t *msg, Buf buffer);
static int  _unpack_checkpoint_resp_msg(checkpoint_resp_msg_t **msg_ptr, 
		Buf buffer);

static void _pack_checkpoint_comp(checkpoint_comp_msg_t *msg, Buf buffer);
static int  _unpack_checkpoint_comp(checkpoint_comp_msg_t **msg_ptr, 
		Buf buffer);

static void _pack_suspend_msg(suspend_msg_t *msg, Buf buffer);
static int  _unpack_suspend_msg(suspend_msg_t **msg_ptr, Buf buffer);

static void _pack_buffer_msg(slurm_msg_t * msg, Buf buffer);

static void _pack_jobacct_data(jobacct_msg_t * msg , Buf buffer );

static int _unpack_jobacct_data(jobacct_msg_t ** msg_ptr , Buf buffer );

static void _pack_kvs_rec(struct kvs_comm *msg_ptr, Buf buffer);
static int  _unpack_kvs_rec(struct kvs_comm **msg_ptr, Buf buffer);
static void _pack_kvs_data(struct kvs_comm_set *msg_ptr, Buf buffer);
static int  _unpack_kvs_data(struct kvs_comm_set **msg_ptr, Buf buffer);
static void _pack_kvs_get(kvs_get_msg_t *msg_ptr, Buf buffer);
static int  _unpack_kvs_get(kvs_get_msg_t **msg_ptr, Buf buffer);

static void _pack_file_bcast(file_bcast_msg_t * msg , Buf buffer );

static int _unpack_file_bcast(file_bcast_msg_t ** msg_ptr , Buf buffer );

/* pack_header
 * packs a slurm protocol header that proceeds every slurm message
 * IN header - the header structure to pack
 * IN/OUT buffer - destination of the pack, contains pointers that are
 *			automatically updated
 */
void
pack_header(header_t * header, Buf buffer)
{
	
	pack16((uint16_t)header->version, buffer);
	pack16((uint16_t)header->flags, buffer);
	pack16((uint16_t) header->msg_type, buffer);
	pack32((uint32_t)header->body_length, buffer);
	pack16((uint16_t)header->forward.cnt, buffer);
	if (header->forward.cnt > 0) {
		_pack_slurm_addr_array(header->forward.addr,
				       header->forward.cnt, buffer);
		packmem(header->forward.name, 
			(header->forward.cnt * MAX_SLURM_NAME), 
			buffer);
		pack32_array(header->forward.node_id, 
			     header->forward.cnt, 
			     buffer);
		pack32((uint32_t)header->forward.timeout, buffer);
	}
	pack16((uint16_t)header->ret_cnt, buffer);	
	if(header->ret_cnt > 0) {
		_pack_ret_list(header->ret_list,
			       header->ret_cnt, buffer);		
	}
	pack32((uint32_t)header->srun_node_id, buffer);
	slurm_pack_slurm_addr(&header->orig_addr, buffer);
}

/* unpack_header
 * unpacks a slurm protocol header that proceeds every slurm message
 * OUT header - the header structure to unpack
 * IN/OUT buffer - source of the unpack data, contains pointers that are
 *			automatically updated
 * RET 0 or error code
 */
int
unpack_header(header_t * header, Buf buffer)
{
	uint16_t uint16_tmp = 0;
	uint32_t uint32_tmp = 0;

	forward_init(&header->forward, NULL);
	header->ret_list = NULL;

	safe_unpack16(&header->version, buffer);
	safe_unpack16(&header->flags, buffer);
	safe_unpack16(&uint16_tmp, buffer);
	header->msg_type = (slurm_msg_type_t) uint16_tmp;
	safe_unpack32(&header->body_length, buffer);
	safe_unpack16(&header->forward.cnt, buffer);
	if (header->forward.cnt > 0) {
		if(_unpack_slurm_addr_array(&(header->forward.addr),
					    &uint16_tmp, buffer))
			goto unpack_error;
		if(uint16_tmp != header->forward.cnt)
			goto unpack_error;
		safe_unpackmem_xmalloc(&(header->forward.name), 
				       &uint16_tmp, 
				       buffer);
		
		safe_unpack32_array(&(header->forward.node_id), 
				    &uint32_tmp,
				    buffer);
		if(uint32_tmp != (uint32_t)header->forward.cnt) {
			goto unpack_error;
		}
		safe_unpack32(&header->forward.timeout, buffer);
	} 
	
	safe_unpack16(&header->ret_cnt, buffer);	
	if(header->ret_cnt > 0) {
		if(_unpack_ret_list(&(header->ret_list),
				    header->ret_cnt, buffer))
			goto unpack_error;
	} else {
		header->ret_list = NULL;
	}
	safe_unpack32(&header->srun_node_id, buffer);
	slurm_unpack_slurm_addr_no_alloc(&header->orig_addr, buffer);
	
	return SLURM_SUCCESS;

      unpack_error:
	error("unpacking header");
	destroy_forward(&header->forward);
	if(header->ret_list)
		list_destroy(header->ret_list);
	return SLURM_ERROR;
}


/* pack_msg
 * packs a generic slurm protocol message body
 * IN msg - the body structure to pack (note: includes message type)
 * IN/OUT buffer - destination of the pack, contains pointers that are 
 *			automatically updated
 * RET 0 or error code
 */
int
pack_msg(slurm_msg_t const *msg, Buf buffer)
{
	switch (msg->msg_type) {
	 case REQUEST_NODE_INFO:
		_pack_node_info_request_msg((node_info_request_msg_t *)
				       msg->data, buffer);
		break;
	 case REQUEST_PARTITION_INFO:
		_pack_part_info_request_msg((part_info_request_msg_t *)
				       msg->data, buffer);
		break;
	 case REQUEST_BUILD_INFO:
	 case REQUEST_ACCTING_INFO:
		 _pack_last_update_msg((last_update_msg_t *)
				       msg->data, buffer);
		 break;
	 case RESPONSE_BUILD_INFO:
		 _pack_slurm_ctl_conf_msg((slurm_ctl_conf_info_msg_t *)
					  msg->data, buffer);
		 break;
	 case RESPONSE_JOB_INFO:
		 _pack_job_info_msg((slurm_msg_t *) msg, buffer);
		 break;
	 case RESPONSE_PARTITION_INFO:
		 _pack_partition_info_msg((slurm_msg_t *) msg, buffer);
		 break;
	 case RESPONSE_NODE_INFO:
		 _pack_node_info_msg((slurm_msg_t *) msg, buffer);
		 break;
	 case MESSAGE_NODE_REGISTRATION_STATUS:
		 _pack_node_registration_status_msg(
			(slurm_node_registration_status_msg_t *) msg->data, 
			buffer);
		 break;
	 case REQUEST_RESOURCE_ALLOCATION:
	 case REQUEST_SUBMIT_BATCH_JOB:
	 case REQUEST_JOB_WILL_RUN:
	 case REQUEST_ALLOCATION_AND_RUN_JOB_STEP:
	 case REQUEST_UPDATE_JOB:
		 _pack_job_desc_msg((job_desc_msg_t *)
				    msg->data, buffer);
		 break;
	 case REQUEST_OLD_JOB_RESOURCE_ALLOCATION:
		 _pack_old_job_desc_msg((old_job_alloc_msg_t *) msg->data,
					buffer);
		 break;
	 case REQUEST_NODE_REGISTRATION_STATUS:
	 case REQUEST_RECONFIGURE:
	 case REQUEST_SHUTDOWN_IMMEDIATE:
	 case REQUEST_PING:
	 case REQUEST_CONTROL:
		 /* Message contains no body/information */
		 break;
	 case REQUEST_SHUTDOWN:
		 _pack_shutdown_msg((shutdown_msg_t *) msg->data, buffer);
		 break;
	 case RESPONSE_SUBMIT_BATCH_JOB:
		 _pack_submit_response_msg((submit_response_msg_t *) 
					   msg->data, buffer);
		 break;
	 case RESPONSE_RESOURCE_ALLOCATION:
	 case RESPONSE_JOB_WILL_RUN:
		 _pack_resource_allocation_response_msg
		     ((resource_allocation_response_msg_t *) msg->data,
		      buffer);
		 break;
	 case RESPONSE_ALLOCATION_AND_RUN_JOB_STEP:
		 _pack_resource_allocation_and_run_response_msg(
			(resource_allocation_and_run_response_msg_t *) 
			msg->data, buffer);
		 break;
	 case REQUEST_UPDATE_NODE:
		 _pack_update_node_msg((update_node_msg_t *) msg->data,
				       buffer);
		 break;
	 case REQUEST_UPDATE_PARTITION:
		 _pack_update_partition_msg((update_part_msg_t *) msg->
					    data, buffer);
		 break;
	 case REQUEST_DELETE_PARTITION:
		 _pack_delete_partition_msg((delete_part_msg_t *) msg->
					    data, buffer);
		 break;
	 case REQUEST_REATTACH_TASKS:
		 _pack_reattach_tasks_request_msg(
			(reattach_tasks_request_msg_t *) msg->data, buffer);
		 break;
	 case RESPONSE_REATTACH_TASKS:
		 _pack_reattach_tasks_response_msg(
			(reattach_tasks_response_msg_t *) msg->data, buffer);
		 break;
	 case REQUEST_LAUNCH_TASKS:
		 _pack_launch_tasks_request_msg(
				(launch_tasks_request_msg_t *) msg->data, 
				buffer);
		 break;
	 case RESPONSE_LAUNCH_TASKS:
		 _pack_launch_tasks_response_msg((launch_tasks_response_msg_t
						  *) msg->data, buffer);
		 break;
	 case REQUEST_SPAWN_TASK:
		 _pack_spawn_task_request_msg(
					(spawn_task_request_msg_t *)
					msg->data, buffer);
		 break;
	 case REQUEST_SIGNAL_TASKS:
	 case REQUEST_TERMINATE_TASKS:
		 _pack_cancel_tasks_msg((kill_tasks_msg_t *) msg->data,
					buffer);
		 break;
	 case REQUEST_JOB_STEP_INFO:
		 _pack_job_step_info_req_msg((job_step_info_request_msg_t
					      *) msg->data, buffer);
		 break;
	 case REQUEST_JOB_INFO:
		_pack_job_info_request_msg((job_info_request_msg_t *) 
				msg->data, buffer);
		break;
	 case REQUEST_CANCEL_JOB_STEP:
		 _pack_job_step_kill_msg((job_step_kill_msg_t *) 
					 msg->data, buffer);
		 break;
	 case REQUEST_COMPLETE_JOB_ALLOCATION:
		 _pack_complete_job_allocation_msg(
			 (complete_job_allocation_msg_t *)msg->data, buffer);
		 break;
	 case REQUEST_COMPLETE_BATCH_SCRIPT:
		 _pack_complete_batch_script_msg(
			 (complete_batch_script_msg_t *)msg->data, buffer);
		 break;
	 case REQUEST_STEP_COMPLETE:
		 _pack_step_complete_msg((step_complete_msg_t *)msg->data,
					 buffer);
		 break;
	 case REQUEST_SIGNAL_JOB:
		 _pack_signal_job_msg((signal_job_msg_t *) msg->data, buffer);
		 break;
	 case REQUEST_KILL_TIMELIMIT:
	 case REQUEST_TERMINATE_JOB:
		 _pack_kill_job_msg((kill_job_msg_t *) msg->data, buffer);
		 break;
	 case MESSAGE_EPILOG_COMPLETE:
		_pack_epilog_comp_msg((epilog_complete_msg_t *) msg->data, buffer);
		break;
	 case REQUEST_UPDATE_JOB_TIME:
		 _pack_update_job_time_msg((job_time_msg_t *)
					     msg->data, buffer);
		 break;
	 case RESPONSE_RECONFIGURE:
	 case RESPONSE_SHUTDOWN:
	 case RESPONSE_CANCEL_JOB_STEP:
		 break;
	 case REQUEST_JOB_ATTACH:
		 break;
	 case RESPONSE_JOB_ATTACH:
		 break;
	 case RESPONSE_JOB_STEP_INFO:
		 _pack_job_step_info_msg((slurm_msg_t *) msg, buffer);
		 break;
	 case REQUEST_JOB_RESOURCE:
		 break;
	 case RESPONSE_JOB_RESOURCE:
		 break;
	 case REQUEST_RUN_JOB_STEP:
		 break;
	 case RESPONSE_RUN_JOB_STEP:
		 break;
	 case MESSAGE_TASK_EXIT:
		 _pack_task_exit_msg((task_exit_msg_t *) msg->data, buffer);
		 break;
	 case REQUEST_BATCH_JOB_LAUNCH:
		 _pack_batch_job_launch_msg((batch_job_launch_msg_t *)
					    msg->data, buffer);
		 break;
	 case MESSAGE_UPLOAD_ACCOUNTING_INFO:
		 break;
	 case RESPONSE_JOB_READY:
	 case RESPONSE_SLURM_RC:
		 _pack_return_code_msg((return_code_msg_t *) msg->data,
				       buffer);
		 break;
	 case RESPONSE_JOB_STEP_CREATE:
		 _pack_job_step_create_response_msg(
				(job_step_create_response_msg_t *)
				msg->data, buffer);
		 break;
	 case REQUEST_JOB_STEP_CREATE:
		 _pack_job_step_create_request_msg(
				(job_step_create_request_msg_t *)
				msg->data, buffer);
		 break;
	 case REQUEST_JOB_ID:
		 _pack_job_id_request_msg(
				(job_id_request_msg_t *)msg->data,
				buffer);
		 break;
	 case RESPONSE_JOB_ID:
		 _pack_job_id_response_msg(
				(job_id_response_msg_t *)msg->data,
				buffer);
		 break;
	 case SRUN_PING:
		_pack_srun_ping_msg((srun_ping_msg_t *)msg->data, buffer);
		break;
	 case SRUN_NODE_FAIL:
		_pack_srun_node_fail_msg((srun_node_fail_msg_t *)msg->data, 
				buffer);
		break;
	 case SRUN_TIMEOUT:
		_pack_srun_timeout_msg((srun_timeout_msg_t *)msg->data, buffer);
		break;
	 case REQUEST_CHECKPOINT:
		_pack_checkpoint_msg((checkpoint_msg_t *)msg->data, buffer);
		break;
	 case REQUEST_CHECKPOINT_COMP:
		_pack_checkpoint_comp((checkpoint_comp_msg_t *)msg->data, 
			buffer);
		break;
	 case RESPONSE_CHECKPOINT:
	 case RESPONSE_CHECKPOINT_COMP:
		_pack_checkpoint_resp_msg((checkpoint_resp_msg_t *)msg->data, 
			buffer);
		break;
	 case REQUEST_SUSPEND:
		_pack_suspend_msg((suspend_msg_t *)msg->data, buffer);
		break;

	 case REQUEST_JOB_READY:
		_pack_job_ready_msg((job_id_msg_t *)msg->data, buffer);
		break;
	 case REQUEST_NODE_SELECT_INFO:
		_pack_node_select_info_req_msg(
			(node_info_select_request_msg_t *) msg->data, buffer);
		break;
	 case RESPONSE_NODE_SELECT_INFO:
		_pack_node_select_info_msg((slurm_msg_t *) msg, buffer);
		break;
	 case MESSAGE_JOBACCT_DATA:
		_pack_jobacct_data((jobacct_msg_t *) msg->data, buffer);
		break;
	 case REQUEST_FILE_BCAST:
		_pack_file_bcast((file_bcast_msg_t *) msg->data, buffer);
		break;
	 case PMI_KVS_PUT_REQ:
	 case PMI_KVS_GET_RESP:
		_pack_kvs_data((struct kvs_comm_set *) msg->data, buffer);
		break;
	 case PMI_KVS_GET_REQ:
		_pack_kvs_get((kvs_get_msg_t *) msg->data, buffer);
		break;
	 case PMI_KVS_PUT_RESP:
		break;	/* no data in message */
	 default:
		 debug("No pack method for msg type %i", msg->msg_type);
		 return EINVAL;
		 break;

	}
	return SLURM_SUCCESS;
}

/* unpack_msg
 * unpacks a generic slurm protocol message body
 * OUT msg - the body structure to unpack (note: includes message type)
 * IN/OUT buffer - source of the unpack, contains pointers that are 
 *			automatically updated
 * RET 0 or error code
 */
int
unpack_msg(slurm_msg_t * msg, Buf buffer)
{
	int rc = SLURM_SUCCESS;
	msg->data = NULL;	/* Initialize to no data for now */
	
	switch (msg->msg_type) {
	 case REQUEST_NODE_INFO:
		rc = _unpack_node_info_request_msg((node_info_request_msg_t **)
				       & (msg->data), buffer);
		break;
	 case REQUEST_PARTITION_INFO:
		rc = _unpack_part_info_request_msg((part_info_request_msg_t **)
				       & (msg->data), buffer);
		break;
	 case REQUEST_BUILD_INFO:
	 case REQUEST_ACCTING_INFO:
		 rc = _unpack_last_update_msg((last_update_msg_t **) &
					      (msg->data), buffer);
		 break;
	 case RESPONSE_BUILD_INFO:
		 rc = _unpack_slurm_ctl_conf_msg((slurm_ctl_conf_info_msg_t
						  **)
						 & (msg->data), buffer);
		 break;
	 case RESPONSE_JOB_INFO:
		 rc = _unpack_job_info_msg((job_info_msg_t **) & (msg->data),
					   buffer);
		 break;
	 case RESPONSE_PARTITION_INFO:
		 rc = _unpack_partition_info_msg((partition_info_msg_t **) &
						 (msg->data), buffer);
		 break;
	 case RESPONSE_NODE_INFO:
		 rc = _unpack_node_info_msg((node_info_msg_t **) &
					    (msg->data), buffer);
		 break;
	 case MESSAGE_NODE_REGISTRATION_STATUS:
		 rc = _unpack_node_registration_status_msg(
				(slurm_node_registration_status_msg_t **)
				& (msg->data), buffer);
		 break;
	 case REQUEST_RESOURCE_ALLOCATION:
	 case REQUEST_SUBMIT_BATCH_JOB:
	 case REQUEST_JOB_WILL_RUN:
	 case REQUEST_ALLOCATION_AND_RUN_JOB_STEP:
	 case REQUEST_UPDATE_JOB:
		 rc = _unpack_job_desc_msg((job_desc_msg_t **) & (msg->data),
					   buffer);
		 break;
	 case REQUEST_OLD_JOB_RESOURCE_ALLOCATION:
		 rc = _unpack_old_job_desc_msg((old_job_alloc_msg_t **) &
					       (msg->data), buffer);
		 break;
	 case REQUEST_NODE_REGISTRATION_STATUS:
	 case REQUEST_RECONFIGURE:
	 case REQUEST_SHUTDOWN_IMMEDIATE:
	 case REQUEST_PING:
	 case REQUEST_CONTROL:
		 /* Message contains no body/information */
		 break;
	 case REQUEST_SHUTDOWN:
		 rc = _unpack_shutdown_msg((shutdown_msg_t **) & (msg->data),
					   buffer);
		 break;
	 case RESPONSE_SUBMIT_BATCH_JOB:
		 rc = _unpack_submit_response_msg((submit_response_msg_t **)
						  & (msg->data), buffer);
		 break;
	 case RESPONSE_RESOURCE_ALLOCATION:
	 case RESPONSE_JOB_WILL_RUN:
		 rc = _unpack_resource_allocation_response_msg(
				(resource_allocation_response_msg_t **)
				& (msg->data), buffer);
		 break;
	 case RESPONSE_ALLOCATION_AND_RUN_JOB_STEP:
		 rc = _unpack_resource_allocation_and_run_response_msg(
				(resource_allocation_and_run_response_msg_t **)
				&(msg->data), buffer);
		 break;
	 case REQUEST_UPDATE_NODE:
		 rc = _unpack_update_node_msg((update_node_msg_t **) &
					      (msg->data), buffer);
		 break;
	 case REQUEST_UPDATE_PARTITION:
		 rc = _unpack_update_partition_msg((update_part_msg_t **) &
						   (msg->data), buffer);
		 break;
	 case REQUEST_DELETE_PARTITION:
		 rc = _unpack_delete_partition_msg((delete_part_msg_t **) &
						   (msg->data), buffer);
		 break;
	 case REQUEST_LAUNCH_TASKS:
		 rc = _unpack_launch_tasks_request_msg(
					(launch_tasks_request_msg_t **)
					& (msg->data), buffer);
		 break;
	 case RESPONSE_LAUNCH_TASKS:
		 rc = _unpack_launch_tasks_response_msg(
					(launch_tasks_response_msg_t **)
					& (msg->data), buffer);
		 break;
	 case REQUEST_SPAWN_TASK:
		 rc = _unpack_spawn_task_request_msg(
					(spawn_task_request_msg_t **)
					& (msg->data), buffer);
		 break;
	 case REQUEST_REATTACH_TASKS:
		 rc = _unpack_reattach_tasks_request_msg(
			(reattach_tasks_request_msg_t **) & msg->data, 
			buffer);
		 break;
	 case RESPONSE_REATTACH_TASKS:
		 rc = _unpack_reattach_tasks_response_msg(
					(reattach_tasks_response_msg_t **) 
					& msg->data, buffer);
		 break;
	 case REQUEST_SIGNAL_TASKS:
	 case REQUEST_TERMINATE_TASKS:
		 rc = _unpack_cancel_tasks_msg((kill_tasks_msg_t **) &
					       (msg->data), buffer);
		 break;
	 case REQUEST_JOB_STEP_INFO:
		 rc = _unpack_job_step_info_req_msg(
				(job_step_info_request_msg_t **)
				& (msg->data), buffer);
		 break;
		/********  job_step_id_t Messages  ********/
	 case REQUEST_JOB_INFO:
		 rc = _unpack_job_info_request_msg((job_info_request_msg_t**)
					      & (msg->data), buffer);
		 break;
	 case REQUEST_CANCEL_JOB_STEP:
		 rc = _unpack_job_step_kill_msg((job_step_kill_msg_t **)
					        & (msg->data), buffer);
		 break;
	 case REQUEST_COMPLETE_JOB_ALLOCATION:
		 rc = _unpack_complete_job_allocation_msg(
			 (complete_job_allocation_msg_t **)&msg->data, buffer);
		 break;
	 case REQUEST_COMPLETE_BATCH_SCRIPT:
		 rc = _unpack_complete_batch_script_msg(
			 (complete_batch_script_msg_t **)&msg->data, buffer);
		 break;
	 case REQUEST_STEP_COMPLETE:
		 rc = _unpack_step_complete_msg((step_complete_msg_t
						 **) & (msg->data),
						buffer);
		 break;
	 case REQUEST_SIGNAL_JOB:
		 rc = _unpack_signal_job_msg((signal_job_msg_t **)&(msg->data),
					     buffer);
		 break;
	 case REQUEST_KILL_TIMELIMIT:
	 case REQUEST_TERMINATE_JOB:
		 rc = _unpack_kill_job_msg((kill_job_msg_t **) & (msg->data), 
					   buffer);
		 break;
	 case MESSAGE_EPILOG_COMPLETE:
		rc = _unpack_epilog_comp_msg((epilog_complete_msg_t **) 
		                           & (msg->data), buffer);
		break;
	 case REQUEST_UPDATE_JOB_TIME:
		 rc = _unpack_update_job_time_msg(
					(job_time_msg_t **)
					& (msg->data), buffer);
		 break;
	 case RESPONSE_RECONFIGURE:
	 case RESPONSE_SHUTDOWN:
	 case RESPONSE_CANCEL_JOB_STEP:
		 break;
	 case REQUEST_JOB_ATTACH:
		 break;
	 case RESPONSE_JOB_ATTACH:
		 break;
	 case RESPONSE_JOB_STEP_INFO:
		 rc = _unpack_job_step_info_response_msg(
					(job_step_info_response_msg_t **)
					& (msg->data), buffer);
		 break;
	 case REQUEST_JOB_RESOURCE:
		 break;
	 case RESPONSE_JOB_RESOURCE:
		 break;
	 case REQUEST_RUN_JOB_STEP:
		 break;
	 case RESPONSE_RUN_JOB_STEP:
		 break;
	 case MESSAGE_TASK_EXIT:
		 rc = _unpack_task_exit_msg((task_exit_msg_t **)
					    & (msg->data), buffer);
		 break;
	 case REQUEST_BATCH_JOB_LAUNCH:
		 rc = _unpack_batch_job_launch_msg((batch_job_launch_msg_t **)
						   & (msg->data), buffer);
		 break;
	 case MESSAGE_UPLOAD_ACCOUNTING_INFO:
		 break;
	 case RESPONSE_JOB_READY:
	 case RESPONSE_SLURM_RC:
		 rc = _unpack_return_code_msg((return_code_msg_t **)
					      & (msg->data), buffer);
		 break;
	 case RESPONSE_JOB_STEP_CREATE:
		 rc = _unpack_job_step_create_response_msg(
				(job_step_create_response_msg_t **)
				& msg->data, buffer);
		 break;
	 case REQUEST_JOB_STEP_CREATE:
		 rc = _unpack_job_step_create_request_msg(
					(job_step_create_request_msg_t **)
					& msg->data, buffer);
		 break;
	 case REQUEST_JOB_ID:
		 rc = _unpack_job_id_request_msg(
				(job_id_request_msg_t **) & msg->data,
				buffer);
		 break;
	 case RESPONSE_JOB_ID:
		 rc = _unpack_job_id_response_msg(
				(job_id_response_msg_t **) & msg->data,
				buffer);
		 break;
	 case SRUN_PING:
		rc = _unpack_srun_ping_msg((srun_ping_msg_t **) & msg->data, 
				buffer);
		break;
	 case SRUN_NODE_FAIL:
		rc = _unpack_srun_node_fail_msg((srun_node_fail_msg_t **)
				& msg->data, buffer);
		break;
	 case SRUN_TIMEOUT:
		rc = _unpack_srun_timeout_msg((srun_timeout_msg_t **)
				& msg->data, buffer);
		break;
	 case REQUEST_CHECKPOINT:
		rc = _unpack_checkpoint_msg((checkpoint_msg_t **)
				& msg->data, buffer);
		break;
	 case REQUEST_CHECKPOINT_COMP:
		rc = _unpack_checkpoint_comp((checkpoint_comp_msg_t **)
				& msg->data, buffer);
		break;
	 case RESPONSE_CHECKPOINT:
	 case RESPONSE_CHECKPOINT_COMP:
		rc = _unpack_checkpoint_resp_msg((checkpoint_resp_msg_t **)
				& msg->data, buffer);
		break;
	 case REQUEST_SUSPEND:
		rc = _unpack_suspend_msg((suspend_msg_t **) &msg->data, 
				buffer);
		break;

	 case REQUEST_JOB_READY:
		rc = _unpack_job_ready_msg((job_id_msg_t **)
				& msg->data, buffer);
		break;
	 case REQUEST_NODE_SELECT_INFO:
		rc = _unpack_node_select_info_req_msg(
				(node_info_select_request_msg_t **) &msg->data,
				buffer);
		break;
	 case RESPONSE_NODE_SELECT_INFO:
		rc = _unpack_node_select_info_msg((node_select_info_msg_t **) &
			(msg->data), buffer);
		break;
	 case MESSAGE_JOBACCT_DATA:
		rc = _unpack_jobacct_data( (jobacct_msg_t **) 
				& msg->data, buffer);
		break;
	 case REQUEST_FILE_BCAST:
		rc = _unpack_file_bcast( (file_bcast_msg_t **)
				& msg->data, buffer);
		break;
	 case PMI_KVS_PUT_REQ:
	 case PMI_KVS_GET_RESP:
		rc = _unpack_kvs_data((struct kvs_comm_set **) &msg->data, 
				buffer);
		break;
	 case PMI_KVS_GET_REQ:
		rc = _unpack_kvs_get((kvs_get_msg_t **) &msg->data, buffer);
		break;
	 case PMI_KVS_PUT_RESP:
		break;	/* no data */
	 default:
		 debug("No unpack method for msg type %i", msg->msg_type);
		 return EINVAL;
		 break;
	}

	if (rc)
		error("Malformed RPC of type %u received", msg->msg_type);
	return rc;
}

static void
_pack_update_node_msg(update_node_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	packstr(msg->node_names, buffer);
	pack16((uint16_t)msg->node_state, buffer);
	packstr(msg->reason, buffer);
}

static int
_unpack_update_node_msg(update_node_msg_t ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	update_node_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(update_node_msg_t));
	*msg = tmp_ptr;

	safe_unpackstr_xmalloc(&tmp_ptr->node_names, &uint16_tmp, buffer);
	safe_unpack16(&tmp_ptr->node_state, buffer);
	safe_unpackstr_xmalloc(&tmp_ptr->reason, &uint16_tmp, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->node_names);
	xfree(tmp_ptr->reason);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_node_registration_status_msg(slurm_node_registration_status_msg_t *
				   msg, Buf buffer)
{
	int i;
	xassert(msg != NULL);

	pack_time(msg->timestamp, buffer);
	pack32((uint32_t)msg->status, buffer);
	packstr(msg->node_name, buffer);
	pack32((uint32_t)msg->cpus, buffer);
	pack32((uint32_t)msg->real_memory_size, buffer);
	pack32((uint32_t)msg->temporary_disk_space, buffer);
	pack32((uint32_t)msg->job_count, buffer);
	for (i = 0; i < msg->job_count; i++) {
		pack32((uint32_t)msg->job_id[i], buffer);
	}
	for (i = 0; i < msg->job_count; i++) {
		pack16((uint16_t)msg->step_id[i], buffer);
	}
	pack16((uint16_t)msg->startup, buffer);
	if (msg->startup)
		switch_g_pack_node_info(msg->switch_nodeinfo, buffer);
}

static int
_unpack_node_registration_status_msg(slurm_node_registration_status_msg_t
				     ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	int i;
	slurm_node_registration_status_msg_t *node_reg_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	node_reg_ptr = xmalloc(sizeof(slurm_node_registration_status_msg_t));
	*msg = node_reg_ptr;

	/* unpack timestamp of snapshot */
	safe_unpack_time(&node_reg_ptr->timestamp, buffer);
	/* load the data values */
	safe_unpack32(&node_reg_ptr->status, buffer);
	safe_unpackstr_xmalloc(&node_reg_ptr->node_name, &uint16_tmp, buffer);
	safe_unpack32(&node_reg_ptr->cpus, buffer);
	safe_unpack32(&node_reg_ptr->real_memory_size, buffer);
	safe_unpack32(&node_reg_ptr->temporary_disk_space, buffer);
	safe_unpack32(&node_reg_ptr->job_count, buffer);
	node_reg_ptr->job_id =
	    xmalloc(sizeof(uint32_t) * node_reg_ptr->job_count);
	for (i = 0; i < node_reg_ptr->job_count; i++) {
		safe_unpack32(&node_reg_ptr->job_id[i], buffer);
	}
	node_reg_ptr->step_id =
	    xmalloc(sizeof(uint16_t) * node_reg_ptr->job_count);
	for (i = 0; i < node_reg_ptr->job_count; i++) {
		safe_unpack16(&node_reg_ptr->step_id[i], buffer);
	}

	safe_unpack16(&node_reg_ptr->startup, buffer);
	if (node_reg_ptr->startup
	&&  (switch_g_alloc_node_info(&node_reg_ptr->switch_nodeinfo)
	||   switch_g_unpack_node_info(node_reg_ptr->switch_nodeinfo, buffer)))
		goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	xfree(node_reg_ptr->node_name);
	xfree(node_reg_ptr->job_id);
	xfree(node_reg_ptr->step_id);
	switch_g_free_node_info(&node_reg_ptr->switch_nodeinfo);
	xfree(node_reg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_resource_allocation_response_msg(resource_allocation_response_msg_t *
				       msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->error_code, buffer);
	pack32((uint32_t)msg->job_id, buffer);
	packstr(msg->node_list, buffer);

	pack16((uint16_t)msg->num_cpu_groups, buffer);
	if (msg->num_cpu_groups) {
		pack32_array(msg->cpus_per_node, msg->num_cpu_groups, buffer);
		pack32_array(msg->cpu_count_reps, msg->num_cpu_groups, buffer);
	}

	pack16((uint16_t)msg->node_cnt, buffer);
	if (msg->node_cnt > 0)
		_pack_slurm_addr_array(msg->node_addr, msg->node_cnt, buffer);

	select_g_pack_jobinfo(msg->select_jobinfo, buffer);
}

static int
_unpack_resource_allocation_response_msg(resource_allocation_response_msg_t
					 ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	uint32_t uint32_tmp;
	resource_allocation_response_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(resource_allocation_response_msg_t));
	*msg = tmp_ptr;

	/* load the data values */
	safe_unpack32(&tmp_ptr->error_code, buffer);
	safe_unpack32(&tmp_ptr->job_id, buffer);
	safe_unpackstr_xmalloc(&tmp_ptr->node_list, &uint16_tmp, buffer);

	safe_unpack16(&tmp_ptr->num_cpu_groups, buffer);
	if (tmp_ptr->num_cpu_groups > 0) {
		safe_unpack32_array((uint32_t **) &
				    (tmp_ptr->cpus_per_node), &uint32_tmp,
				    buffer);
		if (tmp_ptr->num_cpu_groups != uint32_tmp)
			goto unpack_error;
		safe_unpack32_array((uint32_t **) &
				    (tmp_ptr->cpu_count_reps), &uint32_tmp,
				    buffer);
		if (tmp_ptr->num_cpu_groups != uint32_tmp)
			goto unpack_error;
	} else {
		tmp_ptr->cpus_per_node = NULL;
		tmp_ptr->cpu_count_reps = NULL;
	}

	safe_unpack16(&tmp_ptr->node_cnt, buffer);
	if (tmp_ptr->node_cnt > 0) {
		if (_unpack_slurm_addr_array(&(tmp_ptr->node_addr),
					     &uint16_tmp, buffer))
			goto unpack_error;
		if (uint16_tmp != tmp_ptr->node_cnt)
			goto unpack_error;
	} else
		tmp_ptr->node_addr = NULL;

	if (select_g_alloc_jobinfo (&tmp_ptr->select_jobinfo)
	||  select_g_unpack_jobinfo(tmp_ptr->select_jobinfo, buffer))
		goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	select_g_free_jobinfo(&tmp_ptr->select_jobinfo);
	xfree(tmp_ptr->node_list);
	xfree(tmp_ptr->cpus_per_node);
	xfree(tmp_ptr->cpu_count_reps);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
 _pack_resource_allocation_and_run_response_msg
    (resource_allocation_and_run_response_msg_t * msg, Buf buffer) {
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
	packstr(msg->node_list, buffer);
	pack16((uint16_t)msg->num_cpu_groups, buffer);
	if (msg->num_cpu_groups > 0) {
		pack32_array(msg->cpus_per_node, msg->num_cpu_groups, buffer);
		pack32_array(msg->cpu_count_reps, msg->num_cpu_groups, buffer);
	}
	pack32((uint32_t)msg->job_step_id, buffer);

	pack16((uint16_t)msg->node_cnt, buffer);
	if (msg->node_cnt)
		_pack_slurm_addr_array(msg->node_addr, msg->node_cnt, buffer);

	slurm_cred_pack(msg->cred, buffer);
	switch_pack_jobinfo(msg->switch_job, buffer);
}

static int
 _unpack_resource_allocation_and_run_response_msg
    (resource_allocation_and_run_response_msg_t ** msg, Buf buffer) {
	uint16_t uint16_tmp;
	uint32_t uint32_tmp;
	resource_allocation_and_run_response_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(resource_allocation_and_run_response_msg_t));
	*msg = tmp_ptr;

	/* load the data values */
	safe_unpack32(&tmp_ptr->job_id, buffer);
	safe_unpackstr_xmalloc(&tmp_ptr->node_list, &uint16_tmp, buffer);
	safe_unpack16(&tmp_ptr->num_cpu_groups, buffer);

	if (tmp_ptr->num_cpu_groups > 0) {
		tmp_ptr->cpus_per_node = (uint32_t *)
		    xmalloc(sizeof(uint32_t) * tmp_ptr->num_cpu_groups);
		tmp_ptr->cpu_count_reps = (uint32_t *)
		    xmalloc(sizeof(uint32_t) * tmp_ptr->num_cpu_groups);
		safe_unpack32_array((uint32_t **) &
				    (tmp_ptr->cpus_per_node), 
				    &uint32_tmp, buffer);
		if (tmp_ptr->num_cpu_groups != uint32_tmp)
			goto unpack_error;
		safe_unpack32_array((uint32_t **) &
				    (tmp_ptr->cpu_count_reps), 
				    &uint32_tmp, buffer);
		if (tmp_ptr->num_cpu_groups != uint32_tmp)
			goto unpack_error;
	}

	safe_unpack32(&tmp_ptr->job_step_id, buffer);
	safe_unpack16(&tmp_ptr->node_cnt, buffer);
	if (tmp_ptr->node_cnt > 0) {
		if (_unpack_slurm_addr_array(&(tmp_ptr->node_addr),
					     &uint16_tmp, buffer))
			goto unpack_error;
		if (uint16_tmp != tmp_ptr->node_cnt)
			goto unpack_error;
	} else
		tmp_ptr->node_addr = NULL;

	if (!(tmp_ptr->cred = slurm_cred_unpack(buffer)))
		goto unpack_error;

	switch_alloc_jobinfo(&tmp_ptr->switch_job);
	if (switch_unpack_jobinfo(tmp_ptr->switch_job, buffer) < 0) {
		error("switch_unpack_jobinfo: %m");
		switch_free_jobinfo(tmp_ptr->switch_job);
		goto unpack_error;
	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->node_list);
	xfree(tmp_ptr->cpus_per_node);
	xfree(tmp_ptr->cpu_count_reps);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_submit_response_msg(submit_response_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->step_id, buffer);
	pack32((uint32_t)msg->error_code, buffer);
}

static int
_unpack_submit_response_msg(submit_response_msg_t ** msg, Buf buffer)
{
	submit_response_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(submit_response_msg_t));
	*msg = tmp_ptr;

	/* load the data values */
	safe_unpack32(&tmp_ptr->job_id, buffer);
	safe_unpack32(&tmp_ptr->step_id, buffer);
	safe_unpack32(&tmp_ptr->error_code, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static int
_unpack_node_info_msg(node_info_msg_t ** msg, Buf buffer)
{
	int i;
	node_info_t *node = NULL;

	xassert(msg != NULL);
	*msg = xmalloc(sizeof(node_info_msg_t));

	/* load buffer's header (data structure version and time) */
	safe_unpack32(&((*msg)->record_count), buffer);
	safe_unpack_time(&((*msg)->last_update), buffer);

	node = (*msg)->node_array =
	    xmalloc(sizeof(node_info_t) * (*msg)->record_count);

	/* load individual job info */
	for (i = 0; i < (*msg)->record_count; i++) {
		if (_unpack_node_info_members(&node[i], buffer))
			goto unpack_error;

	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(node);
	xfree(*msg);
	*msg = NULL;
	return SLURM_ERROR;
}

static int
_unpack_node_info_members(node_info_t * node, Buf buffer)
{
	uint16_t uint16_tmp;

	xassert(node != NULL);

	safe_unpackstr_xmalloc(&node->name, &uint16_tmp, buffer);
	safe_unpack16(&node->node_state, buffer);
	safe_unpack32(&node->cpus, buffer);
	safe_unpack32(&node->real_memory, buffer);
	safe_unpack32(&node->tmp_disk, buffer);
	safe_unpack32(&node->weight, buffer);
	safe_unpackstr_xmalloc(&node->features, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&node->reason, &uint16_tmp, buffer);

	return SLURM_SUCCESS;

      unpack_error:
	xfree(node->name);
	xfree(node->features);
	xfree(node->reason);
	return SLURM_ERROR;
}

static int _unpack_node_select_info_msg(node_select_info_msg_t ** msg,
		Buf buffer)
{
	xassert(msg != NULL);

	return select_g_unpack_node_info(msg, buffer);
}

static void
_pack_update_partition_msg(update_part_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	packstr(msg->allow_groups, buffer);
	pack16((uint16_t)msg-> default_part, buffer);
	pack32((uint32_t)msg-> max_time,     buffer);
	pack32((uint32_t)msg-> max_nodes,    buffer);
	pack32((uint32_t)msg-> min_nodes,    buffer);
	packstr(msg->name,         buffer);
	packstr(msg->nodes,        buffer);

	pack16((uint16_t)msg-> hidden,      buffer);
	pack16((uint16_t)msg-> root_only,    buffer);
	pack16((uint16_t)msg-> shared,       buffer);
	pack16((uint16_t)msg-> state_up,     buffer);
}

static int
_unpack_update_partition_msg(update_part_msg_t ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	update_part_msg_t *tmp_ptr;

	xassert(msg != NULL);

	/* alloc memory for structure */
	tmp_ptr = xmalloc(sizeof(update_part_msg_t));
	*msg = tmp_ptr;

	safe_unpackstr_xmalloc(&tmp_ptr->allow_groups, &uint16_tmp, buffer);
	safe_unpack16(&tmp_ptr->default_part, buffer);
	safe_unpack32(&tmp_ptr->max_time, buffer);
	safe_unpack32(&tmp_ptr->max_nodes, buffer);
	safe_unpack32(&tmp_ptr->min_nodes, buffer);
	safe_unpackstr_xmalloc(&tmp_ptr->name, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&tmp_ptr->nodes, &uint16_tmp, buffer);

	safe_unpack16(&tmp_ptr->hidden, buffer);
	safe_unpack16(&tmp_ptr->root_only, buffer);
	safe_unpack16(&tmp_ptr->shared, buffer);
	safe_unpack16(&tmp_ptr->state_up, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->name);
	xfree(tmp_ptr->nodes);
	xfree(tmp_ptr->allow_groups);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_delete_partition_msg(delete_part_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	packstr(msg->name,         buffer);
}

static int
_unpack_delete_partition_msg(delete_part_msg_t ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	delete_part_msg_t *tmp_ptr;

	xassert(msg != NULL);

	/* alloc memory for structure */
	tmp_ptr = xmalloc(sizeof(delete_part_msg_t));
	*msg = tmp_ptr;

	safe_unpackstr_xmalloc(&tmp_ptr->name, &uint16_tmp, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->name);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}
static void
_pack_job_step_create_request_msg(job_step_create_request_msg_t
				  * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->user_id, buffer);
	pack32((uint32_t)msg->node_count, buffer);
	pack32((uint32_t)msg->cpu_count, buffer);
	pack32((uint32_t)msg->num_tasks, buffer);

	pack16((uint16_t)msg->relative, buffer);
	pack16((uint16_t)msg->task_dist, buffer);
	pack16((uint16_t)msg->port, buffer);
	packstr(msg->host, buffer);
	packstr(msg->name, buffer);
	packstr(msg->network, buffer);
	packstr(msg->node_list, buffer);
}

static int
_unpack_job_step_create_request_msg(job_step_create_request_msg_t ** msg,
				    Buf buffer)
{
	uint16_t uint16_tmp;
	job_step_create_request_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(job_step_create_request_msg_t));
	*msg = tmp_ptr;

	safe_unpack32(&(tmp_ptr->job_id), buffer);
	safe_unpack32(&(tmp_ptr->user_id), buffer);
	safe_unpack32(&(tmp_ptr->node_count), buffer);
	safe_unpack32(&(tmp_ptr->cpu_count), buffer);
	safe_unpack32(&(tmp_ptr->num_tasks), buffer);

	safe_unpack16(&(tmp_ptr->relative), buffer);
	safe_unpack16(&(tmp_ptr->task_dist), buffer);
	safe_unpack16(&(tmp_ptr->port), buffer);
	safe_unpackstr_xmalloc(&(tmp_ptr->host), &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&(tmp_ptr->name), &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&(tmp_ptr->network), &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&(tmp_ptr->node_list), &uint16_tmp, buffer);

	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->host);
	xfree(tmp_ptr->name);
	xfree(tmp_ptr->network);
	xfree(tmp_ptr->node_list);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_kill_job_msg(kill_job_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id,  buffer);
	pack32((uint32_t)msg->job_uid, buffer);
	packstr(msg->nodes, buffer);
	select_g_pack_jobinfo(msg->select_jobinfo, buffer);
}

static int
_unpack_kill_job_msg(kill_job_msg_t ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	kill_job_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg);
	tmp_ptr = xmalloc(sizeof(kill_job_msg_t));
	*msg = tmp_ptr;

	safe_unpack32(&(tmp_ptr->job_id),  buffer);
	safe_unpack32(&(tmp_ptr->job_uid), buffer);
	safe_unpackstr_xmalloc(&(tmp_ptr->nodes), &uint16_tmp, buffer);
	if (select_g_alloc_jobinfo (&tmp_ptr->select_jobinfo)
	||  select_g_unpack_jobinfo(tmp_ptr->select_jobinfo, buffer))
		goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->nodes);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_signal_job_msg(signal_job_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id,  buffer);
	pack32((uint32_t)msg->signal, buffer);
	debug("_pack_signal_job_msg signal = %d", msg->signal);
}

static int
_unpack_signal_job_msg(signal_job_msg_t ** msg, Buf buffer)
{
	signal_job_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg);
	tmp_ptr = xmalloc(sizeof(signal_job_msg_t));
	*msg = tmp_ptr;

	safe_unpack32(&(tmp_ptr->job_id), buffer);
	safe_unpack32(&(tmp_ptr->signal), buffer);
	debug("_unpack_signal_job_msg signal = %d", tmp_ptr->signal);

	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void 
_pack_epilog_comp_msg(epilog_complete_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->return_code, buffer);
	packstr(msg->node_name, buffer);
	switch_g_pack_node_info(msg->switch_nodeinfo, buffer);
}

static int  
_unpack_epilog_comp_msg(epilog_complete_msg_t ** msg, Buf buffer)
{
	epilog_complete_msg_t *tmp_ptr;
	uint16_t uint16_tmp;

	/* alloc memory for structure */
	xassert(msg);
	tmp_ptr = xmalloc(sizeof(epilog_complete_msg_t));
	*msg = tmp_ptr;

	safe_unpack32(&(tmp_ptr->job_id), buffer);
	safe_unpack32(&(tmp_ptr->return_code), buffer);
	safe_unpackstr_xmalloc(& (tmp_ptr->node_name), &uint16_tmp, buffer);
	if (switch_g_alloc_node_info(&tmp_ptr->switch_nodeinfo)
	||  switch_g_unpack_node_info(tmp_ptr->switch_nodeinfo, buffer))
		goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->node_name);
	switch_g_free_node_info(&tmp_ptr->switch_nodeinfo);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_update_job_time_msg(job_time_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
	pack_time(msg->expiration_time, buffer);
}

static int
_unpack_update_job_time_msg(job_time_msg_t ** msg, Buf buffer)
{
	job_time_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg);
	tmp_ptr = xmalloc(sizeof(job_time_msg_t));
	*msg = tmp_ptr;

	safe_unpack32(&(tmp_ptr->job_id), buffer);
	safe_unpack_time(& (tmp_ptr->expiration_time), buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_job_step_create_response_msg(job_step_create_response_msg_t * msg,
				   Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_step_id, buffer);
	packstr(msg->node_list, buffer);
	slurm_cred_pack(msg->cred, buffer);
	switch_pack_jobinfo(msg->switch_job, buffer);

}

static int
_unpack_job_step_create_response_msg(job_step_create_response_msg_t ** msg,
				     Buf buffer)
{
	uint16_t uint16_tmp;
	job_step_create_response_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(job_step_create_response_msg_t));
	*msg = tmp_ptr;

	safe_unpack32(&tmp_ptr->job_step_id, buffer);
	safe_unpackstr_xmalloc(&tmp_ptr->node_list, &uint16_tmp, buffer);
	if (!(tmp_ptr->cred = slurm_cred_unpack(buffer)))
		goto unpack_error;

	switch_alloc_jobinfo(&tmp_ptr->switch_job);
	if (switch_unpack_jobinfo(tmp_ptr->switch_job, buffer)) {
		error("switch_unpack_jobinfo: %m");
		switch_free_jobinfo(tmp_ptr->switch_job);
		goto unpack_error;
	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr->node_list);
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}


static void
_pack_partition_info_msg(slurm_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	packmem_array(msg->data, msg->data_size, buffer);
}

static int
_unpack_partition_info_msg(partition_info_msg_t ** msg, Buf buffer)
{
	int i;
	partition_info_t *partition = NULL;

	xassert(msg != NULL);
	*msg = xmalloc(sizeof(partition_info_msg_t));

	/* load buffer's header (data structure version and time) */
	safe_unpack32(&((*msg)->record_count), buffer);
	safe_unpack_time(&((*msg)->last_update), buffer);

	partition = (*msg)->partition_array =
	    xmalloc(sizeof(partition_info_t) * (*msg)->record_count);

	/* load individual job info */
	for (i = 0; i < (*msg)->record_count; i++) {
		if (_unpack_partition_info_members(&partition[i], buffer))
			goto unpack_error;
	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(partition);
	xfree(*msg);
	*msg = NULL;
	return SLURM_ERROR;
}


static int
_unpack_partition_info_members(partition_info_t * part, Buf buffer)
{
	uint16_t uint16_tmp;
	char *node_inx_str = NULL;

	safe_unpackstr_xmalloc(&part->name, &uint16_tmp, buffer);
	if (part->name == NULL)
		part->name = xmalloc(1);	/* part->name = "" implicit */
	safe_unpack32(&part->max_time,     buffer);
	safe_unpack32(&part->max_nodes,    buffer);
	safe_unpack32(&part->min_nodes,    buffer);
	safe_unpack32(&part->total_nodes,  buffer);
	safe_unpack32(&part->node_scaling, buffer);
	
	safe_unpack32(&part->total_cpus,   buffer);
	safe_unpack16(&part->default_part, buffer);
	safe_unpack16(&part->hidden,       buffer);
	safe_unpack16(&part->root_only,    buffer);
	safe_unpack16(&part->shared,       buffer);

	safe_unpack16(&part->state_up, buffer);
	safe_unpackstr_xmalloc(&part->allow_groups, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&part->nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&node_inx_str, &uint16_tmp, buffer);
	if (node_inx_str == NULL)
		part->node_inx = bitfmt2int("");
	else {
		part->node_inx = bitfmt2int(node_inx_str);
		xfree(node_inx_str);
		node_inx_str = NULL;
	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(part->name);
	xfree(part->allow_groups);
	xfree(part->nodes);
	xfree(node_inx_str);
	return SLURM_ERROR;
}

/* pack_job_step_info_members
 * pack selected fields of the description of a job into a buffer
 * IN job_id, step_id, user_id, start_time, partition, nodes - job info
 * IN/OUT buffer - destination of the pack, contains pointers that are 
 *			automatically updated
 */
void
pack_job_step_info_members(uint32_t job_id, uint16_t step_id,
			   uint32_t user_id, uint32_t num_tasks,
			   time_t start_time, char *partition, 
			   char *nodes, char *name, char *network,
			   Buf buffer)
{
	pack32((uint32_t)job_id, buffer);
	pack16((uint16_t)step_id, buffer);
	pack32((uint32_t)user_id, buffer);
	pack32((uint32_t)num_tasks, buffer);

	pack_time(start_time, buffer);
	packstr(partition, buffer);
	packstr(nodes, buffer);
	packstr(name, buffer);
	packstr(network, buffer);
}

/* pack_job_step_info
 * packs a slurm job steps info
 * IN step - pointer to the job step info
 * IN/OUT buffer - destination of the pack, contains pointers that are 
 *			automatically updated
 */
void
pack_job_step_info(job_step_info_t * step, Buf buffer)
{
	pack_job_step_info_members(step->job_id,
				   step->step_id,
				   step->user_id,
				   step->num_tasks,
				   step->start_time,
				   step->partition, step->nodes, 
				   step->name, step->network, buffer);
}

/* _unpack_job_step_info_members
 * unpacks a set of slurm job step info for one job step
 * OUT step - pointer to the job step info buffer
 * IN/OUT buffer - source of the unpack, contains pointers that are 
 *			automatically updated
 */
static int
_unpack_job_step_info_members(job_step_info_t * step, Buf buffer)
{
	uint16_t uint16_tmp = 0;

	safe_unpack32(&step->job_id, buffer);
	safe_unpack16(&step->step_id, buffer);
	safe_unpack32(&step->user_id, buffer);
	safe_unpack32(&step->num_tasks, buffer);

	safe_unpack_time(&step->start_time, buffer);
	safe_unpackstr_xmalloc(&step->partition, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&step->nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&step->name, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&step->network, &uint16_tmp, buffer);

	return SLURM_SUCCESS;

      unpack_error:
	xfree(step->partition);
	xfree(step->nodes);
	xfree(step->name);
	xfree(step->network);
	return SLURM_ERROR;
}

static int
_unpack_job_step_info_response_msg(job_step_info_response_msg_t
				   ** msg, Buf buffer)
{
	int i = 0;
	job_step_info_t *step;

	xassert(msg != NULL);
	*msg = xmalloc(sizeof(job_step_info_response_msg_t));

	safe_unpack_time(&(*msg)->last_update, buffer);
	safe_unpack32(&(*msg)->job_step_count, buffer);

	step = (*msg)->job_steps =
	    xmalloc(sizeof(job_step_info_t) * (*msg)->job_step_count);

	for (i = 0; i < (*msg)->job_step_count; i++)
		if (_unpack_job_step_info_members(&step[i], buffer))
			goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	xfree(step);
	xfree(*msg);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_buffer_msg(slurm_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);
	packmem_array(msg->data, msg->data_size, buffer);
}

static int
_unpack_job_info_msg(job_info_msg_t ** msg, Buf buffer)
{
	int i;
	job_info_t *job = NULL;

	xassert(msg != NULL);
	*msg = xmalloc(sizeof(job_info_msg_t));

	/* load buffer's header (data structure version and time) */
	safe_unpack32(&((*msg)->record_count), buffer);
	safe_unpack_time(&((*msg)->last_update), buffer);
	job = (*msg)->job_array =
	    xmalloc(sizeof(job_info_t) * (*msg)->record_count);

	/* load individual job info */
	for (i = 0; i < (*msg)->record_count; i++) {
		if (_unpack_job_info_members(&job[i], buffer))
			goto unpack_error;
	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(job);
	xfree(*msg);
	*msg = NULL;
	return SLURM_ERROR;
}

/* _unpack_job_info_members
 * unpacks a set of slurm job info for one job
 * OUT job - pointer to the job info buffer
 * IN/OUT buffer - source of the unpack, contains pointers that are 
 *			automatically updated
 */
static int
_unpack_job_info_members(job_info_t * job, Buf buffer)
{
	uint16_t uint16_tmp;
	char *node_inx_str;

	safe_unpack32(&job->job_id, buffer);
	safe_unpack32(&job->user_id, buffer);
	safe_unpack32(&job->group_id, buffer);

	safe_unpack16(&job->job_state,  buffer);
	safe_unpack16(&job->batch_flag, buffer);
	safe_unpack32(&job->alloc_sid,  buffer);
	safe_unpack32(&job->time_limit, buffer);

	safe_unpack_time(&job->submit_time, buffer);
	safe_unpack_time(&job->start_time, buffer);
	safe_unpack_time(&job->end_time, buffer);
	safe_unpack_time(&job->suspend_time, buffer);
	safe_unpack_time(&job->pre_sus_time, buffer);
	safe_unpack32(&job->priority, buffer);

	safe_unpackstr_xmalloc(&job->nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job->partition, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job->account, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job->network, &uint16_tmp, buffer);
	safe_unpack32(&job->dependency, buffer);

	safe_unpackstr_xmalloc(&job->name, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job->alloc_node, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&node_inx_str, &uint16_tmp, buffer);
	if (node_inx_str == NULL)
		job->node_inx = bitfmt2int("");
	else {
		job->node_inx = bitfmt2int(node_inx_str);
		xfree(node_inx_str);
	}
	safe_unpack32(&job->num_procs, buffer);
	
	if (select_g_alloc_jobinfo(&job->select_jobinfo) 
	||  select_g_unpack_jobinfo(job->select_jobinfo, buffer))
		goto unpack_error;

	safe_unpackstr_xmalloc(&job->features, &uint16_tmp, buffer);

	safe_unpack32(&job->num_nodes, buffer);
	safe_unpack16(&job->shared, buffer);
	safe_unpack16(&job->contiguous, buffer);
	safe_unpack16(&job->cpus_per_task, buffer);

	safe_unpack32(&job->min_procs, buffer);
	safe_unpack32(&job->min_memory, buffer);
	safe_unpack32(&job->min_tmp_disk, buffer);
	safe_unpack16(&job->wait_reason, buffer);

	safe_unpackstr_xmalloc(&job->req_nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&node_inx_str, &uint16_tmp, buffer);
	if (node_inx_str == NULL)
		job->req_node_inx = bitfmt2int("");
	else {
		job->req_node_inx = bitfmt2int(node_inx_str);
		xfree(node_inx_str);
	}
	safe_unpackstr_xmalloc(&job->exc_nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&node_inx_str, &uint16_tmp, buffer);
	if (node_inx_str == NULL)
		job->exc_node_inx = bitfmt2int("");
	else {
		job->exc_node_inx = bitfmt2int(node_inx_str);
		xfree(node_inx_str);
	}

	return SLURM_SUCCESS;

      unpack_error:
	xfree(job->nodes);
	xfree(job->partition);
	xfree(job->account);
	xfree(job->name);
	xfree(job->alloc_node);
	xfree(job->node_inx);
	select_g_free_jobinfo(&job->select_jobinfo);
	xfree(job->features);
	xfree(job->req_nodes);
	xfree(job->req_node_inx);
	xfree(job->exc_nodes);
	xfree(job->exc_node_inx);
	xfree(job->network);
	return SLURM_ERROR;
}

static void
_pack_slurm_ctl_conf_msg(slurm_ctl_conf_info_msg_t * build_ptr, Buf buffer)
{
	pack_time(build_ptr->last_update, buffer);
	packstr(build_ptr->authtype, buffer);
	packstr(build_ptr->backup_addr, buffer);
	packstr(build_ptr->backup_controller, buffer);
	pack16((uint16_t)build_ptr->cache_groups, buffer);
	packstr(build_ptr->checkpoint_type, buffer);
	packstr(build_ptr->control_addr, buffer);
	packstr(build_ptr->control_machine, buffer);
	packstr(build_ptr->epilog, buffer);
	pack16((uint16_t)build_ptr->fast_schedule, buffer);
	pack32((uint32_t)build_ptr->first_job_id, buffer);
	pack16((uint16_t)build_ptr->inactive_limit, buffer);
	packstr(build_ptr->job_acct_loc, buffer);
	packstr(build_ptr->job_acct_parameters, buffer);
	packstr(build_ptr->job_acct_type, buffer);
	packstr(build_ptr->job_comp_loc, buffer);
	packstr(build_ptr->job_comp_type, buffer);
	pack16((uint16_t)build_ptr->kill_wait, buffer);
	pack16((uint16_t)build_ptr->max_job_cnt, buffer);
	pack16((uint16_t)build_ptr->min_job_age, buffer);
	packstr(build_ptr->mpi_default, buffer);
	packstr(build_ptr->plugindir, buffer);
	packstr(build_ptr->proctrack_type, buffer);
	packstr(build_ptr->prolog, buffer);
        packstr(build_ptr->propagate_rlimits, buffer);
        packstr(build_ptr->propagate_rlimits_except, buffer);
	pack16((uint16_t)build_ptr->ret2service, buffer);
	packstr(build_ptr->schedauth, buffer);
	pack16((uint16_t)build_ptr->schedport, buffer);
	pack16((uint16_t)build_ptr->schedrootfltr, buffer);
	packstr(build_ptr->schedtype, buffer);
	packstr(build_ptr->select_type, buffer);
	pack32((uint32_t)build_ptr->slurm_user_id, buffer);
	packstr(build_ptr->slurm_user_name, buffer);
	pack16((uint16_t)build_ptr->slurmctld_debug, buffer);
	packstr(build_ptr->slurmctld_logfile, buffer);
	packstr(build_ptr->slurmctld_pidfile, buffer);
	pack32((uint32_t)build_ptr->slurmctld_port, buffer);
	pack16((uint16_t)build_ptr->slurmctld_timeout, buffer);
	pack16((uint16_t)build_ptr->slurmd_debug, buffer);
	packstr(build_ptr->slurmd_logfile, buffer);
	packstr(build_ptr->slurmd_pidfile, buffer);
#ifndef MULTIPLE_SLURMD
	pack32((uint32_t)build_ptr->slurmd_port, buffer);
#endif
	packstr(build_ptr->slurmd_spooldir, buffer);
	debug2("Packing string %s", build_ptr->slurmd_spooldir);
	pack16((uint16_t)build_ptr->slurmd_timeout, buffer);
	packstr(build_ptr->slurm_conf, buffer);
	packstr(build_ptr->state_save_location, buffer);
	packstr(build_ptr->switch_type, buffer);
	packstr(build_ptr->task_epilog, buffer);
	packstr(build_ptr->task_prolog, buffer);
	packstr(build_ptr->task_plugin, buffer);
	packstr(build_ptr->tmp_fs, buffer);
	pack16((uint16_t)build_ptr->wait_time, buffer);
	packstr(build_ptr->job_credential_private_key, buffer);
	packstr(build_ptr->job_credential_public_certificate, buffer);
	packstr(build_ptr->srun_prolog, buffer);
	packstr(build_ptr->srun_epilog, buffer);
	packstr(build_ptr->node_prefix, buffer);
	pack16((uint16_t)build_ptr->tree_width, buffer);
}

static int
_unpack_slurm_ctl_conf_msg(slurm_ctl_conf_info_msg_t **
			   build_buffer_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	slurm_ctl_conf_info_msg_t *build_ptr;

	/* alloc memory for structure */
	build_ptr = xmalloc(sizeof(slurm_ctl_conf_t));
	*build_buffer_ptr = build_ptr;

	/* load the data values */
	/* unpack timestamp of snapshot */
	safe_unpack_time(&build_ptr->last_update, buffer);
	safe_unpackstr_xmalloc(&build_ptr->authtype, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->backup_addr, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->backup_controller, &uint16_tmp,
			       buffer);
	safe_unpack16(&build_ptr->cache_groups, buffer);
	safe_unpackstr_xmalloc(&build_ptr->checkpoint_type, &uint16_tmp,
			       buffer);
	safe_unpackstr_xmalloc(&build_ptr->control_addr, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->control_machine, &uint16_tmp,
			       buffer);
	safe_unpackstr_xmalloc(&build_ptr->epilog, &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->fast_schedule, buffer);
	safe_unpack32(&build_ptr->first_job_id, buffer);
	safe_unpack16(&build_ptr->inactive_limit, buffer);
	safe_unpackstr_xmalloc(&build_ptr->job_acct_loc, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->job_acct_parameters, &uint16_tmp,
			       buffer);
	safe_unpackstr_xmalloc(&build_ptr->job_acct_type, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->job_comp_loc, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->job_comp_type, &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->kill_wait, buffer);
	safe_unpack16(&build_ptr->max_job_cnt, buffer);
	safe_unpack16(&build_ptr->min_job_age, buffer);
	safe_unpackstr_xmalloc(&build_ptr->mpi_default, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->plugindir, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->proctrack_type, &uint16_tmp, 
			       buffer);
	safe_unpackstr_xmalloc(&build_ptr->prolog, &uint16_tmp, buffer);
        safe_unpackstr_xmalloc(&build_ptr->propagate_rlimits,
                               &uint16_tmp, buffer);
        safe_unpackstr_xmalloc(&build_ptr->propagate_rlimits_except,
                               &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->ret2service, buffer);
	safe_unpackstr_xmalloc(&build_ptr->schedauth, &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->schedport, buffer);
	safe_unpack16(&build_ptr->schedrootfltr, buffer);
	safe_unpackstr_xmalloc(&build_ptr->schedtype, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->select_type, &uint16_tmp, buffer);
	safe_unpack32(&build_ptr->slurm_user_id, buffer);
	safe_unpackstr_xmalloc(&build_ptr->slurm_user_name,
			       &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->slurmctld_debug, buffer);
	safe_unpackstr_xmalloc(&build_ptr->slurmctld_logfile,
			       &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->slurmctld_pidfile,
			       &uint16_tmp, buffer);
	safe_unpack32(&build_ptr->slurmctld_port, buffer);
	safe_unpack16(&build_ptr->slurmctld_timeout, buffer);
	safe_unpack16(&build_ptr->slurmd_debug, buffer);
	safe_unpackstr_xmalloc(&build_ptr->slurmd_logfile, &uint16_tmp,
			       buffer);
	safe_unpackstr_xmalloc(&build_ptr->slurmd_pidfile, &uint16_tmp,
			       buffer);
#ifndef MULTIPLE_SLURMD
	safe_unpack32(&build_ptr->slurmd_port, buffer);
#endif
	safe_unpackstr_xmalloc(&build_ptr->slurmd_spooldir, &uint16_tmp,
			       buffer);
	safe_unpack16(&build_ptr->slurmd_timeout, buffer);
	safe_unpackstr_xmalloc(&build_ptr->slurm_conf, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->state_save_location,
			       &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->switch_type, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->task_epilog, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->task_prolog, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->task_plugin, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->tmp_fs, &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->wait_time, buffer);
	safe_unpackstr_xmalloc(&build_ptr->job_credential_private_key,
			       &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->
			       job_credential_public_certificate,
			       &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->srun_prolog, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->srun_epilog, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&build_ptr->node_prefix, &uint16_tmp, buffer);
	safe_unpack16(&build_ptr->tree_width, buffer);

	return SLURM_SUCCESS;

      unpack_error:
	xfree(build_ptr->authtype);
	xfree(build_ptr->backup_addr);
	xfree(build_ptr->backup_controller);
	xfree(build_ptr->checkpoint_type);
	xfree(build_ptr->control_addr);
	xfree(build_ptr->control_machine);
	xfree(build_ptr->epilog);
	xfree(build_ptr->job_acct_loc);
	xfree(build_ptr->job_acct_parameters);
	xfree(build_ptr->job_acct_type);
	xfree(build_ptr->job_comp_loc);
	xfree(build_ptr->job_comp_type);
	xfree(build_ptr->job_credential_private_key);
	xfree(build_ptr->job_credential_public_certificate);
	xfree(build_ptr->plugindir);
	xfree(build_ptr->proctrack_type);
	xfree(build_ptr->prolog);
	xfree(build_ptr->propagate_rlimits);
	xfree(build_ptr->propagate_rlimits_except);
	xfree(build_ptr->schedauth);
	xfree(build_ptr->schedtype);
	xfree(build_ptr->select_type);
	xfree(build_ptr->slurm_conf);
	xfree(build_ptr->slurm_user_name);
	xfree(build_ptr->slurmctld_logfile);
	xfree(build_ptr->slurmctld_pidfile);
	xfree(build_ptr->slurmd_logfile);
	xfree(build_ptr->slurmd_pidfile);
	xfree(build_ptr->slurmd_spooldir);
	xfree(build_ptr->state_save_location);
	xfree(build_ptr->switch_type);
	xfree(build_ptr->task_epilog);
	xfree(build_ptr->task_prolog);
	xfree(build_ptr->task_plugin);
	xfree(build_ptr->tmp_fs);
	xfree(build_ptr->srun_prolog);
	xfree(build_ptr->srun_epilog);
	xfree(build_ptr->node_prefix);
	xfree(build_ptr);
	*build_buffer_ptr = NULL;
	return SLURM_ERROR;
}

/* _pack_job_desc_msg
 * packs a job_desc struct 
 * IN job_desc_ptr - pointer to the job descriptor to pack
 * IN/OUT buffer - destination of the pack, contains pointers that are 
 *			automatically updated
 */
static void
_pack_job_desc_msg(job_desc_msg_t * job_desc_ptr, Buf buffer)
{
	select_jobinfo_t jobinfo;

	/* load the data values */
	pack16((uint16_t)job_desc_ptr->contiguous, buffer);
        pack16((uint16_t)job_desc_ptr->exclusive, buffer);
	pack16((uint16_t)job_desc_ptr->kill_on_node_fail, buffer);
	packstr(job_desc_ptr->features, buffer);
	pack32((uint32_t)job_desc_ptr->job_id, buffer);
	packstr(job_desc_ptr->name, buffer);

	packstr(job_desc_ptr->alloc_node, buffer);
	pack32((uint32_t)job_desc_ptr->alloc_sid, buffer);
	pack32((uint32_t)job_desc_ptr->min_procs, buffer);
	pack32((uint32_t)job_desc_ptr->min_memory, buffer);
	pack32((uint32_t)job_desc_ptr->min_tmp_disk, buffer);

	packstr(job_desc_ptr->partition, buffer);
	pack32((uint32_t)job_desc_ptr->priority, buffer);
	pack32((uint32_t)job_desc_ptr->dependency, buffer);
	packstr(job_desc_ptr->account, buffer);
	pack16((uint16_t)job_desc_ptr->nice, buffer);

	packstr(job_desc_ptr->req_nodes, buffer);
	packstr(job_desc_ptr->exc_nodes, buffer);
	packstr_array(job_desc_ptr->environment, job_desc_ptr->env_size,
		      buffer);
	packstr(job_desc_ptr->script, buffer);
	packstr_array(job_desc_ptr->argv, job_desc_ptr->argc, buffer);

	packstr(job_desc_ptr->err, buffer);
	packstr(job_desc_ptr->in, buffer);
	packstr(job_desc_ptr->out, buffer);
	packstr(job_desc_ptr->work_dir, buffer);

	pack16((uint16_t)job_desc_ptr->immediate, buffer);
	pack16((uint16_t)job_desc_ptr->shared, buffer);
	pack16((uint16_t)job_desc_ptr->cpus_per_task, buffer);
	pack16((uint16_t)job_desc_ptr->task_dist, buffer);
	pack32((uint32_t)job_desc_ptr->time_limit, buffer);

	pack32((uint32_t)job_desc_ptr->num_procs, buffer);
	pack32((uint32_t)job_desc_ptr->min_nodes, buffer);
	pack32((uint32_t)job_desc_ptr->max_nodes, buffer);
	pack32((uint32_t)job_desc_ptr->num_tasks, buffer);
	pack32((uint32_t)job_desc_ptr->user_id, buffer);
	pack32((uint32_t)job_desc_ptr->group_id, buffer);

	pack16((uint16_t)job_desc_ptr->port, buffer);
	packstr(job_desc_ptr->host, buffer);
	packstr(job_desc_ptr->network, buffer);
	pack_time(job_desc_ptr->begin_time, buffer);

	pack16((uint16_t)job_desc_ptr->mail_type, buffer);
	packstr(job_desc_ptr->mail_user, buffer);
	select_g_pack_jobinfo(job_desc_ptr->select_jobinfo, buffer);
}

/* _unpack_job_desc_msg
 * unpacks a job_desc struct
 * OUT job_desc_buffer_ptr - place to put pointer to allocated job desc struct
 * IN/OUT buffer - source of the unpack, contains pointers that are 
 *			automatically updated
 */
static int
_unpack_job_desc_msg(job_desc_msg_t ** job_desc_buffer_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	job_desc_msg_t *job_desc_ptr;

	/* alloc memory for structure */
	job_desc_ptr = xmalloc(sizeof(job_desc_msg_t));
	*job_desc_buffer_ptr = job_desc_ptr;

	/* load the data values */
	safe_unpack16(&job_desc_ptr->contiguous, buffer);
        safe_unpack16(&job_desc_ptr->exclusive, buffer);
	safe_unpack16(&job_desc_ptr->kill_on_node_fail, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->features, &uint16_tmp, buffer);
	safe_unpack32(&job_desc_ptr->job_id, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->name, &uint16_tmp, buffer);

	safe_unpackstr_xmalloc(&job_desc_ptr->alloc_node, &uint16_tmp, buffer);
	safe_unpack32(&job_desc_ptr->alloc_sid, buffer);
	safe_unpack32(&job_desc_ptr->min_procs, buffer);
	safe_unpack32(&job_desc_ptr->min_memory, buffer);
	safe_unpack32(&job_desc_ptr->min_tmp_disk, buffer);

	safe_unpackstr_xmalloc(&job_desc_ptr->partition, &uint16_tmp, buffer);
	safe_unpack32(&job_desc_ptr->priority, buffer);
	safe_unpack32(&job_desc_ptr->dependency, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->account, &uint16_tmp, buffer);
	safe_unpack16(&job_desc_ptr->nice, buffer);

	safe_unpackstr_xmalloc(&job_desc_ptr->req_nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->exc_nodes, &uint16_tmp, buffer);
	safe_unpackstr_array(&job_desc_ptr->environment,
			     &job_desc_ptr->env_size, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->script, &uint16_tmp, buffer);
	safe_unpackstr_array(&job_desc_ptr->argv, &job_desc_ptr->argc, buffer);

	safe_unpackstr_xmalloc(&job_desc_ptr->err, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->in, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->out, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->work_dir, &uint16_tmp, buffer);

	safe_unpack16(&job_desc_ptr->immediate, buffer);
	safe_unpack16(&job_desc_ptr->shared, buffer);
	safe_unpack16(&job_desc_ptr->cpus_per_task, buffer);
	safe_unpack16(&job_desc_ptr->task_dist, buffer);
	safe_unpack32(&job_desc_ptr->time_limit, buffer);

	safe_unpack32(&job_desc_ptr->num_procs, buffer);
	safe_unpack32(&job_desc_ptr->min_nodes, buffer);
	safe_unpack32(&job_desc_ptr->max_nodes, buffer);
	safe_unpack32(&job_desc_ptr->num_tasks, buffer);
	safe_unpack32(&job_desc_ptr->user_id, buffer);
	safe_unpack32(&job_desc_ptr->group_id, buffer);

	safe_unpack16(&job_desc_ptr->port, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->host, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->network, &uint16_tmp, buffer);
	safe_unpack_time(&job_desc_ptr->begin_time, buffer);

	safe_unpack16(&job_desc_ptr->mail_type, buffer);
	safe_unpackstr_xmalloc(&job_desc_ptr->mail_user, &uint16_tmp, buffer);

	if (select_g_alloc_jobinfo (&job_desc_ptr->select_jobinfo)
	||  select_g_unpack_jobinfo(job_desc_ptr->select_jobinfo, buffer))
		goto unpack_error;
	
	return SLURM_SUCCESS;

      unpack_error:
	select_g_free_jobinfo(&job_desc_ptr->select_jobinfo);
	xfree(job_desc_ptr->features);
	xfree(job_desc_ptr->name);
	xfree(job_desc_ptr->partition);
	xfree(job_desc_ptr->req_nodes);
	xfree(job_desc_ptr->environment);
	xfree(job_desc_ptr->script);
	xfree(job_desc_ptr->argv);
	xfree(job_desc_ptr->err);
	xfree(job_desc_ptr->in);
	xfree(job_desc_ptr->out);
	xfree(job_desc_ptr->work_dir);
	xfree(job_desc_ptr->host);
	xfree(job_desc_ptr->network);
	xfree(job_desc_ptr->mail_user);
	xfree(job_desc_ptr);
	*job_desc_buffer_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_old_job_desc_msg(old_job_alloc_msg_t * job_desc_ptr, Buf buffer)
{
	/* load the data values */
	pack32((uint32_t)job_desc_ptr->job_id, buffer);
}

static int
_unpack_old_job_desc_msg(old_job_alloc_msg_t **
			 job_desc_buffer_ptr, Buf buffer)
{
	old_job_alloc_msg_t *job_desc_ptr;

	/* alloc memory for structure */
	assert(job_desc_buffer_ptr != NULL);
	job_desc_ptr = xmalloc(sizeof(old_job_alloc_msg_t));
	*job_desc_buffer_ptr = job_desc_ptr;

	/* load the data values */
	safe_unpack32(&job_desc_ptr->job_id, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(job_desc_ptr);
	*job_desc_buffer_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_last_update_msg(last_update_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);
	pack_time(msg->last_update, buffer);
}

static int
_unpack_last_update_msg(last_update_msg_t ** msg, Buf buffer)
{
	last_update_msg_t *last_update_msg;

	xassert(msg != NULL);
	last_update_msg = xmalloc(sizeof(last_update_msg_t));
	*msg = last_update_msg;

	safe_unpack_time(&last_update_msg->last_update, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(last_update_msg);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_return_code_msg(return_code_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);
	pack32((uint32_t)msg->return_code, buffer);
}

static int
_unpack_return_code_msg(return_code_msg_t ** msg, Buf buffer)
{
	return_code_msg_t *return_code_msg;

	xassert(msg != NULL);
	return_code_msg = xmalloc(sizeof(return_code_msg_t));
	*msg = return_code_msg;

	safe_unpack32(&return_code_msg->return_code, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(return_code_msg);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_reattach_tasks_request_msg(reattach_tasks_request_msg_t * msg,
				 Buf buffer)
{
	xassert(msg != NULL);
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_step_id, buffer);
	pack32((uint32_t)msg->srun_node_id, buffer);
	pack16((uint16_t)msg->resp_port, buffer);
	pack16((uint16_t)msg->io_port, buffer);
	packstr(msg->ofname, buffer);
	packstr(msg->efname, buffer);
	packstr(msg->ifname, buffer);
	slurm_cred_pack(msg->cred, buffer);
}

static int
_unpack_reattach_tasks_request_msg(reattach_tasks_request_msg_t ** msg_ptr,
				   Buf buffer)
{
	uint16_t uint16_tmp;
	reattach_tasks_request_msg_t *msg;

	xassert(msg_ptr != NULL);
	msg = xmalloc(sizeof(*msg));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_step_id, buffer);
	safe_unpack32(&msg->srun_node_id, buffer);
	safe_unpack16(&msg->resp_port, buffer);
	safe_unpack16(&msg->io_port, buffer);
	safe_unpackstr_xmalloc(&msg->ofname, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&msg->efname, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&msg->ifname, &uint16_tmp, buffer);

	if (!(msg->cred = slurm_cred_unpack(buffer)))
		goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	slurm_free_reattach_tasks_request_msg(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_reattach_tasks_response_msg(reattach_tasks_response_msg_t * msg,
				  Buf buffer)
{
	xassert(msg != NULL);
	packstr(msg->node_name,   buffer);
	packstr(msg->executable_name, buffer);
	pack32((uint32_t)msg->return_code,  buffer);
	pack32((uint32_t)msg->srun_node_id, buffer);
	pack32((uint32_t)msg->ntasks,       buffer);
	pack32_array(msg->gtids,      msg->ntasks, buffer);
	pack32_array(msg->local_pids, msg->ntasks, buffer);
}

static int
_unpack_reattach_tasks_response_msg(reattach_tasks_response_msg_t ** msg_ptr,
				    Buf buffer)
{
	uint32_t ntasks;
	uint16_t uint16_tmp;
	reattach_tasks_response_msg_t *msg = xmalloc(sizeof(*msg));

	xassert(msg_ptr != NULL);
	*msg_ptr = msg;

	safe_unpackstr_xmalloc(&msg->node_name, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&msg->executable_name, &uint16_tmp, buffer);
	safe_unpack32(&msg->return_code,  buffer);
	safe_unpack32(&msg->srun_node_id, buffer);
	safe_unpack32(&msg->ntasks,       buffer);
	safe_unpack32_array(&msg->gtids,      &ntasks, buffer);
	safe_unpack32_array(&msg->local_pids, &ntasks, buffer);
	if (msg->ntasks != ntasks)
		goto unpack_error;
	return SLURM_SUCCESS;

      unpack_error:
	slurm_free_reattach_tasks_response_msg(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}


static void
_pack_task_exit_msg(task_exit_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);
	pack32((uint32_t)msg->return_code, buffer);
	pack32((uint32_t)msg->num_tasks, buffer);
	pack32_array(msg->task_id_list,
		     msg->num_tasks, buffer);
}

static int
_unpack_task_exit_msg(task_exit_msg_t ** msg_ptr, Buf buffer)
{
	task_exit_msg_t *msg;
	uint32_t uint32_tmp;

	xassert(msg_ptr != NULL);
	msg = xmalloc(sizeof(task_exit_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->return_code, buffer);
	safe_unpack32(&msg->num_tasks, buffer);
	safe_unpack32_array(&msg->task_id_list, &uint32_tmp, buffer);
	if (msg->num_tasks != uint32_tmp)
		goto unpack_error;
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}


static void
_pack_launch_tasks_response_msg(launch_tasks_response_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);
	pack32((uint32_t)msg->return_code, buffer);
	packstr(msg->node_name, buffer);
	pack32((uint32_t)msg->srun_node_id, buffer);
	pack32((uint32_t)msg->count_of_pids, buffer);
	pack32_array(msg->local_pids,
		     msg->count_of_pids, buffer);
}

static int
_unpack_launch_tasks_response_msg(launch_tasks_response_msg_t **
				  msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	uint32_t uint32_tmp;
	launch_tasks_response_msg_t *msg;

	xassert(msg_ptr != NULL);
	msg = xmalloc(sizeof(launch_tasks_response_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->return_code, buffer);
	safe_unpackstr_xmalloc(&msg->node_name, &uint16_tmp, buffer);
	safe_unpack32(&msg->srun_node_id, buffer);
	safe_unpack32(&msg->count_of_pids, buffer);
	safe_unpack32_array(&msg->local_pids, &uint32_tmp, buffer);
	if (msg->count_of_pids != uint32_tmp)
		goto unpack_error;
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg->node_name);
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_launch_tasks_request_msg(launch_tasks_request_msg_t * msg, Buf buffer)
{
	int i=0;
	xassert(msg != NULL);
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_step_id, buffer);
	pack32((uint32_t)msg->nnodes, buffer);
	pack32((uint32_t)msg->nprocs, buffer);
	pack32((uint32_t)msg->uid, buffer);
	pack32((uint32_t)msg->gid, buffer);
	pack32((uint32_t)msg->srun_node_id, buffer);
	slurm_cred_pack(msg->cred, buffer);
	for(i=0; i<msg->nnodes; i++) {
		pack32((uint32_t)msg->tasks_to_launch[i], buffer);
		pack32((uint32_t)msg->cpus_allocated[i], buffer);
		pack16((uint16_t)msg->resp_port[i], buffer);
		pack16((uint16_t)msg->io_port[i], buffer);
		pack32_array(msg->global_task_ids[i], 
			     msg->tasks_to_launch[i], 
			     buffer);	
	}
	slurm_pack_slurm_addr(&msg->orig_addr, buffer);
	packstr_array(msg->env, msg->envc, buffer);
	packstr(msg->cwd, buffer);
	pack32((uint32_t)msg->cpu_bind_type, buffer);
	packstr(msg->cpu_bind, buffer);
	pack32((uint32_t)msg->mem_bind_type, buffer);
	packstr(msg->mem_bind, buffer);
	packstr_array(msg->argv, msg->argc, buffer);
	pack16((uint16_t)msg->task_flags, buffer);
	packstr(msg->ofname, buffer);
	packstr(msg->efname, buffer);
	packstr(msg->ifname, buffer);
	pack8(msg->buffered_stdio, buffer);
	packstr(msg->task_prolog, buffer);
	packstr(msg->task_epilog, buffer);
	pack32((uint32_t)msg->slurmd_debug, buffer);
	switch_pack_jobinfo(msg->switch_job, buffer);
}

static int
_unpack_launch_tasks_request_msg(launch_tasks_request_msg_t **
				 msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	uint32_t uint32_tmp;
	launch_tasks_request_msg_t *msg;
	int i=0;

	xassert(msg_ptr != NULL);
	msg = xmalloc(sizeof(launch_tasks_request_msg_t));
	*msg_ptr = msg;
 
	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_step_id, buffer);
	safe_unpack32(&msg->nnodes, buffer);
	safe_unpack32(&msg->nprocs, buffer);
	safe_unpack32(&msg->uid, buffer);
	safe_unpack32(&msg->gid, buffer);
	safe_unpack32(&msg->srun_node_id, buffer);
	if (!(msg->cred = slurm_cred_unpack(buffer)))
		goto unpack_error;
	msg->tasks_to_launch = xmalloc(sizeof(uint32_t) * msg->nnodes);
	msg->cpus_allocated = xmalloc(sizeof(uint32_t) * msg->nnodes);
	msg->resp_port = xmalloc(sizeof(uint16_t) * msg->nnodes);
	msg->io_port = xmalloc(sizeof(uint16_t) * msg->nnodes);
	msg->global_task_ids = xmalloc(sizeof(uint32_t *) * msg->nnodes);
	for(i=0; i<msg->nnodes; i++) {
		safe_unpack32(&msg->tasks_to_launch[i], buffer);
		safe_unpack32(&msg->cpus_allocated[i], buffer);
		safe_unpack16(&msg->resp_port[i], buffer);
		safe_unpack16(&msg->io_port[i], buffer);
		safe_unpack32_array(&msg->global_task_ids[i], 
				    &uint32_tmp, 
				    buffer);	
		if (msg->tasks_to_launch[i] != uint32_tmp)
			goto unpack_error;

	}
	slurm_unpack_slurm_addr_no_alloc(&msg->orig_addr, buffer);
	safe_unpackstr_array(&msg->env, &msg->envc, buffer);
	safe_unpackstr_xmalloc(&msg->cwd, &uint16_tmp, buffer);
	safe_unpack32(&msg->cpu_bind_type, buffer);
	safe_unpackstr_xmalloc(&msg->cpu_bind, &uint16_tmp, buffer);
	safe_unpack32(&msg->mem_bind_type, buffer);
	safe_unpackstr_xmalloc(&msg->mem_bind, &uint16_tmp, buffer);
	safe_unpackstr_array(&msg->argv, &msg->argc, buffer);
	safe_unpack16(&msg->task_flags, buffer);
	safe_unpackstr_xmalloc(&msg->ofname, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&msg->efname, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&msg->ifname, &uint16_tmp, buffer);
	safe_unpack8(&msg->buffered_stdio, buffer);
	safe_unpackstr_xmalloc(&msg->task_prolog, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&msg->task_epilog, &uint16_tmp, buffer);
	safe_unpack32(&msg->slurmd_debug, buffer);
	
	switch_alloc_jobinfo(&msg->switch_job);
	if (switch_unpack_jobinfo(msg->switch_job, buffer) < 0) {
		error("switch_unpack_jobinfo: %m");
		switch_free_jobinfo(msg->switch_job);
		goto unpack_error;
	}
	return SLURM_SUCCESS;

      unpack_error:
	slurm_free_launch_tasks_request_msg(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_spawn_task_request_msg(spawn_task_request_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_step_id, buffer);
	pack32((uint32_t)msg->nnodes, buffer);
	pack32((uint32_t)msg->nprocs, buffer);
	pack32((uint32_t)msg->uid, buffer);
        pack32((uint32_t)msg->gid, buffer);
	pack32((uint32_t)msg->srun_node_id, buffer);
	slurm_cred_pack(msg->cred, buffer);
	packstr_array(msg->env, msg->envc, buffer);
	packstr(msg->cwd, buffer);
	packstr_array(msg->argv, msg->argc, buffer);
	pack16((uint16_t)msg->io_port, buffer);
	pack16((uint16_t)msg->task_flags, buffer);
	pack16((uint16_t)msg->cpus_allocated, buffer);
	pack32((uint32_t)msg->slurmd_debug, buffer);
	pack32((uint32_t)msg->global_task_id, buffer);
	switch_pack_jobinfo(msg->switch_job, buffer);
}

static int
_unpack_spawn_task_request_msg(spawn_task_request_msg_t **
				 msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	spawn_task_request_msg_t *msg;

	xassert(msg_ptr != NULL);
	msg = xmalloc(sizeof(spawn_task_request_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_step_id, buffer);
	safe_unpack32(&msg->nnodes, buffer);
	safe_unpack32(&msg->nprocs, buffer);
	safe_unpack32(&msg->uid, buffer);
        safe_unpack32(&msg->gid, buffer);
	safe_unpack32(&msg->srun_node_id, buffer);
	if (!(msg->cred = slurm_cred_unpack(buffer)))
		goto unpack_error;
	safe_unpackstr_array(&msg->env, &msg->envc, buffer);
	safe_unpackstr_xmalloc(&msg->cwd, &uint16_tmp, buffer);
	safe_unpackstr_array(&msg->argv, &msg->argc, buffer);
	safe_unpack16(&msg->io_port, buffer);
	safe_unpack16(&msg->task_flags, buffer);
	safe_unpack16(&msg->cpus_allocated, buffer);
	safe_unpack32(&msg->slurmd_debug, buffer);
	safe_unpack32(&msg->global_task_id, buffer);

	switch_alloc_jobinfo(&msg->switch_job);
	if (switch_unpack_jobinfo(msg->switch_job, buffer) < 0) {
		error("switch_unpack_jobinfo: %m");
		switch_free_jobinfo(msg->switch_job);
		goto unpack_error;
	}

	return SLURM_SUCCESS;

      unpack_error:
	slurm_free_spawn_task_request_msg(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_cancel_tasks_msg(kill_tasks_msg_t * msg, Buf buffer)
{
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_step_id, buffer);
	pack32((uint32_t)msg->signal, buffer);
}

static int
_unpack_cancel_tasks_msg(kill_tasks_msg_t ** msg_ptr, Buf buffer)
{
	kill_tasks_msg_t *msg;

	msg = xmalloc(sizeof(kill_tasks_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_step_id, buffer);
	safe_unpack32(&msg->signal, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_shutdown_msg(shutdown_msg_t * msg, Buf buffer)
{
	pack16((uint16_t)msg->core, buffer);
}

static int
_unpack_shutdown_msg(shutdown_msg_t ** msg_ptr, Buf buffer)
{
	shutdown_msg_t *msg;

	msg = xmalloc(sizeof(shutdown_msg_t));
	*msg_ptr = msg;

	safe_unpack16(&msg->core, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

/* _pack_job_step_kill_msg
 * packs a slurm job step signal message
 * IN msg - pointer to the job step signal message
 * IN/OUT buffer - destination of the pack, contains pointers that are 
 *			automatically updated
 */
static void
_pack_job_step_kill_msg(job_step_kill_msg_t * msg, Buf buffer)
{
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_step_id, buffer);
	pack16((uint16_t)msg->signal, buffer);
	pack16((uint16_t)msg->batch_flag, buffer);
}

/* _unpack_job_step_kill_msg
 * unpacks a slurm job step signal message
 * OUT msg_ptr - pointer to the job step signal message buffer
 * IN/OUT buffer - source of the unpack, contains pointers that are 
 *			automatically updated
 */
static int
_unpack_job_step_kill_msg(job_step_kill_msg_t ** msg_ptr, Buf buffer)
{
	job_step_kill_msg_t *msg;

	msg = xmalloc(sizeof(job_step_kill_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_step_id, buffer);
	safe_unpack16(&msg->signal, buffer);
	safe_unpack16(&msg->batch_flag, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_complete_job_allocation_msg(
	complete_job_allocation_msg_t * msg, Buf buffer)
{
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_rc, buffer);
}

static int
_unpack_complete_job_allocation_msg(
	complete_job_allocation_msg_t ** msg_ptr, Buf buffer)
{
	complete_job_allocation_msg_t *msg;
	uint16_t uint16_tmp;

	msg = xmalloc(sizeof(complete_job_allocation_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_rc, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_complete_batch_script_msg(
	complete_batch_script_msg_t * msg, Buf buffer)
{
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_rc, buffer);
	pack32((uint32_t)msg->slurm_rc, buffer);
	packstr(msg->node_name, buffer);
}

static int
_unpack_complete_batch_script_msg(
	complete_batch_script_msg_t ** msg_ptr, Buf buffer)
{
	complete_batch_script_msg_t *msg;
	uint16_t uint16_tmp;

	msg = xmalloc(sizeof(complete_batch_script_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_rc, buffer);
	safe_unpack32(&msg->slurm_rc, buffer);
	safe_unpackstr_xmalloc(&msg->node_name, &uint16_tmp, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_step_complete_msg(step_complete_msg_t * msg, Buf buffer)
{
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->job_step_id, buffer);
	pack32((uint32_t)msg->range_first, buffer);
	pack32((uint32_t)msg->range_last, buffer);
	pack32((uint32_t)msg->step_rc, buffer);
}

static int
_unpack_step_complete_msg(step_complete_msg_t ** msg_ptr, Buf buffer)
{
	step_complete_msg_t *msg;

	msg = xmalloc(sizeof(step_complete_msg_t));
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id, buffer);
	safe_unpack32(&msg->job_step_id, buffer);
	safe_unpack32(&msg->range_first, buffer);
	safe_unpack32(&msg->range_last, buffer);
	safe_unpack32(&msg->step_rc, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void
_pack_job_info_request_msg(job_info_request_msg_t * msg, Buf buffer)
{
	pack_time(msg->last_update, buffer);
	pack16((uint16_t)msg->show_flags, buffer);
}

static int
_unpack_job_info_request_msg(job_info_request_msg_t** msg, 
		Buf buffer)
{
	job_info_request_msg_t*job_info;

	job_info = xmalloc(sizeof(job_step_info_request_msg_t));
	*msg = job_info;

	safe_unpack_time(&job_info->last_update, buffer);
	safe_unpack16(&job_info->show_flags, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(job_info);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_node_select_info_req_msg(node_info_select_request_msg_t *msg, Buf buffer)
{
	pack_time(msg->last_update, buffer);
}

static int
_unpack_node_select_info_req_msg(node_info_select_request_msg_t **msg, 
		Buf buffer)
{
	node_info_select_request_msg_t *node_sel_info;

	node_sel_info = xmalloc(sizeof(node_info_select_request_msg_t));
	*msg = node_sel_info;

	safe_unpack_time(&node_sel_info->last_update, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(node_sel_info);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_job_step_info_req_msg(job_step_info_request_msg_t * msg, Buf buffer)
{
	pack_time(msg->last_update, buffer);
	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->step_id, buffer);
	pack16((uint16_t)msg->show_flags, buffer);
}

static int
_unpack_job_step_info_req_msg(job_step_info_request_msg_t ** msg, Buf buffer)
{
	job_step_info_request_msg_t *job_step_info;

	job_step_info = xmalloc(sizeof(job_step_info_request_msg_t));
	*msg = job_step_info;

	safe_unpack_time(&job_step_info->last_update, buffer);
	safe_unpack32(&job_step_info->job_id, buffer);
	safe_unpack32(&job_step_info->step_id, buffer);
	safe_unpack16(&job_step_info->show_flags, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(job_step_info);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_node_info_request_msg(node_info_request_msg_t * msg, Buf buffer)
{
	pack_time(msg->last_update, buffer);
	pack16((uint16_t)msg->show_flags, buffer);
}

static int
_unpack_node_info_request_msg(node_info_request_msg_t ** msg, Buf buffer)
{
	node_info_request_msg_t* node_info;

	node_info = xmalloc(sizeof(node_info_request_msg_t));
	*msg = node_info;

	safe_unpack_time(&node_info->last_update, buffer);
	safe_unpack16(&node_info->show_flags, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(node_info);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_part_info_request_msg(part_info_request_msg_t * msg, Buf buffer)
{
	pack_time(msg->last_update, buffer);
	pack16((uint16_t)msg->show_flags, buffer);
}

static int
_unpack_part_info_request_msg(part_info_request_msg_t ** msg, Buf buffer)
{
	part_info_request_msg_t* part_info;

	part_info = xmalloc(sizeof(part_info_request_msg_t));
	*msg = part_info;

	safe_unpack_time(&part_info->last_update, buffer);
	safe_unpack16(&part_info->show_flags, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(part_info);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_slurm_addr_array(slurm_addr * slurm_address,
		       uint16_t size_val, Buf buffer)
{
	int i = 0;
	uint16_t nl = htons(size_val);
	pack16((uint16_t)nl, buffer);

	for (i = 0; i < size_val; i++) {
		slurm_pack_slurm_addr(slurm_address + i, buffer);
	}

}

static int
_unpack_slurm_addr_array(slurm_addr ** slurm_address,
			 uint16_t * size_val, Buf buffer)
{
	int i = 0;
	uint16_t nl;

	*slurm_address = NULL;
	safe_unpack16(&nl, buffer);
	*size_val = ntohs(nl);
	*slurm_address = xmalloc((*size_val) * sizeof(slurm_addr));

	for (i = 0; i < *size_val; i++) {
		if (slurm_unpack_slurm_addr_no_alloc((*slurm_address) + i,
						     buffer))
			goto unpack_error;

	}
	return SLURM_SUCCESS;

      unpack_error:
	xfree(*slurm_address);
	*slurm_address = NULL;
	return SLURM_ERROR;
}

static void
_pack_ret_list(List ret_list,
	       uint16_t size_val, Buf buffer)
{
	uint16_t i = 0;
	ListIterator itr;
	ListIterator itr_data;
	ret_types_t *ret_type = NULL;
	ret_data_info_t *ret_data_info = NULL;
	slurm_msg_t msg;

	
	itr = list_iterator_create(ret_list);		
	while((ret_type = list_next(itr)) != NULL) {
		pack32((uint32_t)ret_type->msg_rc, buffer);
		pack32((uint32_t)ret_type->err, buffer);
		pack32((uint32_t)ret_type->type, buffer);
		
		msg.msg_type = ret_type->type;

		i = list_count(ret_type->ret_data_list);
		pack16((uint16_t)i, buffer);
		itr_data = list_iterator_create(ret_type->ret_data_list);
		while((ret_data_info = list_next(itr_data)) != NULL) {
			packstr(ret_data_info->node_name, buffer);
			pack32((uint32_t)ret_data_info->nodeid, buffer);
			msg.data = ret_data_info->data;
			pack_msg(&msg, buffer);
		} 
		list_iterator_destroy(itr_data);		
	}
	list_iterator_destroy(itr);
}

static int
_unpack_ret_list(List *ret_list,
		 uint16_t size_val, Buf buffer)
{
	int i = 0, j = 0;
	uint16_t nl, uint16_tmp;
	ret_types_t *ret_type = NULL;
	ret_data_info_t *ret_data_info = NULL;
	slurm_msg_t msg;
	*ret_list = list_create(destroy_ret_types);
	
	for (i=0; i<size_val; i++) {
		ret_type = xmalloc(sizeof(ret_types_t));
		list_push(*ret_list, ret_type);
		safe_unpack32((uint32_t *)&ret_type->msg_rc, buffer);
		safe_unpack32((uint32_t *)&ret_type->err, buffer);
		safe_unpack32((uint32_t *)&ret_type->type, buffer);
		
		msg.msg_type = ret_type->type;
			
		safe_unpack16(&nl, buffer);
		ret_type->ret_data_list = list_create(destroy_data_info);
		for(j=0; j<nl; j++) {
			ret_data_info = xmalloc(sizeof(ret_data_info_t));
			list_push(ret_type->ret_data_list, ret_data_info);
			safe_unpackstr_xmalloc(&ret_data_info->node_name, 
					       &uint16_tmp, buffer);
			safe_unpack32((uint32_t *)&ret_data_info->nodeid, buffer);
			unpack_msg(&msg, buffer);
			ret_data_info->data = msg.data;
		}
	}
	return SLURM_SUCCESS;

      unpack_error:
	list_destroy(*ret_list);
	*ret_list = NULL;
	return SLURM_ERROR;
}

static void
_pack_batch_job_launch_msg(batch_job_launch_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
	pack32((uint32_t)msg->step_id, buffer);
	pack32((uint32_t)msg->uid, buffer);
	pack32((uint32_t)msg->gid, buffer);
	pack32((uint32_t)msg->nprocs, buffer);

	pack16((uint16_t)msg->num_cpu_groups, buffer);
	pack32_array(msg->cpus_per_node, msg->num_cpu_groups, buffer);
	pack32_array(msg->cpu_count_reps, msg->num_cpu_groups, buffer);

	packstr(msg->nodes, buffer);
	packstr(msg->script, buffer);
	packstr(msg->work_dir, buffer);

	packstr(msg->err, buffer);
	packstr(msg->in, buffer);
	packstr(msg->out, buffer);

	pack16((uint16_t)msg->argc, buffer);
	packstr_array(msg->argv, msg->argc, buffer);

	pack16((uint16_t)msg->envc, buffer);
	packstr_array(msg->environment, msg->envc, buffer);

	select_g_pack_jobinfo(msg->select_jobinfo, buffer);
}

static int
_unpack_batch_job_launch_msg(batch_job_launch_msg_t ** msg, Buf buffer)
{
	uint16_t uint16_tmp;
	uint32_t uint32_tmp;
	batch_job_launch_msg_t *launch_msg_ptr;

	xassert(msg != NULL);
	launch_msg_ptr = xmalloc(sizeof(batch_job_launch_msg_t));
	*msg = launch_msg_ptr;

	safe_unpack32(&launch_msg_ptr->job_id, buffer);
	safe_unpack32(&launch_msg_ptr->step_id, buffer);
	safe_unpack32(&launch_msg_ptr->uid, buffer);
	safe_unpack32(&launch_msg_ptr->gid, buffer);
	safe_unpack32(&launch_msg_ptr->nprocs, buffer);

	safe_unpack16(&launch_msg_ptr->num_cpu_groups, buffer);
	safe_unpack32_array((uint32_t **) &(launch_msg_ptr->cpus_per_node), 
			    &uint32_tmp,
			    buffer);
	if (launch_msg_ptr->num_cpu_groups != uint32_tmp)
		goto unpack_error;
	safe_unpack32_array((uint32_t **) &(launch_msg_ptr->cpu_count_reps), 
			    &uint32_tmp,
			    buffer);
	if (launch_msg_ptr->num_cpu_groups != uint32_tmp)
		goto unpack_error;
	
	safe_unpackstr_xmalloc(&launch_msg_ptr->nodes, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&launch_msg_ptr->script, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&launch_msg_ptr->work_dir, &uint16_tmp,
			       buffer);

	safe_unpackstr_xmalloc(&launch_msg_ptr->err, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&launch_msg_ptr->in, &uint16_tmp, buffer);
	safe_unpackstr_xmalloc(&launch_msg_ptr->out, &uint16_tmp, buffer);

	safe_unpack16(&launch_msg_ptr->argc, buffer);
	safe_unpackstr_array(&launch_msg_ptr->argv,
			     &launch_msg_ptr->argc, buffer);

	safe_unpack16(&launch_msg_ptr->envc, buffer);
	safe_unpackstr_array(&launch_msg_ptr->environment,
			     &launch_msg_ptr->envc, buffer);

	if (select_g_alloc_jobinfo (&launch_msg_ptr->select_jobinfo)
	||  select_g_unpack_jobinfo(launch_msg_ptr->select_jobinfo, buffer))
		goto unpack_error;

	return SLURM_SUCCESS;

      unpack_error:
	xfree(launch_msg_ptr->nodes);
	xfree(launch_msg_ptr->script);
	xfree(launch_msg_ptr->work_dir);
	xfree(launch_msg_ptr->err);
	xfree(launch_msg_ptr->in);
	xfree(launch_msg_ptr->out);
	xfree(launch_msg_ptr->argv);
	xfree(launch_msg_ptr->environment);
	xfree(launch_msg_ptr->cpus_per_node);
	xfree(launch_msg_ptr->cpu_count_reps);
	select_g_free_jobinfo(&launch_msg_ptr->select_jobinfo);
	xfree(launch_msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_job_id_request_msg(job_id_request_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_pid, buffer);
}

static int
_unpack_job_id_request_msg(job_id_request_msg_t ** msg, Buf buffer)
{
	job_id_request_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(job_id_request_msg_t));
	*msg = tmp_ptr;

	/* load the data values */
	safe_unpack32(&tmp_ptr->job_pid, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_job_id_response_msg(job_id_response_msg_t * msg, Buf buffer)
{
	xassert(msg != NULL);

	pack32((uint32_t)msg->job_id, buffer);
}

static int
_unpack_job_id_response_msg(job_id_response_msg_t ** msg, Buf buffer)
{
	job_id_response_msg_t *tmp_ptr;

	/* alloc memory for structure */
	xassert(msg != NULL);
	tmp_ptr = xmalloc(sizeof(job_id_response_msg_t));
	*msg = tmp_ptr;

	/* load the data values */
	safe_unpack32(&tmp_ptr->job_id, buffer);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(tmp_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

static void
_pack_srun_ping_msg(srun_ping_msg_t * msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack32((uint32_t)msg ->job_id  , buffer ) ;
	pack32((uint32_t)msg ->step_id , buffer ) ;
}

static int  
_unpack_srun_ping_msg(srun_ping_msg_t ** msg_ptr, Buf buffer)
{
	srun_ping_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (srun_ping_msg_t) ) ;
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id  , buffer ) ;
	safe_unpack32(&msg->step_id , buffer ) ;
	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree(msg);
	return SLURM_ERROR;
}

static void 
_pack_srun_node_fail_msg(srun_node_fail_msg_t * msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack32((uint32_t)msg ->job_id  , buffer ) ;
	pack32((uint32_t)msg ->step_id , buffer ) ;
	packstr(msg ->nodelist, buffer ) ;
}

static int 
_unpack_srun_node_fail_msg(srun_node_fail_msg_t ** msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	srun_node_fail_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (srun_node_fail_msg_t) ) ;
	*msg_ptr = msg;

	safe_unpack32(&msg->job_id  , buffer ) ;
	safe_unpack32(&msg->step_id , buffer ) ;
	safe_unpackstr_xmalloc ( & msg->nodelist, &uint16_tmp, buffer);

	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree( msg->nodelist );
	xfree( msg );
	return SLURM_ERROR;
}

static void
_pack_job_ready_msg(job_id_msg_t * msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack32((uint32_t)msg->job_id  , buffer ) ;
}

static int
_unpack_job_ready_msg(job_id_msg_t ** msg_ptr, Buf buffer)
{
	job_id_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (job_id_msg_t) );
	*msg_ptr = msg ;

	safe_unpack32(&msg->job_id  , buffer ) ;
	return SLURM_SUCCESS;

      unpack_error:
	*msg_ptr = NULL;
	xfree(msg);
	return SLURM_ERROR;
}

static void
_pack_srun_timeout_msg(srun_timeout_msg_t * msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack32((uint32_t)msg->job_id, buffer ) ;
	pack32((uint32_t)msg->step_id , buffer ) ;
	pack_time ( msg -> timeout, buffer );
}

static int  
_unpack_srun_timeout_msg(srun_timeout_msg_t ** msg_ptr, Buf buffer)
{
	srun_timeout_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (srun_timeout_msg_t) ) ;
	*msg_ptr = msg ;

	safe_unpack32(&msg->job_id, buffer ) ;
	safe_unpack32(&msg->step_id, buffer ) ;
	safe_unpack_time (&msg->timeout, buffer );
	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree(msg);
	return SLURM_ERROR;
}

static void _pack_suspend_msg(suspend_msg_t *msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack16((uint16_t)msg -> op, buffer ) ;
	pack32((uint32_t)msg->job_id,  buffer ) ;
}

static int  _unpack_suspend_msg(suspend_msg_t **msg_ptr, Buf buffer)
{
	suspend_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (suspend_msg_t) );
	*msg_ptr = msg ;

	safe_unpack16(&msg->op ,      buffer ) ;
	safe_unpack32(&msg->job_id  , buffer ) ;
	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree(msg);
	return SLURM_ERROR;
}

static void
_pack_checkpoint_msg(checkpoint_msg_t *msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack16((uint16_t)msg->op,      buffer ) ;
	pack16((uint16_t)msg->data,    buffer ) ;
	pack32((uint32_t)msg->job_id,  buffer ) ;
	pack32((uint32_t)msg->step_id, buffer ) ;
}

static int
_unpack_checkpoint_msg(checkpoint_msg_t **msg_ptr, Buf buffer)
{
	checkpoint_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (checkpoint_msg_t) ) ;
	*msg_ptr = msg ;

	safe_unpack16(&msg->op, buffer ) ;
	safe_unpack16(&msg->data, buffer ) ;
	safe_unpack32(&msg->job_id, buffer ) ;
	safe_unpack32(&msg->step_id, buffer ) ;
	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree(msg);
	return SLURM_ERROR;
}

static void
_pack_checkpoint_comp(checkpoint_comp_msg_t *msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack32((uint32_t)msg -> job_id,  buffer ) ;
	pack32((uint32_t)msg -> step_id, buffer ) ;
	pack32((uint32_t)msg -> error_code, buffer ) ;
	packstr ( msg -> error_msg, buffer ) ;
	pack_time ( msg -> begin_time, buffer ) ;
}

static int
_unpack_checkpoint_comp(checkpoint_comp_msg_t **msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	checkpoint_comp_msg_t * msg;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (checkpoint_comp_msg_t) );
	*msg_ptr = msg ;

	safe_unpack32(& msg -> job_id  , buffer ) ;
	safe_unpack32(& msg -> step_id , buffer ) ;
	safe_unpack32(& msg -> error_code , buffer ) ;
	safe_unpackstr_xmalloc ( & msg -> error_msg, & uint16_tmp , buffer ) ;
	safe_unpack_time ( & msg -> begin_time , buffer ) ;
	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree (msg->error_msg);
	xfree (msg);
	return SLURM_ERROR;
}

static void
_pack_checkpoint_resp_msg(checkpoint_resp_msg_t *msg, Buf buffer)
{
	xassert ( msg != NULL );

	pack_time ( msg -> event_time, buffer ) ;
	pack32((uint32_t)msg -> error_code,  buffer ) ;
	packstr ( msg -> error_msg, buffer ) ;
}

static int
_unpack_checkpoint_resp_msg(checkpoint_resp_msg_t **msg_ptr, Buf buffer)
{
	checkpoint_resp_msg_t * msg;
	uint16_t uint16_tmp;
	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (checkpoint_resp_msg_t) ) ;
	*msg_ptr = msg ;

	safe_unpack_time ( & msg -> event_time, buffer ) ;
	safe_unpack32(& msg -> error_code , buffer ) ;
	safe_unpackstr_xmalloc ( & msg -> error_msg, & uint16_tmp , buffer ) ;
	return SLURM_SUCCESS;

    unpack_error:
	*msg_ptr = NULL;
	xfree(msg);
	return SLURM_ERROR;
}

static void _pack_jobacct_data(jobacct_msg_t * msg , Buf buffer )
{
	xassert ( msg != NULL );
	/* debug("jobacct: packing message"); */

	pack16((uint16_t)msg->len, buffer );
	packmem(msg->data, msg->len, buffer ) ;
}

static int _unpack_jobacct_data(jobacct_msg_t ** msg_ptr , Buf buffer )
{
	uint16_t uint16_tmp;
	jobacct_msg_t *msg ;

	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (jobacct_msg_t) ) ;
	*msg_ptr = msg;

	safe_unpack16(&msg->len, buffer );
	safe_unpackmem_xmalloc(&msg->data, &uint16_tmp , buffer ) ;
	/* debug("jobacct: unpacked message"); */
	return SLURM_SUCCESS;

    unpack_error:
	xfree(msg->data);
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void _pack_file_bcast(file_bcast_msg_t * msg , Buf buffer )
{
	xassert ( msg != NULL );

	pack16 ( msg->block_no, buffer );
	pack16 ( msg->last_block, buffer );
	pack16 ( msg->force, buffer );
	pack16 ( msg->modes, buffer );

	pack32 ( msg->uid, buffer );
	pack32 ( msg->gid, buffer );
	pack32 ( msg->block_len, buffer );

	pack_time ( msg->atime, buffer );
	pack_time ( msg->mtime, buffer );

	packstr( msg->fname, buffer );
	packmem ( msg->data, msg->block_len, buffer ) ;
}

static int _unpack_file_bcast(file_bcast_msg_t ** msg_ptr , Buf buffer )
{
	uint16_t uint16_tmp;
	file_bcast_msg_t *msg ;

	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof (file_bcast_msg_t) ) ;
	*msg_ptr = msg;

	safe_unpack16 ( & msg->block_no, buffer );
	safe_unpack16 ( & msg->last_block, buffer );
	safe_unpack16 ( & msg->force, buffer );
	safe_unpack16 ( & msg->modes, buffer );

	safe_unpack32 ( & msg->uid, buffer );
	safe_unpack32 ( & msg->gid, buffer );
	safe_unpack32 ( & msg->block_len, buffer );

	safe_unpack_time ( & msg->atime, buffer );
	safe_unpack_time ( & msg->mtime, buffer );

	safe_unpackstr_xmalloc ( & msg->fname, &uint16_tmp, buffer );
	safe_unpackmem_xmalloc ( & msg->data, &uint16_tmp , buffer ) ;
	if ( uint16_tmp != msg->block_len )
		goto unpack_error;
	return SLURM_SUCCESS;

    unpack_error:
	xfree( msg -> data );
	xfree( msg -> fname );
	xfree( msg );
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void _pack_kvs_rec(struct kvs_comm *msg_ptr, Buf buffer)
{
	int i;
	xassert(msg_ptr != NULL);

	packstr(msg_ptr->kvs_name, buffer);
	pack16((uint16_t)msg_ptr->kvs_cnt, buffer);
	for (i=0; i<msg_ptr->kvs_cnt; i++) {
		packstr(msg_ptr->kvs_keys[i], buffer);
		packstr(msg_ptr->kvs_values[i], buffer);
	}
}
static int  _unpack_kvs_rec(struct kvs_comm **msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	int i;
	struct kvs_comm *msg;

	msg = xmalloc(sizeof(struct kvs_comm));
	*msg_ptr = msg;
	safe_unpackstr_xmalloc(&msg->kvs_name, &uint16_tmp, buffer);
	safe_unpack16(&msg->kvs_cnt, buffer);
	msg->kvs_keys   = xmalloc(sizeof(char *) * msg->kvs_cnt);
	msg->kvs_values = xmalloc(sizeof(char *) * msg->kvs_cnt);
	for (i=0; i<msg->kvs_cnt; i++) {
		safe_unpackstr_xmalloc(&msg->kvs_keys[i], 
				&uint16_tmp, buffer);
		safe_unpackstr_xmalloc(&msg->kvs_values[i], 
				&uint16_tmp, buffer);
	}
	return SLURM_SUCCESS;

unpack_error:
	return SLURM_ERROR;
}
static void _pack_kvs_data(struct kvs_comm_set *msg_ptr, Buf buffer)
{
	int i;
	xassert(msg_ptr != NULL);

	pack16((uint16_t)msg_ptr->task_id, buffer);
	pack16((uint16_t)msg_ptr->kvs_comm_recs, buffer);
	for (i=0; i<msg_ptr->kvs_comm_recs; i++) 
		_pack_kvs_rec(msg_ptr->kvs_comm_ptr[i], buffer);
}

static int  _unpack_kvs_data(struct kvs_comm_set **msg_ptr, Buf buffer)
{
	struct kvs_comm_set *msg;
	int i, j;

	msg = xmalloc(sizeof(struct kvs_comm_set));
	*msg_ptr = msg;
	safe_unpack16(&msg->task_id, buffer);
	safe_unpack16(&msg->kvs_comm_recs, buffer);
	msg->kvs_comm_ptr = xmalloc(sizeof(struct kvs_comm) * 
			msg->kvs_comm_recs);
	for (i=0; i<msg->kvs_comm_recs; i++) {
		if (_unpack_kvs_rec(&msg->kvs_comm_ptr[i], buffer))
			goto unpack_error;
	}
	return SLURM_SUCCESS;

unpack_error:
	for (i=0; i<msg->kvs_comm_recs; i++) {
		xfree(msg->kvs_comm_ptr[i]->kvs_name);
		for (j=0; j<msg->kvs_comm_ptr[i]->kvs_cnt; j++) {
			xfree(msg->kvs_comm_ptr[i]->kvs_keys[j]);
			xfree(msg->kvs_comm_ptr[i]->kvs_values[j]);
		}
	}
	xfree(msg->kvs_comm_ptr);
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

static void _pack_kvs_get(kvs_get_msg_t *msg_ptr, Buf buffer)
{
	pack16((uint16_t)msg_ptr->task_id, buffer);
	pack16((uint16_t)msg_ptr->size, buffer);
	pack16((uint16_t)msg_ptr->port, buffer);
	packstr(msg_ptr->hostname, buffer);
}

static int  _unpack_kvs_get(kvs_get_msg_t **msg_ptr, Buf buffer)
{
	uint16_t uint16_tmp;
	kvs_get_msg_t *msg;

	msg = xmalloc(sizeof(struct kvs_get_msg));
	*msg_ptr = msg;
	safe_unpack16(&msg->task_id, buffer);
	safe_unpack16(&msg->size, buffer);
	safe_unpack16(&msg->port, buffer);
	safe_unpackstr_xmalloc(&msg->hostname, &uint16_tmp, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}

/* template 
void pack_ ( * msg , Buf buffer )
{
	xassert ( msg != NULL );

	pack16( msg -> , buffer ) ;
	pack32((uint32_t)msg -> , buffer ) ;
	packstr ( msg -> , buffer ) ;
}

int unpack_ ( ** msg_ptr , Buf buffer )
{
	uint16_t uint16_tmp;
	* msg ;

	xassert ( msg_ptr != NULL );

	msg = xmalloc ( sizeof ( ) ) ;
	*msg_ptr = msg;

	safe_unpack16( & msg -> , buffer ) ;
	safe_unpack32(& msg -> , buffer ) ;
	safe_unpackstr_xmalloc ( & msg -> x, & uint16_tmp , buffer ) ;
	return SLURM_SUCCESS;

    unpack_error:
	xfree(msg -> x);
	xfree(msg);
	*msg_ptr = NULL;
	return SLURM_ERROR;
}
*/
