/*****************************************************************************\
 *  proc_msg.c - process incomming messages to slurmctld
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Moe Jette <jette@llnl.gov>, Kevin Tew <tew1@llnl.gov>, et. al.
 *  UCRL-CODE-2002-040.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif				/* WITH_PTHREADS */

#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <slurm/slurm_errno.h>

#include "src/common/daemonize.h"
#include "src/common/fd.h"
#include "src/common/hostlist.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/pack.h"
#include "src/common/read_config.h"
#include "src/common/slurm_auth.h"
#include "src/common/slurm_cred.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/switch.h"
#include "src/common/xstring.h"

#include "src/slurmctld/locks.h"
#include "src/slurmctld/proc_req.h"
#include "src/slurmctld/read_config.h"
#include "src/slurmctld/slurmctld.h"

#define BUF_SIZE	  1024	/* Temporary buffer size */

static void         _fill_ctld_conf(slurm_ctl_conf_t * build_ptr);
static inline bool 	_is_super_user(uid_t uid);
static void         _kill_job_on_msg_fail(uint32_t job_id);
static int          _make_step_cred(struct step_record *step_rec, 
				    slurm_cred_t *slurm_cred);
inline static void  _slurm_rpc_allocate_resources(slurm_msg_t * msg);
inline static void  _slurm_rpc_allocate_and_run(slurm_msg_t * msg);
inline static void  _slurm_rpc_dump_conf(slurm_msg_t * msg);
inline static void  _slurm_rpc_dump_jobs(slurm_msg_t * msg);
inline static void  _slurm_rpc_dump_nodes(slurm_msg_t * msg);
inline static void  _slurm_rpc_dump_partitions(slurm_msg_t * msg);
inline static void  _slurm_rpc_epilog_complete(slurm_msg_t * msg);
inline static void  _slurm_rpc_job_step_kill(slurm_msg_t * msg);
inline static void  _slurm_rpc_job_step_complete(slurm_msg_t * msg);
inline static void  _slurm_rpc_job_step_create(slurm_msg_t * msg);
inline static void  _slurm_rpc_job_step_get_info(slurm_msg_t * msg);
inline static void  _slurm_rpc_job_will_run(slurm_msg_t * msg);
inline static void  _slurm_rpc_node_registration(slurm_msg_t * msg);
inline static void  _slurm_rpc_old_job_alloc(slurm_msg_t * msg);
inline static void  _slurm_rpc_ping(slurm_msg_t * msg);
inline static void  _slurm_rpc_reconfigure_controller(slurm_msg_t * msg);
inline static void  _slurm_rpc_shutdown_controller(slurm_msg_t * msg);
inline static void  _slurm_rpc_shutdown_controller_immediate(slurm_msg_t *
							     msg);
inline static void  _slurm_rpc_submit_batch_job(slurm_msg_t * msg);
inline static void  _slurm_rpc_update_job(slurm_msg_t * msg);
inline static void  _slurm_rpc_update_node(slurm_msg_t * msg);
inline static void  _slurm_rpc_update_partition(slurm_msg_t * msg);
inline static void  _slurm_rpc_delete_partition(slurm_msg_t * msg);
inline static void  _update_cred_key(void);

/*
 * diff_tv_str - build a string showing the time difference between two times
 * IN tv1 - start of event
 * IN tv2 - end of event
 * OUT tv_str - place to put delta time in format "usec=%ld"
 * IN len_tv_str - size of tv_str in bytes
 */
inline void diff_tv_str(struct timeval *tv1,struct timeval *tv2, 
		char *tv_str, int len_tv_str)
{
	long delta_t;
	delta_t  = (tv2->tv_sec  - tv1->tv_sec) * 1000000;
	delta_t +=  tv2->tv_usec - tv1->tv_usec;
	snprintf(tv_str, len_tv_str, "usec=%ld", delta_t);
	if (delta_t > 1000000)
		error("Note very large time: %s",tv_str); 
}

/*
 * slurmctld_req  - Process an individual RPC request
 * IN/OUT msg - the request message, data associated with the message is freed
 */
void slurmctld_req (slurm_msg_t * msg)
{

	switch (msg->msg_type) {
	case REQUEST_RESOURCE_ALLOCATION:
		_slurm_rpc_allocate_resources(msg);
		slurm_free_job_desc_msg(msg->data);
		break;
	case REQUEST_ALLOCATION_AND_RUN_JOB_STEP:
		_slurm_rpc_allocate_and_run(msg);
		slurm_free_job_desc_msg(msg->data);
		break;
	case REQUEST_BUILD_INFO:
		_slurm_rpc_dump_conf(msg);
		slurm_free_last_update_msg(msg->data);
		break;
	case REQUEST_JOB_INFO:
		_slurm_rpc_dump_jobs(msg);
		slurm_free_job_info_request_msg(msg->data);
		break;
	case REQUEST_NODE_INFO:
		_slurm_rpc_dump_nodes(msg);
		slurm_free_last_update_msg(msg->data);
		break;
	case REQUEST_PARTITION_INFO:
		_slurm_rpc_dump_partitions(msg);
		slurm_free_last_update_msg(msg->data);
		break;
	case MESSAGE_EPILOG_COMPLETE:
		_slurm_rpc_epilog_complete(msg);
		slurm_free_epilog_complete_msg(msg->data);
		break;
	case REQUEST_CANCEL_JOB_STEP:
		_slurm_rpc_job_step_kill(msg);
		slurm_free_job_step_kill_msg(msg->data);
		break;
	case REQUEST_COMPLETE_JOB_STEP:
		_slurm_rpc_job_step_complete(msg);
		slurm_free_job_complete_msg(msg->data);
		break;
	case REQUEST_JOB_STEP_CREATE:
		_slurm_rpc_job_step_create(msg);
		slurm_free_job_step_create_request_msg(msg->data);
		break;
	case REQUEST_JOB_STEP_INFO:
		_slurm_rpc_job_step_get_info(msg);
		slurm_free_job_step_info_request_msg(msg->data);
		break;
	case REQUEST_JOB_WILL_RUN:
		_slurm_rpc_job_will_run(msg->data);
		slurm_free_job_desc_msg(msg->data);
		break;
	case MESSAGE_NODE_REGISTRATION_STATUS:
		_slurm_rpc_node_registration(msg);
		slurm_free_node_registration_status_msg(msg->data);
		break;
	case REQUEST_OLD_JOB_RESOURCE_ALLOCATION:
		_slurm_rpc_old_job_alloc(msg);
		slurm_free_old_job_alloc_msg(msg->data);
		break;
	case REQUEST_PING:
		_slurm_rpc_ping(msg);
		/* No body to free */
		break;
	case REQUEST_RECONFIGURE:
		_slurm_rpc_reconfigure_controller(msg);
		/* No body to free */
		break;
	case REQUEST_CONTROL:
		_slurm_rpc_shutdown_controller(msg);
		/* No body to free */
		break;
	case REQUEST_SHUTDOWN:
		_slurm_rpc_shutdown_controller(msg);
		slurm_free_shutdown_msg(msg->data);
		break;
	case REQUEST_SHUTDOWN_IMMEDIATE:
		_slurm_rpc_shutdown_controller_immediate(msg);
		/* No body to free */
		break;
	case REQUEST_SUBMIT_BATCH_JOB:
		_slurm_rpc_submit_batch_job(msg);
		slurm_free_job_desc_msg(msg->data);
		break;
	case REQUEST_UPDATE_JOB:
		_slurm_rpc_update_job(msg);
		slurm_free_job_desc_msg(msg->data);
		break;
	case REQUEST_UPDATE_NODE:
		_slurm_rpc_update_node(msg);
		slurm_free_update_node_msg(msg->data);
		break;
	case REQUEST_UPDATE_PARTITION:
		_slurm_rpc_update_partition(msg);
		slurm_free_update_part_msg(msg->data);
		break;
	case REQUEST_DELETE_PARTITION:
		_slurm_rpc_delete_partition(msg);
		slurm_free_delete_part_msg(msg->data);
		break;
	case REQUEST_NODE_REGISTRATION_STATUS:
		error("slurmctld is talking with itself. SlurmctldPort == SlurmdPort");
		slurm_send_rc_msg(msg, EINVAL);
		break;
	default:
		error("invalid RPC msg_type=%d", msg->msg_type);
		slurm_send_rc_msg(msg, EINVAL);
		break;
	}
}

/*
 * _fill_ctld_conf - make a copy of current slurm configuration
 *	this is done with locks set so the data can change at other times
 * OUT conf_ptr - place to copy configuration to
 */
void _fill_ctld_conf(slurm_ctl_conf_t * conf_ptr)
{
	conf_ptr->last_update         = time(NULL);
	conf_ptr->authtype            = slurmctld_conf.authtype;
	conf_ptr->backup_addr         = slurmctld_conf.backup_addr;
	conf_ptr->backup_controller   = slurmctld_conf.backup_controller;
	conf_ptr->control_addr        = slurmctld_conf.control_addr;
	conf_ptr->control_machine     = slurmctld_conf.control_machine;
	conf_ptr->epilog              = slurmctld_conf.epilog;
	conf_ptr->fast_schedule       = slurmctld_conf.fast_schedule;
	conf_ptr->first_job_id        = slurmctld_conf.first_job_id;
	conf_ptr->hash_base           = slurmctld_conf.hash_base;
	conf_ptr->heartbeat_interval  = slurmctld_conf.heartbeat_interval;
	conf_ptr->inactive_limit      = slurmctld_conf.inactive_limit;
	conf_ptr->job_comp_loc        = slurmctld_conf.job_comp_loc;
	conf_ptr->job_comp_type       = slurmctld_conf.job_comp_type;
	conf_ptr->job_credential_private_key = 
			slurmctld_conf.job_credential_private_key;
	conf_ptr->job_credential_public_certificate = 
			slurmctld_conf.job_credential_public_certificate;
	conf_ptr->kill_wait           = slurmctld_conf.kill_wait;
	conf_ptr->max_job_cnt         = slurmctld_conf.max_job_cnt;
	conf_ptr->min_job_age         = slurmctld_conf.min_job_age;
	conf_ptr->plugindir           = slurmctld_conf.plugindir;
	conf_ptr->prolog              = slurmctld_conf.prolog;
	conf_ptr->ret2service         = slurmctld_conf.ret2service;
	conf_ptr->schedauth           = slurmctld_conf.schedauth;
	conf_ptr->schedport           = slurmctld_conf.schedport;
	conf_ptr->schedtype           = slurmctld_conf.schedtype;
	conf_ptr->slurm_user_id       = slurmctld_conf.slurm_user_id;
	conf_ptr->slurm_user_name     = slurmctld_conf.slurm_user_name;
	conf_ptr->slurmctld_debug     = slurmctld_conf.slurmctld_debug;
	conf_ptr->slurmctld_logfile   = slurmctld_conf.slurmctld_logfile;
	conf_ptr->slurmctld_pidfile   = slurmctld_conf.slurmctld_pidfile;
	conf_ptr->slurmctld_port      = slurmctld_conf.slurmctld_port;
	conf_ptr->slurmctld_timeout   = slurmctld_conf.slurmctld_timeout;
	conf_ptr->slurmd_debug        = slurmctld_conf.slurmd_debug;
	conf_ptr->slurmd_logfile      = slurmctld_conf.slurmd_logfile;
	conf_ptr->slurmd_pidfile      = slurmctld_conf.slurmd_pidfile;
	conf_ptr->slurmd_port         = slurmctld_conf.slurmd_port;
	conf_ptr->slurmd_spooldir     = slurmctld_conf.slurmd_spooldir;
	conf_ptr->slurmd_timeout      = slurmctld_conf.slurmd_timeout;
	conf_ptr->slurm_conf          = slurmctld_conf.slurm_conf;
	conf_ptr->state_save_location = slurmctld_conf.state_save_location;
	conf_ptr->switch_type         = slurmctld_conf.switch_type;
	conf_ptr->tmp_fs              = slurmctld_conf.tmp_fs;
	conf_ptr->wait_time           = slurmctld_conf.wait_time;
	return;
}

/* return true if supplied uid is a super-user: root, self, or SlurmUser */
static inline bool _is_super_user(uid_t uid)
{
	if ( (uid == 0) || 
	     (uid == slurmctld_conf.slurm_user_id) ||
	     (uid == getuid()) )
		return true;
	else
		return false;
}

/* _kill_job_on_msg_fail - The request to create a job record successed, 
 *	but the reply message to srun failed. We kill the job to avoid 
 *	leaving it orphaned */
static void _kill_job_on_msg_fail(uint32_t job_id)
{
	/* Locks: Write job, write node */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };

	error("Job allocate response msg send failure, killing JobId=%u",
		job_id);
	lock_slurmctld(job_write_lock);
	job_complete(job_id, 0, false, 0);
	unlock_slurmctld(job_write_lock);
}

/* create a credential for a given job step, return error code */
static int _make_step_cred(struct step_record *step_rec, 
			   slurm_cred_t *slurm_cred)
{
	slurm_cred_arg_t cred_arg;

	cred_arg.jobid    = step_rec->job_ptr->job_id;
	cred_arg.stepid   = step_rec->step_id;
	cred_arg.uid      = step_rec->job_ptr->user_id;
	cred_arg.hostlist = step_rec->step_node_list;

	if ( (*slurm_cred = slurm_cred_create(slurmctld_config.cred_ctx, 
			&cred_arg)) == NULL) {
		error("slurm_cred_create error");
		return ESLURM_INVALID_JOB_CREDENTIAL;
	}

	return SLURM_SUCCESS;
}

/* _slurm_rpc_allocate_resources:  process RPC to allocate resources for 
 *	a job */
static void _slurm_rpc_allocate_resources(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	slurm_msg_t response_msg;
	DEF_TIMERS;
	job_desc_msg_t *job_desc_msg = (job_desc_msg_t *) msg->data;
	char *node_list_ptr = NULL;
	uint16_t num_cpu_groups = 0;
	uint32_t *cpus_per_node = NULL, *cpu_count_reps = NULL;
	uint32_t job_id = 0;
	resource_allocation_response_msg_t alloc_msg;
	/* Locks: Write job, write node, read partition */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, READ_LOCK };
	uid_t uid;
	uint16_t node_cnt = 0;
	slurm_addr *node_addr = NULL;
	int immediate = job_desc_msg->immediate;
	bool do_unlock = false;

	START_TIMER;
	debug2("Processing RPC: REQUEST_RESOURCE_ALLOCATION");

	/* do RPC call */
	dump_job_desc(job_desc_msg);
	uid = g_slurm_auth_get_uid(msg->cred);
	if ( (uid != job_desc_msg->user_id) && (!_is_super_user(uid)) ) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation, RESOURCE_ALLOCATE from uid=%u",
		      (unsigned int) uid);
	}

	if (error_code == SLURM_SUCCESS) {
		do_unlock = true;
		lock_slurmctld(job_write_lock);
		error_code = job_allocate(job_desc_msg, &job_id,
					  &node_list_ptr, &num_cpu_groups,
					  &cpus_per_node, &cpu_count_reps,
					  immediate, false, true, uid,
					  &node_cnt, &node_addr);
		/* unlock after finished using the job structure data */
		END_TIMER;
	}

	/* return result */
	if ((error_code == SLURM_SUCCESS) ||
	    ((immediate == 0) && 
	     (error_code == ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE))) {
		info("_slurm_rpc_allocate_resources JobId=%u NodeList=%s %s",
			job_id, node_list_ptr, TIME_STR);

		/* send job_ID  and node_name_ptr */
		alloc_msg.cpu_count_reps = xmalloc(sizeof(uint32_t) *
				num_cpu_groups);
		memcpy(alloc_msg.cpu_count_reps, cpu_count_reps,
				(sizeof(uint32_t) * num_cpu_groups));
		alloc_msg.cpus_per_node  = xmalloc(sizeof(uint32_t) *
				num_cpu_groups);
		memcpy(alloc_msg.cpus_per_node, cpus_per_node,
				(sizeof(uint32_t) * num_cpu_groups));
		alloc_msg.error_code     = error_code;
		alloc_msg.job_id         = job_id;
		alloc_msg.node_addr      = xmalloc(sizeof(slurm_addr) *
				node_cnt);
		memcpy(alloc_msg.node_addr, node_addr, 
				(sizeof(slurm_addr) * node_cnt));
		alloc_msg.node_cnt       = node_cnt;
		alloc_msg.node_list      = xstrdup(node_list_ptr);
		alloc_msg.num_cpu_groups = num_cpu_groups;
		unlock_slurmctld(job_write_lock);

		response_msg.msg_type = RESPONSE_RESOURCE_ALLOCATION;
		response_msg.data = &alloc_msg;

		if (slurm_send_node_msg(msg->conn_fd, &response_msg) < 0)
			_kill_job_on_msg_fail(job_id);
		xfree(alloc_msg.cpu_count_reps);
		xfree(alloc_msg.cpus_per_node);
		xfree(alloc_msg.node_addr);
		xfree(alloc_msg.node_list);
		(void) dump_all_job_state();
	} else {	/* allocate error */
		if (do_unlock)
			unlock_slurmctld(job_write_lock);
		info("_slurm_rpc_allocate_resources: %s ", 
		     slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	}
}

/* _slurm_rpc_allocate_and_run: process RPC to allocate resources for a job 
 *	and initiate a job step */
static void _slurm_rpc_allocate_and_run(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	slurm_msg_t response_msg;
	DEF_TIMERS;
	job_desc_msg_t *job_desc_msg = (job_desc_msg_t *) msg->data;
	char *node_list_ptr = NULL;
	uint16_t num_cpu_groups = 0;
	uint32_t *cpus_per_node = NULL, *cpu_count_reps = NULL;
	uint32_t job_id;
	resource_allocation_and_run_response_msg_t alloc_msg;
	struct step_record *step_rec;
	slurm_cred_t slurm_cred;
	job_step_create_request_msg_t req_step_msg;
	/* Locks: Write job, write node, read partition */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, READ_LOCK };
	uid_t uid;
	uint16_t node_cnt;
	slurm_addr *node_addr;
	int immediate = true;   /* implicit job_desc_msg->immediate == true */

	START_TIMER;
	debug2("Processing RPC: REQUEST_ALLOCATE_AND_RUN_JOB_STEP");

	/* do RPC call */
	dump_job_desc(job_desc_msg);
	uid = g_slurm_auth_get_uid(msg->cred);
	if ( (uid != job_desc_msg->user_id) && (!_is_super_user(uid)) ) {
		error("Security violation, ALLOCATE_AND_RUN RPC from uid=%u",
		      (unsigned int) uid);
		slurm_send_rc_msg(msg, ESLURM_USER_ID_MISSING);
		return;
	}

	lock_slurmctld(job_write_lock);
	error_code = job_allocate(job_desc_msg, &job_id,
				  &node_list_ptr, &num_cpu_groups,
				  &cpus_per_node, &cpu_count_reps,
				  immediate, false, true, uid,
				  &node_cnt, &node_addr);

	/* return result */
	if (error_code) {
		unlock_slurmctld(job_write_lock);
		info("_slurm_rpc_allocate_and_run: %s", 
		     slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
		return;
	}

	req_step_msg.job_id     = job_id;
	req_step_msg.user_id    = job_desc_msg->user_id;
	req_step_msg.node_count = INFINITE;
	req_step_msg.cpu_count  = job_desc_msg->num_procs;
	req_step_msg.num_tasks  = job_desc_msg->num_tasks;
	req_step_msg.task_dist  = job_desc_msg->task_dist;
	error_code = step_create(&req_step_msg, &step_rec, true);
	if (error_code == SLURM_SUCCESS) {
		error_code = _make_step_cred(step_rec, &slurm_cred);
		END_TIMER;
	}

	/* note: no need to free step_rec, pointer to global job step record */
	if (error_code) {
		job_complete(job_id, job_desc_msg->user_id, false, 0);
		unlock_slurmctld(job_write_lock);
		info("_slurm_rpc_allocate_and_run creating job step: %s",
			slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {

		info("_slurm_rpc_allocate_and_run JobId=%u NodeList=%s %s", 
			job_id, node_list_ptr, TIME_STR);

		/* send job_ID  and node_name_ptr */
		alloc_msg.job_id         = job_id;
		alloc_msg.node_list      = node_list_ptr;
		alloc_msg.num_cpu_groups = num_cpu_groups;
		alloc_msg.cpus_per_node  = cpus_per_node;
		alloc_msg.cpu_count_reps = cpu_count_reps;
		alloc_msg.job_step_id    = step_rec->step_id;
		alloc_msg.node_cnt       = node_cnt;
		alloc_msg.node_addr      = node_addr;
		alloc_msg.cred           = slurm_cred;
		alloc_msg.switch_job     = g_switch_copy_jobinfo(
						step_rec->switch_job);
		unlock_slurmctld(job_write_lock);
		response_msg.msg_type = RESPONSE_ALLOCATION_AND_RUN_JOB_STEP;
		response_msg.data = &alloc_msg;

		if (slurm_send_node_msg(msg->conn_fd, &response_msg) < 0)
			_kill_job_on_msg_fail(job_id);
		slurm_cred_destroy(slurm_cred);
		g_switch_free_jobinfo(alloc_msg.switch_job);
		(void) dump_all_job_state();	/* Has its own locks */
	}
}

/* _slurm_rpc_dump_conf - process RPC for Slurm configuration information */
static void _slurm_rpc_dump_conf(slurm_msg_t * msg)
{
	DEF_TIMERS;
	slurm_msg_t response_msg;
	last_update_msg_t *last_time_msg = (last_update_msg_t *) msg->data;
	slurm_ctl_conf_info_msg_t config_tbl;
	/* Locks: Read config */
	slurmctld_lock_t config_read_lock = { 
		READ_LOCK, NO_LOCK, NO_LOCK, NO_LOCK };

	START_TIMER;
	debug2("Processing RPC: REQUEST_BUILD_INFO");
	lock_slurmctld(config_read_lock);

	/* check to see if configuration data has changed */
	if ((last_time_msg->last_update - 1) >= slurmctld_conf.last_update) {
		unlock_slurmctld(config_read_lock);
		debug2("_slurm_rpc_dump_conf, no change");
		slurm_send_rc_msg(msg, SLURM_NO_CHANGE_IN_DATA);
	} else {
		_fill_ctld_conf(&config_tbl);
		unlock_slurmctld(config_read_lock);
		END_TIMER;
		debug2("_slurm_rpc_dump_conf %s", TIME_STR);

		/* init response_msg structure */
		response_msg.address = msg->address;
		response_msg.msg_type = RESPONSE_BUILD_INFO;
		response_msg.data = &config_tbl;

		/* send message */
		slurm_send_node_msg(msg->conn_fd, &response_msg);
	}
}

/* _slurm_rpc_dump_jobs - process RPC for job state information */
static void _slurm_rpc_dump_jobs(slurm_msg_t * msg)
{
	DEF_TIMERS;
	char *dump;
	int dump_size;
	slurm_msg_t response_msg;
	job_info_request_msg_t *last_time_msg =
	    (job_info_request_msg_t *) msg->data;
	/* Locks: Read job */
	slurmctld_lock_t job_read_lock = { 
		NO_LOCK, READ_LOCK, NO_LOCK, NO_LOCK };

	START_TIMER;
	debug2("Processing RPC: REQUEST_JOB_INFO");
	lock_slurmctld(job_read_lock);

	if ((last_time_msg->last_update - 1) >= last_job_update) {
		unlock_slurmctld(job_read_lock);
		debug2("_slurm_rpc_dump_jobs, no change");
		slurm_send_rc_msg(msg, SLURM_NO_CHANGE_IN_DATA);
	} else {
		pack_all_jobs(&dump, &dump_size);
		unlock_slurmctld(job_read_lock);
		END_TIMER;
		debug2("_slurm_rpc_dump_jobs, size=%d %s",
		     dump_size, TIME_STR);

		/* init response_msg structure */
		response_msg.address = msg->address;
		response_msg.msg_type = RESPONSE_JOB_INFO;
		response_msg.data = dump;
		response_msg.data_size = dump_size;

		/* send message */
		slurm_send_node_msg(msg->conn_fd, &response_msg);
		xfree(dump);
	}
}

/* _slurm_rpc_dump_nodes - process RPC for node state information */
static void _slurm_rpc_dump_nodes(slurm_msg_t * msg)
{
	DEF_TIMERS;
	char *dump;
	int dump_size;
	slurm_msg_t response_msg;
	last_update_msg_t *last_time_msg = (last_update_msg_t *) msg->data;
	/* Locks: Read node */
	slurmctld_lock_t node_read_lock = { 
		NO_LOCK, NO_LOCK, READ_LOCK, NO_LOCK };

	START_TIMER;
	debug2("Processing RPC: REQUEST_NODE_INFO");
	lock_slurmctld(node_read_lock);

	if ((last_time_msg->last_update - 1) >= last_node_update) {
		unlock_slurmctld(node_read_lock);
		debug2("_slurm_rpc_dump_nodes, no change");
		slurm_send_rc_msg(msg, SLURM_NO_CHANGE_IN_DATA);
	} else {
		pack_all_node(&dump, &dump_size);
		unlock_slurmctld(node_read_lock);
		END_TIMER;
		debug2("_slurm_rpc_dump_nodes, size=%d %s",
		     dump_size, TIME_STR);

		/* init response_msg structure */
		response_msg.address = msg->address;
		response_msg.msg_type = RESPONSE_NODE_INFO;
		response_msg.data = dump;
		response_msg.data_size = dump_size;

		/* send message */
		slurm_send_node_msg(msg->conn_fd, &response_msg);
		xfree(dump);
	}
}

/* _slurm_rpc_dump_partitions - process RPC for partition state information */
static void _slurm_rpc_dump_partitions(slurm_msg_t * msg)
{
	DEF_TIMERS;
	char *dump;
	int dump_size;
	slurm_msg_t response_msg;
	last_update_msg_t *last_time_msg = (last_update_msg_t *) msg->data;
	/* Locks: Read partition */
	slurmctld_lock_t part_read_lock = { 
		NO_LOCK, NO_LOCK, NO_LOCK, READ_LOCK };

	START_TIMER;
	debug2("Processing RPC: REQUEST_PARTITION_INFO");
	lock_slurmctld(part_read_lock);

	if ((last_time_msg->last_update - 1) >= last_part_update) {
		unlock_slurmctld(part_read_lock);
		debug2("_slurm_rpc_dump_partitions, no change");
		slurm_send_rc_msg(msg, SLURM_NO_CHANGE_IN_DATA);
	} else {
		pack_all_part(&dump, &dump_size);
		unlock_slurmctld(part_read_lock);
		END_TIMER;
		debug2("_slurm_rpc_dump_partitions, size=%d %s",
		     dump_size, TIME_STR);

		/* init response_msg structure */
		response_msg.address = msg->address;
		response_msg.msg_type = RESPONSE_PARTITION_INFO;
		response_msg.data = dump;
		response_msg.data_size = dump_size;

		/* send message */
		slurm_send_node_msg(msg->conn_fd, &response_msg);
		xfree(dump);
	}
}

/* _slurm_rpc_epilog_complete - process RPC noting the completion of 
 * the epilog denoting the completion of a job it its entirety */
static void  _slurm_rpc_epilog_complete(slurm_msg_t * msg)
{
	DEF_TIMERS;
	/* Locks: Write job, write node */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };
	uid_t uid;
	epilog_complete_msg_t *epilog_msg = 
		(epilog_complete_msg_t *) msg->data;
	bool run_scheduler = false;

	START_TIMER;
	debug2("Processing RPC: MESSAGE_EPILOG_COMPLETE");
	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error("Security violation, EPILOG_COMPLETE RPC from uid=%u",
		      (unsigned int) uid);
		return;
	}

	lock_slurmctld(job_write_lock);
	if (job_epilog_complete(epilog_msg->job_id, epilog_msg->node_name, 
	                        epilog_msg->return_code))
		run_scheduler = true;
	unlock_slurmctld(job_write_lock);
	END_TIMER;

	if (epilog_msg->return_code)
		error("_slurm_rpc_epilog_complete JobId=%u Node=%s Err=%s %s",
			epilog_msg->job_id, epilog_msg->node_name,
			slurm_strerror(epilog_msg->return_code), TIME_STR);
	else
		debug2("_slurm_rpc_epilog_complete JobId=%u Node=%s %s",
			epilog_msg->job_id, epilog_msg->node_name,
			TIME_STR);

	/* Functions below provide their own locking */
	if (run_scheduler) {
		(void) schedule();
		(void) dump_all_node_state();
		(void) dump_all_job_state();
	}

	/* NOTE: RPC has no response */
}

/* _slurm_rpc_job_step_kill - process RPC to cancel an entire job or 
 * an individual job step */
static void _slurm_rpc_job_step_kill(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	job_step_kill_msg_t *job_step_kill_msg =
	    (job_step_kill_msg_t *) msg->data;
	/* Locks: Write job, write node */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_CANCEL_JOB_STEP");
	uid = g_slurm_auth_get_uid(msg->cred);
	lock_slurmctld(job_write_lock);

	/* do RPC call */
	if (job_step_kill_msg->job_step_id == NO_VAL) {
		error_code = job_signal(job_step_kill_msg->job_id, 
					job_step_kill_msg->signal, 
					job_step_kill_msg->batch_flag, uid);
		unlock_slurmctld(job_write_lock);
		END_TIMER;

		/* return result */
		if (error_code) {
			info("_slurm_rpc_job_step_kill JobId=%u: %s", 
				job_step_kill_msg->job_id, 
				slurm_strerror(error_code));
			slurm_send_rc_msg(msg, error_code);
		} else {
			info("_slurm_rpc_job_step_kill JobId=%u %s",
				job_step_kill_msg->job_id, TIME_STR);
			slurm_send_rc_msg(msg, SLURM_SUCCESS);

			/* Below function provides its own locking */
			(void) dump_all_job_state();
		}
	} else {
		error_code = job_step_signal(job_step_kill_msg->job_id,
					     job_step_kill_msg->job_step_id,
					     job_step_kill_msg->signal,
					     uid);
		unlock_slurmctld(job_write_lock);
		END_TIMER;

		/* return result */
		if (error_code) {
			info("_slurm_rpc_job_step_kill StepId=%u.%u: %s",
				job_step_kill_msg->job_id, 
				job_step_kill_msg->job_step_id, 
				slurm_strerror(error_code));
			slurm_send_rc_msg(msg, error_code);
		} else {
			info("_slurm_rpc_job_step_kill StepId=%u.%u %s",
				job_step_kill_msg->job_id, 
				job_step_kill_msg->job_step_id, TIME_STR);
			slurm_send_rc_msg(msg, SLURM_SUCCESS);

			/* Below function provides its own locking */
			(void) dump_all_job_state();
		}
	}
}

/* _slurm_rpc_job_step_complete - process RPC to note the completion an  
 *	entire job or an individual job step */
static void _slurm_rpc_job_step_complete(slurm_msg_t * msg)
{
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	complete_job_step_msg_t *complete_job_step_msg =
	    (complete_job_step_msg_t *) msg->data;
	/* Locks: Write job, write node */
	slurmctld_lock_t job_write_lock = { NO_LOCK, WRITE_LOCK,
		WRITE_LOCK, NO_LOCK
	};
	uid_t uid;
	bool job_requeue = false;
	bool dump_job = false, dump_node = false;

	/* init */
	START_TIMER;
	debug2("Processing RPC: REQUEST_COMPLETE_JOB_STEP");
	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		/* Don't trust slurm_rc, it is not from slurmd */
		complete_job_step_msg->slurm_rc = SLURM_SUCCESS;
	}
	lock_slurmctld(job_write_lock);

	/* do RPC call */
	/* First set node DOWN if fatal error */
	if (complete_job_step_msg->slurm_rc == ESLURM_ALREADY_DONE) {
		/* race condition on job termination, not a real error */
		info("slurmd error running JobId=%u from node=%s: %s",
		      complete_job_step_msg->job_id,
		      complete_job_step_msg->node_name,
		      slurm_strerror(complete_job_step_msg->slurm_rc));
		complete_job_step_msg->slurm_rc = SLURM_SUCCESS;
	}
	if (complete_job_step_msg->slurm_rc != SLURM_SUCCESS) {
		error("Fatal slurmd error running JobId=%u from node=%s: %s",
		      complete_job_step_msg->job_id,
		      complete_job_step_msg->node_name,
		      slurm_strerror(complete_job_step_msg->slurm_rc));
		if (error_code == SLURM_SUCCESS) {
			update_node_msg_t update_node_msg;
			update_node_msg.node_names =
			    complete_job_step_msg->node_name;
			update_node_msg.node_state = NODE_STATE_DOWN;
			update_node_msg.reason = "step complete failure";
			error_code = update_node(&update_node_msg);
			if (complete_job_step_msg->job_rc != SLURM_SUCCESS)
				job_requeue = true;
			dump_job = true;
			dump_node = true;
		}
	}

	/* Mark job and/or job step complete */
	if (complete_job_step_msg->job_step_id == NO_VAL) {
		error_code = job_complete(complete_job_step_msg->job_id,
					  uid, job_requeue,
					  complete_job_step_msg->job_rc);
		unlock_slurmctld(job_write_lock);
		END_TIMER;

		/* return result */
		if (error_code) {
			info("_slurm_rpc_job_step_complete JobId=%u: %s ",
				complete_job_step_msg->job_id, 
				slurm_strerror(error_code));
			slurm_send_rc_msg(msg, error_code);
		} else {
			debug2("_slurm_rpc_job_step_complete JobId=%u %s", 
				complete_job_step_msg->job_id, TIME_STR);
			slurm_send_rc_msg(msg, SLURM_SUCCESS);
			dump_job = true;
		}
	} else {
		error_code =
		    job_step_complete(complete_job_step_msg->job_id,
				      complete_job_step_msg->job_step_id,
				      uid, job_requeue,
				      complete_job_step_msg->job_rc);
		unlock_slurmctld(job_write_lock);
		END_TIMER;

		/* return result */
		if (error_code) {
			info("_slurm_rpc_job_step_complete StepId=%u.%u: %s",
				complete_job_step_msg->job_id, 
				complete_job_step_msg->job_step_id, 
				slurm_strerror(error_code));
			slurm_send_rc_msg(msg, error_code);
		} else {
			info("_slurm_rpc_job_step_complete StepId=%u.%u %s",
				complete_job_step_msg->job_id, 
				complete_job_step_msg->job_step_id, TIME_STR);
			slurm_send_rc_msg(msg, SLURM_SUCCESS);
			dump_job = true;
		}
	}
	if (dump_job)
		(void) dump_all_job_state();	/* Has own locking */
	if (dump_node)
		(void) dump_all_node_state();	/* Has own locking */
}

/* _slurm_rpc_job_step_create - process RPC to creates/registers a job step 
 *	with the step_mgr */
static void _slurm_rpc_job_step_create(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	slurm_msg_t resp;
	struct step_record *step_rec;
	job_step_create_response_msg_t job_step_resp;
	job_step_create_request_msg_t *req_step_msg =
	    (job_step_create_request_msg_t *) msg->data;
	slurm_cred_t slurm_cred;
	/* Locks: Write jobs, read nodes */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, READ_LOCK, NO_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_JOB_STEP_CREATE");
	dump_step_desc(req_step_msg);
	uid = g_slurm_auth_get_uid(msg->cred);
	if ( (uid != req_step_msg->user_id) && (!_is_super_user(uid)) ) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation, JOB_STEP_CREATE RPC from uid=%u",
			(unsigned int) uid);
	}

	if (error_code == SLURM_SUCCESS) {
		/* issue the RPC */
		lock_slurmctld(job_write_lock);
		error_code = step_create(req_step_msg, &step_rec, false);
	}
	if (error_code == SLURM_SUCCESS)
		error_code = _make_step_cred(step_rec, &slurm_cred);
	END_TIMER;

	/* return result */
	if (error_code) {
		unlock_slurmctld(job_write_lock);
		info("_slurm_rpc_job_step_create: %s", 
			slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		info("_slurm_rpc_job_step_create: StepId=%u.%u %s",
		     step_rec->job_ptr->job_id, step_rec->step_id, TIME_STR);

		job_step_resp.job_step_id = step_rec->step_id;
		job_step_resp.node_list   = xstrdup(step_rec->step_node_list);
		job_step_resp.cred        = slurm_cred;
		job_step_resp.switch_job  = g_switch_copy_jobinfo(
						step_rec->switch_job);
		unlock_slurmctld(job_write_lock);
		resp.address = msg->address;
		resp.msg_type = RESPONSE_JOB_STEP_CREATE;
		resp.data = &job_step_resp;

		slurm_send_node_msg(msg->conn_fd, &resp);
		xfree(job_step_resp.node_list);
		slurm_cred_destroy(slurm_cred);
		g_switch_free_jobinfo(job_step_resp.switch_job);
		(void) dump_all_job_state();	/* Sets own locks */
	}
}

/* _slurm_rpc_job_step_get_info - process request for job step info */
static void _slurm_rpc_job_step_get_info(slurm_msg_t * msg)
{
	DEF_TIMERS;
	void *resp_buffer = NULL;
	int resp_buffer_size = 0;
	int error_code = SLURM_SUCCESS;
	job_step_info_request_msg_t *request =
	    (job_step_info_request_msg_t *) msg->data;
	/* Locks: Read job */
	slurmctld_lock_t job_read_lock = { NO_LOCK, READ_LOCK,
		NO_LOCK, NO_LOCK
	};

	START_TIMER;
	debug2("Processing RPC: REQUEST_JOB_STEP_INFO");
	lock_slurmctld(job_read_lock);

	if ((request->last_update - 1) >= last_job_update) {
		unlock_slurmctld(job_read_lock);
		debug2("_slurm_rpc_job_step_get_info, no change");
		error_code = SLURM_NO_CHANGE_IN_DATA;
	} else {
		Buf buffer;
		buffer = init_buf(BUF_SIZE);
		error_code =
		    pack_ctld_job_step_info_response_msg(request->job_id,
							 request->step_id,
							 buffer);
		unlock_slurmctld(job_read_lock);
		END_TIMER;
		if (error_code) {
			/* job_id:step_id not found or otherwise *\
			\* error message is printed elsewhere    */
			debug2("_slurm_rpc_job_step_get_info: %s",
				slurm_strerror(error_code));
			free_buf(buffer);
		} else {
			resp_buffer_size = get_buf_offset(buffer);
			resp_buffer = xfer_buf_data(buffer);
			debug2("_slurm_rpc_job_step_get_info size=%d %s",
			     resp_buffer_size, TIME_STR);
		}
	}

	if (error_code)
		slurm_send_rc_msg(msg, error_code);
	else {
		slurm_msg_t response_msg;

		response_msg.address = msg->address;
		response_msg.msg_type = RESPONSE_JOB_STEP_INFO;
		response_msg.data = resp_buffer;
		response_msg.data_size = resp_buffer_size;
		slurm_send_node_msg(msg->conn_fd, &response_msg);
		xfree(resp_buffer);
	}
}

/* _slurm_rpc_job_will_run - process RPC to determine if job with given 
 *	configuration can be initiated */
static void _slurm_rpc_job_will_run(slurm_msg_t * msg)
{
	/* init */
	DEF_TIMERS;
	int error_code = SLURM_SUCCESS;
	uint16_t num_cpu_groups = 0;
	uint32_t *cpus_per_node = NULL, *cpu_count_reps = NULL;
	uint32_t job_id;
	job_desc_msg_t *job_desc_msg = (job_desc_msg_t *) msg->data;
	char *node_list_ptr = NULL;
	/* Locks: Write job, read node, read partition */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, READ_LOCK, READ_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_JOB_WILL_RUN");

	/* do RPC call */
	dump_job_desc(job_desc_msg);
	uid = g_slurm_auth_get_uid(msg->cred);
	if ( (uid != job_desc_msg->user_id) && (!_is_super_user(uid)) ) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation, JOB_WILL_RUN RPC from uid=%u",
		      (unsigned int) uid);
	}

	if (error_code == SLURM_SUCCESS) {
		lock_slurmctld(job_write_lock);
		error_code = job_allocate(job_desc_msg, &job_id,
					  &node_list_ptr, &num_cpu_groups,
					  &cpus_per_node, &cpu_count_reps,
					  false, true, true, uid, NULL,
					  NULL);
		unlock_slurmctld(job_write_lock);
		END_TIMER;
	}

	/* return result */
	if (error_code) {
		info("_slurm_rpc_job_will_run: %s", 
			slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		debug2("_slurm_rpc_job_will_run success %s", TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
	}
}

/* _slurm_rpc_node_registration - process RPC to determine if a node's 
 *	actual configuration satisfies the configured specification */
static void _slurm_rpc_node_registration(slurm_msg_t * msg)
{
	/* init */
	DEF_TIMERS;
	int error_code = SLURM_SUCCESS;
	slurm_node_registration_status_msg_t *node_reg_stat_msg =
	    (slurm_node_registration_status_msg_t *) msg->data;
	/* Locks: Write job and node */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: MESSAGE_NODE_REGISTRATION_STATUS");
	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation,  NODE_REGISTER RPC from uid=%u",
		      (unsigned int) uid);
	}
	if (error_code == SLURM_SUCCESS) {
		/* do RPC call */
		lock_slurmctld(job_write_lock);
		validate_jobs_on_node(node_reg_stat_msg->node_name,
				      &node_reg_stat_msg->job_count,
				      node_reg_stat_msg->job_id,
				      node_reg_stat_msg->step_id);
		error_code =
		    validate_node_specs(node_reg_stat_msg->node_name,
					node_reg_stat_msg->cpus,
					node_reg_stat_msg->
					real_memory_size,
					node_reg_stat_msg->
					temporary_disk_space,
					node_reg_stat_msg->job_count,
					node_reg_stat_msg->status);
		unlock_slurmctld(job_write_lock);
		END_TIMER;
	}

	/* return result */
	if (error_code) {
		error("_slurm_rpc_node_registration node=%s: %s",
			node_reg_stat_msg->node_name, 
			slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		debug2("_slurm_rpc_node_registration complete for %s %s",
			node_reg_stat_msg->node_name, TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
	}
}

/* _slurm_rpc_old_job_alloc - process RPC to get details on existing job */
static void _slurm_rpc_old_job_alloc(slurm_msg_t * msg)
{
	int error_code = SLURM_SUCCESS;
	slurm_msg_t response_msg;
	DEF_TIMERS;
	old_job_alloc_msg_t *job_desc_msg =
	    (old_job_alloc_msg_t *) msg->data;
	char *node_list_ptr = NULL;
	uint16_t num_cpu_groups = 0;
	uint32_t *cpus_per_node = NULL, *cpu_count_reps = NULL;
	resource_allocation_response_msg_t alloc_msg;
	/* Locks: Read job, read node */
	slurmctld_lock_t job_read_lock = { 
		NO_LOCK, READ_LOCK, READ_LOCK, NO_LOCK };
	uint16_t node_cnt;
	slurm_addr *node_addr;
	uid_t uid;
	bool do_unlock = false;

	START_TIMER;
	debug2("Processing RPC: REQUEST_OLD_JOB_RESOURCE_ALLOCATION");

	/* do RPC call */
	uid = g_slurm_auth_get_uid(msg->cred);
	if ( (uid != job_desc_msg->uid) && (!_is_super_user(uid)) ) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation, RESOURCE_ALLOCATE from uid=%u",
		      (unsigned int) uid);
	}
	if (error_code == SLURM_SUCCESS) {
		do_unlock = true;
		lock_slurmctld(job_read_lock);
		error_code = old_job_info(job_desc_msg->uid,
					  job_desc_msg->job_id,
					  &node_list_ptr, &num_cpu_groups,
					  &cpus_per_node, &cpu_count_reps,
					  &node_cnt, &node_addr);
		END_TIMER;
	}

	/* return result */
	if (error_code) {
		if (do_unlock)
			unlock_slurmctld(job_read_lock);
		debug2("_slurm_rpc_old_job_alloc: JobId=%u, uid=%u: %s",
			job_desc_msg->job_id, job_desc_msg->uid, 
			slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		debug2("_slurm_rpc_old_job_alloc JobId=%u NodeList=%s %s",
			job_desc_msg->job_id, node_list_ptr, TIME_STR);

		/* send job_ID  and node_name_ptr */
		alloc_msg.cpu_count_reps = xmalloc(sizeof(uint32_t) *
				num_cpu_groups);
		memcpy(alloc_msg.cpu_count_reps, cpu_count_reps,
				(sizeof(uint32_t) * num_cpu_groups));
		alloc_msg.cpus_per_node  = xmalloc(sizeof(uint32_t) *
				num_cpu_groups);
		memcpy(alloc_msg.cpus_per_node, cpus_per_node,
				(sizeof(uint32_t) * num_cpu_groups));
		alloc_msg.error_code     = error_code;
		alloc_msg.job_id         = job_desc_msg->job_id;
		alloc_msg.node_addr      = xmalloc(sizeof(slurm_addr) *
				node_cnt);
		memcpy(alloc_msg.node_addr, node_addr,
				(sizeof(slurm_addr) * node_cnt));
		alloc_msg.node_cnt       = node_cnt;
		alloc_msg.node_list      = xstrdup(node_list_ptr);
		alloc_msg.num_cpu_groups = num_cpu_groups;
		unlock_slurmctld(job_read_lock);

		response_msg.msg_type    = RESPONSE_RESOURCE_ALLOCATION;
		response_msg.data        = &alloc_msg;

		slurm_send_node_msg(msg->conn_fd, &response_msg);
		xfree(alloc_msg.cpu_count_reps);
		xfree(alloc_msg.cpus_per_node);
		xfree(alloc_msg.node_addr);
		xfree(alloc_msg.node_list);
	}
}

/* _slurm_rpc_ping - process ping RPC */
static void _slurm_rpc_ping(slurm_msg_t * msg)
{
	/* We could authenticate here, if desired */

	/* return result */
	slurm_send_rc_msg(msg, SLURM_SUCCESS);
}


/* _slurm_rpc_reconfigure_controller - process RPC to re-initialize 
 *	slurmctld from configuration file */
static void _slurm_rpc_reconfigure_controller(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	/* Locks: Write configuration, job, node and partition */
	slurmctld_lock_t config_write_lock = { 
		WRITE_LOCK, WRITE_LOCK, WRITE_LOCK, WRITE_LOCK };
	/* Locks: Read node */
	slurmctld_lock_t node_read_lock = { 
		NO_LOCK, NO_LOCK, READ_LOCK, NO_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_RECONFIGURE");
	uid = g_slurm_auth_get_uid(msg->cred);\
	if (!_is_super_user(uid)) {
		error("Security violation, RECONFIGURE RPC from uid=%u",
		      (unsigned int) uid);
		error_code = ESLURM_USER_ID_MISSING;
	}

	/* do RPC call */
	if (error_code == SLURM_SUCCESS) {
		lock_slurmctld(config_write_lock);
		error_code = read_slurm_conf(0);
		unlock_slurmctld(config_write_lock);
		if (error_code == SLURM_SUCCESS) {
			lock_slurmctld(node_read_lock);
			msg_to_slurmd(REQUEST_RECONFIGURE);
			unlock_slurmctld(node_read_lock);
		}
	}
	if (error_code == SLURM_SUCCESS) {  /* Stuff to do after unlock */
		_update_cred_key();
		if (slurmctld_config.daemonize && 
		    chdir(slurmctld_conf.state_save_location) < 0) {
			error("chdir to %s error %m",
			      slurmctld_conf.state_save_location);
		}
	}
	END_TIMER;

	/* return result */
	if (error_code) {
		error("_slurm_rpc_reconfigure_controller: %s",
			slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		info("_slurm_rpc_reconfigure_controller: completed %s", 
			TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
		schedule();	/* has its own locks */
		save_all_state();
	}
}

/* _slurm_rpc_shutdown_controller - process RPC to shutdown slurmctld */
static void _slurm_rpc_shutdown_controller(slurm_msg_t * msg)
{
	int error_code = SLURM_SUCCESS, i;
	uint16_t core_arg = 0;
	shutdown_msg_t *shutdown_msg = (shutdown_msg_t *) msg->data;
	uid_t uid;
	/* Locks: Read node */
	slurmctld_lock_t node_read_lock = { 
		NO_LOCK, NO_LOCK, READ_LOCK, NO_LOCK };

	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error("Security violation, SHUTDOWN RPC from uid=%u",
		      (unsigned int) uid);
		error_code = ESLURM_USER_ID_MISSING;
	}
	if (error_code);
	else if (msg->msg_type == REQUEST_CONTROL) {
		info("Performing RPC: REQUEST_CONTROL");
		/* resume backup mode */
		slurmctld_config.resume_backup = true;	
	} else {
		info("Performing RPC: REQUEST_SHUTDOWN");
		core_arg = shutdown_msg->core;
	}

	/* do RPC call */
	if (error_code);
	else if (core_arg)
		info("performing immeditate shutdown without state save");
	else if (slurmctld_config.shutdown_time)
		debug2("shutdown RPC issued when already in progress");
	else {
		if (msg->msg_type == REQUEST_SHUTDOWN) {
			/* This means (msg->msg_type != REQUEST_CONTROL) */
			lock_slurmctld(node_read_lock);
			msg_to_slurmd(REQUEST_SHUTDOWN);
			unlock_slurmctld(node_read_lock);
		}
		if (slurmctld_config.thread_id_sig)	/* signal clean-up */
			pthread_kill(slurmctld_config.thread_id_sig, SIGTERM);
		else {
			error("thread_id_sig undefined, hard shutdown");
			slurmctld_config.shutdown_time = time(NULL);
			/* send REQUEST_SHUTDOWN_IMMEDIATE RPC */
			slurmctld_shutdown();
		}
	}

	if (msg->msg_type == REQUEST_CONTROL) {
		/* Wait for workload to dry up before sending reply.
		 * One thread should remain, this one. */
		for (i = 1; i < CONTROL_TIMEOUT; i++) {
			if (slurmctld_config.server_thread_count <= 1)
				break;
			sleep(1);
		}
		if (slurmctld_config.server_thread_count > 1)
			error("REQUEST_CONTROL reply with %d active threads",
				slurmctld_config.server_thread_count);
		/* save_all_state();	performed by _slurmctld_background */
	}
	slurm_send_rc_msg(msg, error_code);
	if ((error_code == SLURM_SUCCESS) && core_arg &&
	    (slurmctld_config.thread_id_sig))
		pthread_kill(slurmctld_config.thread_id_sig, SIGABRT);
}

/* _slurm_rpc_shutdown_controller_immediate - process RPC to shutdown 
 *	slurmctld */
static void _slurm_rpc_shutdown_controller_immediate(slurm_msg_t * msg)
{
	int error_code = SLURM_SUCCESS;
	uid_t uid;

	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error
		    ("Security violation, SHUTDOWN_IMMEDIATE RPC from uid=%u",
		     (unsigned int) uid);
		error_code = ESLURM_USER_ID_MISSING;
	}

	/* do RPC call */
	/* No op: just used to knock loose accept RPC thread */
	if (error_code == SLURM_SUCCESS)
		debug("Performing RPC: REQUEST_SHUTDOWN_IMMEDIATE");
}

/* _slurm_rpc_submit_batch_job - process RPC to submit a batch job */
static void _slurm_rpc_submit_batch_job(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	uint32_t job_id;
	slurm_msg_t response_msg;
	submit_response_msg_t submit_msg;
	job_desc_msg_t *job_desc_msg = (job_desc_msg_t *) msg->data;
	/* Locks: Write job, read node, read partition */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, READ_LOCK, READ_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_SUBMIT_BATCH_JOB");

	/* do RPC call */
	dump_job_desc(job_desc_msg);
	uid = g_slurm_auth_get_uid(msg->cred);
	if ( (uid != job_desc_msg->user_id) && (!_is_super_user(uid)) ) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation, SUBMIT_JOB from uid=%u",
		      (unsigned int) uid);
	}
	if (error_code == SLURM_SUCCESS) {
		lock_slurmctld(job_write_lock);
		error_code = job_allocate(job_desc_msg, &job_id,
					  (char **) NULL,
					  (uint16_t *) NULL,
					  (uint32_t **) NULL,
					  (uint32_t **) NULL, false, false,
					  false, uid, NULL, NULL);
		unlock_slurmctld(job_write_lock);
		END_TIMER;
	}

	/* return result */
	if ((error_code != SLURM_SUCCESS) &&
	    (error_code != ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE)) {
		info("_slurm_rpc_submit_batch_job: %s",
		     slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		info(
		   "_slurm_rpc_submit_batch_job JobId=%u %s", 
		   job_id, TIME_STR);
		/* send job_ID */
		submit_msg.job_id     = job_id;
		submit_msg.error_code = error_code;
		response_msg.msg_type = RESPONSE_SUBMIT_BATCH_JOB;
		response_msg.data = &submit_msg;
		slurm_send_node_msg(msg->conn_fd, &response_msg);
		schedule();	/* has own locks */
		(void) dump_all_job_state();	/* has own locks */
	}
}

/* _slurm_rpc_update_job - process RPC to update the configuration of a 
 *	job (e.g. priority) */
static void _slurm_rpc_update_job(slurm_msg_t * msg)
{
	/* init */
	int error_code;
	DEF_TIMERS;
	job_desc_msg_t *job_desc_msg = (job_desc_msg_t *) msg->data;
	/* Locks: Write job, read node, read partition */
	slurmctld_lock_t job_write_lock = { 
		NO_LOCK, WRITE_LOCK, READ_LOCK, READ_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_UPDATE_JOB");

	/* do RPC call */
	uid = g_slurm_auth_get_uid(msg->cred);
	lock_slurmctld(job_write_lock);
	error_code = update_job(job_desc_msg, uid);
	unlock_slurmctld(job_write_lock);
	END_TIMER;

	/* return result */
	if (error_code) {
		error("_slurm_rpc_update_job JobId=%u: %s",
		     job_desc_msg->job_id, slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		debug2("_slurm_rpc_update_job complete JobId=%u %s", 
			job_desc_msg->job_id, TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
		/* Below functions provide their own locking */
		schedule();
		(void) dump_all_job_state();
	}
}

/* _slurm_rpc_update_node - process RPC to update the configuration of a 
 *	node (e.g. UP/DOWN) */
static void _slurm_rpc_update_node(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	update_node_msg_t *update_node_msg_ptr =
	    			(update_node_msg_t *) msg->data;
	/* Locks: Write job and write node */
	slurmctld_lock_t node_write_lock = { 
		NO_LOCK, WRITE_LOCK, WRITE_LOCK, NO_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_UPDATE_NODE");
	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error_code = ESLURM_USER_ID_MISSING;
		error("Security violation, UPDATE_NODE RPC from uid=%u",
		      (unsigned int) uid);
	}

	if (error_code == SLURM_SUCCESS) {
		/* do RPC call */
		lock_slurmctld(node_write_lock);
		error_code = update_node(update_node_msg_ptr);
		unlock_slurmctld(node_write_lock);
		END_TIMER;
	}

	/* return result */
	if (error_code) {
		info("_slurm_rpc_update_node for %s: %s",
		      update_node_msg_ptr->node_names,
		      slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		debug2("_slurm_rpc_update_node complete for %s %s", 
			update_node_msg_ptr->node_names, TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);
	}

	/* Below functions provide their own locks */
	if (schedule())
		(void) dump_all_job_state();
	(void) dump_all_node_state();
}

/* _slurm_rpc_update_partition - process RPC to update the configuration 
 *	of a partition (e.g. UP/DOWN) */
static void _slurm_rpc_update_partition(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	update_part_msg_t *part_desc_ptr = (update_part_msg_t *) msg->data;
	/* Locks: Read node, write partition */
	slurmctld_lock_t part_write_lock = { 
		NO_LOCK, NO_LOCK, READ_LOCK, WRITE_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_UPDATE_PARTITION");
	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error_code = ESLURM_USER_ID_MISSING;
		error
		    ("Security violation, UPDATE_PARTITION RPC from uid=%u",
		     (unsigned int) uid);
	}

	if (error_code == SLURM_SUCCESS) {
		/* do RPC call */
		lock_slurmctld(part_write_lock);
		error_code = update_part(part_desc_ptr);
		unlock_slurmctld(part_write_lock);
		END_TIMER;
	}

	/* return result */
	if (error_code) {
		info("_slurm_rpc_update_partition partition=%s: %s",
			part_desc_ptr->name, slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		debug2("_slurm_rpc_update_partition complete for %s %s",
			part_desc_ptr->name, TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);

		/* NOTE: These functions provide their own locks */
		(void) dump_all_part_state();
		if (schedule())
			(void) dump_all_job_state();
	}
}

/* _slurm_rpc_delete_partition - process RPC to delete a partition */
static void _slurm_rpc_delete_partition(slurm_msg_t * msg)
{
	/* init */
	int error_code = SLURM_SUCCESS;
	DEF_TIMERS;
	delete_part_msg_t *part_desc_ptr = (delete_part_msg_t *) msg->data;
	/* Locks: Read job, write partition */
	slurmctld_lock_t part_write_lock = { 
		NO_LOCK, WRITE_LOCK, NO_LOCK, WRITE_LOCK };
	uid_t uid;

	START_TIMER;
	debug2("Processing RPC: REQUEST_DELETE_PARTITION");
	uid = g_slurm_auth_get_uid(msg->cred);
	if (!_is_super_user(uid)) {
		error_code = ESLURM_USER_ID_MISSING;
		error
		    ("Security violation, DELETE_PARTITION RPC from uid=%u",
		     (unsigned int) uid);
	}

	if (error_code == SLURM_SUCCESS) {
		/* do RPC call */
		lock_slurmctld(part_write_lock);
		error_code = delete_partition(part_desc_ptr);
		unlock_slurmctld(part_write_lock);
		END_TIMER;
	}

	/* return result */
	if (error_code) {
		info("_slurm_rpc_delete_partition partition=%s: %s",
			part_desc_ptr->name, slurm_strerror(error_code));
		slurm_send_rc_msg(msg, error_code);
	} else {
		info("_slurm_rpc_delete_partition complete for %s %s",
			part_desc_ptr->name, TIME_STR);
		slurm_send_rc_msg(msg, SLURM_SUCCESS);

		/* NOTE: These functions provide their own locks */
		(void) dump_all_job_state();
		(void) dump_all_part_state();
		if (schedule())
			(void) dump_all_job_state();
	}
}

/* Reset the job credential key based upon configuration parameters */
static void _update_cred_key(void) 
{
	slurm_cred_ctx_key_update(slurmctld_config.cred_ctx, 
				  slurmctld_conf.job_credential_private_key);
}
