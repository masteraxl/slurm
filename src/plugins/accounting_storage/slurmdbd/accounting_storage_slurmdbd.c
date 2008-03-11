/*****************************************************************************\
 *  accounting_storage_slurmdbd.c - accounting interface to slurmdbd.
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#if HAVE_STDINT_H
#  include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#include <stdio.h>
#include <sys/types.h>
#include <pwd.h>

#include <slurm/slurm_errno.h>

#include "src/common/jobacct_common.h"
#include "src/common/read_config.h"
#include "src/common/slurm_accounting_storage.h"
#include "src/common/slurmdbd_defs.h"
#include "src/common/xstring.h"
#include "src/slurmctld/slurmctld.h"

/*
 * These variables are required by the generic plugin interface.  If they
 * are not found in the plugin, the plugin loader will ignore it.
 *
 * plugin_name - a string giving a human-readable description of the
 * plugin.  There is no maximum length, but the symbol must refer to
 * a valid string.
 *
 * plugin_type - a string suggesting the type of the plugin or its
 * applicability to a particular form of data or method of data handling.
 * If the low-level plugin API is used, the contents of this string are
 * unimportant and may be anything.  SLURM uses the higher-level plugin
 * interface which requires this string to be of the form
 *
 *	<application>/<method>
 *
 * where <application> is a description of the intended application of
 * the plugin (e.g., "jobacct" for SLURM job completion logging) and <method>
 * is a description of how this plugin satisfies that application.  SLURM will
 * only load job completion logging plugins if the plugin_type string has a 
 * prefix of "jobacct/".
 *
 * plugin_version - an unsigned 32-bit integer giving the version number
 * of the plugin.  If major and minor revisions are desired, the major
 * version number may be multiplied by a suitable magnitude constant such
 * as 100 or 1000.  Various SLURM versions will likely require a certain
 * minimum versions for their plugins as the job accounting API 
 * matures.
 */
const char plugin_name[] = "Accounting storage SLURMDBD plugin";
const char plugin_type[] = "accounting_storage/slurmdbd";
const uint32_t plugin_version = 100;

static char *cluster_name       = NULL;
static char *slurmdbd_auth_info = NULL;
static List local_association_list = NULL;

/*
 * init() is called when the plugin is loaded, before any other functions
 * are called.  Put global initialization here.
 */
extern int init ( void )
{
	static int first = 1;

	if (first) {
		/* since this can be loaded from many different places
		   only tell us once. */
		if (!(cluster_name = slurm_get_cluster_name()))
			fatal("%s requires ClusterName in slurm.conf",
			      plugin_name);
		
		slurmdbd_auth_info = slurm_get_accounting_storage_pass();
		if(!slurmdbd_auth_info)
			
			verbose("%s loaded AuthInfo=%s",
				plugin_name, slurmdbd_auth_info);
		slurm_open_slurmdbd_conn(slurmdbd_auth_info);

		first = 0;
	} else {
		debug4("%s loaded", plugin_name);
	}

	return SLURM_SUCCESS;
}

extern int fini ( void )
{
	xfree(cluster_name);
	xfree(slurmdbd_auth_info);

	if(local_association_list)
		list_destroy(local_association_list);

	slurm_close_slurmdbd_conn();
	return SLURM_SUCCESS;
}

extern void *acct_storage_p_get_connection()
{
	return NULL;
}

extern int acct_storage_p_close_connection(void *db_conn)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_add_users(void *db_conn,
				    List user_list)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_add_coord(void *db_conn,
				    char *acct, acct_user_cond_t *user_q)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_add_accts(void *db_conn,
				    List acct_list)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_add_clusters(void *db_conn,
				       List cluster_list)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_add_associations(void *db_conn,
					   List association_list)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_get_assoc_id(void *db_conn,
				       acct_association_rec_t *assoc)
{
	ListIterator itr = NULL;
	acct_association_rec_t * found_assoc = NULL;
	acct_association_rec_t * ret_assoc = NULL;

	if(!local_association_list) {
		acct_association_cond_t assoc_q;
		char *cluster_name = NULL;
		memset(&assoc_q, 0, sizeof(acct_association_cond_t));
		assoc_q.cluster_list = list_create(slurm_destroy_char);
		cluster_name = slurm_get_cluster_name();
		if(!cluster_name) {
			error("acct_storage_p_get_assoc_id: "
			      "no cluster name here going to get "
			      "all associations.");
		} else 
			list_append(assoc_q.cluster_list, cluster_name);
		local_association_list =
			acct_storage_g_get_associations(db_conn, &assoc_q);
		list_destroy(assoc_q.cluster_list);
		if(!local_association_list) {
			error("acct_storage_p_get_assoc_id: "
			      "no list was made.");
			return SLURM_SUCCESS;
		}
	}
	if((!assoc->cluster && !assoc->acct) && !assoc->id) {
		error("acct_storage_p_get_assoc_id: "
		      "You need to supply a cluster and account name to get "
		      "an association.");
		return SLURM_ERROR;
	}

	itr = list_iterator_create(local_association_list);
	while((found_assoc = list_next(itr))) {
		if(assoc->id) {
			if(assoc->id == found_assoc->id) {
				ret_assoc = found_assoc;
				break;
			}
			continue;
		} else {
			if((!found_assoc->acct 
			    || strcasecmp(assoc->acct,
					  found_assoc->acct))
			   || (!assoc->cluster 
			       || strcasecmp(assoc->cluster,
					     found_assoc->cluster))
			   || (assoc->user 
			       && (!found_assoc->user 
				   || strcasecmp(assoc->user,
						 found_assoc->user)))
			   || (!assoc->user && found_assoc->user 
			       && strcasecmp("none",
					     found_assoc->user)))
				continue;
			if(assoc->partition
			   && (!assoc->partition 
			       || strcasecmp(assoc->partition, 
					     found_assoc->partition))) {
				ret_assoc = found_assoc;
				continue;
			}
		}
		ret_assoc = found_assoc;
		break;
	}
	list_iterator_destroy(itr);

	if(!ret_assoc)
		return SLURM_ERROR;

	assoc->id = ret_assoc->id;
	if(!assoc->user)
		assoc->user = ret_assoc->user;
	if(!assoc->acct)
		assoc->acct = ret_assoc->acct;
	if(!assoc->cluster)
		assoc->cluster = ret_assoc->cluster;
	if(!assoc->partition)
		assoc->partition = ret_assoc->partition;

	return SLURM_SUCCESS;
}

extern int acct_storage_p_validate_assoc_id(void *db_conn,
					    uint32_t assoc_id)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_modify_users(void *db_conn,
				       acct_user_cond_t *user_q,
				       acct_user_rec_t *user)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_modify_user_admin_level(void *db_conn,
						  acct_user_cond_t *user_q)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_modify_accts(void *db_conn,
				       acct_account_cond_t *acct_q,
				       acct_account_rec_t *acct)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_modify_clusters(void *db_conn,
					  acct_cluster_cond_t *cluster_q,
					  acct_cluster_rec_t *cluster)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_modify_associations(void *db_conn,
					      acct_association_cond_t *assoc_q,
					      acct_association_rec_t *assoc)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_remove_users(void *db_conn,
				       acct_user_cond_t *user_q)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_remove_coord(void *db_conn,
				       char *acct, acct_user_cond_t *user_q)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_remove_accts(void *db_conn,
				       acct_account_cond_t *acct_q)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_remove_clusters(void *db_conn,
					  acct_account_cond_t *cluster_q)
{
	return SLURM_SUCCESS;
}

extern int acct_storage_p_remove_associations(void *db_conn,
					      acct_association_cond_t *assoc_q)
{
	return SLURM_SUCCESS;
}

extern List acct_storage_p_get_users(void *db_conn,
				     acct_user_cond_t *user_q)
{
	return NULL;
}

extern List acct_storage_p_get_accts(void *db_conn,
				     acct_account_cond_t *acct_q)
{
	return NULL;
}

extern List acct_storage_p_get_clusters(void *db_conn,
					acct_account_cond_t *cluster_q)
{
	return NULL;
}

extern List acct_storage_p_get_associations(void *db_conn,
					    acct_association_cond_t *assoc_q)
{
	
	return NULL;
}

extern int acct_storage_p_get_hourly_usage(void *db_conn,
					   acct_association_rec_t *acct_assoc,
					   time_t start, time_t end)
{
	int rc = SLURM_SUCCESS;

	return rc;
}

extern int acct_storage_p_get_daily_usage(void *db_conn,
					  acct_association_rec_t *acct_assoc,
					  time_t start, time_t end)
{
	int rc = SLURM_SUCCESS;

	return rc;
}

extern int acct_storage_p_get_monthly_usage(void *db_conn,
					    acct_association_rec_t *acct_assoc,
					    time_t start, time_t end)
{
	int rc = SLURM_SUCCESS;
	return rc;
}

extern int clusteracct_storage_p_node_down(void *db_conn,
					   char *cluster,
					   struct node_record *node_ptr,
					   time_t event_time, char *reason)
{
	slurmdbd_msg_t msg;
	dbd_node_state_msg_t req;

	req.cluster_name = cluster;
	req.hostlist   = node_ptr->name;
	req.new_state  = DBD_NODE_STATE_DOWN;
	req.event_time = event_time;
	req.reason     = reason;
	msg.msg_type   = DBD_NODE_STATE;
	msg.data       = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}
extern int clusteracct_storage_p_node_up(void *db_conn,
					 char *cluster,
					 struct node_record *node_ptr,
					 time_t event_time)
{
	slurmdbd_msg_t msg;
	dbd_node_state_msg_t req;

	req.cluster_name = cluster;
	req.hostlist   = node_ptr->name;
	req.new_state  = DBD_NODE_STATE_UP;
	req.event_time = event_time;
	req.reason     = NULL;
	msg.msg_type   = DBD_NODE_STATE;
	msg.data       = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}
extern int clusteracct_storage_p_cluster_procs(void *db_conn,
					       char *cluster,
					       uint32_t procs,
					       time_t event_time)
{
	slurmdbd_msg_t msg;
	dbd_cluster_procs_msg_t req;

	req.cluster_name = cluster;
	req.proc_count   = procs;
	req.event_time   = event_time;
	msg.msg_type     = DBD_CLUSTER_PROCS;
	msg.data         = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}

extern int clusteracct_storage_p_get_hourly_usage(void *db_conn,
				    
						  acct_cluster_rec_t *cluster_rec, time_t start, 
						  time_t end, void *params)
{

	return SLURM_SUCCESS;
}

extern int clusteracct_storage_p_get_daily_usage(void *db_conn,
				    
						 acct_cluster_rec_t *cluster_rec, time_t start, 
						 time_t end, void *params)
{
	
	return SLURM_SUCCESS;
}

extern int clusteracct_storage_p_get_monthly_usage(void *db_conn,
				    
						   acct_cluster_rec_t *cluster_rec, time_t start, 
						   time_t end, void *params)
{
	
	return SLURM_SUCCESS;
}

/* 
 * load into the storage the start of a job
 */
extern int jobacct_storage_p_job_start(void *db_conn,
				       struct job_record *job_ptr)
{
	slurmdbd_msg_t msg, msg_rc;
	dbd_job_start_msg_t req;
	dbd_job_start_rc_msg_t *resp;
	char *block_id = NULL;
	int rc = SLURM_SUCCESS;

	if (!job_ptr->details || !job_ptr->details->submit_time) {
		error("jobacct_storage_p_job_start: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	req.alloc_cpus    = job_ptr->total_procs;
	req.assoc_id      = job_ptr->assoc_id;
#ifdef HAVE_BG
	select_g_get_jobinfo(job_ptr->select_jobinfo, 
			     SELECT_DATA_BLOCK_ID, 
			     &block_id);
#endif
	req.block_id      = block_id;
	xfree(block_id);
	if (job_ptr->details) 
		req.eligible_time = job_ptr->details->begin_time;
	req.gid           = job_ptr->group_id;
	req.job_id        = job_ptr->job_id;
	req.job_state     = job_ptr->job_state & (~JOB_COMPLETING);
	req.name          = job_ptr->name;
	req.nodes         = job_ptr->nodes;
	req.partition     = job_ptr->partition;
	req.req_cpus      = job_ptr->num_procs;
	req.priority      = job_ptr->priority;
	req.start_time    = job_ptr->start_time;
	if (job_ptr->details)
		req.submit_time   = job_ptr->details->submit_time;

	msg.msg_type      = DBD_JOB_START;
	msg.data          = &req;
	rc = slurm_send_recv_slurmdbd_msg(&msg, &msg_rc);
	if (rc != SLURM_SUCCESS) {
		if (slurm_send_slurmdbd_msg(&msg) < 0)
			return SLURM_ERROR;
	} else if (msg_rc.msg_type != DBD_JOB_START_RC) {
		error("slurmdbd: response type not DBD_GOT_JOBS: %u", 
		      msg_rc.msg_type);
	} else {
		resp = (dbd_job_start_rc_msg_t *) msg_rc.data;
		job_ptr->db_index = resp->db_index;
		slurm_dbd_free_job_start_rc_msg(resp);
	}
	
	return rc;
}

/* 
 * load into the storage the end of a job
 */
extern int jobacct_storage_p_job_complete(void *db_conn,
					  struct job_record *job_ptr)
{
	slurmdbd_msg_t msg;
	dbd_job_comp_msg_t req;

	if (!job_ptr->db_index 
	    && (!job_ptr->details || !job_ptr->details->submit_time)) {
		error("jobacct_storage_p_job_complete: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	req.assoc_id    = job_ptr->assoc_id;
	req.db_index    = job_ptr->db_index;
	req.end_time    = job_ptr->end_time;
	req.exit_code   = job_ptr->exit_code;
	req.job_id      = job_ptr->job_id;
	req.job_state   = job_ptr->job_state & (~JOB_COMPLETING);
	req.nodes       = job_ptr->nodes;
	req.start_time  = job_ptr->start_time;
	if (job_ptr->details)
		req.submit_time   = job_ptr->details->submit_time;

	msg.msg_type    = DBD_JOB_COMPLETE;
	msg.data        = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}

/* 
 * load into the storage the start of a job step
 */
extern int jobacct_storage_p_step_start(void *db_conn,
					struct step_record *step_ptr)
{
	uint32_t cpus = 0;
	char node_list[BUFFER_SIZE];
	slurmdbd_msg_t msg;
	dbd_step_start_msg_t req;

#ifdef HAVE_BG
	char *ionodes = NULL;

	cpus = step_ptr->job_ptr->num_procs;
	select_g_get_jobinfo(step_ptr->job_ptr->select_jobinfo, 
			     SELECT_DATA_IONODES, 
			     &ionodes);
	if (ionodes) {
		snprintf(node_list, BUFFER_SIZE, 
			 "%s[%s]", step_ptr->job_ptr->nodes, ionodes);
		xfree(ionodes);
	} else {
		snprintf(node_list, BUFFER_SIZE, "%s",
			 step_ptr->job_ptr->nodes);
	}
	
#else
	if (!step_ptr->step_layout || !step_ptr->step_layout->task_cnt) {
		cpus = step_ptr->job_ptr->total_procs;
		snprintf(node_list, BUFFER_SIZE, "%s",
			 step_ptr->job_ptr->nodes);
	} else {
		cpus = step_ptr->step_layout->task_cnt;
		snprintf(node_list, BUFFER_SIZE, "%s", 
			 step_ptr->step_layout->node_list);
	}
#endif

	if (!step_ptr->job_ptr->db_index 
	    && (!step_ptr->job_ptr->details
		|| !step_ptr->job_ptr->details->submit_time)) {
		error("jobacct_storage_p_step_start: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	req.assoc_id    = step_ptr->job_ptr->assoc_id;
	req.db_index    = step_ptr->job_ptr->db_index;
	req.job_id      = step_ptr->job_ptr->job_id;
	req.name        = step_ptr->name;
	req.nodes       = node_list;
	req.start_time  = step_ptr->start_time;
	if (step_ptr->job_ptr->details)
		req.job_submit_time   = step_ptr->job_ptr->details->submit_time;
	req.step_id     = step_ptr->step_id;
	req.total_procs = cpus;

	msg.msg_type    = DBD_STEP_START;
	msg.data        = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}

/* 
 * load into the storage the end of a job step
 */
extern int jobacct_storage_p_step_complete(void *db_conn,
					   struct step_record *step_ptr)
{
	uint32_t cpus = 0;
	char node_list[BUFFER_SIZE];
	slurmdbd_msg_t msg;
	dbd_step_comp_msg_t req;

#ifdef HAVE_BG
	char *ionodes = NULL;

	cpus = step_ptr->job_ptr->num_procs;
	select_g_get_jobinfo(step_ptr->job_ptr->select_jobinfo, 
			     SELECT_DATA_IONODES, 
			     &ionodes);
	if (ionodes) {
		snprintf(node_list, BUFFER_SIZE, 
			 "%s[%s]", step_ptr->job_ptr->nodes, ionodes);
		xfree(ionodes);
	} else {
		snprintf(node_list, BUFFER_SIZE, "%s",
			 step_ptr->job_ptr->nodes);
	}
	
#else
	if (!step_ptr->step_layout || !step_ptr->step_layout->task_cnt) {
		cpus = step_ptr->job_ptr->total_procs;
		snprintf(node_list, BUFFER_SIZE, "%s", step_ptr->job_ptr->nodes);
	} else {
		cpus = step_ptr->step_layout->task_cnt;
		snprintf(node_list, BUFFER_SIZE, "%s", 
			 step_ptr->step_layout->node_list);
	}
#endif

	if (!step_ptr->job_ptr->db_index 
	    && (!step_ptr->job_ptr->details
		|| !step_ptr->job_ptr->details->submit_time)) {
		error("jobacct_storage_p_step_complete: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	req.assoc_id    = step_ptr->job_ptr->assoc_id;
	req.db_index    = step_ptr->job_ptr->db_index;
	req.end_time    = time(NULL);	/* called at step completion */
	req.jobacct     = step_ptr->jobacct;
	req.job_id      = step_ptr->job_ptr->job_id;
	req.req_uid     = step_ptr->job_ptr->requid;
	req.start_time  = step_ptr->start_time;
	if (step_ptr->job_ptr->details)
		req.job_submit_time   = step_ptr->job_ptr->details->submit_time;
	req.step_id     = step_ptr->step_id;
	req.total_procs = cpus;

	msg.msg_type    = DBD_STEP_COMPLETE;
	msg.data        = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}

/* 
 * load into the storage a suspention of a job
 */
extern int jobacct_storage_p_suspend(void *db_conn,
				     struct job_record *job_ptr)
{
	slurmdbd_msg_t msg;
	dbd_job_suspend_msg_t req;

	req.assoc_id     = 0;	/* FIXME */
	req.job_id       = job_ptr->job_id;
	req.job_state    = job_ptr->job_state & (~JOB_COMPLETING);
	if (job_ptr->details)
		req.submit_time   = job_ptr->details->submit_time;
	req.suspend_time = job_ptr->suspend_time;
	msg.msg_type     = DBD_JOB_SUSPEND;
	msg.data         = &req;

	if (slurm_send_slurmdbd_msg(&msg) < 0)
		return SLURM_ERROR;

	return SLURM_SUCCESS;
}

/* 
 * get info from the storage 
 * returns List of job_rec_t *
 * note List needs to be freed when called
 */
extern List jobacct_storage_p_get_jobs(void *db_conn,
				       List selected_steps,
				       List selected_parts,
				       sacct_parameters_t *params)
{
	slurmdbd_msg_t req, resp;
	dbd_get_jobs_msg_t get_msg;
	dbd_list_msg_t *got_msg;
	int rc;
	List job_list = NULL;
	struct passwd *pw = NULL;

	get_msg.selected_steps = selected_steps;
	get_msg.selected_parts = selected_parts;
	get_msg.cluster_name = params->opt_cluster;
	get_msg.gid = params->opt_gid;
	
	if (params->opt_uid >=0 && (pw=getpwuid(params->opt_uid)))
		get_msg.user = pw->pw_name;
	else
		get_msg.user = NULL;

	req.msg_type = DBD_GET_JOBS;
	req.data = &get_msg;
	rc = slurm_send_recv_slurmdbd_msg(&req, &resp);

	if (rc != SLURM_SUCCESS)
		error("slurmdbd: DBD_GET_JOBS failure: %m");
	else if (resp.msg_type != DBD_GOT_JOBS) {
		error("slurmdbd: response type not DBD_GOT_JOBS: %u", 
		      resp.msg_type);
	} else {
		got_msg = (dbd_list_msg_t *) resp.data;
		job_list = got_msg->ret_list;
		got_msg->ret_list = NULL;
		slurm_dbd_free_list_msg(got_msg);
	}

	return job_list;
}

/* 
 * Expire old info from the storage
 * Not applicable for any database
 */
extern void jobacct_storage_p_archive(void *db_conn,
				      List selected_parts,
				      void *params)
{
	return;
}
