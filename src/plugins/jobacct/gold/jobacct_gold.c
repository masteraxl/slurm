/*****************************************************************************\
 *  jobacct_gold.c - jobacct interface to gold.
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
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
#include "gold_interface.h"

#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>

#include "src/common/xmalloc.h"
#include "src/common/list.h"
#include "src/common/xstring.h"
#include "src/common/uid.h"
#include <src/common/parse_time.h>

#include "src/slurmctld/slurmctld.h"
#include "src/slurmd/slurmd/slurmd.h"
#include "src/common/slurm_jobacct.h"
#include "src/common/slurm_protocol_api.h"


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
const char plugin_name[] = "Job accounting GOLD plugin";
const char plugin_type[] = "jobacct/gold";
const uint32_t plugin_version = 100;

/* for this first draft we are only supporting one cluster per slurm
 * 1.3 will probably do better than this.
 */

static char *cluster_name = NULL;

/* _check_for_job 
 * IN jobid - job id to check for 
 * IN submit - timestamp for submit time of job
 * RET 0 for not found 1 for found
 */

static int _check_for_job(uint32_t jobid, time_t submit) 
{
	gold_request_t *gold_request = create_gold_request(GOLD_OBJECT_JOB,
							   GOLD_ACTION_QUERY);
	gold_response_t *gold_response = NULL;
	char tmp_buff[50];
	int rc = 0;

	if(!gold_request) 
		return rc;

	gold_request_add_selection(gold_request, "JobId");

	snprintf(tmp_buff, sizeof(tmp_buff), "%u", jobid);
	gold_request_add_condition(gold_request, "JobId", tmp_buff);

	snprintf(tmp_buff, sizeof(tmp_buff), "%u", (int)submit);
	gold_request_add_condition(gold_request, "SubmitTime", tmp_buff);

	gold_response = get_gold_response(gold_request);
	destroy_gold_request(gold_request);

	if(!gold_response) {
		error("_check_for_job: no response received");
		return 0;
	}

	if(gold_response->entry_cnt > 0) 
		rc = 1;
	destroy_gold_response(gold_response);
	
	return rc;
}

static char *_get_account_id(char *user, char *project, char *machine)
{
	gold_request_t *gold_request = create_gold_request(GOLD_OBJECT_ACCOUNT,
							   GOLD_ACTION_QUERY);
	gold_response_t *gold_response = NULL;
	char *gold_account_id = NULL;
	gold_response_entry_t *resp_entry = NULL;
	gold_name_value_t *name_val = NULL;

	gold_request_add_selection(gold_request, "Id");
	gold_request_add_condition(gold_request, "User", user);
	if(project)
		gold_request_add_condition(gold_request, "Project", project);
	gold_request_add_condition(gold_request, "Machine", machine);
		
	gold_response = get_gold_response(gold_request);
	destroy_gold_request(gold_request);

	if(!gold_response) {
		error("_get_account_id: no response received");
		return NULL;
	}

	if(gold_response->entry_cnt > 0) {
		resp_entry = list_pop(gold_response->entries);
		name_val = list_pop(resp_entry->name_val);

		gold_account_id = xstrdup(name_val->value);

		destroy_gold_name_value(name_val);
		destroy_gold_response_entry(resp_entry);
	} else {
		error("no account found returning 0");
		gold_account_id = xstrdup("0");
	}

	destroy_gold_response(gold_response);

	return gold_account_id;
}

static int _add_edit_job(struct job_record *job_ptr, gold_object_t action)
{
	gold_request_t *gold_request = create_gold_request(GOLD_OBJECT_JOB,
							   action);
	gold_response_t *gold_response = NULL;
	char tmp_buff[50];
	int rc = SLURM_ERROR;
	char *gold_account_id = NULL;
	char *user = uid_to_string((uid_t)job_ptr->user_id);
	char *jname = NULL;
	int ncpus=0, tmp = 0, i = 0;
	char *account = NULL;
	char *nodes = "(null)";

	if(!gold_request) 
		return rc;

	if ((tmp = strlen(job_ptr->name))) {
		jname = xmalloc(++tmp);
		for (i=0; i<tmp; i++) {
			if (isspace(job_ptr->name[i]))
				jname[i]='_';
			else
				jname[i]=job_ptr->name[i];
		}
	} else
		jname = xstrdup("allocation");
	
	if (job_ptr->account && job_ptr->account[0])
		account = job_ptr->account;
	
	if (job_ptr->nodes && job_ptr->nodes[0])
		nodes = job_ptr->nodes;
	
	for (i=0; i < job_ptr->num_cpu_groups; i++) {
		ncpus += (job_ptr->cpus_per_node[i])
			* (job_ptr->cpu_count_reps[i]);
		//info("got %d from %d * %d", ncpus, job_ptr->cpus_per_node[i],
		//   job_ptr->cpu_count_reps[i]);
	}
//info("total procs is  %d", job_ptr->details->total_procs);
	if(action == GOLD_ACTION_CREATE) {
		snprintf(tmp_buff, sizeof(tmp_buff), "%u", job_ptr->job_id);
		gold_request_add_assignment(gold_request, "JobId", tmp_buff);
		
		gold_account_id = _get_account_id(user, account, 
						  cluster_name);

		gold_request_add_assignment(gold_request, "GoldAccountId",
					    gold_account_id);
		xfree(gold_account_id);

		snprintf(tmp_buff, sizeof(tmp_buff), "%u",
			 (int)job_ptr->details->submit_time);
		gold_request_add_assignment(gold_request, "SubmitTime",
					    tmp_buff);

	} else if (action == GOLD_ACTION_MODIFY) {
		snprintf(tmp_buff, sizeof(tmp_buff), "%u", job_ptr->job_id);
		gold_request_add_condition(gold_request, "JobId", tmp_buff);
		
		snprintf(tmp_buff, sizeof(tmp_buff), "%u",
			 (int)job_ptr->details->submit_time);
		gold_request_add_condition(gold_request, "SubmitTime",
					   tmp_buff);
								
		snprintf(tmp_buff, sizeof(tmp_buff), "%u",
			 (int)job_ptr->end_time);
		gold_request_add_assignment(gold_request, "EndTime",
					    tmp_buff);		
		
		snprintf(tmp_buff, sizeof(tmp_buff), "%u",
			 (int)job_ptr->exit_code);
		gold_request_add_assignment(gold_request, "ExitCode",
					    tmp_buff);

	} else {
		destroy_gold_request(gold_request);
		error("_add_edit_job: bad action given %d", action);		
		return rc;
	}

	gold_request_add_assignment(gold_request, "Partition",
				    job_ptr->partition);
	
	snprintf(tmp_buff, sizeof(tmp_buff), "%u", job_ptr->num_procs);
	gold_request_add_assignment(gold_request, "RequestedCPUS",
				    tmp_buff);
	snprintf(tmp_buff, sizeof(tmp_buff), "%u", ncpus);
	gold_request_add_assignment(gold_request, "AllocatedCPUS",
				    tmp_buff);

	gold_request_add_assignment(gold_request, "NodeList",
				    nodes);

	gold_request_add_assignment(gold_request, "JobName",
				    jname);
	xfree(jname);

/* 	gold_request_add_assignment(gold_request, "CPUSecondsReserved", */
/* 	     		            ); */


	snprintf(tmp_buff, sizeof(tmp_buff), "%u",
		 (int)job_ptr->details->begin_time);
	gold_request_add_assignment(gold_request, "EligibleTime",
				    tmp_buff);

	snprintf(tmp_buff, sizeof(tmp_buff), "%u",
		 (int)job_ptr->start_time);
	gold_request_add_assignment(gold_request, "StartTime",
				    tmp_buff);
		
	snprintf(tmp_buff, sizeof(tmp_buff), "%u",
		 job_ptr->job_state & (~JOB_COMPLETING));
	gold_request_add_assignment(gold_request, "State",
				    tmp_buff);

	

	gold_response = get_gold_response(gold_request);	
	destroy_gold_request(gold_request);

	if(!gold_response) {
		error("_add_edit_job: no response received");
		return rc;
	}

	if(!gold_response->rc) 
		rc = SLURM_SUCCESS;
	else {
		error("gold_response has non-zero rc(%d): %s",
		      gold_response->rc,
		      gold_response->message);
	}
	destroy_gold_response(gold_response);

	return rc;
}

/*
 * init() is called when the plugin is loaded, before any other functions
 * are called.  Put global initialization here.
 */
extern int init ( void )
{
	verbose("%s loaded", plugin_name);
	return SLURM_SUCCESS;
}

extern int fini ( void )
{
	return SLURM_SUCCESS;
}

/*
 * The following routines are called by slurmctld
 */

/*
 * The following routines are called by slurmd
 */
int jobacct_p_init_struct(struct jobacctinfo *jobacct, 
			  jobacct_id_t *jobacct_id)
{
	return SLURM_SUCCESS;
}

struct jobacctinfo *jobacct_p_alloc(jobacct_id_t *jobacct_id)
{
	return NULL;
}

void jobacct_p_free(struct jobacctinfo *jobacct)
{
	return;
}

int jobacct_p_setinfo(struct jobacctinfo *jobacct, 
		      enum jobacct_data_type type, void *data)
{
	return SLURM_SUCCESS;
	
}

int jobacct_p_getinfo(struct jobacctinfo *jobacct, 
		      enum jobacct_data_type type, void *data)
{
	return SLURM_SUCCESS;
}

void jobacct_p_aggregate(struct jobacctinfo *dest, struct jobacctinfo *from)
{
	return;
}

void jobacct_p_2_sacct(sacct_t *sacct, struct jobacctinfo *jobacct)
{
	return;
}

void jobacct_p_pack(struct jobacctinfo *jobacct, Buf buffer)
{
	return;
}

int jobacct_p_unpack(struct jobacctinfo **jobacct, Buf buffer)
{
	return SLURM_SUCCESS;
}


int jobacct_p_init_slurmctld(char *gold_info)
{
	char *total = "localhost:/etc/gold/auth_key:localhost:7112";
	int found = 0;
	int i=0, j=0;
	char *host = NULL;
	char *keyfile = NULL;
	uint16_t port = 0;

	debug2("jobacct_init() called");
	if(gold_info) 
		total = gold_info;
	i = 0;
	while(total[j]) {
		if(total[j] == ':') {
			switch(found) {
			case 0: // cluster_name name
			        cluster_name = xstrndup(total+i, j-i);
				break;
			case 1: // keyfile name
				keyfile = xstrndup(total+i, j-i);
				break;
			case 2: // host name
				host = xstrndup(total+i, j-i);
				break;
			case 3: // port
				port = atoi(total+i);
				break;
			}
			found++;
			i = j+1;	
		}
		j++;
	}
	if(!port) 
		port = atoi(total+i);

	if(!cluster_name)
		fatal("JobAcctLogfile should be in the format of "
		      "cluster_name:gold_auth_key_file_path:"
		      "goldd_host:goldd_port "
		      "bad cluster_name");
	if (!keyfile || *keyfile != '/')
		fatal("JobAcctLogfile should be in the format of "
		      "cluster_name:gold_auth_key_file_path:"
		      "goldd_host:goldd_port "
		      "bad key file");
	if(!host)
		fatal("JobAcctLogfile should be in the format of "
		      "cluster_name:gold_auth_key_file_path:"
		      "goldd_host:goldd_port "
		      "bad host");
	if(!port) 
		fatal("JobAcctLogfile should be in the format of "
		      "cluster_name:gold_auth_key_file_path:"
		      "goldd_host:goldd_port "
		      "bad port");
	
	debug2("connecting from %s to gold with keyfile='%s' for %s(%d)",
	       cluster_name, keyfile, host, port);

	init_gold(cluster_name, keyfile, host, port);
		
	xfree(keyfile);
	xfree(host);

	return SLURM_SUCCESS;
}

int jobacct_p_fini_slurmctld()
{
	xfree(cluster_name);
	fini_gold();
	return SLURM_SUCCESS;
}

int jobacct_p_job_start_slurmctld(struct job_record *job_ptr)
{
	gold_object_t action = GOLD_ACTION_CREATE;
	
	if(_check_for_job(job_ptr->job_id, job_ptr->details->submit_time)) {
		error("It looks like this job is already in GOLD.  "
		      "This shouldn't happen, we are going to overwrite "
		      "old info.");
		action = GOLD_ACTION_MODIFY;
	}

	return _add_edit_job(job_ptr, action);
}

int jobacct_p_job_complete_slurmctld(struct job_record *job_ptr) 
{
	gold_object_t action = GOLD_ACTION_MODIFY;
	
	if(!_check_for_job(job_ptr->job_id, job_ptr->details->submit_time)) {
		error("Couldn't find this job entry.  "
		      "This shouldn't happen, we are going to create one.");
		action = GOLD_ACTION_CREATE;
	}

	return _add_edit_job(job_ptr, action);
}

int jobacct_p_step_start_slurmctld(struct step_record *step)
{
	gold_object_t action = GOLD_ACTION_MODIFY;
	
	if(!_check_for_job(step->job_ptr->job_id,
			   step->job_ptr->details->submit_time)) {
		error("Couldn't find this job entry.  "
		      "This shouldn't happen, we are going to create one.");
		action = GOLD_ACTION_CREATE;
	}

	return _add_edit_job(step->job_ptr, action);

}

int jobacct_p_step_complete_slurmctld(struct step_record *step)
{
	return SLURM_SUCCESS;	
}

int jobacct_p_suspend_slurmctld(struct job_record *job_ptr)
{
	return SLURM_SUCCESS;
}

int jobacct_p_startpoll(int frequency)
{
	info("jobacct GOLD plugin loaded");
	debug3("slurmd_jobacct_init() called");
	
	return SLURM_SUCCESS;
}

int jobacct_p_endpoll()
{
	return SLURM_SUCCESS;
}

int jobacct_p_set_proctrack_container_id(uint32_t id)
{
	return SLURM_SUCCESS;
}

int jobacct_p_add_task(pid_t pid, jobacct_id_t *jobacct_id)
{
	return SLURM_SUCCESS;
}

struct jobacctinfo *jobacct_p_stat_task(pid_t pid)
{
	return NULL;
}

struct jobacctinfo *jobacct_p_remove_task(pid_t pid)
{
	return NULL;
}

void jobacct_p_suspend_poll()
{
	return;
}

void jobacct_p_resume_poll()
{
	return;
}

#if (0)
/* If defined and FastSchedule=0 in slurm.conf, then report the CPU count that a 
 * node registers with rather than the CPU count defined for the node in slurm.conf */
#define SLURM_NODE_ACCT_REGISTER 1
#endif
#define DEBUG 1

extern void jobacct_p_node_down_slurmctld(struct node_record *node_ptr)
{
	char tmp[32];
	time_t now = time(NULL);
	uint16_t cpus;

	if (slurmctld_conf.fast_schedule)
		cpus = node_ptr->config_ptr->cpus;
	else
		cpus = node_ptr->cpus;

	slurm_make_time_str(&now, tmp, sizeof(tmp));
#if _DEBUG
	info("Node_acct_down: %s at %s with %u cpus due to %s", 
	     node_ptr->name, tmp, cpus, node_ptr->reason);
#endif
/* FIXME: WRITE TO DATABASE HERE */
}

extern void jobacct_p_node_all_down_slurmctld(char *reason)
{
	char *state_file, tmp[32];
	struct stat stat_buf;
	struct node_record *node_ptr;
	int i;

	state_file = xstrdup (slurmctld_conf.state_save_location);
	xstrcat (state_file, "/node_state");
	if (stat(state_file, &stat_buf)) {
		error("node_acct_all_down: could not stat(%s) to record "
		      "node down time", state_file);
		xfree(state_file);
		return;
	}
	xfree(state_file);

	slurm_make_time_str(&stat_buf.st_mtime, tmp, sizeof(tmp));
	node_ptr = node_record_table_ptr;
	for (i = 0; i < node_record_count; i++, node_ptr++) {
		if (node_ptr->name == '\0')
			continue;
		jobacct_p_node_down_slurmctld(node_ptr);
	}
}

extern void jobacct_p_node_up_slurmctld(struct node_record *node_ptr)
{
	char tmp[32];
	time_t now = time(NULL);

	slurm_make_time_str(&now, tmp, sizeof(tmp));
#if _DEBUG
	info("Node_acct_up: %s at %s", node_ptr->name, tmp);
#endif
	/* FIXME: WRITE TO DATABASE HERE */
}
extern void jobacct_p_cluster_procs(char *cluster_name, uint32_t procs)
{
	static uint32_t last_procs = 0;
	char tmp[32];
	time_t now = time(NULL);
	gold_request_t *gold_request = NULL;
	gold_response_t *gold_response = NULL;
	char tmp_buff[50];
	int rc = SLURM_ERROR;

	if (procs == last_procs)
		return;
	last_procs = procs;

	/* Record the processor count */
	slurm_make_time_str(&now, tmp, sizeof(tmp));
#if _DEBUG
	info("Node_acct_procs: %s has %u total CPUs at %s", 
	     cluster_name, procs, tmp);
#endif
	
	gold_request = create_gold_request(GOLD_OBJECT_EVENT,
					   GOLD_ACTION_MODIFY);
	if(!gold_request) 
		return;
	
	gold_request_add_condition(gold_request, "Machine", cluster_name);
	gold_request_add_condition(gold_request, "EndTime", "0");

	snprintf(tmp_buff, sizeof(tmp_buff), "%u", (int)now);
	gold_request_add_assignment(gold_request, "EndTime", tmp_buff);		
			
	gold_response = get_gold_response(gold_request);	
	destroy_gold_request(gold_request);

	if(!gold_response) {
		error("jobacct_p_cluster_procs: no response received");
		return;
	}

	if(!gold_response->rc) 
		rc = SLURM_SUCCESS;
	else {
		error("gold_response has non-zero rc(%d): %s",
		      gold_response->rc,
		      gold_response->message);
	}
	destroy_gold_response(gold_response);

	/* now add the new one */
	gold_request = create_gold_request(GOLD_OBJECT_EVENT,
					   GOLD_ACTION_CREATE);
	if(!gold_request) 
		return;
	
	gold_request_add_assignment(gold_request, "Machine", cluster_name);
	snprintf(tmp_buff, sizeof(tmp_buff), "%u", (int)now);
	gold_request_add_assignment(gold_request, "StartTime", tmp_buff);
	snprintf(tmp_buff, sizeof(tmp_buff), "%u", (int)procs);
	gold_request_add_assignment(gold_request, "CPUCount", tmp_buff);
			
	gold_response = get_gold_response(gold_request);	
	destroy_gold_request(gold_request);

	if(!gold_response) {
		error("jobacct_p_cluster_procs: no response received");
		return;
	}

	if(!gold_response->rc) 
		rc = SLURM_SUCCESS;
	else {
		error("gold_response has non-zero rc(%d): %s",
		      gold_response->rc,
		      gold_response->message);
	}
	destroy_gold_response(gold_response);
}
extern void jobacct_p_cluster_ready()
{
	uint32_t procs = 0;
	struct node_record *node_ptr;
	int i, j;
	char *cluster_name = NULL;

	node_ptr = node_record_table_ptr;
	for (i = 0; i < node_record_count; i++, node_ptr++) {
		if (node_ptr->name == '\0')
			continue;
		if (cluster_name == NULL) {
			cluster_name = xstrdup(node_ptr->name);
			for (j = 0; cluster_name[j]; j++) {
				if (isdigit(cluster_name[j])) {
					cluster_name[j] = '\0';
					break;
				}
			}
		}
#ifdef SLURM_NODE_ACCT_REGISTER
		if (slurmctld_conf.fast_schedule)
			procs += node_ptr->config_ptr->cpus;
		else
			procs += node_ptr->cpus;
#else
		procs += node_ptr->config_ptr->cpus;
#endif
	}

	jobacct_p_cluster_procs(cluster_name, procs);
	xfree(cluster_name);
}
