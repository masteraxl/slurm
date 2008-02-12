/*****************************************************************************\
 *  sacct_stat.c - stat slurmd for percise job information
 *
 *  $Id: options.c 7541 2006-03-18 01:44:58Z da $
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>.
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

#include "sacct.h"
#include <pthread.h>
#include "src/common/forward.h"
#include "src/common/slurm_auth.h"

jobacct_step_rec_t step;
	
int thr_finished = 0;
	
void *_stat_thread(void *args);
int _sacct_query(slurm_step_layout_t *step_layout, uint32_t job_id, 
		 uint32_t step_id);
int _process_results();


int _sacct_query(slurm_step_layout_t *step_layout, uint32_t job_id,
		 uint32_t step_id)
{
	slurm_msg_t msg;
	stat_jobacct_msg_t r;
	stat_jobacct_msg_t *jobacct_msg = NULL;
	ListIterator itr;
	List ret_list = NULL;
	sacct_t temp_sacct;
	ret_data_info_t *ret_data_info = NULL;
	int rc = SLURM_SUCCESS;
	int ntasks = 0;

	debug("getting the stat of job %d on %d nodes", 
	      job_id, step_layout->node_cnt);

	memset(&temp_sacct, 0, sizeof(sacct_t));
	temp_sacct.min_cpu = (float)NO_VAL;
	memset(&step.sacct, 0, sizeof(sacct_t));
	step.sacct.min_cpu = (float)NO_VAL;

	step.header.jobnum = job_id;
	step.header.partition = NULL;
#ifdef HAVE_BG
	step.header.blockid = NULL;
#endif
	step.stepnum = step_id;
	step.nodes = step_layout->node_list;
	step.stepname = NULL;
	step.status = JOB_RUNNING;
	step.ntasks = 0;
	slurm_msg_t_init(&msg);
	/* Common message contents */
	r.job_id      = job_id;
	r.step_id     = step_id;
	r.jobacct     = jobacct_gather_g_create(NULL);
	msg.msg_type        = MESSAGE_STAT_JOBACCT;
	msg.data            = &r;
	
	
	ret_list = slurm_send_recv_msgs(step_layout->node_list, &msg, 0);
	if (!ret_list) {
		error("got an error no list returned");
		goto cleanup;
	}
	
	itr = list_iterator_create(ret_list);		
	while((ret_data_info = list_next(itr))) {
		switch (ret_data_info->type) {
		case MESSAGE_STAT_JOBACCT:
			jobacct_msg = (stat_jobacct_msg_t *)
				ret_data_info->data;
			if(jobacct_msg) {
				debug2("got it back for job %d", 
				       jobacct_msg->job_id);
				jobacct_gather_g_2_sacct(
					&temp_sacct, 
					jobacct_msg->jobacct);
				ntasks += jobacct_msg->num_tasks;
				aggregate_sacct(&step.sacct, &temp_sacct);
			}
			break;
		case RESPONSE_SLURM_RC:
			rc = slurm_get_return_code(ret_data_info->type, 
						   ret_data_info->data);
			error("there was an error with the request rc = %s", 
			      slurm_strerror(rc));
			break;
		default:
			rc = slurm_get_return_code(ret_data_info->type, 
						   ret_data_info->data);
			error("unknown return given %d rc = %s", 
			      ret_data_info->type, slurm_strerror(rc));
			break;
		}
	}
	list_iterator_destroy(itr);
	list_destroy(ret_list);

	step.ntasks += ntasks;		
cleanup:
	
	if(step.ntasks) {
		step.sacct.ave_rss *= 1024;
		step.sacct.max_rss *= 1024;
		step.sacct.ave_vsize *= 1024;
		step.sacct.max_vsize *= 1024;

		step.sacct.ave_cpu /= step.ntasks;
		step.sacct.ave_cpu /= 100;
		step.sacct.min_cpu /= 100;
		step.sacct.ave_rss /= step.ntasks;
		step.sacct.ave_vsize /= step.ntasks;
		step.sacct.ave_pages /= step.ntasks;
	}
	jobacct_gather_g_destroy(r.jobacct);	
	return SLURM_SUCCESS;
}

int _process_results()
{
	print_fields(JOBSTEP, &step);
	return SLURM_SUCCESS;
}

int sacct_stat(uint32_t jobid, uint32_t stepid)
{
	slurm_msg_t req_msg;
	slurm_msg_t resp_msg;
	job_step_id_msg_t req;
	slurm_step_layout_t *step_layout = NULL;
	int rc = SLURM_SUCCESS;

	slurm_msg_t_init(&req_msg);
	slurm_msg_t_init(&resp_msg);
	debug("requesting info for job %u.%u", jobid, stepid);
	req.job_id = jobid;
	req.step_id = stepid;
	req_msg.msg_type = REQUEST_STEP_LAYOUT;
	req_msg.data     = &req;
	
	if (slurm_send_recv_controller_msg(&req_msg, &resp_msg) < 0) {
		return SLURM_ERROR;
	}
		
	switch (resp_msg.msg_type) {
	case RESPONSE_STEP_LAYOUT:
		step_layout = (slurm_step_layout_t *)resp_msg.data;
		break;
	case RESPONSE_SLURM_RC:
		rc = ((return_code_msg_t *) resp_msg.data)->return_code;
		slurm_free_return_code_msg(resp_msg.data);	
		printf("problem getting job: %s\n", slurm_strerror(rc));
		slurm_seterrno_ret(rc);
		break;
	default:
		slurm_seterrno_ret(SLURM_UNEXPECTED_MSG_ERROR);
		break;
	}
		
	if(!step_layout) {
		error("didn't get the job record rc = %s", slurm_strerror(rc));
		return rc;
	}

	_sacct_query(step_layout, jobid, stepid);
	
	_process_results();
	
	slurm_step_layout_destroy(step_layout);	
	
	return rc;
}
