/*****************************************************************************\
 *  job_step_info.c - get/print the job step state information of slurm
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>, 
 *             Joey Ekstrom <ekstrom1@llnl.gov>,  et. al.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <slurm/slurm.h>

#include "src/api/job_info.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"

#define HUGE_BUF 4096
/*
 * slurm_print_job_step_info_msg - output information about all Slurm 
 *	job steps based upon message as loaded using slurm_get_job_steps
 * IN out - file to write to
 * IN job_step_info_msg_ptr - job step information message pointer
 * IN one_liner - print as a single line if true
 */
void 
slurm_print_job_step_info_msg ( FILE* out, 
		job_step_info_response_msg_t * job_step_info_msg_ptr, 
		int one_liner )
{
	int i;
	job_step_info_t *job_step_ptr = job_step_info_msg_ptr->job_steps ;
	char time_str[32];

	slurm_make_time_str ((time_t *)&job_step_info_msg_ptr->last_update, 
			time_str, sizeof(time_str));
	fprintf( out, "Job step data as of %s, record count %d\n",
		time_str, job_step_info_msg_ptr->job_step_count);

	for (i = 0; i < job_step_info_msg_ptr-> job_step_count; i++) 
	{
		slurm_print_job_step_info ( out, & job_step_ptr[i], 
					    one_liner ) ;
	}
}

/*
 * slurm_print_job_step_info - output information about a specific Slurm 
 *	job step based upon message as loaded using slurm_get_job_steps
 * IN out - file to write to
 * IN job_ptr - an individual job step information record pointer
 * IN one_liner - print as a single line if true
 */
void
slurm_print_job_step_info ( FILE* out, job_step_info_t * job_step_ptr, 
			    int one_liner )
{
	char *print_this = slurm_sprint_job_step_info(job_step_ptr, one_liner);
	fprintf ( out, "%s", print_this);
	xfree(print_this);
}

/*
 * slurm_sprint_job_step_info - output information about a specific Slurm 
 *	job step based upon message as loaded using slurm_get_job_steps
 * IN job_ptr - an individual job step information record pointer
 * IN one_liner - print as a single line if true
 * RET out - char * containing formatted output (must be freed after call)
 *           NULL is returned on failure.
 */
char *
slurm_sprint_job_step_info ( job_step_info_t * job_step_ptr, 
			    int one_liner )
{
	char time_str[32];
	int j;
	char tmp2[100];
	char *out = xmalloc(HUGE_BUF);
	int len = 0;
	int lentmp2 = 0;
	int tmplen = 0;

	/****** Line 1 ******/
	slurm_make_time_str ((time_t *)&job_step_ptr->start_time, time_str,
		sizeof(time_str));
	sprintf ( tmp2, "StepId=%u.%u UserId=%u Tasks=%u StartTime=%s", 
		job_step_ptr->job_id, job_step_ptr->step_id, 
		job_step_ptr->user_id, job_step_ptr->num_tasks, time_str);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out, "%s", tmp2);
	len += lentmp2;

	if (one_liner)
		sprintf ( tmp2, " ");
	else
		sprintf ( tmp2, "\n   ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;


	/****** Line 2 ******/
	sprintf ( tmp2, "Partition=%s Nodes=%s Name=%s Network=%s\n\n", 
		job_step_ptr->partition, job_step_ptr->nodes,
		job_step_ptr->name, job_step_ptr->network);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;

	return out;
}

/*
 * slurm_get_job_steps - issue RPC to get specific slurm job step   
 *	configuration information if changed since update_time.
 *	a job_id value of zero implies all jobs, a step_id value of 
 *	zero implies all steps
 * IN update_time - time of current configuration data
 * IN job_id - get information for specific job id, zero for all jobs
 * IN step_id - get information for specific job step id, zero for all 
 *	job steps
 * IN job_info_msg_pptr - place to store a job configuration pointer
 * IN show_flags - job step filtering options
 * RET 0 on success, otherwise return -1 and set errno to indicate the error
 * NOTE: free the response using slurm_free_job_step_info_response_msg
 */
int
slurm_get_job_steps (time_t update_time, uint32_t job_id, uint32_t step_id, 
		     job_step_info_response_msg_t **resp, uint16_t show_flags)
{
	int rc;
	slurm_msg_t req_msg;
	slurm_msg_t resp_msg;
	job_step_info_request_msg_t req;

	req.last_update  = update_time;
	req.job_id	= job_id;
	req.step_id	= step_id;
	req.show_flags	= show_flags;
	req_msg.msg_type = REQUEST_JOB_STEP_INFO;
	req_msg.data	= &req;

	if (slurm_send_recv_controller_msg(&req_msg, &resp_msg) < 0)
		return SLURM_ERROR;

	switch (resp_msg.msg_type) {
	case RESPONSE_JOB_STEP_INFO:
		*resp = (job_step_info_response_msg_t *) resp_msg.data;
		break;
	case RESPONSE_SLURM_RC:
		rc = ((return_code_msg_t *) resp_msg.data)->return_code;
		slurm_free_return_code_msg(resp_msg.data);	
		if (rc) 
			slurm_seterrno_ret(rc);
		*resp = NULL;
		break;
	default:
		slurm_seterrno_ret(SLURM_UNEXPECTED_MSG_ERROR);
		break;
	}

	return SLURM_PROTOCOL_SUCCESS;
}

