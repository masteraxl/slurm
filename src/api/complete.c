/*****************************************************************************\
 *  complete.c - note the completion a slurm job or job step
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Moe Jette <jette1@llnl.gov>.
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <slurm/slurm.h>

#include "src/common/read_config.h"
#include "src/common/slurm_protocol_api.h"

/*
 * slurm_complete_job - note the completion of a job and all of its steps 
 * IN job_id - the job's id
 * IN job_return_code - the highest exit code of any task of the job
 * IN system_return_code - any slurm/system exit code
 * RET 0 on success or slurm error code
 */
int 
slurm_complete_job ( uint32_t job_id, uint32_t job_return_code,
                     uint32_t system_return_code )
{
	return slurm_complete_job_step ( job_id, NO_VAL, job_return_code, 
	                                 system_return_code);
}

/*
 * slurm_complete_job_step - note the completion of a specific job step 
 * IN job_id - the job's id
 * IN step_id - the job step's id or NO_VAL for all of the job's steps
 * IN job_return_code - the highest exit code of any task of the job
 * IN system_return_code - any slurm/system exit code
 * RET 0 on success or slurm error code
 */
int 
slurm_complete_job_step ( uint32_t job_id, uint32_t step_id, 
                          uint32_t job_return_code, 
                          uint32_t system_return_code )
{
	int rc;
	slurm_msg_t req_msg;
	complete_job_step_msg_t req;
	char host[128];

	(void) getnodename (host, sizeof(host)) ;

	req.job_id      = job_id;
	req.job_step_id	= step_id;
	req.job_rc      = job_return_code;
	req.slurm_rc    = system_return_code;
	req.node_name   = host;

	req_msg.msg_type= REQUEST_COMPLETE_JOB_STEP;
	req_msg.data	= &req;

	if (slurm_send_recv_controller_rc_msg(&req_msg, &rc) < 0)
	       return SLURM_ERROR;	
	
	if (rc)
		slurm_seterrno_ret(rc);

	return SLURM_PROTOCOL_SUCCESS;
}
