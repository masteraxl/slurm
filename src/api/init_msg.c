/*****************************************************************************\
 *  init_msg.c - initialize RPC messages contents
 *****************************************************************************
 *  Copyright (C) 2002-2007 The Regents of the University of California.
 *  Copyright (C) 2008-2009 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>.
 *  LLNL-CODE-402394.
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

#include <errno.h>
#include <stdio.h>

#include <slurm/slurm.h>

#include "src/common/slurm_protocol_api.h"
#include "src/common/forward.h"

/*
 * slurm_init_job_desc_msg - initialize job descriptor with 
 *	default values 
 * IN/OUT job_desc_msg - user defined job descriptor
 */
void slurm_init_job_desc_msg(job_desc_msg_t * job_desc_msg)
{
	job_desc_msg->account		= NULL;
	job_desc_msg->acctg_freq	= (uint16_t) NO_VAL;
	job_desc_msg->alloc_node	= NULL;
	job_desc_msg->alloc_resp_port	= 0;
	job_desc_msg->alloc_sid		= NO_VAL;
	job_desc_msg->argc		= 0;
	job_desc_msg->argv		= ((char **) NULL);
	job_desc_msg->begin_time	= 0;
	job_desc_msg->blrtsimage	= NULL;
	job_desc_msg->comment		= NULL;
	job_desc_msg->conn_type		= (uint16_t) NO_VAL;
	job_desc_msg->contiguous	= (uint16_t) NO_VAL;
	job_desc_msg->cpu_bind		= NULL;
	job_desc_msg->cpu_bind_type	= (uint16_t) NO_VAL;
	job_desc_msg->cpus_per_task	= (uint16_t) NO_VAL;
	job_desc_msg->dependency	= NULL;
	job_desc_msg->environment	= ((char **) NULL);
	job_desc_msg->env_size		= 0;
	job_desc_msg->err		= NULL;
	job_desc_msg->exc_nodes		= NULL;
	job_desc_msg->features		= NULL;
#if SYSTEM_DIMENSIONS
{
	int i;
	for (i=0; i<SYSTEM_DIMENSIONS; i++)
		job_desc_msg->geometry[i] = (uint16_t) NO_VAL;
}
#endif
	job_desc_msg->group_id		= NO_VAL;
	job_desc_msg->immediate		= 0;
	job_desc_msg->in		= NULL;
	job_desc_msg->job_id		= NO_VAL;
	job_desc_msg->job_min_cores	= (uint16_t) NO_VAL;
	job_desc_msg->job_min_procs	= (uint16_t) NO_VAL;
	job_desc_msg->job_min_sockets	= (uint16_t) NO_VAL;
	job_desc_msg->job_min_threads	= (uint16_t) NO_VAL;
	job_desc_msg->job_min_memory    = NO_VAL;
	job_desc_msg->job_min_tmp_disk  = NO_VAL;
	job_desc_msg->kill_on_node_fail = (uint16_t) NO_VAL;
	job_desc_msg->licenses		= NULL;
	job_desc_msg->linuximage	= NULL;
	job_desc_msg->mail_type		= 0;
	job_desc_msg->mail_user		= NULL;
	job_desc_msg->max_cores		= (uint16_t) NO_VAL;
	job_desc_msg->max_nodes		= NO_VAL;
	job_desc_msg->max_sockets	= (uint16_t) NO_VAL;
	job_desc_msg->max_threads	= (uint16_t) NO_VAL;
	job_desc_msg->mem_bind		= NULL;
	job_desc_msg->mem_bind_type	= (uint16_t) NO_VAL;
	job_desc_msg->min_cores		= (uint16_t) NO_VAL;
	job_desc_msg->min_nodes		= NO_VAL;
	job_desc_msg->min_sockets	= (uint16_t) NO_VAL;
	job_desc_msg->min_threads	= (uint16_t) NO_VAL;
	job_desc_msg->mloaderimage	= NULL;
	job_desc_msg->name		= NULL;
	job_desc_msg->network		= NULL;
	job_desc_msg->nice		= NICE_OFFSET;
	job_desc_msg->ntasks_per_core	= (uint16_t) NO_VAL;
	job_desc_msg->ntasks_per_node	= (uint16_t) NO_VAL;
	job_desc_msg->ntasks_per_socket	= (uint16_t) NO_VAL;
	job_desc_msg->num_procs		= NO_VAL;
	job_desc_msg->num_tasks		= NO_VAL;
	job_desc_msg->open_mode		= 0;	/* system default */
	job_desc_msg->other_port	= 0;
	job_desc_msg->out		= NULL;
	job_desc_msg->overcommit	= (uint8_t) NO_VAL;
	job_desc_msg->partition		= NULL;
	job_desc_msg->plane_size	= (uint16_t) NO_VAL;
	job_desc_msg->priority		= NO_VAL;
	job_desc_msg->ramdiskimage	= NULL;
	job_desc_msg->reboot		= (uint16_t) NO_VAL;
	job_desc_msg->resp_host		= NULL;
	job_desc_msg->req_nodes		= NULL;
	job_desc_msg->requeue		= (uint16_t) NO_VAL;
	job_desc_msg->reservation	= NULL;
	job_desc_msg->rotate		= (uint16_t) NO_VAL;
	job_desc_msg->script		= NULL;
	job_desc_msg->select_jobinfo	= NULL;
	job_desc_msg->shared		= (uint16_t) NO_VAL;
	job_desc_msg->task_dist		= (uint16_t) NO_VAL;
	job_desc_msg->time_limit	= NO_VAL;
	job_desc_msg->user_id		= NO_VAL;
	job_desc_msg->wckey		= NULL;
	job_desc_msg->work_dir		= NULL;
}

/*
 * slurm_init_part_desc_msg - initialize partition descriptor with 
 *	default values 
 * IN/OUT update_part_msg - user defined partition descriptor
 */
void slurm_init_part_desc_msg (update_part_msg_t * update_part_msg)
{
	update_part_msg->name 		= NULL;
	update_part_msg->nodes 		= NULL;
	update_part_msg->allow_groups 	= NULL;
	update_part_msg->max_time 	= (uint32_t) NO_VAL;
	update_part_msg->max_nodes 	= NO_VAL;
	update_part_msg->min_nodes 	= NO_VAL;
	update_part_msg->hidden 	= (uint16_t) NO_VAL;
	update_part_msg->default_part 	= (uint16_t) NO_VAL;
	update_part_msg->root_only 	= (uint16_t) NO_VAL;
	update_part_msg->max_share 	= (uint16_t) NO_VAL;
	update_part_msg->priority 	= (uint16_t) NO_VAL;
	update_part_msg->state_up 	= (uint16_t) NO_VAL;
}

/*
 * slurm_init_resv_desc_msg - initialize reservation descriptor with 
 *	default values 
 * OUT job_desc_msg - user defined partition descriptor
 */
void slurm_init_resv_desc_msg (reserve_request_msg_t * update_resv_msg)
{
	update_resv_msg->name		= NULL;
	update_resv_msg->start_time	= (time_t) NO_VAL;
	update_resv_msg->end_time	= (time_t) NO_VAL;
	update_resv_msg->duration	= NO_VAL;
	update_resv_msg->type		= 0;
	update_resv_msg->node_cnt	= NO_VAL;
	update_resv_msg->node_list	= NULL;
	update_resv_msg->features	= NULL;
	update_resv_msg->partition	= NULL;
	update_resv_msg->users		= NULL;
	update_resv_msg->accounts	= NULL;
}
