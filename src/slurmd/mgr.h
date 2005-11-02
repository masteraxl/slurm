/*****************************************************************************\
 * src/slurmd/mgr.c - job management functions for slurmd
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <mgrondona@llnl.gov>.
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
#ifndef _MGR_H
#define _MGR_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include "src/common/slurm_protocol_defs.h"

#include "src/slurmd/slurmd_job.h"
#include "src/slurmd/slurmd.h"

/* Spawn a task / job step on this node
 */
int mgr_spawn_task(spawn_task_request_msg_t *msg, slurm_addr *client,
		   slurm_addr *self);

/* Launch a job step on this node
 */
int mgr_launch_tasks(launch_tasks_request_msg_t *msg, slurm_addr *client,
		     slurm_addr *self);

/* 
 * Launch batch script on this node
 */
int mgr_launch_batch_job(batch_job_launch_msg_t *msg, slurm_addr *client);

#endif
