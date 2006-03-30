/*****************************************************************************\
 * src/slurmd/slurmstepd/mgr.c - job management functions for slurmstepd
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <mgrondona@llnl.gov>.
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
#ifndef _MGR_H
#define _MGR_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include "src/common/slurm_protocol_defs.h"

#include "src/slurmd/slurmd/slurmd.h"
#include "src/slurmd/slurmstepd/slurmstepd_job.h"

/*
 * Initialize a slurmd_job_t structure for a spawn task
 */
slurmd_job_t *mgr_spawn_task_setup(spawn_task_request_msg_t *msg,
				   slurm_addr *client, slurm_addr *self);

/*
 * Initialize a slurmd_job_t structure for a launch tasks
 */
slurmd_job_t *mgr_launch_tasks_setup(launch_tasks_request_msg_t *msg,
				     slurm_addr *client, slurm_addr *self);

/* 
 * Initialize a slurmd_job_t structure for a batch job
 */
slurmd_job_t *mgr_launch_batch_job_setup(batch_job_launch_msg_t *msg,
					 slurm_addr *client);

/*
 * Finalize a batch job.
 */
void mgr_launch_batch_job_cleanup(slurmd_job_t *job, int rc);

/*
 * Launch and manage the tasks in a job step.
 */
int job_manager(slurmd_job_t *job);

/*
 * with step completion add totals together.
 */
void aggregate_job_data(struct rusage rusage, int psize, int vsize);

/*
 * Register passwd entries so that we do not need to call initgroups(2)
 * frequently.
 */
extern void init_initgroups(int);


#endif
