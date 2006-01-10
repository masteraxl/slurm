/*****************************************************************************\
 * src/srun/allocate.h - node allocation functions for srun
 * $Id$
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

#ifndef _HAVE_ALLOCATE_H
#define _HAVE_ALLOCATE_H

#include <slurm/slurm.h>

#include "src/srun/srun_job.h"

/* 
 * Allocate nodes from the slurm controller -- retrying the attempt
 * if the controller appears to be down, and optionally waiting for
 * resources if none are currently available (see opt.immediate)
 *
 * Returns a pointer to a resource_allocation_response_msg which must
 * be freed with slurm_free_resource_allocation_response_msg()
 */
resource_allocation_response_msg_t * allocate_nodes(void);

/*
 * Test if an allocation would occur now given the job request.
 * Do not actually allocate resources
 */
int allocate_test(void);

/*
 * Create a job_desc_msg_t object, filled in from the current srun options
 * (see opt.h), if script != NULL then this is a batch job.
 * The resulting memory must be freed with  job_desc_msg_destroy()
 */
job_desc_msg_t * job_desc_msg_create_from_opts (char *script);

/* 
 * Destroy (free memory from) a job_desc_msg_t object allocated with
 * job_desc_msg_create()
 */
void job_desc_msg_destroy (job_desc_msg_t *j);

/*
 * Check for SLURM_JOBID environment variable, and if it is a valid
 * jobid, return a pseudo allocation response pointer.
 *
 * Returns NULL if SLURM_JOBID is not present or is invalid.
 */
resource_allocation_response_msg_t * existing_allocation(void);

/*
 * Return the jobid number stored in SLURM_JOBID env var
 *
 * Returns 0 if SLURM_JOBID is not set in current environment, or
 * is invalid.
 */
uint32_t jobid_from_env(void);

/*
 * Create a job step given the job information stored in 'j'
 * After returning, 'j' is filled in with information for job step.
 *
 * Returns -1 if job step creation failure, 0 otherwise
 */
int create_job_step(srun_job_t *j,
		    resource_allocation_response_msg_t *alloc_resp);

/* set the job for debugging purpose */
void set_allocate_job(srun_job_t *job);

#endif /* !_HAVE_ALLOCATE_H */
