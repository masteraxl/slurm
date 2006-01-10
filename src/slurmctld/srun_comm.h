/*****************************************************************************\
 *  srun_comm.h - definitions srun communications
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette@llnl.gov> et. al.
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

#ifndef _HAVE_SRUN_COMM_H
#define _HAVE_SRUN_COMM_H

#include <sys/types.h>
#include <time.h>

/*
 * srun_allocate - notify srun of a resource allocation
 * IN job_id - id of the job allocated resource
 */
extern void srun_allocate (uint32_t job_id);

/*
 * srun_node_fail - notify srun of a node's failure
 * IN job_id    - id of job to notify
 * IN node_name - name of failed node
 */
extern void srun_node_fail (uint32_t job_id, char *node_name);

/* srun_ping - ping all srun commands that have not been heard from recently */
extern void srun_ping (void);

/*
 * srun_response - note that srun has responded
 * IN job_id  - id of job responding
 * IN step_id - id of step responding or NO_VAL if not a step
 */
extern void srun_response(uint32_t job_id, uint32_t step_id);

/*
 * srun_timeout - notify srun of a job's timeout
 * IN job_id  - if of job to notify
 * IN timeout - when job is scheduled to be killed
 */
extern void srun_timeout (uint32_t job_id, time_t timeout);

#endif /* !_HAVE_SRUN_COMM_H */
