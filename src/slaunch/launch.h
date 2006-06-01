/*****************************************************************************\
 * src/srun/launch.h - header for srun launch thread
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

#ifndef _HAVE_LAUNCH_H
#define _HAVE_LAUNCH_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif

#include "src/common/macros.h"
#include "src/common/slurm_protocol_api.h"

#include "src/srun/opt.h"
#include "src/srun/srun_job.h"

typedef struct launch_thr {
	pthread_t	thread;
	pthread_attr_t  attr;
	char            *host;	       /* name of host on which to run       */
	int             ntasks;	       /* number of tasks to initiate on host*/
	int 		*taskid;       /* list of global task ids            */
	int 		i;	       /* temporary index into array	     */
} launch_thr_t;

int    launch_thr_create(srun_job_t *job);
void * launch(void *arg);

#endif /* !_HAVE_LAUNCH_H */
