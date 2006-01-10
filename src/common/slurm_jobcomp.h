/*****************************************************************************\
 *  slurm_jobcomp.h - implementation-independent job completion logging 
 *  API definitions
 *****************************************************************************
 *  Copyright (C) 2003 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.com> et. al.
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

#ifndef __SLURM_JOBCOMP_H__
#define __SLURM_JOBCOMP_H__

#if HAVE_CONFIG_H
#  include "config.h"
#endif
#if HAVE_STDINT_H
#  include <stdint.h>           /* for uint16_t, uint32_t definitions */
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>         /* for uint16_t, uint32_t definitions */
#endif
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "src/slurmctld/slurmctld.h"

typedef struct slurm_jobcomp_context * slurm_jobcomp_context_t;

/* initialization of job completion logging */
extern int g_slurm_jobcomp_init(char *jobcomp_loc);

/* terminate pthreads and free, general clean-up for termination */
extern int g_slurm_jobcomp_fini(void);

/* write record of a job's completion */
extern int g_slurm_jobcomp_write(struct job_record *job_ptr);

/* return error code */
extern int g_slurm_jobcomp_errno(void);

/* convert job completion logger specific error code to a string */
extern char *g_slurm_jobcomp_strerror(int errnum);

#endif /*__SLURM_JOBCOMP_H__*/

