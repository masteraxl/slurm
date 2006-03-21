/*****************************************************************************\
 * src/slurmd/slurmstepd/slurmstepd.h - slurmstepd general header file
 * $Id$
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher J. Morrone <morrone2@llnl.gov>.
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

#ifndef _SLURMSTEPD_H
#define _SLURMSTEPD_H

#include "src/common/bitstring.h"

extern int slurmstepd_blocked_signals[];

typedef struct {
	pthread_cond_t cond;
	pthread_mutex_t lock;
	int rank;
	int parent_rank;
	slurm_addr parent_addr;
	int children;
	bitstr_t *bits;
} step_complete_t;

extern step_complete_t step_complete;

#endif /* !_SLURMSTEPD_H */
