/*****************************************************************************\
 **  mpich1.h - Library routines for initiating jobs on with mpich.
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
 *  UCRL-CODE-226842.
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

#ifndef _HAVE_MPI_MPICH1_H
#define _HAVE_MPI_MPICH1_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include "src/common/slurm_xlator.h"
#include "src/common/mpi.h"
#include "src/common/env.h"

#define MPI_SMP       1
#define MPI_MVAPICH   2

#define MPI_ARCH      MPI_SMP

#if (MPI_ARCH == MPI_MVAPICH)
typedef struct mpich1_state mpich1_state_t;
extern mpich1_state_t *mpich1_thr_create(const mpi_plugin_client_info_t *job,
extern int mpich1_thr_destroy(mpich1_state_t *state);
#endif

#endif	/* !_HAVE_MPI_MPICH1_H */
