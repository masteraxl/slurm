/*****************************************************************************\
 *  src/slurmd/step_msg_api.h - slurmd_step message API
 *  $Id: $
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher Morrone <morrone2@llnl.gov>
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

#ifndef _STEP_MSG_API_H
#define _STEP_MSG_API_H

#include <inttypes.h>

typedef enum {
	REQUEST_SIGNAL = 0,
	REQUEST_TERMINATE,
	REQUEST_STATUS
} step_msg_t;

typedef struct step_location {
	uint32_t jobid;
	uint32_t stepid;
	char *nodename;
	char *directory;
} step_loc_t;

int step_request_status(step_loc_t step);
int step_request_signal(step_loc_t step, int signal);

#endif /* _STEP_MSG_API_H */
