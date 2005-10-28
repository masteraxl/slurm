/*****************************************************************************\
 * src/slurmd/run_script.h - code shared between slurmd and slurmd_step
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

#ifndef _RUN_SCRIPT_H
#define _RUN_SCRIPT_H

#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>

/*
 * Run a prolog or epilog script
 * name IN: class of program (prolog, epilog, etc.),
 *	if prefix is "user" then also set uid
 * path IN: pathname of program to run
 * jobid, uidIN: info on associated job
 * max_wait IN: maximum time to wait in seconds, -1 for no limit
 * env IN: environment variables to use on exec, sets minimal environment 
 *	if NULL
 * RET 0 on success, -1 on failure.
 */
int run_script(const char *name, const char *path, uint32_t jobid, 
	       uid_t uid, int max_wait, char **env);

#endif /* _RUN_SCRIPT_H */
