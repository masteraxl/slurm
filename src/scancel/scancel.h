/*****************************************************************************\
 *  scancel.h - definitions for scancel data structures and functions
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette<jette1@llnl.gov>, et. al.
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

#ifndef _HAVE_SCANCEL_H
#define _HAVE_SCANCEL_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

typedef struct scancel_options {
	bool batch;		/* --batch, -b			*/
	bool interactive;	/* --interactive, -i		*/
	char *job_name;		/* --name=n, -nn		*/
	char *partition;	/* --partition=n, -pn		*/
	uint16_t signal;	/* --signal=n, -sn		*/
	enum job_states state;	/* --state=n, -tn		*/
	uid_t user_id;		/* --user=n, -un		*/
	char *user_name;	/* --user=n, -un		*/
	int verbose;		/* --verbose, -v		*/

	uint16_t job_cnt;	/* count of job_id's specified	*/
	uint32_t *job_id;	/* list of job_id's		*/
	uint32_t *step_id;	/* list of job step id's	*/
} opt_t;

opt_t opt;

/* process options:
 * 1. set defaults
 * 2. update options with env vars
 * 3. update options with commandline args
 * 4. perform some verification that options are reasonable
 */
int initialize_and_process_args(int argc, char *argv[]);

#endif	/* _HAVE_SCANCEL_H */
