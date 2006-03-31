/*****************************************************************************\
 *  jobacct.h - process and record information about process accountablity
 *
 *  $Id: jobacct.h 7620 2006-03-29 17:42:21Z da $
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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
#ifndef _HAVE_JOBACCT_H
#define _HAVE_JOBACCT_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <slurm/slurm_errno.h>
#include <sys/stat.h>
#include "src/common/xstring.h"
#include "slurmctld.h"

int jobacct_init(char *job_acct_log);
int jobacct_job_start(struct job_record *job_ptr);
int jobacct_step_start(struct step_record *step);
int jobacct_step_complete(struct step_record *step);
int jobacct_job_complete(struct job_record *job_ptr);


#endif /* _HAVE_JOBACCT_H */
