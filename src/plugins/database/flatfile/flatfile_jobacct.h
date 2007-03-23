/*****************************************************************************\
 *  flatfile_jobacct.h - functions the flatfile jobacct database.
 *****************************************************************************
 *
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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
 *
 *  This file is patterned after jobcomp_linux.c, written by Morris Jette and
 *  Copyright (C) 2002 The Regents of the University of California.
\*****************************************************************************/

#ifndef _HAVE_FLATFILE_JOBACCT_H
#define _HAVE_FLATFILE_JOBACCT_H

#include "src/plugins/jobacct/common/jobacct_common.h"
#include "src/slurmctld/slurmctld.h"

extern int flatfile_jobacct_init();
extern int flatfile_jobacct_fini();
extern int flatfile_jobacct_job_start(struct job_record *job_ptr);
extern int flatfile_jobacct_job_complete(struct job_record *job_ptr);
extern int flatfile_jobacct_step_start(struct step_record *step_ptr);
extern int flatfile_jobacct_step_complete(struct step_record *step_ptr);
extern int flatfile_jobacct_suspend(struct job_record *job_ptr);
extern List flatfile_jobacct_getdata(List selected_steps, List selected_parts,
				     void *params);
extern void flatfile_jobacct_do_expire(List selected_parts, void *params);

#endif
