/*****************************************************************************\
 *  jobacct_common.h - common functions for almost all jobacct plugins.
 *****************************************************************************
 *
 *  Copyright (C) 2005 Hewlett-Packard Development Company, L.P.
 *  Written by Danny Auble, <da@llnl.gov>
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
 *
 *  This file is patterned after jobcomp_linux.c, written by Morris Jette and
 *  Copyright (C) 2002 The Regents of the University of California.
\*****************************************************************************/

#ifndef _HAVE_JOBACCT_COMMON_H
#define _HAVE_JOBACCT_COMMON_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#if HAVE_STDINT_H
#  include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#include <dirent.h>
#include <sys/stat.h>

#include "src/common/slurm_jobacct.h"
#include "src/common/xmalloc.h"
#include "src/common/list.h"
#include "src/common/xstring.h"
#include "src/common/node_select.h"

#include <ctype.h>

#define BUFFER_SIZE 4096

struct jobacctinfo {
	pid_t pid;
	struct rusage rusage; /* returned by wait3 */
	uint32_t max_vsize; /* max size of virtual memory */
	uint16_t max_vsize_task; /* contains which task number it was on */
	uint32_t tot_vsize; /* total virtual memory 
			       (used to figure out ave later) */
	uint32_t max_rss; /* max Resident Set Size */
	uint16_t max_rss_task; /* contains which task it was on */
	uint32_t tot_rss; /* total rss 
			     (used to figure out ave later) */
	uint32_t max_pages; /* max pages */
	uint16_t max_pages_task; /* contains which task it was on */
	uint32_t tot_pages; /* total pages
			     (used to figure out ave later) */ 
	uint32_t min_cpu; /* min cpu time */
	uint16_t min_cpu_task; /* contains which task it was on */
	uint32_t tot_cpu; /* total cpu time 
				 (used to figure out ave later) */
};

/* Define jobacctinfo_t below to avoid including extraneous slurm headers */
#ifndef __jobacctinfo_t_defined
#  define  __jobacctinfo_t_defined
   typedef struct jobacctinfo *jobacctinfo_t;     /* opaque data type */
#endif


/* in jobacct_common.c */
extern int common_init_struct(struct jobacctinfo *jobacct, uint16_t tid);
extern struct jobacctinfo *common_alloc_jobacct();
extern void common_free_jobacct(void *object);
extern int common_setinfo(struct jobacctinfo *jobacct, 
			  enum jobacct_data_type type, void *data);
extern int common_getinfo(struct jobacctinfo *jobacct, 
			  enum jobacct_data_type type, void *data);
extern void common_aggregate(struct jobacctinfo *dest, 
			     struct jobacctinfo *from);
extern void common_2_sacct(sacct_t *sacct, struct jobacctinfo *jobacct);
extern void common_pack(struct jobacctinfo *jobacct, Buf buffer);
extern int common_unpack(struct jobacctinfo **jobacct, Buf buffer);

/*in common_slurmctld.c */
extern int common_init_slurmctld(char *job_acct_log);
extern int common_fini_slurmctld();
extern int common_job_start_slurmctld(struct job_record *job_ptr);
extern int common_job_complete_slurmctld(struct job_record *job_ptr);
extern int common_step_start_slurmctld(struct step_record *step);
extern int common_step_complete_slurmctld(struct step_record *step);
extern int common_suspend_slurmctld(struct job_record *job_ptr);

/*in common slurmstepd.c */
extern int common_endpoll();
extern int common_add_task(pid_t pid, uint16_t tid);
extern struct jobacctinfo *common_stat_task(pid_t pid);
extern struct jobacctinfo *common_remove_task(pid_t pid);
extern void common_suspendpoll();

extern bool fini;
extern bool suspended;
extern List task_list;
extern pthread_mutex_t jobacct_lock;

#endif
