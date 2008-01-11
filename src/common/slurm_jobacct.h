/*****************************************************************************\
 *  slurm_jobacct.h - implementation-independent job completion logging 
 *  API definitions
 *****************************************************************************
 *  Copyright (C) 2003 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette@llnl.com> et. al.
 *  UCRL-CODE-226842.
 *  
 *  Copyright (C) 2005 Hewlett-Packard Development Company, L.P.
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

/*****************************************************************************\
 *  Modification history
 *
 *  19 Jan 2005 by Andy Riebs <andy.riebs@hp.com>
 *       This file is derived from the file slurm_JOBACCT.c, written by
 *       Morris Jette, et al.
\*****************************************************************************/


#ifndef __SLURM_JOBACCT_H__
#define __SLURM_JOBACCT_H__

#if HAVE_STDINT_H
#  include <stdint.h>           /* for uint16_t, uint32_t definitions */
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>         /* for uint16_t, uint32_t definitions */
#endif
#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "src/slurmd/slurmstepd/slurmstepd_job.h"
#include "src/slurmctld/slurmctld.h"
#include "src/sacct/sacct_stat.h"

/* common */
extern int jobacct_init(void); /* load the plugin */
extern int jobacct_g_init_struct(jobacctinfo_t *jobacct, 
				 jobacct_id_t *jobacct_id);
/* must free jobacctinfo_t if not NULL */
extern jobacctinfo_t *jobacct_g_alloc(jobacct_id_t *jobacct_id);
extern void jobacct_g_free(jobacctinfo_t *jobacct);
extern int jobacct_g_setinfo(jobacctinfo_t *jobacct, 
			     enum jobacct_data_type type, void *data);
extern int jobacct_g_getinfo(jobacctinfo_t *jobacct, 
			     enum jobacct_data_type type, void *data);
extern void jobacct_g_aggregate(jobacctinfo_t *dest, jobacctinfo_t *from);
extern void jobacct_g_2_sacct(sacct_t *sacct, jobacctinfo_t *jobacct);
extern void jobacct_g_pack(jobacctinfo_t *jobacct, Buf buffer);
extern int jobacct_g_unpack(jobacctinfo_t **jobacct, Buf buffer);

/*functions used in slurmctld */
extern int jobacct_g_init_slurmctld(char *job_acct_log);
extern int jobacct_g_fini_slurmctld();
extern int jobacct_g_job_start_slurmctld(struct job_record *job_ptr);
extern int jobacct_g_job_complete_slurmctld(struct job_record *job_ptr); 
extern int jobacct_g_step_start_slurmctld(struct step_record *step);
extern int jobacct_g_step_complete_slurmctld(struct step_record *step);
extern int jobacct_g_suspend_slurmctld(struct job_record *job_ptr);

/*functions used in slurmstepd */
extern int jobacct_g_startpoll(int frequency);
extern int jobacct_g_endpoll();
extern int jobacct_g_set_proctrack_container_id(uint32_t id);
extern int jobacct_g_add_task(pid_t pid, jobacct_id_t *jobacct_id);
/* must free jobacctinfo_t if not NULL */
extern jobacctinfo_t *jobacct_g_stat_task(pid_t pid);
/* must free jobacctinfo_t if not NULL */
extern jobacctinfo_t *jobacct_g_remove_task(pid_t pid);
extern void jobacct_g_suspend_poll();
extern void jobacct_g_resume_poll();

/* functions only to be used by the gold plugin in 1.2 since we had to
 * make some mods to get the plugin to work with node states also
 */
extern void jobacct_g_node_all_down_slurmctld(char *reason);
extern void jobacct_g_node_down_slurmctld(struct node_record *node_ptr);
extern void jobacct_g_node_up_slurmctld(struct node_record *node_ptr);
extern void jobacct_g_cluster_procs(char *cluster_name, uint32_t procs);
extern void jobacct_g_cluster_ready();



#endif /*__SLURM_JOBACCT_H__*/

