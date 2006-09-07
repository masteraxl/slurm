/*****************************************************************************\
 *  slurm_resource_info.c - Functions to determine number of available resources 
 *  $Id: slurm_resource_info.c,v 1.7 2006/08/31 20:11:21 palermo Exp $
 *****************************************************************************
 *  Copyright (C) 2006 Hewlett-Packard Development Company, L.P.
 *  Written by Susanne M. Balle, <susanne.balle@hp.com>
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
#if HAVE_CONFIG_H
#  include "config.h"
#endif

#if HAVE_STRING_H
#  include <string.h>
#endif

#include <sys/types.h>
#include "src/common/log.h"
#include <slurm/slurm.h>
#include "src/common/slurm_resource_info.h"

/*
 * slurm_get_avail_procs - Get the number of "available" cpus on a node
 *	given this number given the number of cpus_per_task and
 *	maximum sockets, cores, threads.  Note that the value of
 *	cpus is the lowest-level logical processor (LLLP).
 * IN mxsockets      - Job requested max sockets
 * IN mxcores        - Job requested max cores
 * IN mxthreads      - Job requested max threads
 * IN cpuspertask    - Job requested cpus per task
 * IN/OUT cpus       - Available cpu count
 * IN/OUT sockets    - Available socket count
 * IN/OUT cores      - Available core count
 * IN/OUT threads    - Available thread count
 * IN alloc_sockets  - Allocated socket count to other jobs
 * IN alloc_lps      - Allocated cpu count to other jobs
 * IN cr_type        - Consumable Resource type
 *
 * Note: used in both the select/{linear,cons_res} plugins.
 */
int slurm_get_avail_procs(const int mxsockets,
				 const int mxcores,
				 const int mxthreads,
				 const int cpuspertask,
				 int *cpus, 
				 int *sockets, 
				 int *cores, 
				 int *threads,
				 const int alloc_sockets,
				 const int alloc_lps,
				 const select_type_plugin_info_t cr_type)
{
	int avail_cpus = 0, max_cpus = 0;
	int max_sockets   = mxsockets;
	int max_cores     = mxcores;
	int max_threads   = mxthreads;
	int cpus_per_task = cpuspertask;

        /* pick defaults for any unspecified items */
	if (cpus_per_task <= 0)
		cpus_per_task = 1;
	if (max_sockets <= 0)
		max_sockets = INT_MAX;
	if (max_cores <= 0)
		max_cores = INT_MAX;
	if (max_threads <= 0)
		max_threads = INT_MAX;

	if (*threads <= 0)
	    	*threads = 1;
	if (*cores <= 0)
	    	*cores = 1;
	if (*sockets <= 0)
	    	*sockets = *cpus / *cores / *threads;
#if(0)
	info("SMB User_ sockets %d cores %d threads %d ", max_sockets, max_cores, max_threads);
	info("SMB HW_   sockets %d cores %d threads %d ", *sockets, *cores, *threads);
	info("SMB cr_type %d Allocated sockets %d lps %d ", cr_type, alloc_sockets, alloc_lps);
#endif
	if ((*threads <= 0) || (*cores <= 0) || (*sockets <= 0))
		fatal(" ((threads <= 0) || (cores <= 0) || (sockets <= 0))");
		
	switch(cr_type) { 
	case CR_SOCKET:
		*sockets -= alloc_sockets; /* sockets count */
		if (*sockets < 0) 
			fatal(" cons_res: *sockets < 0");
		
		*cpus     -= alloc_lps;
		if (*cpus < 0) 
			fatal(" cons_res: *cpus < 0");
		
		avail_cpus = (*cpus / cpus_per_task) * cpus_per_task;
		
		/*** honor socket/core/thread maximums ***/
		*sockets = MIN(*sockets, max_sockets);  /* socket count      */
		*cores   = MIN(*cores,   max_cores);    /* cores per socket  */
		*threads = MIN(*threads, max_threads);  /* threads per cores */
		
		max_cpus = *sockets * *cores * *threads;
		max_cpus *= cpus_per_task;
		
		avail_cpus = MIN(avail_cpus, max_cpus);
		break;
	case CR_CORE:
		/* No yet implemented */
		break;
	case CR_DEFAULT:
		/* no notion of socket, core, threads. Only one level
                   of logical processors */
		*cpus     -= alloc_lps;
		if (*cpus < 0) 
			fatal(" cons_res: *cpus < 0");
		
		avail_cpus = (*cpus / cpus_per_task) * cpus_per_task;
		break;
	default:
		/*** round down based on cpus_per_task ***/
		avail_cpus = (*cpus / cpus_per_task) * cpus_per_task;
		
		/*** honor socket/core/thread maximums ***/
		*sockets = MIN(*sockets, max_sockets);
		*cores   = MIN(*cores,   max_cores);
		*threads = MIN(*threads, max_threads);
		
		max_cpus = (*sockets * *cores * *threads);
		max_cpus = max_cpus * cpus_per_task;
		
		avail_cpus = MIN(avail_cpus, max_cpus);
		break;
	}

	return(avail_cpus);
}

/*
 * slurm_sprint_cpu_bind_type
 *
 * Given a cpu_bind_type, report all flag settings in str
 * IN  - cpu_bind_type
 * OUT - str
 */
void slurm_sprint_cpu_bind_type(char *str, cpu_bind_type_t cpu_bind_type)
{
    	if (!str)
		return;

	str[0] = '\0';

	if (cpu_bind_type & CPU_BIND_TO_THREADS)
		strcat(str, "threads,");
	if (cpu_bind_type & CPU_BIND_TO_CORES)
		strcat(str, "cores,");
	if (cpu_bind_type & CPU_BIND_TO_SOCKETS)
		strcat(str, "sockets,");
	if (cpu_bind_type & CPU_BIND_VERBOSE)
		strcat(str, "verbose,");
	if (cpu_bind_type & CPU_BIND_NONE)
		strcat(str, "none,");
	if (cpu_bind_type & CPU_BIND_RANK)
		strcat(str, "rank,");
	if (cpu_bind_type & CPU_BIND_MAP)
		strcat(str, "mapcpu,");
	if (cpu_bind_type & CPU_BIND_MASK)
		strcat(str, "maskcpu,");

	if (*str) {
		str[strlen(str)-1] = '\0';	/* remove trailing ',' */
	} else {
	    	strcat(str, "(null type)");	/* no bits set */
	}
}

/*
 * slurm_sprint_mem_bind_type
 *
 * Given a mem_bind_type, report all flag settings in str
 * IN  - mem_bind_type
 * OUT - str
 */
void slurm_sprint_mem_bind_type(char *str, mem_bind_type_t mem_bind_type)
{
    	if (!str)
		return;

	str[0] = '\0';

	if (mem_bind_type & MEM_BIND_VERBOSE)
		strcat(str, "verbose,");
	if (mem_bind_type & MEM_BIND_NONE)
		strcat(str, "none,");
	if (mem_bind_type & MEM_BIND_RANK)
		strcat(str, "rank,");
	if (mem_bind_type & MEM_BIND_LOCAL)
		strcat(str, "local,");
	if (mem_bind_type & MEM_BIND_MAP)
		strcat(str, "mapmem,");
	if (mem_bind_type & MEM_BIND_MASK)
		strcat(str, "maskmem,");

	if (*str) {
		str[strlen(str)-1] = '\0';	/* remove trailing ',' */
	} else {
	    	strcat(str, "(null type)");	/* no bits set */
	}
}
