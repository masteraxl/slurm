/*****************************************************************************\
 *  hold_wrapper.c - Hold all newly arriving jobs if there is a file
 *  "/etc/slurm.hold", otherwise use SLURM's internal scheduler.
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> et. al.
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

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <slurm/slurm_errno.h>

#include "src/common/plugin.h"
#include "src/common/log.h"

const char		plugin_name[]	= "SLURM Hold Scheduler plugin";
const char		plugin_type[]	= "sched/hold";
const uint32_t		plugin_version	= 90;

/* A plugin-global errno. */
static int plugin_errno = SLURM_SUCCESS;

/**************************************************************************/
/*  TAG(                              init                              ) */
/**************************************************************************/
int init( void )
{
	verbose( "Hold scheduler plugin loaded" );
	return SLURM_SUCCESS;
}

/**************************************************************************/
/*  TAG(                              fini                              ) */
/**************************************************************************/
void fini( void )
{
	/* Empty. */
}


/***************************************************************************/
/*  TAG(                   slurm_sched_plugin_schedule                   ) */
/***************************************************************************/
int
slurm_sched_plugin_schedule( void )
{
	verbose( "Hold plugin: schedule() is a NO-OP" );
	return SLURM_SUCCESS;
}


/**************************************************************************/
/* TAG(                   slurm_sched_plugin_initial_priority           ) */ 
/**************************************************************************/
u_int32_t
slurm_sched_plugin_initial_priority( u_int32_t max_prio )
{
	struct stat buf;

	if (stat("/etc/slurm.hold", &buf) == 0)
		return 0;	/* hold all new jobs */

	if (max_prio >= 2)
		return (max_prio - 1);
	else
		return 1;
}

/**************************************************************************/
/* TAG(              slurm_sched_plugin_job_is_pending                  ) */
/**************************************************************************/
void slurm_sched_plugin_job_is_pending( void )
{
	/* Empty. */
}

/**************************************************************************/
/* TAG(              slurm_sched_get_errno                              ) */
/**************************************************************************/
int slurm_sched_get_errno( void )
{
	return plugin_errno;
}

/**************************************************************************/
/* TAG(              slurm_sched_strerror                               ) */
/**************************************************************************/
char *slurm_sched_strerror( int errnum )
{
	return NULL;
}

