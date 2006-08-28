/*****************************************************************************\
 *  sched_wiki.c - Wiki plugin for Moab and Maui schedulers.
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#include <stdio.h>
#include <slurm/slurm_errno.h>

#include "src/common/plugin.h"
#include "src/common/log.h"
#include "./msg.h"

const char		plugin_name[]	= "Wiki (Maui and Moab) Scheduler plugin";
const char		plugin_type[]	= "sched/wiki2";
const uint32_t		plugin_version	= 90;

/* A plugin-global errno. */
static int plugin_errno = SLURM_SUCCESS;

/**************************************************************************/
/*  TAG(                              init                              ) */
/**************************************************************************/
extern int init( void )
{
	verbose( "Wiki scheduler plugin loaded" );
	return spawn_msg_thread();
}

/**************************************************************************/
/*  TAG(                              fini                              ) */
/**************************************************************************/
extern void fini( void )
{
	term_msg_thread();
}


/***************************************************************************/
/*  TAG(                   slurm_sched_plugin_schedule                   ) */
/***************************************************************************/
extern int slurm_sched_plugin_schedule( void )
{
	/* No action required */
	return SLURM_SUCCESS;
}


/**************************************************************************/
/* TAG(                   slurm_sched_plugin_initial_priority           ) */ 
/**************************************************************************/
extern uint32_t slurm_sched_plugin_initial_priority( uint32_t last_prio )
{
	if (init_prio_mode == PRIO_DECREMENT) {
		if (last_prio >= 2)
			return (last_prio - 1);
		else
			return 1;
	} else 
		return 0;
}

/**************************************************************************/
/* TAG(              slurm_sched_plugin_job_is_pending                  ) */
/**************************************************************************/
void slurm_sched_plugin_job_is_pending( void )
{
	/* No action required */
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

