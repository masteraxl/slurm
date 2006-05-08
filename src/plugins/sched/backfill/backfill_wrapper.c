/*****************************************************************************\
 *  backfill_wrapper.c - plugin for SLURM backfill scheduler.
 *  Operates like FIFO, but backfill scheduler daemon will explicitly modify
 *  the priority of jobs as needed to achieve backfill scheduling.
 *****************************************************************************
 *  Copyright (C) 2003 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Jay Windley <jwindley@lnxi.com>, Morris Jette <jette1@llnl.gov>
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

#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <slurm/slurm_errno.h>

#include "src/common/plugin.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "backfill.h"

const char		plugin_name[]	= "SLURM Backfill Scheduler plugin";
const char		plugin_type[]	= "sched/backfill";
const uint32_t		plugin_version	= 90;

/* A plugin-global errno. */
static int plugin_errno = SLURM_SUCCESS;

static pthread_t backfill_thread = 0;
static pthread_mutex_t thread_flag_mutex = PTHREAD_MUTEX_INITIALIZER;

/**************************************************************************/
/*  TAG(                              init                              ) */
/**************************************************************************/
int init( void )
{
#ifdef HAVE_BG
	/* Backfill scheduling on Blue Gene is possible, 
	 * but difficult and would require substantial 
	 * software development to accomplish. 
	 * It would need to consider each job's geometry, 
	 * ability to rotate, node-use (coprocessor or virtual)
	 * and conn-type (mesh, torus or nav). */
	fatal("Backfill scheduler incompatable with Blue Gene");
#else
	pthread_attr_t attr;

	verbose( "Backfill scheduler plugin loaded" );

	pthread_mutex_lock( &thread_flag_mutex );
	if ( backfill_thread ) {
		debug2( "Backfill thread already running, not starting another" );
		pthread_mutex_unlock( &thread_flag_mutex );
		return SLURM_ERROR;
	}

	slurm_attr_init( &attr );
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_DETACHED );
	if (pthread_create( &backfill_thread, &attr, backfill_agent, NULL))
		error("Unable to start backfill thread: %m");
	pthread_mutex_unlock( &thread_flag_mutex );
	slurm_attr_destroy( &attr );
#endif
	return SLURM_SUCCESS;
}

/**************************************************************************/
/*  TAG(                              fini                              ) */
/**************************************************************************/
static void _cancel_thread (pthread_t thread_id)
{
	int i;

	for (i=0; i<4; i++) {
		if (pthread_cancel(thread_id))
			return;
		usleep(1000);
	}
	error("Could not kill backfill sched pthread");
}

void fini( void )
{
	pthread_mutex_lock( &thread_flag_mutex );
	if ( backfill_thread ) {
		verbose( "Backfill scheduler plugin shutting down" );
		_cancel_thread( backfill_thread );
		backfill_thread = false;
	}
	pthread_mutex_unlock( &thread_flag_mutex );
}


/***************************************************************************/
/*  TAG(                   slurm_sched_plugin_schedule                   ) */
/***************************************************************************/
int
slurm_sched_plugin_schedule( void )
{
	return SLURM_SUCCESS;
}


/**************************************************************************/
/* TAG(                   slurm_sched_plugin_initial_priority           ) */
/**************************************************************************/
u_int32_t
slurm_sched_plugin_initial_priority( u_int32_t max_prio )
{
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
	run_backfill();
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

