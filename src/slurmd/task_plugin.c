/*****************************************************************************\
 *  task_plugin.h - task launch plugin stub.
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#include <pthread.h>

#include "src/common/plugin.h"
#include "src/common/plugrack.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/slurmd/slurmd_job.h"

typedef struct slurmd_task_ops {
	int		(*pre_launch)		( slurmd_job_t *job );
	int		(*post_term)		( slurmd_job_t *job );
} slurmd_task_ops_t;


typedef struct slurmd_task_context {
	char			*task_type;
	plugrack_t		plugin_list;
	plugin_handle_t		cur_plugin;
	slurmd_task_ops_t	ops;
} slurmd_task_context_t;

static slurmd_task_context_t	*g_task_context = NULL;
static pthread_mutex_t		g_task_context_lock = PTHREAD_MUTEX_INITIALIZER;


static slurmd_task_ops_t *
_slurmd_task_get_ops( slurmd_task_context_t *c )
{
	/*
	 * Must be synchronized with slurmd_task_ops_t above.
	 */
	static const char *syms[] = {
		"task_pre_launch",
		"task_post_term",
	};
	int n_syms = sizeof( syms ) / sizeof( char * );

	/* Get plugin list. */
	if ( c->plugin_list == NULL ) {
		char *plugin_dir;
		c->plugin_list = plugrack_create();
		if ( c->plugin_list == NULL ) {
			error( "cannot create plugin manager" );
			return NULL;
		}
		plugrack_set_major_type( c->plugin_list, "task" );
		plugrack_set_paranoia( c->plugin_list,
				       PLUGRACK_PARANOIA_NONE,
				       0 );
		plugin_dir = slurm_get_plugin_dir();
		plugrack_read_dir( c->plugin_list, plugin_dir );
		xfree(plugin_dir);
	}

	c->cur_plugin = plugrack_use_by_type( c->plugin_list, c->task_type );
	if ( c->cur_plugin == PLUGIN_INVALID_HANDLE ) {
		error( "cannot find task plugin for %s", c->task_type );
		return NULL;
	}

	/* Dereference the API. */
	if ( plugin_get_syms( c->cur_plugin, n_syms, syms,
			(void **) &c->ops ) < n_syms ) {
		error( "incomplete task plugin detected" );
		return NULL;
	}

	return &c->ops;
}


static slurmd_task_context_t *
_slurmd_task_context_create( const char *task_plugin_type )
{
	slurmd_task_context_t *c;

	if ( task_plugin_type == NULL ) {
		debug3( "task_plugin_type is NULL" );
		return NULL;
	}

	c = xmalloc( sizeof( slurmd_task_context_t ) );
	c->task_type	= xstrdup( task_plugin_type );
	c->plugin_list	= NULL;
	c->cur_plugin	= PLUGIN_INVALID_HANDLE;

	return c;
}


static int
_slurmd_task_context_destroy( slurmd_task_context_t *c )
{
	/*
	 * Must check return code here because plugins might still
	 * be loaded and active.
	 */
	if ( c->plugin_list ) {
		if ( plugrack_destroy( c->plugin_list ) != SLURM_SUCCESS ) {
			return SLURM_ERROR;
		}
	}

	xfree( c->task_type );
	xfree( c );

	return SLURM_SUCCESS;
}


/*
 * Initialize the task plugin.
 *
 * RET - slurm error code
 */
extern int slurmd_task_init( void )
{
	int retval = SLURM_SUCCESS;
	char *task_plugin_type = NULL;
	
	slurm_mutex_lock( &g_task_context_lock );

	if ( g_task_context )
		goto done;

	task_plugin_type = slurm_get_task_plugin();
	g_task_context = _slurmd_task_context_create( task_plugin_type );
	if ( g_task_context == NULL ) {
		error( "cannot create task context for %s",
			 task_plugin_type );
		retval = SLURM_ERROR;
		goto done;
	}

	if ( _slurmd_task_get_ops( g_task_context ) == NULL ) {
		error( "cannot resolve task plugin operations" );
		_slurmd_task_context_destroy( g_task_context );
		g_task_context = NULL;
		retval = SLURM_ERROR;
	}

 done:
	slurm_mutex_unlock( &g_task_context_lock );
	xfree(task_plugin_type);
	return retval;
}

/*
 * Terminate the task plugin, free memory.
 *
 * RET - slurm error code
 */
extern int slurmd_task_fini( void )
{
	int rc;

	if (!g_task_context)
		return SLURM_SUCCESS;

	rc = _slurmd_task_context_destroy(g_task_context);
	g_task_context = NULL;
	return rc;
}

/*
 * Note that a task launch is about to occur.
 *
 * RET - slurm error code
 */
extern int pre_launch( slurmd_job_t *job )
{
	if ( slurmd_task_init() )
		return SLURM_ERROR;

	return (*(g_task_context->ops.pre_launch))(job);
}

/*
 * Note that a task has terminated.
 *
 * RET - slurm error code
 */
extern int post_term( slurmd_job_t *job )
{
	if ( slurmd_task_init() )
		return SLURM_ERROR;

	return (*(g_task_context->ops.post_term))(job);
}

