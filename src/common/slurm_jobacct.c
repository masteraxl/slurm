/*****************************************************************************\
 *  slurm_jobacct.c - implementation-independent job accounting logging 
 *  functions
 *****************************************************************************
 *  Copyright (C) 2003 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Jay Windley <jwindley@lnxi.com>, Morris Jette <jette1@llnl.com>
 *  UCRL-CODE-217948.
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
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

/*****************************************************************************\
 *  Modification history
 *
 *  19 Jan 2005 by Andy Riebs <andy.riebs@hp.com>
 *  	 This file is derived from the file slurm_jobcomp.c, written by
 *  	 Morris Jette, et al.
\*****************************************************************************/

#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "src/common/macros.h"
#include "src/common/plugin.h"
#include "src/common/plugrack.h"
#include "src/common/slurm_jobacct.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/slurmd/slurmstepd/slurmstepd_job.h"

/*
 * WARNING:  Do not change the order of these fields or add additional
 * fields at the beginning of the structure.  If you do, job accounting
 * plugins will stop working.  If you need to add fields, add them 
 * at the end of the structure.
 */
typedef struct slurm_jobacct_ops {
	int (*jobacct_init)		(int frequency);
	int (*jobacct_fini)		(slurmd_job_t *job);
	int (*jobacct_suspend)	        ();
} slurm_jobacct_ops_t;

/*
 * These strings must be in the same order as the fields declared
 * for slurm_jobacct_ops_t.
 */
static const char *syms[] = {
	"jobacct_p_init",
	"jobacct_p_fini",
	"jobacct_p_suspend",
};

/*
 * A global job accounting context.  "Global" in the sense that there's
 * only one, with static bindings.  We don't export it.
 */

struct slurm_jobacct_context {
	char *			jobacct_type;
	plugrack_t		plugin_list;
	plugin_handle_t		cur_plugin;
	int			jobacct_errno;
	slurm_jobacct_ops_t	ops;
};

static slurm_jobacct_context_t g_context = NULL;
static pthread_mutex_t      context_lock = PTHREAD_MUTEX_INITIALIZER;

static int _plugin_init(void);
static slurm_jobacct_context_t
_slurm_jobacct_context_create( const char *jobacct_type)
{
	slurm_jobacct_context_t c;

	if ( jobacct_type == NULL ) {
		error( "_slurm_jobacct_context_create: no jobacct type" );
		return NULL;
	}

	c = xmalloc( sizeof( struct slurm_jobacct_context ) );

	c->jobacct_errno = SLURM_SUCCESS;

	/* Copy the job completion job completion type. */
	c->jobacct_type = xstrdup( jobacct_type );
	if ( c->jobacct_type == NULL ) {
		error( "can't make local copy of jobacct type" );
		xfree( c );
		return NULL; 
	}

	/* Plugin rack is demand-loaded on first reference. */
	c->plugin_list = NULL; 
	c->cur_plugin = PLUGIN_INVALID_HANDLE; 

	return c;
}

static int
_slurm_jobacct_context_destroy( slurm_jobacct_context_t c )
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

	xfree( c->jobacct_type );
	xfree( c );

	return SLURM_SUCCESS;
}

/*
 * Resolve the operations from the plugin.
 */
static slurm_jobacct_ops_t *
_slurm_jobacct_get_ops( slurm_jobacct_context_t c )
{
        int n_syms = sizeof( syms ) / sizeof( char * );

        /* Get the plugin list, if needed. */
        if ( c->plugin_list == NULL ) {
		char *plugin_dir;
                c->plugin_list = plugrack_create();
                if ( c->plugin_list == NULL ) {
                        error( "Unable to create a plugin manager" );
                        return NULL;
                }

                plugrack_set_major_type( c->plugin_list, "jobacct" );
                plugrack_set_paranoia( c->plugin_list, 
				       PLUGRACK_PARANOIA_NONE, 
				       0 );
		plugin_dir = slurm_get_plugin_dir();
                plugrack_read_dir( c->plugin_list, plugin_dir );
		xfree(plugin_dir);
        }
  
        /* Find the correct plugin. */
        c->cur_plugin = 
		plugrack_use_by_type( c->plugin_list, c->jobacct_type );
        if ( c->cur_plugin == PLUGIN_INVALID_HANDLE ) {
                error( "can't find a plugin for type %s", c->jobacct_type );
                return NULL;
        }  

        /* Dereference the API. */
        if ( plugin_get_syms( c->cur_plugin,
                              n_syms,
                              syms,
                              (void **) &c->ops ) < n_syms ) {
                error( "incomplete jobacct plugin detected" );
                return NULL;
        }

        return &c->ops;
}

static int
_plugin_init(void)
{
	char	*jobacct_type;
	int	retval=SLURM_SUCCESS;

	if ( g_context )
		_slurm_jobacct_context_destroy(g_context);

	jobacct_type = slurm_get_jobacct_type();
	g_context = _slurm_jobacct_context_create( jobacct_type );
	if ( g_context == NULL ) {
		error( "cannot create a context for %s", jobacct_type );
		xfree(jobacct_type);
		retval = SLURM_ERROR;
		goto done;
	}
	xfree(jobacct_type);

	if ( _slurm_jobacct_get_ops( g_context ) == NULL ) {
		error( "cannot resolve job accounting plugin operations" );
		_slurm_jobacct_context_destroy( g_context );
		g_context = NULL;
		retval = SLURM_ERROR;
	}

  done:
	return(retval);
}


extern int jobacct_g_init(int frequency)
{
	int retval = SLURM_SUCCESS;

	slurm_mutex_lock( &context_lock );
	_plugin_init(); 
	if ( g_context )
		retval = (*(g_context->ops.jobacct_init))(frequency);
	slurm_mutex_unlock( &context_lock );
	return retval;
}

extern int jobacct_g_fini(slurmd_job_t *job)
{
	return (*(g_context->ops.jobacct_fini))(job);
}

extern int jobacct_g_suspend()
{
	return (*(g_context->ops.jobacct_suspend))();
}


