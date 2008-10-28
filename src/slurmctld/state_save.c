/*****************************************************************************\
 *  state_save.c - Keep saved slurmctld state current 
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
 *  LLNL-CODE-402394.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif                          /* WITH_PTHREADS */

#include "src/common/macros.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/trigger_mgr.h"

static pthread_mutex_t state_save_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  state_save_cond = PTHREAD_COND_INITIALIZER;
static int save_jobs = 0, save_nodes = 0, save_parts = 0, save_triggers = 0;
static bool run_save_thread = true;

/* Queue saving of job state information */
extern void schedule_job_save(void)
{
	slurm_mutex_lock(&state_save_lock);
	save_jobs++;
	slurm_mutex_unlock(&state_save_lock);
	pthread_cond_broadcast(&state_save_cond);
}

/* Queue saving of node state information */
extern void schedule_node_save(void)
{
	slurm_mutex_lock(&state_save_lock);
	save_nodes++;
	slurm_mutex_unlock(&state_save_lock);
	pthread_cond_broadcast(&state_save_cond);
}

/* Queue saving of partition state information */
extern void schedule_part_save(void)
{
	slurm_mutex_lock(&state_save_lock);
	save_parts++;
	slurm_mutex_unlock(&state_save_lock);
	pthread_cond_broadcast(&state_save_cond);
}

/* Queue saving of trigger state information */
extern void schedule_trigger_save(void)
{
	slurm_mutex_lock(&state_save_lock);
	save_triggers++;
	slurm_mutex_unlock(&state_save_lock);
	pthread_cond_broadcast(&state_save_cond);
}

/* shutdown the slurmctld_state_save thread */
extern void shutdown_state_save(void)
{
	slurm_mutex_lock(&state_save_lock);
	run_save_thread = false;
	slurm_mutex_unlock(&state_save_lock);
	pthread_cond_broadcast(&state_save_cond);
}

/*
 * Run as pthread to keep saving slurmctld state information as needed,
 * Use schedule_job_save(),  schedule_node_save(), and schedule_part_save()
 * to queue state save of each data structure 
 * no_data IN - unused
 * RET - NULL
 */
extern void *slurmctld_state_save(void *no_data)
{
	bool run_save;

	while (1) {
		/* wait for work to perform */
		slurm_mutex_lock(&state_save_lock);
		while (1) {
			if (save_jobs + save_nodes + save_parts + 
			    save_triggers)
				break;		/* do the work */
			else if (!run_save_thread) {
				run_save_thread = true;
				slurm_mutex_unlock(&state_save_lock);
				return NULL;	/* shutdown */
			} else 			/* wait for more work */
				pthread_cond_wait(&state_save_cond, 
				                  &state_save_lock);
		}

		/* save job info if necessary */
		run_save = false;
		/* slurm_mutex_lock(&state_save_lock); done above */
		if (save_jobs) {
			run_save = true;
			save_jobs = 0;
		}
		slurm_mutex_unlock(&state_save_lock);
		if (run_save)
			(void)dump_all_job_state();

		/* save node info if necessary */
		run_save = false;
		slurm_mutex_lock(&state_save_lock);
		if (save_nodes) {
			run_save = true;
			save_nodes = 0;
		}
		slurm_mutex_unlock(&state_save_lock);
		if (run_save)
			(void)dump_all_node_state();

		/* save partition info if necessary */
		run_save = false;
		slurm_mutex_lock(&state_save_lock);
		if (save_parts) {
			run_save = true;
			save_parts = 0;
		}
		slurm_mutex_unlock(&state_save_lock);
		if (run_save)
			(void)dump_all_part_state();

		/* save trigger info if necessary */
		run_save = false;
		slurm_mutex_lock(&state_save_lock);
		if (save_triggers) {
			run_save = true;
			save_triggers = 0;
		}
		slurm_mutex_unlock(&state_save_lock);
		if (run_save)
			(void)trigger_state_save();
	}
}

