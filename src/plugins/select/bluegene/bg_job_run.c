/*****************************************************************************\
 *  bg_job_run.c - blue gene job execution (e.g. initiation and termination)
 *  functions.
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004-2006 The Regents of the University of California.
 *  Copyright (C) 2011 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>, Danny Auble <da@llnl.gov>
 *
 *  This file is part of SLURM, a resource management program.
 *  For details, see <https://computing.llnl.gov/linux/slurm/>.
 *  Please also read the included file: DISCLAIMER.
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
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  endif
#  if WITH_PTHREADS
#    include <pthread.h>
#  endif
#endif

#include <signal.h>
#include <unistd.h>

#include <slurm/slurm_errno.h>

#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/common/node_select.h"
#include "src/common/uid.h"
#include "src/common/xstring.h"
#include "src/slurmctld/proc_req.h"
#include "bg_core.h"

bool deleting_old_blocks_flag = 0;

enum update_op {START_OP, TERM_OP, SYNC_OP};

typedef struct {
	char *bg_block_id;
	char *blrtsimage;       /* BlrtsImage for this block */
	uint16_t conn_type[SYSTEM_DIMENSIONS]; /* needed to boot small
				   blocks into HTC mode or not */
	struct job_record *job_ptr;	/* pointer to job running on
					 * block or NULL if no job */
	char *linuximage;       /* LinuxImage for this block */
	char *mloaderimage;     /* mloaderImage for this block */
	enum update_op op;	/* start | terminate | sync */
	char *ramdiskimage;     /* RamDiskImage for this block */
	uint16_t reboot;	/* reboot block before starting job */
} bg_action_t;

static void	_destroy_bg_action(void *x);
static int	_excise_block(List block_list,
			      char *bg_block_id,
			      char *nodes);
static List	_get_all_allocated_blocks(void);
static void *	_block_agent(void *args);
static void	_block_op(bg_action_t *bg_action_ptr);
static void	_start_agent(bg_action_t *bg_action_ptr);
static void	_sync_agent(bg_action_t *bg_action_ptr);

/* block_state_mutex should be locked before
 * calling this function.  This should only be called in _start_agent.
 * RET 1 if exists 0 if not, and job is requeued.
 */
static int _make_sure_block_still_exists(bg_action_t *bg_action_ptr,
					 bg_record_t *bg_record)
{
	/* check to make sure this block still exists since
	 * something could had happened and the block is no
	 * longer in existance */
	if ((bg_record->magic != BLOCK_MAGIC)
	    || !block_ptr_exist_in_list(bg_lists->main, bg_record)) {
		slurm_mutex_unlock(&block_state_mutex);
		debug("The block %s disappeared while starting "
		      "job %u requeueing if possible.",
		      bg_action_ptr->bg_block_id,
		      bg_action_ptr->job_ptr->job_id);
		bg_requeue_job(bg_action_ptr->job_ptr->job_id, 1);
		return 0;
	}
	return 1;
}


/* Delete a bg_action_t record */
static void _destroy_bg_action(void *x)
{
	bg_action_t *bg_action_ptr = (bg_action_t *) x;

	if (bg_action_ptr) {
		xfree(bg_action_ptr->blrtsimage);
		xfree(bg_action_ptr->linuximage);
		xfree(bg_action_ptr->mloaderimage);
		xfree(bg_action_ptr->ramdiskimage);
		xfree(bg_action_ptr->bg_block_id);
		xfree(bg_action_ptr);
	}
}

/* Update block user and reboot as needed */
static void _sync_agent(bg_action_t *bg_action_ptr)
{
	bg_record_t * bg_record = NULL;

	slurm_mutex_lock(&block_state_mutex);
	bg_record = find_bg_record_in_list(bg_lists->main,
					   bg_action_ptr->bg_block_id);
	if (!bg_record) {
		slurm_mutex_unlock(&block_state_mutex);
		error("No block %s", bg_action_ptr->bg_block_id);
		bg_requeue_job(bg_action_ptr->job_ptr->job_id, 1);
		return;
	}

	last_bg_update = time(NULL);
	bg_action_ptr->job_ptr->total_cpus =
		bg_action_ptr->job_ptr->details->min_cpus = bg_record->cpu_cnt;
	bg_record->job_running = bg_action_ptr->job_ptr->job_id;
	bg_record->job_ptr = bg_action_ptr->job_ptr;

	if (!block_ptr_exist_in_list(bg_lists->job_running, bg_record)) {
		list_push(bg_lists->job_running, bg_record);
		num_unused_cpus -= bg_record->cpu_cnt;
	}
	if (!block_ptr_exist_in_list(bg_lists->booted, bg_record))
		list_push(bg_lists->booted, bg_record);

	if (bg_record->state == BG_BLOCK_INITED) {
		if (bg_record->job_ptr) {
			bg_record->job_ptr->job_state &= (~JOB_CONFIGURING);
			last_job_update = time(NULL);
		}
		if (bg_record->user_uid != bg_action_ptr->job_ptr->user_id) {
			int set_user_rc = SLURM_SUCCESS;

			debug("User isn't correct for job %d on %s, "
			      "fixing...",
			      bg_action_ptr->job_ptr->job_id,
			      bg_action_ptr->bg_block_id);
			xfree(bg_record->target_name);
			bg_record->target_name =
				uid_to_string(bg_action_ptr->job_ptr->user_id);
			set_user_rc = set_block_user(bg_record);
			slurm_mutex_unlock(&block_state_mutex);

			if (set_user_rc == SLURM_ERROR)
				(void) slurm_fail_job(bg_record->job_running);
		} else
			slurm_mutex_unlock(&block_state_mutex);

	} else {
		if (bg_record->state != BG_BLOCK_BOOTING) {
			error("Block %s isn't ready and isn't "
			      "being configured! Starting job again.",
			      bg_action_ptr->bg_block_id);
		} else {
			debug("Block %s is booting, job ok",
			      bg_action_ptr->bg_block_id);
		}
		slurm_mutex_unlock(&block_state_mutex);
		_start_agent(bg_action_ptr);
	}
}

/* Perform job initiation work */
static void _start_agent(bg_action_t *bg_action_ptr)
{
	int rc, set_user_rc = SLURM_SUCCESS;
	bg_record_t *bg_record = NULL;
	bg_record_t *found_record = NULL;
	ListIterator itr;
	List delete_list = NULL;
	int requeue_job = 0;

	slurm_mutex_lock(&block_state_mutex);
	bg_record = find_bg_record_in_list(bg_lists->main,
					   bg_action_ptr->bg_block_id);

	if (!bg_record) {
		slurm_mutex_unlock(&block_state_mutex);
		error("block %s not found in bg_lists->main",
		      bg_action_ptr->bg_block_id);
		bg_requeue_job(bg_action_ptr->job_ptr->job_id, 1);
		return;
	}

	if (bg_record->job_running <= NO_JOB_RUNNING) {
		// bg_reset_block(bg_record); should already happened
		slurm_mutex_unlock(&block_state_mutex);
		debug("job %u finished during the queueing job "
		      "(everything is ok)",
		      bg_action_ptr->job_ptr->job_id);
		return;
	}
	if (bg_record->state == BG_BLOCK_TERM) {
		debug("Block is in Deallocating state, waiting for free.");
		bg_free_block(bg_record, 1, 1);
		/* no reason to reboot here since we are already
		   deallocating */
		bg_action_ptr->reboot = 0;
		/* Since bg_free_block will unlock block_state_mutex
		   we need to make sure the block we want is still
		   around.  Failure will unlock this so no need to
		   unlock before return.
		*/
		if (!_make_sure_block_still_exists(bg_action_ptr, bg_record))
			return;
	}

	delete_list = list_create(NULL);
	itr = list_iterator_create(bg_lists->main);
	while ((found_record = list_next(itr))) {
		if ((!found_record) || (bg_record == found_record))
			continue;

		if (!blocks_overlap(bg_record, found_record)) {
			debug2("block %s isn't part of %s",
			       found_record->bg_block_id,
			       bg_record->bg_block_id);
			continue;
		}

		if (found_record->job_ptr) {
			error("Trying to start job %u on block %s, "
			      "but there is a job %u running on an overlapping "
			      "block %s it will not end until %ld.  "
			      "This should never happen.",
			      bg_action_ptr->job_ptr->job_id,
			      bg_record->bg_block_id,
			      found_record->job_ptr->job_id,
			      found_record->bg_block_id,
			      found_record->job_ptr->end_time);
			requeue_job = 1;
			break;
		}

		debug2("need to make sure %s is free, it's part of %s",
		       found_record->bg_block_id,
		       bg_record->bg_block_id);
		list_push(delete_list, found_record);
	}
	list_iterator_destroy(itr);

	if (requeue_job) {
		list_destroy(delete_list);

		bg_reset_block(bg_record);

		slurm_mutex_unlock(&block_state_mutex);
		bg_requeue_job(bg_action_ptr->job_ptr->job_id, 0);
		return;
	}

	slurm_mutex_unlock(&block_state_mutex);

	rc = free_block_list(bg_action_ptr->job_ptr->job_id, delete_list, 0, 1);
	list_destroy(delete_list);
	if (rc != SLURM_SUCCESS) {
		error("Problem with deallocating blocks to run job %u "
		      "on block %s", bg_action_ptr->job_ptr->job_id,
		      bg_action_ptr->bg_block_id);
		if (IS_JOB_CONFIGURING(bg_action_ptr->job_ptr))
			bg_requeue_job(bg_action_ptr->job_ptr->job_id, 0);
		return;
	}

	slurm_mutex_lock(&block_state_mutex);
	/* Failure will unlock block_state_mutex so no need to unlock before
	   return. */
	if (!_make_sure_block_still_exists(bg_action_ptr, bg_record))
		return;

	if (bg_record->job_running <= NO_JOB_RUNNING) {
		// bg_reset_block(bg_record); should already happened
		slurm_mutex_unlock(&block_state_mutex);
		debug("job %u already finished before boot",
		      bg_action_ptr->job_ptr->job_id);
		return;
	}

	rc = 0;
#ifdef HAVE_BGL
	if (bg_action_ptr->blrtsimage
	   && strcasecmp(bg_action_ptr->blrtsimage, bg_record->blrtsimage)) {
		debug3("changing BlrtsImage from %s to %s",
		       bg_record->blrtsimage, bg_action_ptr->blrtsimage);
		xfree(bg_record->blrtsimage);
		bg_record->blrtsimage = xstrdup(bg_action_ptr->blrtsimage);
		rc = 1;
	}
#elif defined HAVE_BGP
	if ((bg_action_ptr->conn_type[0] >= SELECT_SMALL)
	   && (bg_action_ptr->conn_type[0] != bg_record->conn_type[0])) {
		debug3("changing small block mode from %s to %s",
		       conn_type_string(bg_record->conn_type[0]),
		       conn_type_string(bg_action_ptr->conn_type[0]));
		rc = 1;
# ifndef HAVE_BG_FILES
		/* since we don't check state on an emulated system we
		 * have to change it here
		 */
		bg_record->conn_type[0] = bg_action_ptr->conn_type[0];
# endif
	}
#endif

#ifdef HAVE_BG_L_P
	if (bg_action_ptr->linuximage
	   && strcasecmp(bg_action_ptr->linuximage, bg_record->linuximage)) {
# ifdef HAVE_BGL
		debug3("changing LinuxImage from %s to %s",
		       bg_record->linuximage, bg_action_ptr->linuximage);
# else
		debug3("changing CnloadImage from %s to %s",
		       bg_record->linuximage, bg_action_ptr->linuximage);
# endif
		xfree(bg_record->linuximage);
		bg_record->linuximage = xstrdup(bg_action_ptr->linuximage);
		rc = 1;
	}
	if (bg_action_ptr->ramdiskimage
	   && strcasecmp(bg_action_ptr->ramdiskimage,
			 bg_record->ramdiskimage)) {
# ifdef HAVE_BGL
		debug3("changing RamDiskImage from %s to %s",
		       bg_record->ramdiskimage, bg_action_ptr->ramdiskimage);
# else
		debug3("changing IoloadImage from %s to %s",
		       bg_record->ramdiskimage, bg_action_ptr->ramdiskimage);
# endif
		xfree(bg_record->ramdiskimage);
		bg_record->ramdiskimage = xstrdup(bg_action_ptr->ramdiskimage);
		rc = 1;
	}
#endif
	if (bg_action_ptr->mloaderimage
	   && strcasecmp(bg_action_ptr->mloaderimage,
			 bg_record->mloaderimage)) {
		debug3("changing MloaderImage from %s to %s",
		       bg_record->mloaderimage, bg_action_ptr->mloaderimage);
		xfree(bg_record->mloaderimage);
		bg_record->mloaderimage = xstrdup(bg_action_ptr->mloaderimage);
		rc = 1;
	}

	if (rc) {
		bg_record->modifying = 1;

		bg_free_block(bg_record, 1, 1);

		/* Since bg_free_block will unlock block_state_mutex
		   we need to make sure the block we want is still
		   around.  Failure will unlock block_state_mutex so
		   no need to unlock before return.
		*/
		if (!_make_sure_block_still_exists(bg_action_ptr, bg_record))
			return;

#if defined HAVE_BG_FILES && defined HAVE_BG_L_P
#ifdef HAVE_BGL
		if ((rc = bridge_block_modify(bg_record->bg_block_id,
					      RM_MODIFY_BlrtsImg,
					      bg_record->blrtsimage))
		    != SLURM_SUCCESS)
			error("bridge_block_modify(RM_MODIFY_BlrtsImg): %s",
			      bg_err_str(rc));

		if ((rc = bridge_block_modify(bg_record->bg_block_id,
					      RM_MODIFY_LinuxImg,
					      bg_record->linuximage))
		    != SLURM_SUCCESS)
			error("bridge_block_modify(RM_MODIFY_LinuxImg): %s",
			      bg_err_str(rc));

		if ((rc = bridge_block_modify(bg_record->bg_block_id,
					      RM_MODIFY_RamdiskImg,
					      bg_record->ramdiskimage))
		    != SLURM_SUCCESS)
			error("bridge_block_modify(RM_MODIFY_RamdiskImg): %s",
			      bg_err_str(rc));

#elif defined HAVE_BGP
		if ((rc = bridge_block_modify(bg_record->bg_block_id,
					      RM_MODIFY_CnloadImg,
					      bg_record->linuximage))
		    != SLURM_SUCCESS)
			error("bridge_block_modify(RM_MODIFY_CnloadImg): %s",
			      bg_err_str(rc));

		if ((rc = bridge_block_modify(bg_record->bg_block_id,
					      RM_MODIFY_IoloadImg,
					      bg_record->ramdiskimage))
		    != SLURM_SUCCESS)
			error("bridge_block_modify(RM_MODIFY_IoloadImg): %s",
			      bg_err_str(rc));

		if (bg_action_ptr->conn_type[0] > SELECT_SMALL) {
			char *conn_type = NULL;
			switch(bg_action_ptr->conn_type[0]) {
			case SELECT_HTC_S:
				conn_type = "s";
				break;
			case SELECT_HTC_D:
				conn_type = "d";
				break;
			case SELECT_HTC_V:
				conn_type = "v";
				break;
			case SELECT_HTC_L:
				conn_type = "l";
				break;
			default:
				break;
			}
			/* the option has to be set before the pool can be
			   set */
			if ((rc = bridge_block_modify(
				     bg_record->bg_block_id,
				     RM_MODIFY_Options,
				     conn_type)) != SLURM_SUCCESS)
				error("bridge_set_data(RM_MODIFY_Options): %s",
				      bg_err_str(rc));
		}
#endif
		if ((rc = bridge_block_modify(bg_record->bg_block_id,
					      RM_MODIFY_MloaderImg,
					      bg_record->mloaderimage))
		    != SLURM_SUCCESS)
			error("bridge_block_modify(RM_MODIFY_MloaderImg): %s",
			      bg_err_str(rc));

#endif
		bg_record->modifying = 0;
	} else if (bg_action_ptr->reboot) {
		bg_record->modifying = 1;

		bg_free_block(bg_record, 1, 1);

		/* Since bg_free_block will unlock block_state_mutex
		   we need to make sure the block we want is still
		   around.  Failure will unlock block_state_mutex so
		   no need to unlock before return.
		*/
		if (!_make_sure_block_still_exists(bg_action_ptr, bg_record))
			return;

		bg_record->modifying = 0;
	}

	if (bg_record->state == BG_BLOCK_FREE) {
		if ((rc = boot_block(bg_record)) != SLURM_SUCCESS) {
			/* Since boot_block could unlock block_state_mutex
			   on error we need to make sure the block we
			   want is still around.  Failure will unlock
			   block_state_mutex so no need to unlock
			   before return.
			*/
			if (!_make_sure_block_still_exists(bg_action_ptr,
							   bg_record))
				return;
			bg_reset_block(bg_record);
			slurm_mutex_unlock(&block_state_mutex);
			bg_requeue_job(bg_action_ptr->job_ptr->job_id, 1);
			return;
		}
	} else if (bg_record->state == BG_BLOCK_BOOTING) {
#ifdef HAVE_BG_FILES
		bg_record->boot_state = 1;
#else
		if (!block_ptr_exist_in_list(bg_lists->booted, bg_record))
			list_push(bg_lists->booted, bg_record);
		bg_record->state = BG_BLOCK_INITED;
		last_bg_update = time(NULL);
#endif
	}


	if (bg_record->job_running <= NO_JOB_RUNNING) {
		slurm_mutex_unlock(&block_state_mutex);
		debug("job %u finished during the start of the boot "
		      "(everything is ok)",
		      bg_action_ptr->job_ptr->job_id);
		return;
	}

	/* Don't reset boot_count, it will be reset when state
	   changes, and needs to outlast a job allocation.
	*/
	/* bg_record->boot_count = 0; */
	xfree(bg_record->target_name);
	bg_record->target_name = uid_to_string(bg_action_ptr->job_ptr->user_id);
	debug("setting the target_name for Block %s to %s",
	      bg_record->bg_block_id, bg_record->target_name);

	if (bg_record->state == BG_BLOCK_INITED) {
		debug("block %s is ready.", bg_record->bg_block_id);
		set_user_rc = set_block_user(bg_record);
		if (bg_action_ptr->job_ptr) {
			bg_action_ptr->job_ptr->job_state &= (~JOB_CONFIGURING);
			last_job_update = time(NULL);
		}
	}
	slurm_mutex_unlock(&block_state_mutex);

	if (set_user_rc == SLURM_ERROR) {
		sleep(2);
		/* wait for the slurmd to begin
		   the batch script, slurm_fail_job()
		   is a no-op if issued prior
		   to the script initiation do clean up just
		   incase the fail job isn't ran */
		(void) slurm_fail_job(bg_record->job_running);
		slurm_mutex_lock(&block_state_mutex);
		if (remove_from_bg_list(bg_lists->job_running, bg_record)
		    == SLURM_SUCCESS)
			num_unused_cpus += bg_record->cpu_cnt;

		slurm_mutex_unlock(&block_state_mutex);
	}
}

static void *_block_agent(void *args)
{
	bg_action_t *bg_action_ptr = (bg_action_t *)args;

	if (bg_action_ptr->op == START_OP)
		_start_agent(bg_action_ptr);
	else if (bg_action_ptr->op == TERM_OP)
		bridge_block_post_job(bg_action_ptr->bg_block_id);
	else if (bg_action_ptr->op == SYNC_OP)
		_sync_agent(bg_action_ptr);
	_destroy_bg_action(bg_action_ptr);

	return NULL;
}

/* Perform an operation upon a BG block (block) for starting or
 * terminating a job */
static void _block_op(bg_action_t *bg_action_ptr)
{
	pthread_attr_t attr_agent;
	pthread_t thread_agent;
	int retries;

	/* spawn an agent */
	slurm_attr_init(&attr_agent);
	if (pthread_attr_setdetachstate(&attr_agent, PTHREAD_CREATE_DETACHED))
		error("pthread_attr_setdetachstate error %m");

	retries = 0;
	while (pthread_create(&thread_agent, &attr_agent,
			      _block_agent, bg_action_ptr)) {
		error("pthread_create error %m");
		if (++retries > MAX_PTHREAD_RETRIES)
			fatal("Can't create pthread");
		usleep(1000);	/* sleep and retry */
	}
	slurm_attr_destroy(&attr_agent);
}


/* get a list of all BG blocks with users */
static List _get_all_allocated_blocks(void)
{
	List ret_list = list_create(destroy_bg_record);
	ListIterator itr;
	bg_record_t *bg_record = NULL;
	bg_record_t *rm_record = NULL;

	if (!ret_list)
		fatal("malloc error");

	if (bg_lists->main) {
		itr = list_iterator_create(bg_lists->main);
		while ((bg_record = (bg_record_t *) list_next(itr))) {
			if ((bg_record->user_name == NULL)
			    ||  (bg_record->user_name[0] == '\0')
			    ||  (bg_record->bg_block_id == NULL)
			    ||  (bg_record->bg_block_id[0] == '0'))
				continue;
			rm_record = xmalloc(sizeof(bg_record_t));
			rm_record->magic = BLOCK_MAGIC;
			rm_record->bg_block_id =
				xstrdup(bg_record->bg_block_id);
			rm_record->mp_str = xstrdup(bg_record->mp_str);

			list_append(ret_list, rm_record);
		}
		list_iterator_destroy(itr);
	} else {
		error("_get_all_allocated_blocks: no bg_lists->main");
	}

	return ret_list;
}

/* remove a BG block from the given list */
static int _excise_block(List block_list, char *bg_block_id,
			 char *nodes)
{
	int rc = SLURM_SUCCESS;
	ListIterator iter;
	bg_record_t *bg_record = NULL;

	if (block_list) {
		iter = list_iterator_create(block_list);
		xassert(iter);
		while ((bg_record = list_next(iter))) {
			rc = SLURM_ERROR;
			if (strcmp(bg_record->bg_block_id, bg_block_id))
				continue;
			if (strcmp(bg_record->mp_str, nodes)) {
				/* changed bgblock */
				error("bg_block_id:%s old_nodes:%s "
				      "new_nodes:%s",
				      bg_block_id, nodes, bg_record->mp_str);
				break;
			}

			/* exact match of name and node list */
			debug("synced Block %s", bg_block_id);
			list_delete_item(iter);
			rc = SLURM_SUCCESS;
			break;
		}
		list_iterator_destroy(iter);
	} else {
		error("_excise_block: No block_list");
		rc = SLURM_ERROR;
	}
	return rc;
}

/*
 * Perform any work required to terminate a jobs on a block.
 * bg_block_id IN - block name
 * RET - SLURM_SUCCESS or an error code
 *
 * NOTE: The job is killed before the function returns. This can take
 * many seconds. Do not call from slurmctld  or any other entity that
 * can not wait.
 */
int term_jobs_on_block(char *bg_block_id)
{
	int rc = SLURM_SUCCESS;
	bg_action_t *bg_action_ptr;

	bg_action_ptr = xmalloc(sizeof(bg_action_t));
	bg_action_ptr->op = TERM_OP;
	bg_action_ptr->bg_block_id = xstrdup(bg_block_id);
	_block_op(bg_action_ptr);

	return rc;
}

/*
 * Perform any setup required to initiate a job
 * job_ptr IN - pointer to the job being initiated
 * RET - SLURM_SUCCESS or an error code
 *
 * NOTE: This happens in parallel with srun and slurmd spawning
 * the job. A prolog script is expected to defer initiation of
 * the job script until the BG block is available for use.
 */
extern int start_job(struct job_record *job_ptr)
{
	int rc = SLURM_SUCCESS;
	bg_record_t *bg_record = NULL;

	bg_action_t *bg_action_ptr = NULL;

	bg_action_ptr = xmalloc(sizeof(bg_action_t));
	bg_action_ptr->op = START_OP;
	bg_action_ptr->job_ptr = job_ptr;

	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_BLOCK_ID,
			   &(bg_action_ptr->bg_block_id));
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_REBOOT,
			   &(bg_action_ptr->reboot));
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_CONN_TYPE,
			   &(bg_action_ptr->conn_type));
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_MLOADER_IMAGE,
			   &(bg_action_ptr->mloaderimage));
#ifdef HAVE_BG_L_P
# ifdef HAVE_BGL
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_BLRTS_IMAGE,
			   &(bg_action_ptr->blrtsimage));
	if (!bg_action_ptr->blrtsimage) {
		bg_action_ptr->blrtsimage =
			xstrdup(bg_conf->default_blrtsimage);
		set_select_jobinfo(job_ptr->select_jobinfo->data,
				   SELECT_JOBDATA_BLRTS_IMAGE,
				   bg_action_ptr->blrtsimage);
	}
# elif HAVE_BGP
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_CONN_TYPE,
			   &(bg_action_ptr->conn_type));
# endif
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_LINUX_IMAGE,
			   &(bg_action_ptr->linuximage));
	if (!bg_action_ptr->linuximage) {
		bg_action_ptr->linuximage =
			xstrdup(bg_conf->default_linuximage);
		set_select_jobinfo(job_ptr->select_jobinfo->data,
				   SELECT_JOBDATA_LINUX_IMAGE,
				   bg_action_ptr->linuximage);
	}

	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_RAMDISK_IMAGE,
			   &(bg_action_ptr->ramdiskimage));
	if (!bg_action_ptr->ramdiskimage) {
		bg_action_ptr->ramdiskimage =
			xstrdup(bg_conf->default_ramdiskimage);
		set_select_jobinfo(job_ptr->select_jobinfo->data,
				   SELECT_JOBDATA_RAMDISK_IMAGE,
				   bg_action_ptr->ramdiskimage);
	}

#endif
	if (!bg_action_ptr->mloaderimage) {
		bg_action_ptr->mloaderimage =
			xstrdup(bg_conf->default_mloaderimage);
		set_select_jobinfo(job_ptr->select_jobinfo->data,
				   SELECT_JOBDATA_MLOADER_IMAGE,
				   bg_action_ptr->mloaderimage);
	}

	slurm_mutex_lock(&block_state_mutex);
	bg_record = find_bg_record_in_list(bg_lists->main,
					   bg_action_ptr->bg_block_id);
	if (!bg_record) {
		slurm_mutex_unlock(&block_state_mutex);
		error("bg_record %s doesn't exist, requested for job (%d)",
		      bg_action_ptr->bg_block_id, job_ptr->job_id);
		_destroy_bg_action(bg_action_ptr);
		return SLURM_ERROR;
	}

	last_bg_update = time(NULL);
	job_ptr->total_cpus = job_ptr->details->min_cpus = bg_record->cpu_cnt;
	bg_record->job_running = bg_action_ptr->job_ptr->job_id;
	bg_record->job_ptr = bg_action_ptr->job_ptr;
	if (!block_ptr_exist_in_list(bg_lists->job_running, bg_record)) {
		list_push(bg_lists->job_running, bg_record);
		num_unused_cpus -= bg_record->cpu_cnt;
	}
	if (!block_ptr_exist_in_list(bg_lists->booted, bg_record))
		list_push(bg_lists->booted, bg_record);
	slurm_mutex_unlock(&block_state_mutex);

	info("Queue start of job %u in BG block %s",
	     job_ptr->job_id,
	     bg_action_ptr->bg_block_id);
	_block_op(bg_action_ptr);
	return rc;
}


/*
 * Perform any work required to terminate a job
 * job_ptr IN - pointer to the job being terminated
 * RET - SLURM_SUCCESS or an error code
 *
 * NOTE: This happens in parallel with srun and slurmd terminating
 * the job. Insure that this function, mpirun and the epilog can
 * all deal with termination race conditions.
 */
int term_job(struct job_record *job_ptr)
{
	int rc = SLURM_SUCCESS;
	bg_action_t *bg_action_ptr = NULL;

	bg_action_ptr = xmalloc(sizeof(bg_action_t));
	bg_action_ptr->op = TERM_OP;
	bg_action_ptr->job_ptr = job_ptr;
	get_select_jobinfo(job_ptr->select_jobinfo->data,
			   SELECT_JOBDATA_BLOCK_ID,
			   &(bg_action_ptr->bg_block_id));
	info("Queue termination of job %u in BG block %s",
	     job_ptr->job_id, bg_action_ptr->bg_block_id);
	_block_op(bg_action_ptr);

	return rc;
}

/*
 * Synchronize BG block state to that of currently active jobs.
 * This can recover from slurmctld crashes when block usership
 * changes were queued
 */
extern int sync_jobs(List job_list)
{
	ListIterator job_iterator;
	struct job_record  *job_ptr = NULL;
	bg_action_t *bg_action_ptr = NULL;
	List block_list = NULL;
	static bool run_already = false;

	/* Execute only on initial startup. We don't support bgblock
	 * creation on demand today, so there is no need to re-sync data. */
	if (run_already)
		return SLURM_SUCCESS;
	run_already = true;

	/* Insure that all running jobs own the specified block */
	block_list = _get_all_allocated_blocks();
	if (job_list) {
		job_iterator = list_iterator_create(job_list);
		while ((job_ptr = (struct job_record *)
			list_next(job_iterator))) {
			bool good_block = true;
			if (!IS_JOB_RUNNING(job_ptr))
				continue;

			bg_action_ptr = xmalloc(sizeof(bg_action_t));
			bg_action_ptr->op = SYNC_OP;
			bg_action_ptr->job_ptr = job_ptr;

			get_select_jobinfo(job_ptr->select_jobinfo->data,
					   SELECT_JOBDATA_BLOCK_ID,
					   &(bg_action_ptr->bg_block_id));
#ifdef HAVE_BG_L_P
# ifdef HAVE_BGL
			get_select_jobinfo(job_ptr->select_jobinfo->data,
					   SELECT_JOBDATA_BLRTS_IMAGE,
					   &(bg_action_ptr->blrtsimage));
# else
			get_select_jobinfo(job_ptr->select_jobinfo->data,
					   SELECT_JOBDATA_CONN_TYPE,
					   &(bg_action_ptr->conn_type));
# endif
			get_select_jobinfo(job_ptr->select_jobinfo->data,
					   SELECT_JOBDATA_LINUX_IMAGE,
					   &(bg_action_ptr->linuximage));
			get_select_jobinfo(job_ptr->select_jobinfo->data,
					   SELECT_JOBDATA_RAMDISK_IMAGE,
					   &(bg_action_ptr->ramdiskimage));
#endif
			get_select_jobinfo(job_ptr->select_jobinfo->data,
					   SELECT_JOBDATA_MLOADER_IMAGE,
					   &(bg_action_ptr->mloaderimage));

			if (bg_action_ptr->bg_block_id == NULL) {
				error("Running job %u has bgblock==NULL",
				      job_ptr->job_id);
				good_block = false;
			} else if (job_ptr->nodes == NULL) {
				error("Running job %u has nodes==NULL",
				      job_ptr->job_id);
				good_block = false;
			} else if (_excise_block(block_list,
						 bg_action_ptr->bg_block_id,
						 job_ptr->nodes)
				   != SLURM_SUCCESS) {
				error("Kill job %u belongs to defunct "
				      "bgblock %s",
				      job_ptr->job_id,
				      bg_action_ptr->bg_block_id);
				good_block = false;
			}
			if (!good_block) {
				job_ptr->job_state = JOB_FAILED
					| JOB_COMPLETING;
				job_ptr->end_time = time(NULL);
				last_job_update = time(NULL);
				_destroy_bg_action(bg_action_ptr);
				continue;
			}

			debug3("Queue sync of job %u in BG block %s "
			       "ending at %ld",
			       job_ptr->job_id,
			       bg_action_ptr->bg_block_id,
			       job_ptr->end_time);
			_block_op(bg_action_ptr);
		}
		list_iterator_destroy(job_iterator);
	} else {
		error("sync_jobs: no job_list");
		list_destroy(block_list);
		return SLURM_ERROR;
	}
	/* Insure that all other blocks are free of users */
	if (block_list) {
		bridge_reset_block_list(block_list);
		list_destroy(block_list);
	} else {
		/* this should never happen,
		 * vestigial logic */
		error("sync_jobs: no block_list");
		return SLURM_ERROR;
	}
	return SLURM_SUCCESS;
}

/*
 * Boot a block. Block state expected to be FREE upon entry.
 * NOTE: This function does not wait for the boot to complete.
 * the slurm prolog script needs to perform the waiting.
 * NOTE: block_state_mutex needs to be locked before entering.
 */
extern int boot_block(bg_record_t *bg_record)
{
	int rc;
	if (bg_record->magic != BLOCK_MAGIC) {
		error("boot_block: magic was bad");
		return SLURM_ERROR;
	}

	info("Booting block %s", bg_record->bg_block_id);
	if ((rc = bridge_block_boot(bg_record)) != SLURM_SUCCESS) {
		/* error("bridge_create_block(%s): %s", */
		/*       bg_record->bg_block_id, bg_err_str(rc)); */
		/* if (rc == INCOMPATIBLE_STATE) { */
		/* 	char reason[200]; */
		/* 	snprintf(reason, sizeof(reason), */
		/* 		 "boot_block: " */
		/* 		 "Block %s is in an incompatible state.  " */
		/* 		 "This usually means hardware is allocated " */
		/* 		 "by another block (maybe outside of SLURM).", */
		/* 		 bg_record->bg_block_id); */
		/* 	bg_record->boot_state = 0; */
		/* 	bg_record->boot_count = 0; */
		/* 	slurm_mutex_unlock(&block_state_mutex); */
		/* 	requeue_and_error(bg_record, reason); */
		/* 	slurm_mutex_lock(&block_state_mutex); */
		/* } */
		return SLURM_ERROR;
	}

	return SLURM_SUCCESS;
}
