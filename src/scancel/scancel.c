/*****************************************************************************\
 *  scancel - cancel specified job(s) and/or job step(s)
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#if HAVE_INTTYPES_H
#  include <inttypes.h>
#else  /* !HAVE_INTTYPES_H */
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#endif  /* HAVE_INTTYPES_H */

#include <slurm/slurm.h>

#include "src/common/log.h"
#include "src/common/xstring.h"
#include "src/common/xmalloc.h"
#include "src/scancel/scancel.h"

#define MAX_CANCEL_RETRY 10

static void _cancel_jobs (void);
static void _cancel_job_id (uint32_t job_id, uint16_t signal);
static void _cancel_step_id (uint32_t job_id, uint32_t step_id, 
			     uint16_t signal);
static int  _confirmation (int i, uint32_t step_id);
static void _filter_job_records (void);
static void _load_job_records (void);

static job_info_msg_t * job_buffer_ptr = NULL;

int
main (int argc, char *argv[]) 
{
	log_options_t log_opts = LOG_OPTS_STDERR_ONLY ;

	log_init (xbasename(argv[0]), log_opts, SYSLOG_FACILITY_DAEMON, NULL);
	initialize_and_process_args(argc, argv);
	if (opt.verbose) {
		log_opts.stderr_level += opt.verbose;
		log_alter (log_opts, SYSLOG_FACILITY_DAEMON, NULL);
	} 

	if ((opt.interactive) ||
	    (opt.job_name) ||
	    (opt.partition) ||
	    (opt.state != JOB_END) ||
	    (opt.user_name)) {
		_load_job_records ();
		_filter_job_records ();
	}

	_cancel_jobs ();

	exit (0);
}


/* _load_job_records - load all job information for filtering and verification */
static void 
_load_job_records (void)
{
	int error_code;

	error_code = slurm_load_jobs ((time_t) NULL, &job_buffer_ptr, 1);

	if (error_code) {
		slurm_perror ("slurm_load_jobs error");
		exit (1);
	}
}


/* _filter_job_records - filtering job information per user specification */
static void 
_filter_job_records (void)
{
	int i, j;
	job_info_t *job_ptr = NULL;

	job_ptr = job_buffer_ptr->job_array ;
	for (i = 0; i < job_buffer_ptr->record_count; i++) {
		if (job_ptr[i].job_id == 0) 
			continue;

		if ((job_ptr[i].job_state != JOB_PENDING)
		&&  (job_ptr[i].job_state != JOB_RUNNING)
		&&  (job_ptr[i].job_state != JOB_SUSPENDED)) {
			job_ptr[i].job_id = 0;
			continue;
		}

		if ((opt.job_name != NULL) &&
		    (strcmp(job_ptr[i].name,opt.job_name) != 0)) {
			job_ptr[i].job_id = 0;
			continue;
		}

		if ((opt.partition != NULL) &&
		    (strcmp(job_ptr[i].partition,opt.partition) != 0)) {
			job_ptr[i].job_id = 0;
			continue;
		}

		if ((opt.state != JOB_END) &&
		    (job_ptr[i].job_state != opt.state)) {
			job_ptr[i].job_id = 0;
			continue;
		}

		if ((opt.user_name != NULL) &&
		    (job_ptr[i].user_id != opt.user_id)) {
			job_ptr[i].job_id = 0;
			continue;
		}

		if (opt.job_cnt == 0)
			continue;
		for (j = 0; j < opt.job_cnt; j++) {
			if (job_ptr[i].job_id == opt.job_id[j])
				break;
		}
		if (j >= opt.job_cnt) { /* not found */
			job_ptr[i].job_id = 0;
			continue;
		}
	}
}


/* _cancel_jobs - filter then cancel jobs or job steps per request */
static void
_cancel_jobs (void)
{
	int i, j;
	job_info_t *job_ptr = NULL;

	if (opt.job_cnt && opt.interactive) {	/* confirm cancel */
		job_ptr = job_buffer_ptr->job_array ;
		for (j = 0; j < opt.job_cnt; j++ ) {
			for (i = 0; i < job_buffer_ptr->record_count; i++) {
				if (job_ptr[i].job_id != opt.job_id[j]) 
					continue;
				if (opt.interactive && 
				    (_confirmation(i, opt.step_id[j]) == 0))
					break;
				if (opt.step_id[j] == SLURM_BATCH_SCRIPT)
					_cancel_job_id (opt.job_id[j], 
							opt.signal);
				else
					_cancel_step_id (opt.job_id[j], 
					                opt.step_id[j],
							opt.signal);
				break;
			}
			if (i >= job_buffer_ptr->record_count)
				fprintf (stderr, "Job %u not found\n", 
				         opt.job_id[j]);
		}

	} else if (opt.job_cnt) {	/* delete specific jobs */
		for (j = 0; j < opt.job_cnt; j++ ) {
			if (opt.step_id[j] == SLURM_BATCH_SCRIPT)
				_cancel_job_id (opt.job_id[j], 
						opt.signal);
			else
				_cancel_step_id (opt.job_id[j], 
				                opt.step_id[j], 
						opt.signal);
		}

	} else {		/* delete all jobs per filtering */
		job_ptr = job_buffer_ptr->job_array ;
		for (i = 0; i < job_buffer_ptr->record_count; i++) {
			if (job_ptr[i].job_id == 0) 
				continue;
			if (opt.interactive && 
			    (_confirmation(i, SLURM_BATCH_SCRIPT) == 0))
				continue;
			_cancel_job_id (job_ptr[i].job_id, opt.signal);
		}
	}
}

static void
_cancel_job_id (uint32_t job_id, uint16_t signal)
{
	int error_code = SLURM_SUCCESS, i;

	for (i=0; i<MAX_CANCEL_RETRY; i++) {
		if (signal == (uint16_t)-1) {
			verbose("Signal %u to job %u", SIGKILL, job_id);
			error_code = slurm_kill_job (job_id, SIGKILL,
						     (uint16_t)opt.batch);
		} else {
			verbose("Signal %u to job %u", signal, job_id);
			if (opt.batch)
				error_code = slurm_signal_job_step(job_id,
							   SLURM_BATCH_SCRIPT,
							   signal);
			else
				error_code = slurm_signal_job (job_id, signal);
		}
		if ((error_code == 0) || 
		    (errno != ESLURM_TRANSITION_STATE_NO_UPDATE))
			break;
		verbose("Job is in transistional state, retrying");
		sleep ( 5 + i );
	}
	if (error_code) {
		error_code = slurm_get_errno();
		if ((opt.verbose >= 0) || 
		    ((error_code != ESLURM_ALREADY_DONE) &&
		     (error_code != ESLURM_INVALID_JOB_ID)))
			error("Kill job error on job id %u: %s", 
				job_id, slurm_strerror(slurm_get_errno()));
	}
}

static void
_cancel_step_id (uint32_t job_id, uint32_t step_id, uint16_t signal)
{
	int error_code = SLURM_SUCCESS, i;

	for (i=0; i<MAX_CANCEL_RETRY; i++) {
		if (signal == (uint16_t)-1) {
			verbose("Signal %u to step %u.u",
				SIGKILL, job_id, step_id);
			error_code = slurm_kill_job_step (job_id, step_id,
							  SIGKILL);
		} else {
			verbose("Signal %u to step %u.u",
				signal, job_id, step_id);
			error_code = slurm_signal_job_step(job_id, step_id,
							   signal);
		}
		if ((error_code == 0) || 
		    (errno != ESLURM_TRANSITION_STATE_NO_UPDATE))
			break;
		verbose("Job is in transistional state, retrying");
		sleep ( 5 + i );
	}
	if (error_code) {
		error_code = slurm_get_errno();
		if ((opt.verbose >= 0) || (error_code != ESLURM_ALREADY_DONE ))
			error("Kill job error on job id %u.%u: %s", 
		 		job_id, step_id, slurm_strerror(slurm_get_errno()));
	}
}

/* _confirmation - Confirm job cancel request interactively */
static int 
_confirmation (int i, uint32_t step_id)
{
	char in_line[128];
	job_info_t *job_ptr = NULL;

	job_ptr = job_buffer_ptr->job_array ;
	while (1) {
		if (step_id == SLURM_BATCH_SCRIPT) {
			printf ("Cancel job_id=%u name=%s partition=%s [y/n]? ", 
			        job_ptr[i].job_id, job_ptr[i].name, 
				job_ptr[i].partition);
		} else {
			printf ("Cancel step_id=%u.%u name=%s partition=%s [y/n]? ", 
			        job_ptr[i].job_id, step_id, job_ptr[i].name, 
				job_ptr[i].partition);
		}

		fgets (in_line, sizeof (in_line), stdin);
		if ((in_line[0] == 'y') || (in_line[0] == 'Y'))
			return 1;
		if ((in_line[0] == 'n') || (in_line[0] == 'N'))
			return 0;
	}

}
