/*****************************************************************************\
 *  get_jobs.c - Process Wiki get job info request
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

#include <grp.h>
#include <sys/types.h>

#include "./msg.h"
#include "src/common/list.h"
#include "src/common/uid.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/slurmctld.h"

static char *	_dump_all_jobs(int *job_cnt, int state_info);
static char *	_dump_job(struct job_record *job_ptr, int state_info);
static char *	_get_group_name(gid_t gid);
static uint32_t	_get_job_end_time(struct job_record *job_ptr);
static uint32_t	_get_job_min_disk(struct job_record *job_ptr);
static uint32_t	_get_job_min_mem(struct job_record *job_ptr);
static uint32_t	_get_job_min_nodes(struct job_record *job_ptr);
static char *	_get_job_state(struct job_record *job_ptr);
static uint32_t	_get_job_submit_time(struct job_record *job_ptr);
static uint32_t	_get_job_suspend_time(struct job_record *job_ptr);
static uint32_t	_get_job_tasks(struct job_record *job_ptr);
static uint32_t	_get_job_time_limit(struct job_record *job_ptr);

#define SLURM_INFO_ALL		0
#define SLURM_INFO_VOLITILE	1
#define SLURM_INFO_STATE	2

/*
 * get_jobs - get information on specific job(s) changed since some time
 * cmd_ptr IN - CMD=GETJOBS ARG=[<UPDATETIME>:<JOBID>[:<JOBID>]...]
 *                              [<UPDATETIME>:ALL]
 * RET 0 on success, -1 on failure
 *
 * Response format
 * ARG=<cnt>#<JOBID>;UPDATE_TIME=<uts>;STATE=<state>;UCLIMIT=<time_limit>;
 *                    TASKS=<cpus>;QUEUETIME=<submit_time>;STARTTIME=<time>;
 *                    UNAME=<user>;GNAME=<group>;RCLASS=<part>;
 *                    NODES=<node_cnt>;RMEM=<mem_size>;RDISK=<disk_space>;
 *                    COMMENT=<comment>;[COMPLETETIME=<end_time>;]
 *         [#<JOBID>;...];
 */
/* RET 0 on success, -1 on failure */
extern int	get_jobs(char *cmd_ptr, int *err_code, char **err_msg)
{
	char *arg_ptr, *tmp_char, *tmp_buf, *buf = NULL;
	time_t update_time;
	/* Locks: read job, partition */
	slurmctld_lock_t job_read_lock = {
		NO_LOCK, READ_LOCK, NO_LOCK, READ_LOCK };
	int job_rec_cnt = 0, buf_size = 0, state_info;

	arg_ptr = strstr(cmd_ptr, "ARG=");
	if (arg_ptr == NULL) {
		*err_code = -300;
		*err_msg = "GETJOBS lacks ARG";
		error("wiki: GETJOBS lacks ARG");
		return -1;
	}
	update_time = (time_t) strtoul(arg_ptr+4, &tmp_char, 10);
	if (tmp_char[0] != ':') {
		*err_code = -300;
		*err_msg = "Invalid ARG value";
		error("wiki: GETJOBS has invalid ARG value");
		return -1;
	}
	tmp_char++;
	lock_slurmctld(job_read_lock);
	if (update_time == 0)
		state_info = SLURM_INFO_ALL;
	else if (update_time > last_job_update)
		state_info = SLURM_INFO_STATE;
	else
		state_info = SLURM_INFO_VOLITILE;

	if (strncmp(tmp_char, "ALL", 3) == 0) {
		/* report all jobs */
		buf = _dump_all_jobs(&job_rec_cnt, state_info);
	} else {
		struct job_record *job_ptr;
		char *job_name, *tmp2_char;
		uint32_t job_id;

		job_name = strtok_r(tmp_char, ":", &tmp2_char);
		while (job_name) {
			job_id = (uint32_t) strtoul(job_name, NULL, 10);
			job_ptr = find_job_record(job_id);
			tmp_buf = _dump_job(job_ptr, state_info);
			if (job_rec_cnt > 0)
				xstrcat(buf, "#");
			xstrcat(buf, tmp_buf);
			xfree(tmp_buf);
			job_rec_cnt++;
			job_name = strtok_r(NULL, ":", &tmp2_char);
		}
	}
	unlock_slurmctld(job_read_lock);

	/* Prepend ("ARG=%d", job_rec_cnt) to reply message */
	if (buf)
		buf_size = strlen(buf);
	tmp_buf = xmalloc(buf_size + 32);
	sprintf(tmp_buf, "SC=0 ARG=%d#%s", job_rec_cnt, buf);
	xfree(buf);
	*err_code = 0;
	*err_msg = tmp_buf;
	return 0;
}

static char *   _dump_all_jobs(int *job_cnt, int state_info)
{
	int cnt = 0;
	struct job_record *job_ptr;
	ListIterator job_iterator;
	char *tmp_buf, *buf = NULL;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		tmp_buf = _dump_job(job_ptr, state_info);
		if (cnt > 0)
			xstrcat(buf, "#");
		xstrcat(buf, tmp_buf);
		xfree(tmp_buf);
		cnt++;
	}
	*job_cnt = cnt;
	return buf;
}

static char *	_dump_job(struct job_record *job_ptr, int state_info)
{
	char tmp[512], *buf = NULL;
	uint32_t end_time, suspend_time;

	if (!job_ptr)
		return NULL;

	/* SLURM_INFO_STATE or SLURM_INFO_VOLITILE or SLURM_INFO_ALL */
	snprintf(tmp, sizeof(tmp), "%u:STATE=%s;",
		job_ptr->job_id, _get_job_state(job_ptr));
	xstrcat(buf, tmp);

	if (state_info == SLURM_INFO_STATE)
		return buf;

	/* SLURM_INFO_VOLITILE or SLURM_INFO_ALL */
	if ((job_ptr->job_state == JOB_PENDING)
	&&  (job_ptr->details)
	&&  (job_ptr->details->req_nodes)
	&&  (job_ptr->details->req_nodes[0])) {
		snprintf(tmp, sizeof(tmp),
			"HOSTLIST=%s;",
			bitmap2wiki_node_name(
			job_ptr->details->req_node_bitmap));
		xstrcat(buf, tmp);
	}

	snprintf(tmp, sizeof(tmp), 
		"UPDATETIME=%u;WCLIMIT=%u;",
		(uint32_t) job_ptr->time_last_active,
		(uint32_t) _get_job_time_limit(job_ptr));
	xstrcat(buf, tmp);

	snprintf(tmp, sizeof(tmp),
		"TASKS=%u;QUEUETIME=%u;STARTTIME=%u;",
		_get_job_tasks(job_ptr),
		_get_job_submit_time(job_ptr),
		(uint32_t) job_ptr->start_time);
	xstrcat(buf, tmp);

	snprintf(tmp, sizeof(tmp),
		"RCLASS=%s;NODES=%u;",
		job_ptr->partition,
		_get_job_min_nodes(job_ptr));
	xstrcat(buf, tmp);

	snprintf(tmp, sizeof(tmp),
		"RMEM=%u;RDISK=%u;",
		_get_job_min_mem(job_ptr),
		_get_job_min_disk(job_ptr));
	xstrcat(buf, tmp);

	if (job_ptr->comment && job_ptr->comment[0]) {
		snprintf(tmp, sizeof(tmp),
			"COMMENT=%s;", job_ptr->comment);
		xstrcat(buf, tmp);
	}

	end_time = _get_job_end_time(job_ptr);
	if (end_time) {
		snprintf(tmp, sizeof(tmp),
			"COMPLETETIME=%u;", end_time);
		xstrcat(buf, tmp);
	}

	suspend_time = _get_job_suspend_time(job_ptr);
	if (suspend_time) {
		snprintf(tmp, sizeof(tmp),
			"SUSPENDTIME=%u;", suspend_time);
		xstrcat(buf, tmp);
	}

	if (state_info == SLURM_INFO_VOLITILE)
		return buf;

	/* SLURM_INFO_ALL only */
	snprintf(tmp, sizeof(tmp),
		"UNAME=%s;GNAME=%s;",
		uid_to_string((uid_t) job_ptr->user_id),
		_get_group_name(job_ptr->group_id));
	xstrcat(buf, tmp);

	return buf;
}

static uint32_t _get_job_min_mem(struct job_record *job_ptr)
{
	if (job_ptr->details)
		return job_ptr->details->job_min_memory;
	return (uint32_t) 0;
}

static uint32_t _get_job_min_disk(struct job_record *job_ptr)
	
{
	if (job_ptr->details)
		return job_ptr->details->job_min_tmp_disk;
	return (uint32_t) 0;
}

static uint32_t	_get_job_min_nodes(struct job_record *job_ptr)
{
	if (job_ptr->details)
		return job_ptr->details->min_nodes;
	return (uint32_t) 1;
}

static char *	_get_group_name(gid_t gid)
{
	struct group *grp;

	grp = getgrgid(gid);
	if (grp)
		return grp->gr_name;
	return "nobody";
}

static uint32_t _get_job_submit_time(struct job_record *job_ptr)
{
	if (job_ptr->details)
		return (uint32_t) job_ptr->details->submit_time;
	return (uint32_t) 0;
}

static uint32_t _get_job_tasks(struct job_record *job_ptr)
{
	if (job_ptr->num_procs)
		return job_ptr->num_procs;
	return (uint32_t) 1;
}

static uint32_t	_get_job_time_limit(struct job_record *job_ptr)
{
	uint32_t limit = job_ptr->time_limit;

	if ((limit == NO_VAL) || (limit == INFINITE))
		return (uint32_t) 0;
	else
		return (limit * 60);	/* seconds, not minutes */
}

/* NOTE: if job has already completed, we append "EXITCODE=#" to 
 * the state name */
static char *	_get_job_state(struct job_record *job_ptr)
{
	uint16_t state = job_ptr->job_state;
	uint16_t base_state = state & (~JOB_COMPLETING);
	char *state_str;
	static char return_msg[128];

	if (base_state == JOB_PENDING)
		return "Idle";
	if (base_state == JOB_RUNNING)
		return "Running";

	if (state & JOB_COMPLETING) {
		/* Give 60 seconds to clear out, then 
		 * then consider job done. Let Moab 
		 * deal with inconsistency between 
		 * job state (DONE) and node state
		 * (some IDLE and others still 
		 * BUSY). */
		int age = (int) difftime(time(NULL), 
			job_ptr->end_time);
		if (age < 60)
			return "Running";
	}

	if (base_state == JOB_COMPLETE)
		state_str = "Completed";
	else if (base_state == JOB_SUSPENDED)
		state_str = "Suspended";
	else /* JOB_CANCELLED, JOB_FAILED, JOB_TIMEOUT, JOB_NODE_FAIL */
		state_str = "Removed";
	snprintf(return_msg, sizeof(return_msg), "%s;EXITCODE=%u",
		state_str, job_ptr->exit_code);
	return return_msg;
}

static uint32_t	_get_job_end_time(struct job_record *job_ptr)
{
	if (IS_JOB_FINISHED(job_ptr))
		return (uint32_t) job_ptr->end_time;
	return (uint32_t) 0;
}

/* returns how long job has been suspended, in seconds */
static uint32_t	_get_job_suspend_time(struct job_record *job_ptr)
{
	if (job_ptr->job_state == JOB_SUSPENDED) {
		time_t now = time(NULL);
		return (uint32_t) difftime(now, 
				job_ptr->suspend_time);
	}
	return (uint32_t) 0;
}
