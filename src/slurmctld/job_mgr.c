/*****************************************************************************\
 *  job_mgr.c - manage the job information of slurm
 *	Note: there is a global job list (job_list), time stamp 
 *	(last_job_update), and hash table (job_hash)
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <slurm/slurm_errno.h>

#include "src/api/job_info.h"
#include "src/common/bitstring.h"
#include "src/common/hostlist.h"
#include "src/common/node_select.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_jobcomp.h"
#include "src/common/switch.h"
#include "src/common/xassert.h"
#include "src/common/xstring.h"
#include "src/common/forward.h"
#include "src/common/slurm_jobacct.h"

#include "src/slurmctld/agent.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/node_scheduler.h"
#include "src/slurmctld/proc_req.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/sched_plugin.h"
#include "src/slurmctld/srun_comm.h"

#define BUFFER_SIZE 1024
#define DETAILS_FLAG 0xdddd
#define HUGE_BUF_SIZE (1024*16)
#define MAX_RETRIES  10
#define SLURM_CREATE_JOB_FLAG_NO_ALLOCATE_0 0
#define STEP_FLAG 0xbbbb
#define TOP_PRIORITY 0xffff0000	/* large, but leave headroom for higher */

#define JOB_HASH_INX(_job_id)	(_job_id % hash_table_size)

#define JOB_STATE_VERSION      "VER003"

/* Global variables */
List   job_list = NULL;		/* job_record list */
time_t last_job_update;		/* time of last update to job records */

/* Local variables */
static uint32_t maximum_prio = TOP_PRIORITY;
static int      hash_table_size = 0;
static int      job_count = 0;		/* job's in the system */
static uint32_t job_id_sequence = 0;	/* first job_id to assign new job */
static struct   job_record **job_hash = NULL;

/* Local functions */
static void _add_job_hash(struct job_record *job_ptr);
static int  _copy_job_desc_to_file(job_desc_msg_t * job_desc,
				   uint32_t job_id);
static int  _copy_job_desc_to_job_record(job_desc_msg_t * job_desc,
					 struct job_record **job_ptr,
					 struct part_record *part_ptr,
					 bitstr_t ** exc_bitmap,
					 bitstr_t ** req_bitmap);
static char *_copy_nodelist_no_dup(char *node_list);
static void _del_batch_list_rec(void *x);
static void _delete_job_desc_files(uint32_t job_id);
static void _dump_job_details(struct job_details *detail_ptr,
				    Buf buffer);
static void _dump_job_state(struct job_record *dump_job_ptr, Buf buffer);
static void _dump_job_step_state(struct step_record *step_ptr, Buf buffer);
static void _excise_node_from_job(struct job_record *job_ptr, 
				  struct node_record *node_ptr);
static int  _find_batch_dir(void *x, void *key);
static void _get_batch_job_dir_ids(List batch_dirs);
static void _job_timed_out(struct job_record *job_ptr);
static int  _job_create(job_desc_msg_t * job_specs, int allocate, int will_run,
		        struct job_record **job_rec_ptr, uid_t submit_uid);
static void _list_delete_job(void *job_entry);
static int  _list_find_job_id(void *job_entry, void *key);
static int  _list_find_job_old(void *job_entry, void *key);
static int  _load_job_details(struct job_record *job_ptr, Buf buffer);
static int  _load_job_state(Buf buffer);
static int  _load_step_state(struct job_record *job_ptr, Buf buffer);
static void _pack_default_job_details(struct job_details *detail_ptr, Buf buffer);
static void _pack_pending_job_details(struct job_details *detail_ptr, Buf buffer);
static int  _purge_job_record(uint32_t job_id);
static void _purge_lost_batch_jobs(int node_inx, time_t now);
static void _read_data_array_from_file(char *file_name, char ***data,
				       uint16_t * size);
static void _read_data_from_file(char *file_name, char **data);
static void _remove_defunct_batch_dirs(List batch_dirs);
static int  _reset_detail_bitmaps(struct job_record *job_ptr);
static void _reset_step_bitmaps(struct job_record *job_ptr);
static int  _resume_job_nodes(struct job_record *job_ptr);
static void _set_job_id(struct job_record *job_ptr);
static void _set_job_prio(struct job_record *job_ptr);
static void _signal_batch_job(struct job_record *job_ptr, uint16_t signal);
static void _signal_job(struct job_record *job_ptr, int signal);
static void _suspend_job(struct job_record *job_ptr, uint16_t op);
static int  _suspend_job_nodes(struct job_record *job_ptr);
static bool _top_priority(struct job_record *job_ptr);
static int  _validate_job_create_req(job_desc_msg_t * job_desc);
static int  _validate_job_desc(job_desc_msg_t * job_desc_msg, int allocate,
				uid_t submit_uid);
static void _validate_job_files(List batch_dirs);
static int  _write_data_to_file(char *file_name, char *data);
static int  _write_data_array_to_file(char *file_name, char **data,
				     uint16_t size);
static void _xmit_new_end_time(struct job_record *job_ptr);

/* 
 * create_job_record - create an empty job_record including job_details.
 *	load its values with defaults (zeros, nulls, and magic cookie)
 * IN/OUT error_code - set to zero if no error, errno otherwise
 * RET pointer to the record or NULL if error
 * global: job_list - global job list
 *	job_count - number of jobs in the system
 *	last_job_update - time of last job table update
 * NOTE: allocates memory that should be xfreed with _list_delete_job
 */
struct job_record *create_job_record(int *error_code)
{
	struct job_record  *job_ptr;
	struct job_details *detail_ptr;

	if (job_count >= slurmctld_conf.max_job_cnt) {
		error("create_job_record: job_count exceeds limit");
		*error_code = EAGAIN;
		return NULL;
	}

	job_count++;
	*error_code = 0;
	last_job_update = time(NULL);

	job_ptr    = (struct job_record *) xmalloc(sizeof(struct job_record));
	detail_ptr = (struct job_details *)xmalloc(sizeof(struct job_details));

	xassert (job_ptr->magic = JOB_MAGIC); /* sets value */
	job_ptr->details = detail_ptr;
	job_ptr->step_list = list_create(NULL);
	if (job_ptr->step_list == NULL)
		fatal("memory allocation failure");

	xassert (detail_ptr->magic = DETAILS_MAGIC); /* set value */
	detail_ptr->submit_time = time(NULL);

	if (list_append(job_list, job_ptr) == 0)
		fatal("list_append memory allocation failure");

	return job_ptr;
}


/* 
 * delete_job_details - delete a job's detail record and clear it's pointer
 *	this information can be deleted as soon as the job is allocated  
 *	resources and running (could need to restart batch job)
 * IN job_entry - pointer to job_record to clear the record of
 */
void delete_job_details(struct job_record *job_entry)
{
	int i;

	if (job_entry->details == NULL)
		return;

	_delete_job_desc_files(job_entry->job_id);
	xassert (job_entry->details->magic == DETAILS_MAGIC);
	for (i=0; i<job_entry->details->argc; i++)
		xfree(job_entry->details->argv[i]);
	xfree(job_entry->details->argv);
	xfree(job_entry->details->req_nodes);
	xfree(job_entry->details->exc_nodes);
	FREE_NULL_BITMAP(job_entry->details->req_node_bitmap);
	FREE_NULL_BITMAP(job_entry->details->exc_node_bitmap);
	xfree(job_entry->details->features);
	xfree(job_entry->details->err);
	xfree(job_entry->details->in);
	xfree(job_entry->details->out);
	xfree(job_entry->details->work_dir);
	xfree(job_entry->details);
}

/* _delete_job_desc_files - delete job descriptor related files */
static void _delete_job_desc_files(uint32_t job_id)
{
	char *dir_name, job_dir[20], *file_name;
	struct stat sbuf;

	dir_name = xstrdup(slurmctld_conf.state_save_location);

	sprintf(job_dir, "/job.%d", job_id);
	xstrcat(dir_name, job_dir);

	file_name = xstrdup(dir_name);
	xstrcat(file_name, "/environment");
	(void) unlink(file_name);
	xfree(file_name);

	file_name = xstrdup(dir_name);
	xstrcat(file_name, "/script");
	(void) unlink(file_name);
	xfree(file_name);

	if (stat(dir_name, &sbuf) == 0)	/* remove job directory as needed */
		(void) rmdir(dir_name);
	xfree(dir_name);
}

/* dump_all_job_state - save the state of all jobs to file for checkpoint
 * RET 0 or error code */
int dump_all_job_state(void)
{
	static int high_buffer_size = (1024 * 1024);
	int error_code = 0, log_fd;
	char *old_file, *new_file, *reg_file;
	/* Locks: Read config and job */
	slurmctld_lock_t job_read_lock =
	    { READ_LOCK, READ_LOCK, NO_LOCK, NO_LOCK };
	ListIterator job_iterator;
	struct job_record *job_ptr;
	Buf buffer = init_buf(high_buffer_size);
	DEF_TIMERS;

	START_TIMER;

        /*
         * write header: The version of the "job_state" file format.
         * Putting a version in the header comes in handy for cases where
         * we need to modify the format of the "job_state" file.
         */
	packstr( JOB_STATE_VERSION, buffer);

	/* write header: time */
	pack_time(time(NULL), buffer);

        /*
         * write header: job id
         * This is needed so that the job id remains persistent even after
         * slurmctld is restarted.
         */
	pack32( job_id_sequence, buffer);

	debug3("Writing job id %u to header record of job_state file",
	       job_id_sequence);

	/* write individual job records */
	lock_slurmctld(job_read_lock);
	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		xassert (job_ptr->magic == JOB_MAGIC);
		_dump_job_state(job_ptr, buffer);
	}
	/* Maintain config lock until we get the state_save_location *\
	\* unlock_slurmctld(job_read_lock);         - see below      */
	list_iterator_destroy(job_iterator);

	/* write the buffer to file */
	old_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(old_file, "/job_state.old");
	reg_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(reg_file, "/job_state");
	new_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(new_file, "/job_state.new");
	unlock_slurmctld(job_read_lock);

	lock_state_files();
	log_fd = creat(new_file, 0600);
	if (log_fd == 0) {
		error("Can't save state, create file %s error %m",
		      new_file);
		error_code = errno;
	} else {
		int pos = 0, nwrite = get_buf_offset(buffer), amount;
		char *data = (char *)get_buf_data(buffer);
		high_buffer_size = MAX(nwrite, high_buffer_size);
		while (nwrite > 0) {
			amount = write(log_fd, &data[pos], nwrite);
			if ((amount < 0) && (errno != EINTR)) {
				error("Error writing file %s, %m", new_file);
				error_code = errno;
				break;
			}
			nwrite -= amount;
			pos    += amount;
		}
		fsync(log_fd);
		close(log_fd);
	}
	if (error_code)
		(void) unlink(new_file);
	else {			/* file shuffle */
		(void) unlink(old_file);
		(void) link(reg_file, old_file);
		(void) unlink(reg_file);
		(void) link(new_file, reg_file);
		(void) unlink(new_file);
	}
	xfree(old_file);
	xfree(reg_file);
	xfree(new_file);
	unlock_state_files();

	free_buf(buffer);
	END_TIMER;
	debug3("dump_all_job_state %s", TIME_STR);
	return error_code;
}

/*
 * load_all_job_state - load the job state from file, recover from last 
 *	checkpoint. Execute this after loading the configuration file data.
 * RET 0 or error code
 */
int load_all_job_state(void)
{
	int data_allocated, data_read = 0, error_code = 0;
	uint32_t data_size = 0;
	int state_fd, job_cnt = 0;
	char *data = NULL, *state_file;
	Buf buffer;
	time_t buf_time;
	uint32_t saved_job_id;
	char *ver_str = NULL;
	uint16_t ver_str_len;

	/* read the file */
	state_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(state_file, "/job_state");
	lock_state_files();
	state_fd = open(state_file, O_RDONLY);
	if (state_fd < 0) {
		info("No job state file (%s) to recover", state_file);
		error_code = ENOENT;
	} else {
		data_allocated = HUGE_BUF_SIZE;
		data = xmalloc(data_allocated);
		while (1) {
			data_read = read(state_fd, &data[data_size],
					HUGE_BUF_SIZE);
			if (data_read < 0) {
				if (errno == EINTR)
					continue;
				else {
					error("Read error on %s: %m", 
						state_file);
					break;
				}
			} else if (data_read == 0)	/* eof */
				break;
			data_size      += data_read;
			data_allocated += data_read;
			xrealloc(data, data_allocated);
		}
		close(state_fd);
	}
	xfree(state_file);
	unlock_state_files();

	if (job_id_sequence == 0)
		job_id_sequence = slurmctld_conf.first_job_id;

	buffer = create_buf(data, data_size);

        /*
         * The old header of the "job_state" file simply contained a
         * timestamp, while the new header contains a "VERXXX" at the
         * beginning (VER001, VER002, etc), a timestamp, and the last
         * job id. To determine if we're looking at an old header or
         * new header, we first check if the file begins with "VER".
         *
         * Each field is preceeded by two bytes which contains the field
         * size.  Since we are bypassing the "pack" functions in order
         * see if the header contains a "VERXXX" string, we need to make
         * sure that there is enough data in the buffer to compare against.
         */
	if (size_buf(buffer) >= sizeof(uint16_t) + strlen(JOB_STATE_VERSION))
	{
	        char *ptr = get_buf_data(buffer);

	        if (memcmp( &ptr[sizeof(uint16_t)], JOB_STATE_VERSION, 3) == 0)
		{
		        safe_unpackstr_xmalloc( &ver_str, &ver_str_len, buffer);
		        debug3("Version string in job_state header is %s",
				ver_str);
		}
	}

	safe_unpack_time(&buf_time, buffer);

        /*
         * If the header has the version string then it also has the job id.
         */
	if (ver_str != NULL) {
	        safe_unpack32( &saved_job_id, buffer);
	        debug3("Job id in job_state header is %d", saved_job_id);
	}

	while (remaining_buf(buffer) > 0) {
		error_code = _load_job_state(buffer);
		if (error_code != SLURM_SUCCESS)
			goto unpack_error;
		job_cnt++;
	}

        /*
         * If the header has the version string then it also has the job id.
	 * Use MAX of preserved value or configuration parameter 
	 * FirstJobId (set above).
         */
	if (ver_str != NULL) {
		job_id_sequence = MAX(saved_job_id, job_id_sequence);
		debug3("Set job_id_sequence to %u", job_id_sequence);
	}

	free_buf(buffer);
	xfree(ver_str);
	info("Recovered state of %d jobs", job_cnt);
	return error_code;

      unpack_error:
	error("Incomplete job data checkpoint file");
	info("State of %d jobs recovered", job_cnt);
	free_buf(buffer);
	xfree(ver_str);
	return SLURM_FAILURE;
}

/*
 * _dump_job_state - dump the state of a specific job, its details, and 
 *	steps to a buffer
 * IN dump_job_ptr - pointer to job for which information is requested
 * IN/OUT buffer - location to store data, pointers automatically advanced
 */
static void _dump_job_state(struct job_record *dump_job_ptr, Buf buffer)
{
	struct job_details *detail_ptr;
	ListIterator step_iterator;
	struct step_record *step_ptr;

	/* Dump basic job info */
	pack32(dump_job_ptr->job_id, buffer);
	pack32(dump_job_ptr->user_id, buffer);
	pack32(dump_job_ptr->group_id, buffer);
	pack32(dump_job_ptr->time_limit, buffer);
	pack32(dump_job_ptr->priority, buffer);
	pack32(dump_job_ptr->alloc_sid, buffer);
	pack32(dump_job_ptr->dependency, buffer);
	pack32(dump_job_ptr->num_procs, buffer);
#if 0
//FIXME: Update in slurm v1.2
Update JOB_STATE_VERSION
	pack32(dump_job_ptr->exit_code, buffer);
#endif

	pack_time(dump_job_ptr->start_time, buffer);
	pack_time(dump_job_ptr->end_time, buffer);
	pack_time(dump_job_ptr->suspend_time, buffer);
        pack_time(dump_job_ptr->pre_sus_time, buffer);

	pack16((uint16_t) dump_job_ptr->job_state, buffer);
	pack16((uint16_t)dump_job_ptr->next_step_id, buffer);
	pack16((uint16_t)dump_job_ptr->kill_on_node_fail, buffer);
	pack16((uint16_t)dump_job_ptr->kill_on_step_done, buffer);
	pack16((uint16_t)dump_job_ptr->batch_flag, buffer);
	pack16((uint16_t)dump_job_ptr->alloc_resp_port, buffer);
	pack16((uint16_t)dump_job_ptr->other_port, buffer);
	pack16((uint16_t)dump_job_ptr->mail_type, buffer);

	packstr(dump_job_ptr->alloc_resp_host, buffer);
	packstr(dump_job_ptr->other_host, buffer);
	packstr(dump_job_ptr->nodes, buffer);
	packstr(dump_job_ptr->partition, buffer);
	packstr(dump_job_ptr->name, buffer);
	packstr(dump_job_ptr->alloc_node, buffer);
	packstr(dump_job_ptr->account, buffer);
	packstr(dump_job_ptr->comment, buffer);
	packstr(dump_job_ptr->network, buffer);
	packstr(dump_job_ptr->mail_user, buffer);

	select_g_pack_jobinfo(dump_job_ptr->select_jobinfo,
		buffer);

	/* Dump job details, if available */
	detail_ptr = dump_job_ptr->details;
	if (detail_ptr) {
		xassert (detail_ptr->magic == DETAILS_MAGIC);
		pack16((uint16_t) DETAILS_FLAG, buffer);
		_dump_job_details(detail_ptr, buffer);
	} else
		pack16((uint16_t) 0, buffer);	/* no details flag */

	/* Dump job steps */
	step_iterator = list_iterator_create(dump_job_ptr->step_list);
	while ((step_ptr = (struct step_record *) 
				list_next(step_iterator))) {
		pack16((uint16_t) STEP_FLAG, buffer);
		_dump_job_step_state(step_ptr, buffer);
	}
	list_iterator_destroy(step_iterator);
	pack16((uint16_t) 0, buffer);	/* no step flag */
}

/* Unpack a job's state information from a buffer */
static int _load_job_state(Buf buffer)
{
	uint32_t job_id, user_id, group_id, time_limit, priority, alloc_sid;
	uint32_t dependency, num_procs;
	time_t start_time, end_time, suspend_time, pre_sus_time;
	uint16_t job_state, next_step_id, details, batch_flag, step_flag;
	uint16_t kill_on_node_fail, kill_on_step_done, name_len;
	uint16_t alloc_resp_port, other_port, mail_type;
	char *nodes = NULL, *partition = NULL, *name = NULL;
	char *alloc_node = NULL, *alloc_resp_host = NULL, *other_host = NULL;
	char *account = NULL, *network = NULL, *mail_user = NULL;
	char *comment = NULL;
	struct job_record *job_ptr;
	struct part_record *part_ptr;
	int error_code;
	select_jobinfo_t select_jobinfo = NULL;

	safe_unpack32(&job_id, buffer);
	safe_unpack32(&user_id, buffer);
	safe_unpack32(&group_id, buffer);
	safe_unpack32(&time_limit, buffer);
	safe_unpack32(&priority, buffer);
	safe_unpack32(&alloc_sid, buffer);
	safe_unpack32(&dependency, buffer);
	safe_unpack32(&num_procs, buffer);
#if 0
//FIXME: Update in slurm v1.2
//Also expose for get_job_info
uint32_t exit_code;
	safe_unpack32(&exit_code, buffer);
job_ptr->exit_code = exit_code;
#endif

	safe_unpack_time(&start_time, buffer);
	safe_unpack_time(&end_time, buffer);
	safe_unpack_time(&suspend_time, buffer);
	safe_unpack_time(&pre_sus_time, buffer);

	safe_unpack16(&job_state, buffer);
	safe_unpack16(&next_step_id, buffer);
	safe_unpack16(&kill_on_node_fail, buffer);
	safe_unpack16(&kill_on_step_done, buffer);
	safe_unpack16(&batch_flag, buffer);
	safe_unpack16(&alloc_resp_port, buffer);
	safe_unpack16(&other_port, buffer);
	safe_unpack16(&mail_type, buffer);

	safe_unpackstr_xmalloc(&alloc_resp_host, &name_len, buffer);
	safe_unpackstr_xmalloc(&other_host, &name_len, buffer);
	safe_unpackstr_xmalloc(&nodes, &name_len, buffer);
	safe_unpackstr_xmalloc(&partition, &name_len, buffer);
	safe_unpackstr_xmalloc(&name, &name_len, buffer);
	safe_unpackstr_xmalloc(&alloc_node, &name_len, buffer);
	safe_unpackstr_xmalloc(&account, &name_len, buffer);
	safe_unpackstr_xmalloc(&comment, &name_len, buffer);
	safe_unpackstr_xmalloc(&network, &name_len, buffer);
	safe_unpackstr_xmalloc(&mail_user, &name_len, buffer);

	if (select_g_alloc_jobinfo(&select_jobinfo)
	||  select_g_unpack_jobinfo(select_jobinfo, buffer))
		goto unpack_error;

	/* validity test as possible */
	if (((job_state & (~JOB_COMPLETING)) >= JOB_END) || 
	    (batch_flag > 2)) {
		error("Invalid data for job %u: job_state=%u batch_flag=%u",
		      job_id, job_state, batch_flag);
		goto unpack_error;
	}
	if (kill_on_step_done > KILL_ON_STEP_DONE) {
		error("Invalid data for job %u: kill_on_step_done=%u",
		      job_id, kill_on_step_done);
		goto unpack_error;
	}
	if (kill_on_node_fail > 1) {
		error("Invalid data for job %u: kill_on_node_fail=%u",
		      job_id, kill_on_node_fail);
		goto unpack_error;
	}
	part_ptr = find_part_record (partition);
	if (part_ptr == NULL) {
		verbose("Invalid partition (%s) for job_id %u", 
			partition, job_id);
		/* not a fatal error, partition could have been removed, 
		 * reset_job_bitmaps() will clean-up this job */
	}

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		job_ptr = create_job_record(&error_code);
		if (error_code) {
			error("Create job entry failed for job_id %u",
			      job_id);
			goto unpack_error;
		}
		job_ptr->job_id = job_id;
		_add_job_hash(job_ptr);
	}

	if ((maximum_prio >= priority) && (priority > 1))
		maximum_prio = priority;
	if (job_id_sequence <= job_id)
		job_id_sequence = job_id + 1;

	safe_unpack16(&details, buffer);
	if ((details == DETAILS_FLAG) && 
	    (_load_job_details(job_ptr, buffer))) {
		job_ptr->job_state = JOB_FAILED;
		job_ptr->exit_code = 1;
		job_ptr->end_time = time(NULL);
		goto unpack_error;
	}

	job_ptr->user_id      = user_id;
	job_ptr->group_id     = group_id;
	job_ptr->time_limit   = time_limit;
	job_ptr->priority     = priority;
	job_ptr->alloc_sid    = alloc_sid;
	job_ptr->start_time   = start_time;
	job_ptr->end_time     = end_time;
	job_ptr->suspend_time = suspend_time;
	job_ptr->pre_sus_time = pre_sus_time;
	job_ptr->job_state    = job_state;
	job_ptr->next_step_id = next_step_id;
	job_ptr->dependency   = dependency;
	job_ptr->num_procs    = num_procs;
	job_ptr->time_last_active = time(NULL);
	strncpy(job_ptr->name, name, MAX_JOBNAME_LEN);
	xfree(name);
	xfree(job_ptr->nodes);
	job_ptr->nodes  = nodes;
	nodes           = NULL;	/* reused, nothing left to free */
	xfree(job_ptr->alloc_node);
	job_ptr->alloc_node = alloc_node;
	alloc_node          = NULL;	/* reused, nothing left to free */
	strncpy(job_ptr->partition, partition, MAX_SLURM_NAME);
	xfree(partition);
	job_ptr->account = account;
	account          = NULL;  /* reused, nothing left to free */
	job_ptr->comment = comment;
	comment          = NULL;  /* reused, nothing left to free */
	job_ptr->network = network;
	network          = NULL;  /* reused, nothing left to free */
	job_ptr->part_ptr = part_ptr;
	job_ptr->kill_on_node_fail = kill_on_node_fail;
	job_ptr->kill_on_step_done = kill_on_step_done;
	job_ptr->batch_flag        = batch_flag;
	job_ptr->alloc_resp_port   = alloc_resp_port;
	job_ptr->alloc_resp_host   = alloc_resp_host;
	job_ptr->other_port        = other_port;
	job_ptr->other_host        = other_host;
	job_ptr->mail_type         = mail_type;
	job_ptr->mail_user         = mail_user;
	mail_user = NULL;	/* reused, nothing left to free */
	job_ptr->select_jobinfo = select_jobinfo;

	build_node_details(job_ptr);	/* set: num_cpu_groups, cpus_per_node, 
					 *	cpu_count_reps, node_cnt, and
					 *	node_addr */
	info("recovered job id %u", job_id);

	safe_unpack16(&step_flag, buffer);
	while (step_flag == STEP_FLAG) {
		if ((error_code = _load_step_state(job_ptr, buffer)))
			goto unpack_error;
		safe_unpack16(&step_flag, buffer);
	}

	return SLURM_SUCCESS;

      unpack_error:
	error("Incomplete job record");
	xfree(alloc_resp_host);
	xfree(other_host);
	xfree(nodes);
	xfree(partition);
	xfree(name);
	xfree(alloc_node);
	xfree(account);
	xfree(comment);
	xfree(mail_user);
	select_g_free_jobinfo(&select_jobinfo);
	return SLURM_FAILURE;
}

/*
 * _dump_job_details - dump the state of a specific job details to 
 *	a buffer
 * IN detail_ptr - pointer to job details for which information is requested
 * IN/OUT buffer - location to store data, pointers automatically advanced
 */
void _dump_job_details(struct job_details *detail_ptr, Buf buffer)
{
	pack32((uint32_t) detail_ptr->min_nodes, buffer);
	pack32((uint32_t) detail_ptr->max_nodes, buffer);
	pack32((uint32_t) detail_ptr->min_sockets, buffer);
	pack32((uint32_t) detail_ptr->max_sockets, buffer);
	pack32((uint32_t) detail_ptr->min_cores, buffer);
	pack32((uint32_t) detail_ptr->max_cores, buffer);
	pack32((uint32_t) detail_ptr->min_threads, buffer);
	pack32((uint32_t) detail_ptr->max_threads, buffer);
	pack32((uint32_t) detail_ptr->total_procs, buffer);
	pack32((uint32_t) detail_ptr->num_tasks, buffer);

	pack16((uint16_t) detail_ptr->shared, buffer);
	pack16((uint16_t) detail_ptr->contiguous, buffer);
	pack16((uint16_t) detail_ptr->cpus_per_task, buffer);
	pack16((uint16_t) detail_ptr->ntasks_per_node, buffer);
	pack16((uint16_t) detail_ptr->ntasks_per_socket, buffer);
	pack16((uint16_t) detail_ptr->ntasks_per_core, buffer);
	pack16((uint16_t) detail_ptr->no_requeue, buffer);
	pack16((uint16_t) detail_ptr->overcommit, buffer);

	pack32((uint32_t) detail_ptr->job_min_procs, buffer);
	pack32((uint32_t) detail_ptr->job_min_sockets, buffer);
	pack32((uint32_t) detail_ptr->job_min_cores, buffer);
	pack32((uint32_t) detail_ptr->job_min_threads, buffer);
	pack32((uint32_t) detail_ptr->job_min_memory, buffer);
	pack32((uint32_t) detail_ptr->job_max_memory, buffer);
	pack32((uint32_t) detail_ptr->job_min_tmp_disk, buffer);
	pack_time(detail_ptr->begin_time, buffer);
	pack_time(detail_ptr->submit_time, buffer);

	packstr(detail_ptr->req_nodes, buffer);
	packstr(detail_ptr->exc_nodes, buffer);
	packstr(detail_ptr->features,  buffer);

	packstr(detail_ptr->err,       buffer);
	packstr(detail_ptr->in,        buffer);
	packstr(detail_ptr->out,       buffer);
	packstr(detail_ptr->work_dir,  buffer);

	packstr_array(detail_ptr->argv, detail_ptr->argc, buffer);
}

/* _load_job_details - Unpack a job details information from buffer */
static int _load_job_details(struct job_record *job_ptr, Buf buffer)
{
	char *req_nodes = NULL, *exc_nodes = NULL, *features = NULL;
	char *err = NULL, *in = NULL, *out = NULL, *work_dir = NULL;
	char **argv = (char **) NULL;
	uint32_t min_nodes, max_nodes, min_sockets, max_sockets;
	uint32_t min_cores, max_cores, min_threads, max_threads;
	uint32_t job_min_procs, job_min_sockets, job_min_cores, job_min_threads;
	uint32_t job_min_memory, job_max_memory, job_min_tmp_disk;
	uint32_t num_tasks;
	uint16_t argc = 0, shared, contiguous;
	uint16_t cpus_per_task, name_len, no_requeue, overcommit;
	uint16_t ntasks_per_node, ntasks_per_socket, ntasks_per_core;
	uint32_t total_procs;
	time_t begin_time, submit_time;
	int i;

	/* unpack the job's details from the buffer */
	safe_unpack32(&min_nodes, buffer);
	safe_unpack32(&max_nodes, buffer);
	safe_unpack32(&min_sockets, buffer);
	safe_unpack32(&max_sockets, buffer);
	safe_unpack32(&min_cores, buffer);
	safe_unpack32(&max_cores, buffer);
	safe_unpack32(&min_threads, buffer);
	safe_unpack32(&max_threads, buffer);
	safe_unpack32(&total_procs, buffer);
	safe_unpack32(&num_tasks, buffer);

	safe_unpack16(&shared, buffer);
	safe_unpack16(&contiguous, buffer);
	safe_unpack16(&cpus_per_task, buffer);
	safe_unpack16(&ntasks_per_node, buffer);
	safe_unpack16(&ntasks_per_socket, buffer);
	safe_unpack16(&ntasks_per_core, buffer);
	safe_unpack16(&no_requeue, buffer);
	safe_unpack16(&overcommit, buffer);

	safe_unpack32(&job_min_procs, buffer);
	safe_unpack32(&job_min_sockets, buffer);
	safe_unpack32(&job_min_cores, buffer);
	safe_unpack32(&job_min_threads, buffer);
	safe_unpack32(&job_min_memory, buffer);
	safe_unpack32(&job_max_memory, buffer);
	safe_unpack32(&job_min_tmp_disk, buffer);
	safe_unpack_time(&begin_time, buffer);
	safe_unpack_time(&submit_time, buffer);

	safe_unpackstr_xmalloc(&req_nodes, &name_len, buffer);
	safe_unpackstr_xmalloc(&exc_nodes, &name_len, buffer);
	safe_unpackstr_xmalloc(&features,  &name_len, buffer);

	safe_unpackstr_xmalloc(&err, &name_len, buffer);
	safe_unpackstr_xmalloc(&in,  &name_len, buffer);
	safe_unpackstr_xmalloc(&out, &name_len, buffer);
	safe_unpackstr_xmalloc(&work_dir, &name_len, buffer);

	safe_unpackstr_array(&argv, &argc, buffer);

	/* validity test as possible */
	if (contiguous > 1) {
		error("Invalid data for job %u: contiguous=%u",
			job_ptr->job_id, contiguous);
		goto unpack_error;
	}
	if ((no_requeue > 1) || (overcommit > 1)) {
		error("Invalid data for job %u: no_requeue=%u overcommit=%u",
			no_requeue, overcommit);
		goto unpack_error;
	}



	/* free any left-over detail data */
	xfree(job_ptr->details->req_nodes);
	xfree(job_ptr->details->exc_nodes);
	xfree(job_ptr->details->features);
	xfree(job_ptr->details->err);
	xfree(job_ptr->details->in);
	xfree(job_ptr->details->out);
	xfree(job_ptr->details->work_dir);
	for (i=0; i<job_ptr->details->argc; i++)
		xfree(job_ptr->details->argv[i]);
	xfree(job_ptr->details->argv);

	/* now put the details into the job record */
	job_ptr->details->min_nodes = min_nodes;
	job_ptr->details->max_nodes = max_nodes;
	job_ptr->details->min_sockets = min_sockets;
	job_ptr->details->max_sockets = max_sockets;
	job_ptr->details->min_cores = min_cores;
	job_ptr->details->max_cores = max_cores;
	job_ptr->details->min_threads = min_threads;
	job_ptr->details->max_threads = max_threads;
	job_ptr->details->total_procs = total_procs;
	job_ptr->details->num_tasks = num_tasks;
	job_ptr->details->shared = shared;
	job_ptr->details->contiguous = contiguous;
	job_ptr->details->cpus_per_task = cpus_per_task;
	job_ptr->details->ntasks_per_node = ntasks_per_node;
	job_ptr->details->ntasks_per_socket = ntasks_per_socket;
	job_ptr->details->ntasks_per_core = ntasks_per_core;
	job_ptr->details->job_min_procs = job_min_procs;
	job_ptr->details->job_min_sockets = job_min_sockets;
	job_ptr->details->job_min_cores = job_min_cores;
	job_ptr->details->job_min_threads = job_min_threads;
	job_ptr->details->job_min_memory = job_min_memory;
	job_ptr->details->job_max_memory = job_max_memory;
	job_ptr->details->job_min_tmp_disk = job_min_tmp_disk;
	job_ptr->details->no_requeue = no_requeue;
	job_ptr->details->overcommit = overcommit;
	job_ptr->details->begin_time = begin_time;
	job_ptr->details->submit_time = submit_time;
	job_ptr->details->req_nodes = req_nodes;
	job_ptr->details->exc_nodes = exc_nodes;
	job_ptr->details->features = features;
	job_ptr->details->err = err;
	job_ptr->details->in = in;
	job_ptr->details->out = out;
	job_ptr->details->work_dir = work_dir;
	job_ptr->details->argc = argc;
	job_ptr->details->argv = argv;

	return SLURM_SUCCESS;

      unpack_error:
	xfree(req_nodes);
	xfree(exc_nodes);
	xfree(features);
	xfree(err);
	xfree(in);
	xfree(out);
	xfree(work_dir);
/*	for (i=0; i<argc; i++) 
		xfree(argv[i]);  Don't trust this on unpack error */
	xfree(argv);
	return SLURM_FAILURE;
}


/*
 * _dump_job_step_state - dump the state of a specific job step to a buffer
 * IN detail_ptr - pointer to job step for which information is requested
 * IN/OUT buffer - location to store data, pointers automatically advanced
 */
static void _dump_job_step_state(struct step_record *step_ptr, Buf buffer)
{
	pack16((uint16_t) step_ptr->step_id, buffer);
	pack16((uint16_t) step_ptr->cyclic_alloc, buffer);
	pack16((uint16_t)step_ptr->port, buffer);
	pack32(step_ptr->exit_code, buffer);
	if (step_ptr->exit_code != NO_VAL) {
		pack_bit_fmt(step_ptr->exit_node_bitmap, buffer);
		pack16((uint16_t) _bitstr_bits(step_ptr->exit_node_bitmap), 
			buffer);
	}

	pack_time(step_ptr->start_time, buffer);
	packstr(step_ptr->host,  buffer);
	packstr(step_ptr->name, buffer);
	packstr(step_ptr->network, buffer);
	pack16((uint16_t)step_ptr->batch_step, buffer);
	if (!step_ptr->batch_step) {
		pack_slurm_step_layout(step_ptr->step_layout, buffer);
		switch_pack_jobinfo(step_ptr->switch_job, buffer);
	}
	checkpoint_pack_jobinfo(step_ptr->check_job, buffer);
}

/* Unpack job step state information from a buffer */
static int _load_step_state(struct job_record *job_ptr, Buf buffer)
{
	struct step_record *step_ptr = NULL;
	uint16_t step_id, cyclic_alloc, name_len, port, batch_step, bit_cnt;
	uint32_t exit_code;
	time_t start_time;
	char *host = NULL;
	char *name = NULL, *network = NULL, *bit_fmt = NULL;
	switch_jobinfo_t switch_tmp = NULL;
	check_jobinfo_t check_tmp = NULL;
	slurm_step_layout_t *step_layout = NULL;
	
	safe_unpack16(&step_id, buffer);
	safe_unpack16(&cyclic_alloc, buffer);
	safe_unpack16(&port, buffer);
	safe_unpack32(&exit_code, buffer);
	if (exit_code != NO_VAL) {
		safe_unpackstr_xmalloc(&bit_fmt, &name_len, buffer);
		safe_unpack16(&bit_cnt, buffer);
	}
	
	safe_unpack_time(&start_time, buffer);
	safe_unpackstr_xmalloc(&host, &name_len, buffer);
	safe_unpackstr_xmalloc(&name, &name_len, buffer);
	safe_unpackstr_xmalloc(&network, &name_len, buffer);
	safe_unpack16(&batch_step, buffer);
	if (!batch_step) {
		if (unpack_slurm_step_layout(&step_layout, buffer))
			goto unpack_error;
		switch_alloc_jobinfo(&switch_tmp);
        	if (switch_unpack_jobinfo(switch_tmp, buffer))
                	goto unpack_error;
	}
	checkpoint_alloc_jobinfo(&check_tmp);
        if (checkpoint_unpack_jobinfo(check_tmp, buffer))
                goto unpack_error;

	/* validity test as possible */
	if (cyclic_alloc > 1) {
		error("Invalid data for job %u.%u: cyclic_alloc=%u",
		      job_ptr->job_id, step_id, cyclic_alloc);
		goto unpack_error;
	}

	step_ptr = find_step_record(job_ptr, step_id);
	if (step_ptr == NULL)
		step_ptr = create_step_record(job_ptr);
	if (step_ptr == NULL)
		goto unpack_error;

	/* set new values */
	step_ptr->step_id      = step_id;
	step_ptr->cyclic_alloc = cyclic_alloc;
	step_ptr->name         = name;
	step_ptr->network      = network;
	step_ptr->port         = port;
	step_ptr->host         = host;
	step_ptr->batch_step   = batch_step;
	host                   = NULL;  /* re-used, nothing left to free */
	step_ptr->start_time   = start_time;

	slurm_step_layout_destroy(step_ptr->step_layout);
	step_ptr->step_layout = step_layout;
	
	step_ptr->switch_job   = switch_tmp;
	step_ptr->check_job    = check_tmp;

	step_ptr->exit_code    = exit_code;
	if (bit_fmt) {
		/* NOTE: This is only recovered if a job step completion
		 * is actively in progress at step save time. Otherwise
		 * the bitmap is NULL. */ 
		step_ptr->exit_node_bitmap = bit_alloc(bit_cnt);
		if (step_ptr->exit_node_bitmap == NULL)
			fatal("bit_alloc: %m");
		if (bit_unfmt(step_ptr->exit_node_bitmap, bit_fmt)) {
			error("error recovering exit_node_bitmap from %s",
				bit_fmt);
		}
		xfree(bit_fmt);
	}

	switch_g_job_step_allocated(switch_tmp, 
				    step_ptr->step_layout->node_list);
	info("recovered job step %u.%u", job_ptr->job_id, step_id);
	return SLURM_SUCCESS;

      unpack_error:
	xfree(host);
	xfree(name);
	xfree(network);
	xfree(bit_fmt);
	if (switch_tmp)
		switch_free_jobinfo(switch_tmp);
	slurm_step_layout_destroy(step_layout);
	return SLURM_FAILURE;
}

/* _add_job_hash - add a job hash entry for given job record, job_id must  
 *	already be set
 * IN job_ptr - pointer to job record
 * Globals: hash table updated
 */
void _add_job_hash(struct job_record *job_ptr)
{
	int inx;

	inx = JOB_HASH_INX(job_ptr->job_id);
	job_ptr->job_next = job_hash[inx];
	job_hash[inx] = job_ptr;
}


/* 
 * find_job_record - return a pointer to the job record with the given job_id
 * IN job_id - requested job's id
 * RET pointer to the job's record, NULL on error
 * global: job_list - global job list pointer
 *	job_hash - hash table into job records
 */
struct job_record *find_job_record(uint32_t job_id)
{
	struct job_record *job_ptr;

	job_ptr = job_hash[JOB_HASH_INX(job_id)];
	while (job_ptr) {
		if (job_ptr->job_id == job_id)
			return job_ptr;
		job_ptr = job_ptr->job_next;
	}

	return NULL;
}

/*
 * kill_job_by_part_name - Given a partition name, deallocate resource for 
 *	its jobs and kill them. All jobs associated with this partition 
 *	will have their partition pointer cleared.
 * IN part_name - name of a partition
 * RET number of jobs associated with this partition
 */
extern int kill_job_by_part_name(char *part_name)
{
	ListIterator job_iterator;
	struct job_record  *job_ptr;
	struct part_record *part_ptr;
	int job_count = 0;

	part_ptr = find_part_record (part_name);
	if (part_ptr == NULL)	/* No such partition */
		return 0;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		bool suspended = false;
		if (job_ptr->part_ptr != part_ptr)
			continue;
		job_ptr->part_ptr = NULL;

		if (job_ptr->job_state == JOB_SUSPENDED)
			suspended = true;
		if ((job_ptr->job_state == JOB_RUNNING) || suspended) {
			job_count++;
			info("Killing job_id %u on defunct partition %s",
			      job_ptr->job_id, part_name);
			job_ptr->job_state = JOB_NODE_FAIL | JOB_COMPLETING;
			job_ptr->exit_code = MAX(job_ptr->exit_code, 1);
			if (suspended)
				job_ptr->end_time = job_ptr->suspend_time;
			else
				job_ptr->end_time = time(NULL);
			job_completion_logger(job_ptr);
			deallocate_nodes(job_ptr, false, suspended);
		}

	}
	list_iterator_destroy(job_iterator);

	if (job_count)
		last_job_update = time(NULL);
	return job_count;
}

/*
 * kill_running_job_by_node_name - Given a node name, deallocate RUNNING 
 *	or COMPLETING jobs from the node or kill them 
 * IN node_name - name of a node
 * IN step_test - if true, only kill the job if a step is running on the node
 * RET number of killed jobs
 */
extern int kill_running_job_by_node_name(char *node_name, bool step_test)
{
	ListIterator job_iterator;
	struct job_record *job_ptr;
	struct node_record *node_ptr;
	int bit_position;
	int job_count = 0;

	node_ptr = find_node_record(node_name);
	if (node_ptr == NULL)	/* No such node */
		return 0;
	bit_position = node_ptr - node_record_table_ptr;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		bool suspended = false;
		if ((job_ptr->node_bitmap == NULL) ||
		    (!bit_test(job_ptr->node_bitmap, bit_position)))
			continue;	/* job not on this node */
		if (job_ptr->job_state == JOB_SUSPENDED)
			suspended = true;
		if (job_ptr->job_state & JOB_COMPLETING) {
			job_count++;
			bit_clear(job_ptr->node_bitmap, bit_position);
			if (job_ptr->node_cnt)
				(job_ptr->node_cnt)--;
			else
				error("node_cnt underflow on JobId=%u", 
			   	      job_ptr->job_id);
			if (job_ptr->node_cnt == 0) {
				job_ptr->job_state &= (~JOB_COMPLETING);
				slurm_sched_schedule();
			}
			if (node_ptr->comp_job_cnt)
				(node_ptr->comp_job_cnt)--;
			else
				error("Node %s comp_job_cnt underflow, "
					"JobId=%u", 
					node_ptr->name, job_ptr->job_id);
		} else if ((job_ptr->job_state == JOB_RUNNING) || suspended) {
			if (step_test && 
			    (step_on_node(job_ptr, node_ptr) == 0))
				continue;

			job_count++;
			srun_node_fail(job_ptr->job_id, node_name);
			if ((job_ptr->details == NULL) ||
			    (job_ptr->kill_on_node_fail) ||
			    (job_ptr->node_cnt <= 1)) {
				error("Killing job_id %u on failed node %s",
				      job_ptr->job_id, node_name);
				job_ptr->job_state = JOB_NODE_FAIL | 
						     JOB_COMPLETING;
				job_ptr->exit_code = MAX(job_ptr->exit_code, 1);
				if (suspended)
					job_ptr->end_time = job_ptr->suspend_time;
				else
					job_ptr->end_time = time(NULL);
				job_completion_logger(job_ptr);
				deallocate_nodes(job_ptr, false, suspended);
			} else {
				error("Removing failed node %s from job_id %u",
				      node_name, job_ptr->job_id);
				_excise_node_from_job(job_ptr, node_ptr);
			}
		}

	}
	list_iterator_destroy(job_iterator);
	if (job_count)
		last_job_update = time(NULL);

	return job_count;
}

/* Remove one node from a job's allocation */
static void _excise_node_from_job(struct job_record *job_ptr, 
				  struct node_record *node_ptr)
{
	make_node_idle(node_ptr, job_ptr); /* updates bitmap */
	job_ptr->nodes = bitmap2node_name(job_ptr->node_bitmap);
	xfree(job_ptr->cpus_per_node);
	xfree(job_ptr->cpu_count_reps);
	xfree(job_ptr->node_addr);

	/* build_node_details rebuilds everything from node_bitmap */
	build_node_details(job_ptr);
}


/*
 * dump_job_desc - dump the incoming job submit request message
 * IN job_specs - job specification from RPC
 */
void dump_job_desc(job_desc_msg_t * job_specs)
{
	long job_id;
	long job_min_procs, job_min_sockets, job_min_cores, job_min_threads;
	long job_min_memory, job_max_memory, job_min_tmp_disk, num_procs;
	long time_limit, priority, contiguous;
	long kill_on_node_fail, shared, immediate, dependency;
	long cpus_per_task, no_requeue, num_tasks, overcommit;
	long ntasks_per_node, ntasks_per_socket, ntasks_per_core;
	char buf[100];

	if (job_specs == NULL)
		return;

	job_id = (job_specs->job_id != NO_VAL) ? 
			(long) job_specs->job_id : -1L;
	debug3("JobDesc: user_id=%u job_id=%ld partition=%s name=%s",
	       job_specs->user_id, job_id,
	       job_specs->partition, job_specs->name);

	num_procs = (job_specs->num_procs != NO_VAL) ? 
			(long) job_specs->num_procs : -1L;
	debug3("   num_procs=%ld", num_procs);

	debug3("   -N min-[max]: %d-[%d]:%d-[%d]:%d-[%d]:%d-[%d]",
		job_specs->min_nodes,   job_specs->max_nodes,
		job_specs->min_sockets, job_specs->max_sockets,
		job_specs->min_cores,   job_specs->max_cores,
		job_specs->min_threads, job_specs->max_threads);

	job_min_procs    = (job_specs->job_min_procs != NO_VAL) ? 
			(long) job_specs->job_min_procs : -1L;
	job_min_sockets  = (job_specs->job_min_sockets != NO_VAL) ? 
			(long) job_specs->job_min_sockets : -1L;
	job_min_cores    = (job_specs->job_min_cores != NO_VAL) ? 
			(long) job_specs->job_min_cores : -1L;
	job_min_threads  = (job_specs->job_min_threads != NO_VAL) ? 
			(long) job_specs->job_min_threads : -1L;
	debug3("   job_min_procs=%ld job_min_sockets=%ld",
	       job_min_procs, job_min_sockets);
	debug3("   job_min_cores=%ld job_min_threads=%ld",
	       job_min_cores, job_min_threads);

	job_min_memory   = (job_specs->job_min_memory != NO_VAL) ? 
			(long) job_specs->job_min_memory : -1L;
	job_max_memory   = (job_specs->job_max_memory != NO_VAL) ? 
			(long) job_specs->job_max_memory : -1L;
	job_min_tmp_disk = (job_specs->job_min_tmp_disk != NO_VAL) ? 
			(long) job_specs->job_min_tmp_disk : -1L;
	debug3("   job_min_memory=%ld job_max_memory=%ld job_min_tmp_disk=%ld",
	       job_min_memory, job_max_memory, job_min_tmp_disk);
	immediate = (job_specs->immediate == 0) ? 0L : 1L;
	debug3("   immediate=%ld features=%s",
		immediate, job_specs->features);

	debug3("   req_nodes=%s exc_nodes=%s", 
	       job_specs->req_nodes, job_specs->exc_nodes);

	time_limit = (job_specs->time_limit != NO_VAL) ? 
			(long) job_specs->time_limit : -1L;
	priority   = (job_specs->priority != NO_VAL) ? 
			(long) job_specs->priority : -1L;
	contiguous = (job_specs->contiguous != (uint16_t) NO_VAL) ? 
			(long) job_specs->contiguous : -1L;
	shared = (job_specs->shared != (uint16_t) NO_VAL) ? 
			(long) job_specs->shared : -1L;
	debug3("   time_limit=%ld priority=%ld contiguous=%ld shared=%ld",
	       time_limit, priority, contiguous, shared);

	kill_on_node_fail = (job_specs->kill_on_node_fail != 
			     (uint16_t) NO_VAL) ? 
			(long) job_specs->kill_on_node_fail : -1L;
	if (job_specs->script)	/* log has problem with string len & null */
		debug3("   kill_on_node_fail=%ld script=%.40s...",
			kill_on_node_fail, job_specs->script);
	else
		debug3("   kill_on_node_fail=%ld script=%s",
			kill_on_node_fail, job_specs->script);

	if (job_specs->argc == 1)
		debug3("   argv=\"%s\"", 
			job_specs->argv[0]);
	else if (job_specs->argc == 2)
		debug3("   argv=%s,%s",
		       job_specs->argv[0],
		       job_specs->argv[1]);
	else if (job_specs->argc > 2)
		debug3("   argv=%s,%s,%s,...",
		       job_specs->argv[0],
		       job_specs->argv[1],
		       job_specs->argv[2]);

	if (job_specs->env_size == 1)
		debug3("   environment=\"%s\"", 
			job_specs->environment[0]);
	else if (job_specs->env_size == 2)
		debug3("   environment=%s,%s",
		       job_specs->environment[0],
		       job_specs->environment[1]);
	else if (job_specs->env_size > 2)
		debug3("   environment=%s,%s,%s,...",
		       job_specs->environment[0],
		       job_specs->environment[1],
		       job_specs->environment[2]);

	debug3("   in=%s out=%s err=%s",
	       job_specs->in, job_specs->out, job_specs->err);

	debug3("   work_dir=%s alloc_node:sid=%s:%u",
	       job_specs->work_dir,
	       job_specs->alloc_node, job_specs->alloc_sid);

	dependency = (job_specs->dependency != NO_VAL) ?
                        (long) job_specs->dependency : -1L;
	debug3("   alloc_resp_hostname=%s alloc_resp_port=%u",
	       job_specs->alloc_resp_hostname, job_specs->alloc_resp_port);
	debug3("   other_hostname=%s other_port=%u",
	       job_specs->other_hostname, job_specs->other_port);
	debug3("   dependency=%ld account=%s comment=%s",
	       dependency, job_specs->account, job_specs->comment);

	num_tasks = (job_specs->num_tasks != (uint16_t) NO_VAL) ?
			(long) job_specs->num_tasks : -1L;
	overcommit = (job_specs->overcommit != (uint16_t) NO_VAL) ?
			(long) job_specs->overcommit : -1L;
	debug3("   mail_type=%u mail_user=%s nice=%d num_tasks=%d overcommit=%d",
		job_specs->mail_type, job_specs->mail_user,
		(int)job_specs->nice - NICE_OFFSET, num_tasks, overcommit);

	slurm_make_time_str(&job_specs->begin_time, buf, sizeof(buf));
	cpus_per_task = (job_specs->cpus_per_task != (uint16_t) NO_VAL) ?
			(long) job_specs->cpus_per_task : -1L;
	no_requeue = (job_specs->no_requeue != (uint16_t) NO_VAL) ?
			(long) job_specs->no_requeue : -1L;
	debug3("   network=%s begin=%s cpus_per_task=%ld no_requeue=%ld", 
		job_specs->network, buf, cpus_per_task, no_requeue);

	ntasks_per_node = (job_specs->ntasks_per_node != (uint16_t) NO_VAL) ?
			(long) job_specs->ntasks_per_node : -1L;
	ntasks_per_socket = (job_specs->ntasks_per_socket != (uint16_t) NO_VAL) ?
			(long) job_specs->ntasks_per_socket : -1L;
	ntasks_per_core = (job_specs->ntasks_per_core != (uint16_t) NO_VAL) ?
			(long) job_specs->ntasks_per_core : -1L;
	debug3("   ntasks_per_node=%ld ntasks_per_socket=%ld ntasks_per_core=%ld", 
		ntasks_per_node, ntasks_per_socket, ntasks_per_core);

	select_g_sprint_jobinfo(job_specs->select_jobinfo, 
		buf, sizeof(buf), SELECT_PRINT_MIXED);
	if (buf[0] != '\0')
		debug3("   %s", buf);
}


/* 
 * init_job_conf - initialize the job configuration tables and values. 
 *	this should be called after creating node information, but 
 *	before creating any job entries. Pre-existing job entries are 
 *	left unchanged. 
 *	NOTE: The job hash table size does not change after initial creation.
 * RET 0 if no error, otherwise an error code
 * global: last_job_update - time of last job table update
 *	job_list - pointer to global job list
 */
int init_job_conf(void)
{
	if (job_list == NULL) {
		job_count = 0;
		job_list = list_create(_list_delete_job);
		if (job_list == NULL)
			fatal ("Memory allocation failure");
	}

	last_job_update = time(NULL);
	return SLURM_SUCCESS;
}

/*
 * rehash_jobs - Create or rebuild the job hash table.
 * NOTE: run lock_slurmctld before entry: Read config, write job
 */
extern void rehash_jobs(void)
{
	if (job_hash == NULL) {
		hash_table_size = slurmctld_conf.max_job_cnt;
		job_hash = (struct job_record **) xmalloc(hash_table_size *
					sizeof(struct job_record *));
	} else if (hash_table_size < (slurmctld_conf.max_job_cnt / 2)) {
		/* If the MaxJobCount grows by too much, the hash table will 
		 * be ineffective without rebuilding. We don't presently bother
		 * to rebuild the hash table, but cut MaxJobCount back as 
		 * needed. */ 
		error ("MaxJobCount reset too high, restart slurmctld");
		slurmctld_conf.max_job_cnt = hash_table_size;
	}
}

/*
 * job_allocate - create job_records for the suppied job specification and 
 *	allocate nodes for it.
 * IN job_specs - job specifications
 * IN immediate - if set then either initiate the job immediately or fail
 * IN will_run - don't initiate the job if set, just test if it could run 
 *	now or later
 * IN allocate - resource allocation request if set, not a full job
 * IN submit_uid -uid of user issuing the request
 * OUT job_pptr - set to pointer to job record
 * RET 0 or an error code. If the job would only be able to execute with 
 *	some change in partition configuration then 
 *	ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE is returned
 * NOTE: If allocating nodes lx[0-7] to a job and those nodes have cpu counts  
 *	of 4, 4, 4, 4, 8, 8, 4, 4 then num_cpu_groups=3, cpus_per_node={4,8,4}
 *	and cpu_count_reps={4,2,2}
 * globals: job_list - pointer to global job list 
 *	list_part - global list of partition info
 *	default_part_loc - pointer to default partition
 * NOTE: lock_slurmctld on entry: Read config Write job, Write node, Read part
 */
extern int job_allocate(job_desc_msg_t * job_specs, int immediate, int will_run, 
		int allocate, uid_t submit_uid, struct job_record **job_pptr)
{
	int error_code;
	bool no_alloc, top_prio, test_only, too_fragmented, independent;
	struct job_record *job_ptr;
	error_code = _job_create(job_specs, allocate, will_run,
				 &job_ptr, submit_uid);
	*job_pptr = job_ptr;
	
	if (error_code) {
		if (immediate && job_ptr) {
			job_ptr->job_state = JOB_FAILED;
			job_ptr->exit_code = 1;
			job_ptr->start_time = job_ptr->end_time = time(NULL);
			job_completion_logger(job_ptr);
		}
		return error_code;
	}
	xassert(job_ptr);

	independent = job_independent(job_ptr);

	/* Avoid resource fragmentation if important */
	if (independent && switch_no_frag() && 
	    (submit_uid || (job_specs->req_nodes == NULL)) && 
	    job_is_completing())
		too_fragmented = true;	/* Don't pick nodes for job now */
		/* FIXME: Ideally we only want to refuse the request if the 
		 * required node list is insufficient to satisfy the job's
		 * processor or node count requirements, but the overhead is
		 * rather high to do that right here. We let requests from
		 * user root proceed if a node list is specified, for 
		 * meta-schedulers (e.g. LCRM). */
	else
		too_fragmented = false;

	if (independent && (!too_fragmented))
		top_prio = _top_priority(job_ptr);
	else
		top_prio = true;	/* don't bother testing, 
					 * it is not runable anyway */
	if (immediate && (too_fragmented || (!top_prio) || (!independent))) {
		job_ptr->job_state  = JOB_FAILED;
		job_ptr->exit_code  = 1;
		job_ptr->start_time = job_ptr->end_time = time(NULL);
		job_completion_logger(job_ptr);
		if (!independent)
			return ESLURM_DEPENDENCY;
		else if (too_fragmented)
			return ESLURM_FRAGMENTATION;
		else
			return ESLURM_NOT_TOP_PRIORITY;
	}

	test_only = will_run || (allocate == 0);
	if (!test_only)
		last_job_update = time(NULL);

	no_alloc = test_only || too_fragmented || 
			(!top_prio) || (!independent);
	error_code = select_nodes(job_ptr, no_alloc, NULL);

	if ((error_code == ESLURM_NODES_BUSY)
	||  (error_code == ESLURM_JOB_HELD)
	||  (error_code == ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE)) {
		/* Not fatal error, but job can't be scheduled right now */
		if (immediate) {
			job_ptr->job_state  = JOB_FAILED;
			job_ptr->exit_code  = 1;
			job_ptr->start_time = job_ptr->end_time = time(NULL);
			job_completion_logger(job_ptr);
		} else		/* job remains queued */
			if (error_code == ESLURM_NODES_BUSY) {
				error_code = SLURM_SUCCESS;
			}
		return error_code;
	}

	if (error_code) {	/* fundamental flaw in job request */
		job_ptr->job_state  = JOB_FAILED;
		job_ptr->exit_code  = 1;
		job_ptr->start_time = job_ptr->end_time = time(NULL);
		job_completion_logger(job_ptr);
		return error_code;
	}

	if (will_run) {		/* job would run, flag job destruction */
		job_ptr->job_state  = JOB_FAILED;
		job_ptr->exit_code  = 1;
		job_ptr->start_time = job_ptr->end_time = time(NULL);
	} 
	return SLURM_SUCCESS;
}

/*
 * job_fail - terminate a job due to initiation failure
 * IN job_id - id of the job to be killed
 * RET 0 on success, otherwise ESLURM error code
 */
extern int job_fail(uint32_t job_id)
{
	struct job_record *job_ptr;
	time_t now = time(NULL);
	bool suspended = false;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		error("job_fail: invalid job id %u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}

	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;
	if (job_ptr->job_state == JOB_SUSPENDED)
		suspended = true;
	if ((job_ptr->job_state == JOB_RUNNING) || suspended) {
		/* No need to signal steps, deallocate kills them */
		job_ptr->time_last_active       = now;
		if (suspended)
			job_ptr->end_time       = job_ptr->suspend_time;
		else
			job_ptr->end_time       = now;
		last_job_update                 = now;
		job_ptr->job_state = JOB_FAILED | JOB_COMPLETING;
		job_ptr->exit_code = 1;
		deallocate_nodes(job_ptr, false, suspended);
		job_completion_logger(job_ptr);
		return SLURM_SUCCESS;
	}
	/* All other states */
	verbose("job_fail: job %u can't be killed from state=%s",
		job_id, job_state_string(job_ptr->job_state));
	return ESLURM_TRANSITION_STATE_NO_UPDATE;

}

/* 
 * job_signal - signal the specified job
 * IN job_id - id of the job to be signaled
 * IN signal - signal to send, SIGKILL == cancel the job
 * IN batch_flag - signal batch shell only if set
 * IN uid - uid of requesting user
 * RET 0 on success, otherwise ESLURM error code 
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 */
extern int job_signal(uint32_t job_id, uint16_t signal, uint16_t batch_flag, 
		uid_t uid)
{
	struct job_record *job_ptr;
	time_t now = time(NULL);
	bool super_user;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		info("job_signal: invalid job id %u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}

	super_user = ((uid == 0) || (uid == getuid()));
	if ((job_ptr->user_id != uid) && (!super_user)) {
		error("Security violation, JOB_CANCEL RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}
	if ((!super_user) && job_ptr->part_ptr
	&&  (job_ptr->part_ptr->root_only)) {
		info("Attempt to cancel job in RootOnly partition from uid %d",
			uid);
		return ESLURM_USER_ID_MISSING;
	}

	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;

	/* save user ID of the one who requested the job be cancelled */
	if(signal == SIGKILL)
		job_ptr->requid = uid;
	
	if ((job_ptr->job_state == JOB_PENDING) &&
	    (signal == SIGKILL)) {
		last_job_update		= now;
		job_ptr->job_state	= JOB_CANCELLED;
		job_ptr->start_time	= now;
		job_ptr->end_time	= now;
		job_completion_logger(job_ptr);
		delete_job_details(job_ptr);
		verbose("job_signal of pending job %u successful", job_id);
		return SLURM_SUCCESS;
	}

	if ((job_ptr->job_state == JOB_SUSPENDED)
	&&  (signal == SIGKILL)) {
		last_job_update         = now;
		job_ptr->end_time       = job_ptr->suspend_time;
		job_ptr->job_state      = JOB_CANCELLED | JOB_COMPLETING;
		deallocate_nodes(job_ptr, false, true);
		job_completion_logger(job_ptr);
		verbose("job_signal %u of suspended job %u successful",
			signal, job_id);
		return SLURM_SUCCESS;
	}
	
	if (job_ptr->job_state == JOB_RUNNING) {
		if (signal == SIGKILL) {
			/* No need to signal steps, deallocate kills them */
			job_ptr->time_last_active	= now;
			job_ptr->end_time		= now;
			last_job_update			= now;
			job_ptr->job_state = JOB_CANCELLED | JOB_COMPLETING;
			deallocate_nodes(job_ptr, false, false);
			job_completion_logger(job_ptr);
		} else if (batch_flag) {
			if (job_ptr->batch_flag)
				_signal_batch_job(job_ptr, signal);
			else
				return ESLURM_JOB_SCRIPT_MISSING;
		} else {
			_signal_job(job_ptr, signal);
		}
		verbose("job_signal %u of running job %u successful", 
			signal, job_id);
		return SLURM_SUCCESS;
	}

	verbose("job_signal: job %u can't be sent signal %u from state=%s",
		job_id, signal, job_state_string(job_ptr->job_state));
	return ESLURM_TRANSITION_STATE_NO_UPDATE;
}

static void
_signal_batch_job(struct job_record *job_ptr, uint16_t signal)
{
	bitoff_t i;
	kill_tasks_msg_t *kill_tasks_msg = NULL;
	agent_arg_t *agent_args = NULL;

	xassert(job_ptr);
	i = bit_ffs(job_ptr->node_bitmap);
	if (i < 0) {
		error("_signal_batch_job JobId=%u lacks assigned nodes");
		return;
	}

	agent_args = xmalloc(sizeof(agent_arg_t));
	agent_args->msg_type	= REQUEST_SIGNAL_TASKS;
	agent_args->retry	= 1;
	agent_args->node_count  = 1;
	agent_args->hostlist	= 
		hostlist_create(node_record_table_ptr[i].name);	
	kill_tasks_msg = xmalloc(sizeof(kill_tasks_msg_t));
	kill_tasks_msg->job_id      = job_ptr->job_id;
	kill_tasks_msg->job_step_id = NO_VAL;
	kill_tasks_msg->signal      = signal;

	agent_args->msg_args = kill_tasks_msg;
	agent_args->node_count = 1; /* slurm/477 be sure to update node_count */
	agent_queue_request(agent_args);
	return;
}

/* 
 * job_complete - note the normal termination the specified job
 * IN job_id - id of the job which completed
 * IN uid - user id of user issuing the RPC
 * IN requeue - job should be run again if possible
 * IN job_return_code - job's return code, if set then set state to FAILED
 * RET - 0 on success, otherwise ESLURM error code 
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 */
extern int job_complete(uint32_t job_id, uid_t uid, bool requeue,
	     uint32_t job_return_code)
{
	struct job_record *job_ptr;
	time_t now = time(NULL);
	uint32_t job_comp_flag = 0;
	bool suspended = false;
	info("completing job %u", job_id);
	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		info("job_complete: invalid JobId=%u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}

	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;

	if ((job_ptr->user_id != uid) && (uid != 0) && (uid != getuid())) {
		error("Security violation, JOB_COMPLETE RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}
	if (job_ptr->job_state & JOB_COMPLETING)
		return SLURM_SUCCESS;	/* avoid replay */

	if (job_ptr->job_state == JOB_RUNNING)
		job_comp_flag = JOB_COMPLETING;
	if (job_ptr->job_state == JOB_SUSPENDED) {
		job_comp_flag = JOB_COMPLETING;
		suspended = true;
	}

	if (requeue && (job_ptr->batch_flag > 1)) {
		/* Failed one requeue, just kill it */
		requeue = 0;
		if (job_return_code == 0)
			job_return_code = 1;
		info("Batch job launch failure, JobId=%u", job_ptr->job_id);
	}

	if (requeue && job_ptr->details && job_ptr->batch_flag) {
		job_ptr->batch_flag++;	/* only one retry */
		job_ptr->job_state = JOB_PENDING | job_comp_flag;
		info("Non-responding node, requeue JobId=%u", job_ptr->job_id);
	} else if (job_ptr->job_state == JOB_PENDING) {
		job_ptr->job_state  = JOB_CANCELLED;
		job_ptr->start_time = now;
		job_ptr->end_time   = now;
		job_completion_logger(job_ptr);
	} else {
		if (job_return_code == NO_VAL) {
			job_ptr->job_state = JOB_CANCELLED| job_comp_flag;
			if (job_ptr->requid == -1)
 				job_ptr->requid = uid;
		} else if (WEXITSTATUS(job_return_code)) {
			job_ptr->job_state = JOB_FAILED   | job_comp_flag;
			job_ptr->exit_code = job_return_code;
		}
		else if (job_comp_flag &&		/* job was running */
			 (job_ptr->end_time < now)) {	/* over time limit */
			job_ptr->job_state = JOB_TIMEOUT  | job_comp_flag;
			job_ptr->exit_code = MAX(job_ptr->exit_code, 1);
		} else
			job_ptr->job_state = JOB_COMPLETE | job_comp_flag;
		if (suspended)
			job_ptr->end_time = job_ptr->suspend_time;
		else
			job_ptr->end_time = now;
		job_completion_logger(job_ptr);
	}

	last_job_update = now;
	if (job_comp_flag) 	/* job was running */
		deallocate_nodes(job_ptr, false, suspended);
	info("job_complete for JobId=%u successful", job_id);

	return SLURM_SUCCESS;
}

/*
 * _job_create - create a job table record for the supplied specifications.
 *	this performs only basic tests for request validity (access to 
 *	partition, nodes count in partition, and sufficient processors in 
 *	partition).
 * input: job_specs - job specifications
 * IN allocate - resource allocation request if set rather than job submit
 * IN will_run - job is not to be created, test of validity only
 * OUT job_pptr - pointer to the job (NULL on error)
 * RET 0 on success, otherwise ESLURM error code. If the job would only be
 *	able to execute with some change in partition configuration then
 *	ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE is returned
 * globals: job_list - pointer to global job list 
 *	list_part - global list of partition info
 *	default_part_loc - pointer to default partition 
 *	job_hash - hash table into job records
 */

static int _job_create(job_desc_msg_t * job_desc, int allocate, int will_run,
		       struct job_record **job_pptr, uid_t submit_uid)
{
	int error_code = SLURM_SUCCESS, i;
	struct job_details *detail_ptr;
	enum job_wait_reason fail_reason;
	struct part_record *part_ptr;
	bitstr_t *req_bitmap = NULL, *exc_bitmap = NULL;
	bool super_user = false;
	struct job_record *job_ptr;
	uint32_t total_nodes;
#if SYSTEM_DIMENSIONS
	uint16_t geo[SYSTEM_DIMENSIONS];
#endif

	debug2("before alteration asking for nodes %d-%d procs %d", 
		     job_desc->min_nodes, job_desc->max_nodes,
		     job_desc->num_procs);
	
	select_g_alter_node_cnt(SELECT_SET_NODE_CNT, job_desc);
	select_g_get_jobinfo(job_desc->select_jobinfo,
			     SELECT_DATA_MAX_PROCS, &i);
	
	debug2("after alteration asking for nodes %d-%d procs %d-%d", 
		     job_desc->min_nodes, job_desc->max_nodes,
		     job_desc->num_procs, i);
	
	*job_pptr = (struct job_record *) NULL;
	if ((error_code = _validate_job_desc(job_desc, allocate, submit_uid)))
		return error_code;

	/* find selected partition */
	if (job_desc->partition) {
		part_ptr = list_find_first(part_list, &list_find_part,
					   job_desc->partition);
		if (part_ptr == NULL) {
			info("_job_create: invalid partition specified: %s", 
			     job_desc->partition);
			error_code = ESLURM_INVALID_PARTITION_NAME;
			return error_code;
		}
	} else {
		if (default_part_loc == NULL) {
			error("_job_create: default partition not set.");
			error_code = ESLURM_DEFAULT_PARTITION_NOT_SET;
			return error_code;
		}
		part_ptr = default_part_loc;
	}

	/* can this user access this partition */
	if ((part_ptr->root_only) && (submit_uid != 0)) {
		info("_job_create: uid %u access to partition %s denied, %s",
		     (unsigned int) submit_uid, part_ptr->name, "not root");
		error_code = ESLURM_ACCESS_DENIED;
		return error_code;
	}
	if (validate_group(part_ptr, job_desc->user_id) == 0) {
		info("_job_create: uid %u access to partition %s denied, %s",
			(unsigned int) job_desc->user_id, part_ptr->name, 
			"bad group");
		error_code = ESLURM_JOB_MISSING_REQUIRED_PARTITION_GROUP;
		return error_code;
	}

	/* check if select partition has sufficient resources to satisfy
	 * the request */

	/* insure that selected nodes are in this partition */
	if (job_desc->req_nodes) {
		error_code = node_name2bitmap(job_desc->req_nodes, false,  
					      &req_bitmap);
		if (error_code) {
			error_code = ESLURM_INVALID_NODE_NAME;
			goto cleanup;
		}
		if (job_desc->contiguous)
			bit_fill_gaps(req_bitmap);
		if (bit_super_set(req_bitmap, part_ptr->node_bitmap) != 1) {
			char *tmp = bitmap2node_name(req_bitmap);
			info("_job_create: requested nodes %s not in "
				"partition %s", tmp, part_ptr->name);
			xfree(tmp);
			error_code = ESLURM_REQUESTED_NODES_NOT_IN_PARTITION;
			goto cleanup;
		}
		i = bit_set_count(req_bitmap);
		if (i > job_desc->min_nodes)
			job_desc->min_nodes = i;
		if (i > job_desc->num_procs)
			job_desc->num_procs = i;
	}
	if (job_desc->exc_nodes) {
		error_code = node_name2bitmap(job_desc->exc_nodes, false,
					      &exc_bitmap);
		if (error_code) {
			error_code = ESLURM_INVALID_NODE_NAME;
			goto cleanup;
		}
	}
	if (exc_bitmap && req_bitmap) {
		bitstr_t *tmp_bitmap = NULL;
		bitoff_t first_set;
		tmp_bitmap = bit_copy(exc_bitmap);
		if (tmp_bitmap == NULL)
			fatal("bit_copy malloc failure");
		bit_and(tmp_bitmap, req_bitmap);
		first_set = bit_ffs(tmp_bitmap);
		FREE_NULL_BITMAP(tmp_bitmap);
		if (first_set != -1) {
			info("Job's required and excluded node lists overlap");
			error_code = ESLURM_INVALID_NODE_NAME;
			goto cleanup;
		}
	}

	if (job_desc->min_nodes == NO_VAL)
		job_desc->min_nodes = 1;

#if SYSTEM_DIMENSIONS
	select_g_get_jobinfo(job_desc->select_jobinfo,
			     SELECT_DATA_GEOMETRY,
			     &geo);
	if ((geo[0] != (uint16_t) NO_VAL) && (geo[0] != 0)) {
		int i, tot = 1;
		for (i=0; i<SYSTEM_DIMENSIONS; i++) {
			tot *= geo[i];
		}
		if (job_desc->min_nodes > tot) {
			info("MinNodes(%d) > GeometryNodes(%d)", 
				job_desc->min_nodes, tot);
			error_code = ESLURM_TOO_MANY_REQUESTED_CPUS;
			goto cleanup;
		}
		job_desc->min_nodes = tot;
	}
#endif

	if (job_desc->max_nodes == NO_VAL)
		job_desc->max_nodes = 0;
	if ((part_ptr->state_up)
	&&  (job_desc->num_procs > part_ptr->total_cpus)) {
		info("Job requested too many cpus (%d) of partition %s(%d)", 
		     job_desc->num_procs, part_ptr->name, 
		     part_ptr->total_cpus);
		error_code = ESLURM_TOO_MANY_REQUESTED_CPUS;
		goto cleanup;
	}
	total_nodes = part_ptr->total_nodes;
	select_g_alter_node_cnt(SELECT_APPLY_NODE_MIN_OFFSET,
			&total_nodes);
	if ((part_ptr->state_up) &&  (job_desc->min_nodes > total_nodes)) {
		info("Job requested too many nodes (%d) of partition %s(%d)", 
		     job_desc->min_nodes, part_ptr->name, 
		     part_ptr->total_nodes);
		error_code = ESLURM_TOO_MANY_REQUESTED_NODES;
		goto cleanup;
	}
	if (job_desc->max_nodes && 
	    (job_desc->max_nodes < job_desc->min_nodes)) {
		info("Job's max_nodes < min_nodes");
		error_code = ESLURM_TOO_MANY_REQUESTED_NODES;
		goto cleanup;
	}


	if ((error_code =_validate_job_create_req(job_desc)))
		goto cleanup;
	
	if ((error_code = _copy_job_desc_to_job_record(job_desc,
						       job_pptr,
						       part_ptr,
						       &req_bitmap,
						       &exc_bitmap))) {
		error_code = ESLURM_ERROR_ON_DESC_TO_RECORD_COPY;
		goto cleanup;
	}
	
	job_ptr = *job_pptr;
	if (job_ptr->dependency == job_ptr->job_id) {
		info("User specified self as dependent job");
		error_code = ESLURM_DEPENDENCY;
		goto cleanup;
	}

	if (job_desc->script
	&&  (!will_run)) {	/* don't bother with copy if just a test */
		if ((error_code = _copy_job_desc_to_file(job_desc,
							 job_ptr->job_id))) {
			job_ptr->job_state = JOB_FAILED;
			job_ptr->exit_code = 1;
			job_ptr->start_time = job_ptr->end_time = time(NULL);
			error_code = ESLURM_WRITING_TO_FILE;
			goto cleanup;
		}
		job_ptr->batch_flag = 1;
	} else
		job_ptr->batch_flag = 0;

	/* Insure that requested partition is valid right now, 
	 * otherwise leave job queued and provide warning code */
	detail_ptr = job_ptr->details;
	fail_reason= WAIT_NO_REASON;
	if ((job_desc->user_id == 0) ||
	    (job_desc->user_id == slurmctld_conf.slurm_user_id))
		super_user = true;
	if ((!super_user) && 
	    (job_desc->min_nodes > part_ptr->max_nodes)) {
		info("Job %u requested too many nodes (%d) of "
			"partition %s(%d)", 
			job_ptr->job_id, job_desc->min_nodes, 
			part_ptr->name, part_ptr->max_nodes);
		fail_reason = WAIT_PART_NODE_LIMIT;
	} else if ((!super_user) &&
	           (job_desc->max_nodes != 0) &&    /* no max_nodes for job */
		   (job_desc->max_nodes < part_ptr->min_nodes)) {
		info("Job %u requested too few nodes (%d) of partition %s(%d)",
			job_ptr->job_id, job_desc->max_nodes, 
			part_ptr->name, part_ptr->min_nodes);
		fail_reason = WAIT_PART_NODE_LIMIT;
	} else if (part_ptr->state_up == 0) {
		info("Job %u requested down partition %s", 
			job_ptr->job_id, part_ptr->name);
		fail_reason = WAIT_PART_STATE;
	}
	if (fail_reason != WAIT_NO_REASON) {
		error_code = ESLURM_REQUESTED_PART_CONFIG_UNAVAILABLE;
		job_ptr->priority = 1;      /* Move to end of queue */
		 if (detail_ptr)
			detail_ptr->wait_reason = fail_reason;
	}
	jobacct_g_job_start_slurmctld(job_ptr);
	
      cleanup:
	FREE_NULL_BITMAP(req_bitmap);
	FREE_NULL_BITMAP(exc_bitmap);
	return error_code;
}

/* Perform some size checks on strings we store to prevent
 * malicious user filling slurmctld's memory
 * RET 0 or error code */
static int _validate_job_create_req(job_desc_msg_t * job_desc)
{
	if (job_desc->err && (strlen(job_desc->err) > BUFFER_SIZE)) {
		info("_validate_job_create_req: strlen(err) too big (%d)",
		     strlen(job_desc->err));
		return ESLURM_PATHNAME_TOO_LONG;
	}
	if (job_desc->in && (strlen(job_desc->in) > BUFFER_SIZE)) {
		info("_validate_job_create_req: strlen(in) too big (%d)",
		     strlen(job_desc->in));
		return  ESLURM_PATHNAME_TOO_LONG;
	}
	if (job_desc->out && (strlen(job_desc->out) > BUFFER_SIZE)) {
		info("_validate_job_create_req: strlen(out) too big (%d)",
		     strlen(job_desc->out));
		return  ESLURM_PATHNAME_TOO_LONG;
	}
	if (job_desc->work_dir && (strlen(job_desc->work_dir) > BUFFER_SIZE)) {
		info("_validate_job_create_req: strlen(work_dir) too big (%d)",
		     strlen(job_desc->work_dir));
		return  ESLURM_PATHNAME_TOO_LONG;
	}
	return SLURM_SUCCESS;
}

/* _copy_job_desc_to_file - copy the job script and environment from the RPC  
 *	structure into a file */
static int
_copy_job_desc_to_file(job_desc_msg_t * job_desc, uint32_t job_id)
{
	int error_code = 0;
	char *dir_name, job_dir[20], *file_name;

	/* Create state_save_location directory */
	dir_name = xstrdup(slurmctld_conf.state_save_location);

	/* Create job_id specific directory */
	sprintf(job_dir, "/job.%d", job_id);
	xstrcat(dir_name, job_dir);
	if (mkdir(dir_name, 0700)) {
		error("mkdir(%s) error %m", dir_name);
		xfree(dir_name);
		return ESLURM_WRITING_TO_FILE;
	}

	/* Create environment file, and write data to it */
	file_name = xstrdup(dir_name);
	xstrcat(file_name, "/environment");
	error_code = _write_data_array_to_file(file_name,
					       job_desc->environment,
					       job_desc->env_size);
	xfree(file_name);

	if (error_code == 0) {
		/* Create script file */
		file_name = xstrdup(dir_name);
		xstrcat(file_name, "/script");
		error_code =
		    _write_data_to_file(file_name, job_desc->script);
		xfree(file_name);
	}

	xfree(dir_name);
	return error_code;
}

/*
 * Create file with specified name and write the supplied data array to it
 * IN file_name - file to create and write to
 * IN data - array of pointers to strings (e.g. env)
 * IN size - number of elements in data
 */
static int
_write_data_array_to_file(char *file_name, char **data, uint16_t size)
{
	int fd, i, pos, nwrite, amount;

	fd = creat(file_name, 0600);
	if (fd < 0) {
		error("Error creating file %s, %m", file_name);
		return ESLURM_WRITING_TO_FILE;
	}

	amount = write(fd, &size, sizeof(uint16_t));
	if (amount < sizeof(uint16_t)) {
		error("Error writing file %s, %m", file_name);
		close(fd);
		return ESLURM_WRITING_TO_FILE;
	}

	if (data == NULL)
		return SLURM_SUCCESS;

	for (i = 0; i < size; i++) {
		nwrite = strlen(data[i]) + 1;
		pos = 0;
		while (nwrite > 0) {
			amount = write(fd, &data[i][pos], nwrite);
			if ((amount < 0) && (errno != EINTR)) {
				error("Error writing file %s, %m",
				      file_name);
				close(fd);
				return ESLURM_WRITING_TO_FILE;
			}
			nwrite -= amount;
			pos    += amount;
		}
	}

	close(fd);
	return SLURM_SUCCESS;
}

/*
 * Create file with specified name and write the supplied data array to it
 * IN file_name - file to create and write to
 * IN data - pointer to string
 */
static int _write_data_to_file(char *file_name, char *data)
{
	int fd, pos, nwrite, amount;

	if (data == NULL) {
		(void) unlink(file_name);
		return SLURM_SUCCESS;
	}

	fd = creat(file_name, 0700);
	if (fd < 0) {
		error("Error creating file %s, %m", file_name);
		return ESLURM_WRITING_TO_FILE;
	}

	nwrite = strlen(data) + 1;
	pos = 0;
	while (nwrite > 0) {
		amount = write(fd, &data[pos], nwrite);
		if ((amount < 0) && (errno != EINTR)) {
			error("Error writing file %s, %m", file_name);
			close(fd);
			return ESLURM_WRITING_TO_FILE;
		}
		nwrite -= amount;
		pos    += amount;
	}
	close(fd);
	return SLURM_SUCCESS;
}

/*
 * get_job_env - return the environment variables and their count for a 
 *	given job
 * IN job_ptr - pointer to job for which data is required
 * OUT env_size - number of elements to read
 * RET point to array of string pointers containing environment variables
 * NOTE: READ lock_slurmctld config before entry
 */
char **get_job_env(struct job_record *job_ptr, uint16_t * env_size)
{
	char job_dir[30], *file_name, **environment = NULL;

	file_name = xstrdup(slurmctld_conf.state_save_location);
	sprintf(job_dir, "/job.%d/environment", job_ptr->job_id);
	xstrcat(file_name, job_dir);

	_read_data_array_from_file(file_name, &environment, env_size);

	xfree(file_name);
	return environment;
}

/* 
 * get_job_script - return the script for a given job
 * IN job_ptr - pointer to job for which data is required
 * RET point to string containing job script
 * NOTE: READ lock_slurmctld config before entry
 */
char *get_job_script(struct job_record *job_ptr)
{
	char job_dir[30], *file_name, *script = NULL;

	file_name = xstrdup(slurmctld_conf.state_save_location);
	sprintf(job_dir, "/job.%d/script", job_ptr->job_id);
	xstrcat(file_name, job_dir);

	_read_data_from_file(file_name, &script);

	xfree(file_name);
	return script;
}

/*
 * Read a collection of strings from a file
 * IN file_name - file to read from
 * OUT data - pointer to array of pointers to strings (e.g. env),
 *	must be xfreed when no longer needed
 * OUT size - number of elements in data
 */
static void
_read_data_array_from_file(char *file_name, char ***data, uint16_t * size)
{
	int fd, pos, buf_size, amount, i;
	char *buffer, **array_ptr;
	uint16_t rec_cnt;

	xassert(file_name);
	xassert(data);
	xassert(size);
	*data = NULL;
	*size = 0;

	fd = open(file_name, 0);
	if (fd < 0) {
		error("Error opening file %s, %m", file_name);
		return;
	}

	amount = read(fd, &rec_cnt, sizeof(uint16_t));
	if (amount < sizeof(uint16_t)) {
		if (amount != 0)	/* incomplete write */
			error("Error reading file %s, %m", file_name);
		else 
			verbose("File %s has zero size", file_name); 
		close(fd);
		return;
	}

	if (rec_cnt == 0) {
		*data = NULL;
		*size = 0;
		return;
	}

	pos = 0;
	buf_size = HUGE_BUF_SIZE;
	buffer = xmalloc(buf_size);
	while (1) {
		amount = read(fd, &buffer[pos], buf_size);
		if (amount < 0) {
			error("Error reading file %s, %m", file_name);
			xfree(buffer);
			close(fd);
			return;
		}
		if (amount < buf_size)	/* end of file */
			break;
		pos += amount;
		xrealloc(buffer, (pos + buf_size));
	}
	close(fd);

	/* We have all the data, now let's compute the pointers */
	pos = 0;
	array_ptr = xmalloc(rec_cnt * sizeof(char *));
	for (i = 0; i < rec_cnt; i++) {
		array_ptr[i] = &buffer[pos];
		pos += strlen(&buffer[pos]) + 1;
		if ((pos > buf_size) && ((i + 1) < rec_cnt)) {
			error("Bad environment file %s", file_name);
			break;
		}
	}

	*size = rec_cnt;
	*data = array_ptr;
	return;
}

/*
 * Read a string from a file
 * IN file_name - file to read from
 * OUT data - pointer to  string 
 *	must be xfreed when no longer needed
 */
void _read_data_from_file(char *file_name, char **data)
{
	int fd, pos, buf_size, amount;
	char *buffer;

	xassert(file_name);
	xassert(data);
	*data = NULL;

	fd = open(file_name, 0);
	if (fd < 0) {
		error("Error opening file %s, %m", file_name);
		return;
	}

	pos = 0;
	buf_size = HUGE_BUF_SIZE;
	buffer = xmalloc(buf_size);
	while (1) {
		amount = read(fd, &buffer[pos], buf_size);
		if (amount < 0) {
			error("Error reading file %s, %m", file_name);
			xfree(buffer);
			close(fd);
			return;
		}
		if (amount < buf_size)	/* end of file */
			break;
		pos += amount;
		xrealloc(buffer, (pos + buf_size));
	}

	*data = buffer;
	close(fd);
	return;
}

/* _copy_job_desc_to_job_record - copy the job descriptor from the RPC  
 *	structure into the actual slurmctld job record */
static int
_copy_job_desc_to_job_record(job_desc_msg_t * job_desc,
			     struct job_record **job_rec_ptr,
			     struct part_record *part_ptr,
			     bitstr_t ** req_bitmap,
			     bitstr_t ** exc_bitmap)
{
	int error_code;
	struct job_details *detail_ptr;
	struct job_record *job_ptr;

	job_ptr = create_job_record(&error_code);
	if (error_code)
		return error_code;

	strncpy(job_ptr->partition, part_ptr->name, MAX_SLURM_NAME);
	job_ptr->part_ptr = part_ptr;
	if (job_desc->job_id != NO_VAL)		/* already confirmed unique */
		job_ptr->job_id = job_desc->job_id;
	else
		_set_job_id(job_ptr);
	_add_job_hash(job_ptr);

	if (job_desc->name) {
		strncpy(job_ptr->name, job_desc->name, MAX_JOBNAME_LEN);
	}
	job_ptr->user_id    = (uid_t) job_desc->user_id;
	job_ptr->group_id   = (gid_t) job_desc->group_id;
	job_ptr->job_state  = JOB_PENDING;
	job_ptr->time_limit = job_desc->time_limit;
	job_ptr->alloc_sid  = job_desc->alloc_sid;
	job_ptr->alloc_node = xstrdup(job_desc->alloc_node);
	job_ptr->account    = xstrdup(job_desc->account);
	job_ptr->network    = xstrdup(job_desc->network);
	job_ptr->comment    = xstrdup(job_desc->comment);
	if (job_desc->dependency != NO_VAL) /* leave as zero */
		job_ptr->dependency = job_desc->dependency;

	if (job_desc->priority != NO_VAL) /* already confirmed submit_uid==0 */
		job_ptr->priority = job_desc->priority;
	else {
		_set_job_prio(job_ptr);
		job_ptr->priority -= ((int)job_desc->nice - NICE_OFFSET);
	}

	if (job_desc->kill_on_node_fail != (uint16_t) NO_VAL)
		job_ptr->kill_on_node_fail = job_desc->kill_on_node_fail;

	job_ptr->alloc_resp_port = job_desc->alloc_resp_port;
	job_ptr->alloc_resp_host = xstrdup(job_desc->alloc_resp_hostname);
	job_ptr->other_port = job_desc->other_port;
	job_ptr->other_host = xstrdup(job_desc->other_hostname);
	job_ptr->time_last_active = time(NULL);
	job_ptr->num_procs = job_desc->num_procs;
        job_ptr->cr_enabled = 0;

	job_ptr->mail_type = job_desc->mail_type;
	job_ptr->mail_user = xstrdup(job_desc->mail_user);

	detail_ptr = job_ptr->details;
	detail_ptr->argc = job_desc->argc;
	detail_ptr->argv = job_desc->argv;
	job_desc->argv   = (char **) NULL; /* nothing left */
	job_desc->argc   = 0;		   /* nothing left */
	detail_ptr->min_nodes = job_desc->min_nodes;
	detail_ptr->max_nodes = job_desc->max_nodes;
	detail_ptr->min_sockets = job_desc->min_sockets;
	detail_ptr->max_sockets = job_desc->max_sockets;
	detail_ptr->min_cores = job_desc->min_cores;
	detail_ptr->max_cores = job_desc->max_cores;
	detail_ptr->min_threads = job_desc->min_threads;
	detail_ptr->max_threads = job_desc->max_threads;
	if (job_desc->req_nodes) {
		detail_ptr->req_nodes = 
				_copy_nodelist_no_dup(job_desc->req_nodes);
		detail_ptr->req_node_bitmap = *req_bitmap;
		*req_bitmap = NULL;	/* Reused nothing left to free */
	}
	if (job_desc->exc_nodes) {
		detail_ptr->exc_nodes = 
				_copy_nodelist_no_dup(job_desc->exc_nodes);
		detail_ptr->exc_node_bitmap = *exc_bitmap;
		*exc_bitmap = NULL;	/* Reused nothing left to free */
	}
	if (job_desc->features)
		detail_ptr->features = xstrdup(job_desc->features);
	detail_ptr->shared = job_desc->shared;
	if (job_desc->contiguous != (uint16_t) NO_VAL)
		detail_ptr->contiguous = job_desc->contiguous;
        if (job_desc->task_dist != (uint32_t) NO_VAL)
                detail_ptr->task_dist = job_desc->task_dist;
        if (job_desc->plane_size != (uint32_t) NO_VAL)
                detail_ptr->plane_size = job_desc->plane_size;
	if (job_desc->cpus_per_task != (uint16_t) NO_VAL)
		detail_ptr->cpus_per_task = job_desc->cpus_per_task;
	if (job_desc->ntasks_per_node != (uint16_t) NO_VAL)
		detail_ptr->ntasks_per_node = job_desc->ntasks_per_node;
	if (job_desc->ntasks_per_socket != (uint16_t) NO_VAL)
		detail_ptr->ntasks_per_socket = job_desc->ntasks_per_socket;
	if (job_desc->ntasks_per_core != (uint16_t) NO_VAL)
		detail_ptr->ntasks_per_core = job_desc->ntasks_per_core;
	if (job_desc->no_requeue != (uint16_t) NO_VAL)
		detail_ptr->no_requeue = job_desc->no_requeue;
	if (job_desc->job_min_procs != NO_VAL)
		detail_ptr->job_min_procs = job_desc->job_min_procs;
	detail_ptr->job_min_procs = MAX(detail_ptr->job_min_procs,
			detail_ptr->cpus_per_task);
	if (job_desc->job_min_sockets != NO_VAL)
		detail_ptr->job_min_sockets = job_desc->job_min_sockets;
	if (job_desc->job_min_cores != NO_VAL)
		detail_ptr->job_min_cores = job_desc->job_min_cores;
	if (job_desc->job_min_threads != NO_VAL)
		detail_ptr->job_min_threads = job_desc->job_min_threads;
	if (job_desc->job_min_memory != NO_VAL)
		detail_ptr->job_min_memory = job_desc->job_min_memory;
	if (job_desc->job_max_memory != NO_VAL)
		detail_ptr->job_max_memory = job_desc->job_max_memory;
	if (job_desc->job_min_tmp_disk != NO_VAL)
		detail_ptr->job_min_tmp_disk = job_desc->job_min_tmp_disk;
	if (job_desc->num_tasks != NO_VAL)
		detail_ptr->num_tasks = job_desc->num_tasks;
	if (job_desc->err)
		detail_ptr->err = xstrdup(job_desc->err);
	if (job_desc->in)
		detail_ptr->in = xstrdup(job_desc->in);
	if (job_desc->out)
		detail_ptr->out = xstrdup(job_desc->out);
	if (job_desc->work_dir)
		detail_ptr->work_dir = xstrdup(job_desc->work_dir);
	if (job_desc->overcommit != (uint16_t) NO_VAL)
		detail_ptr->overcommit = job_desc->overcommit;
	detail_ptr->begin_time = job_desc->begin_time;
	job_ptr->select_jobinfo = 
		select_g_copy_jobinfo(job_desc->select_jobinfo);	

	*job_rec_ptr = job_ptr;
	return SLURM_SUCCESS;
}

/*
 * _copy_nodelist_no_dup - Take a node_list string and convert it to an 
 *	expression without duplicate names. For example, we want to convert 
 *	a users request for nodes "lx1,lx2,lx1,lx3" to "lx[1-3]"
 * node_list IN - string describing a list of nodes
 * RET a compact node expression, must be xfreed by the user
 */
static char *_copy_nodelist_no_dup(char *node_list)
{
	char buf[8192];

	hostlist_t hl = hostlist_create(node_list);
	if (hl == NULL)
		return NULL;
	hostlist_uniq(hl);
	hostlist_ranged_string(hl, 8192, buf);
	hostlist_destroy(hl);

	return xstrdup(buf);
}

/* 
 * job_time_limit - terminate jobs which have exceeded their time limit
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 * NOTE: READ lock_slurmctld config before entry
 */
void job_time_limit(void)
{
	ListIterator job_iterator;
	struct job_record *job_ptr;
	time_t now = time(NULL);
	time_t old = now - slurmctld_conf.inactive_limit;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr =
		(struct job_record *) list_next(job_iterator))) {
		xassert (job_ptr->magic == JOB_MAGIC);
		if (job_ptr->job_state != JOB_RUNNING)
			continue;

		/* Consider a job active if it has any active steps */
		if (job_ptr->step_list
		&&  (list_count(job_ptr->step_list) > 0))
			job_ptr->time_last_active = now;

		if (slurmctld_conf.inactive_limit
		&&  (job_ptr->time_last_active <= old)
		&&  (job_ptr->part_ptr)
		&&  (job_ptr->part_ptr->root_only == 0)) {
			/* job inactive, kill it */
			info("Inactivity time limit reached for JobId=%u",
				job_ptr->job_id);
			_job_timed_out(job_ptr);
			continue;
		}
		if ((job_ptr->time_limit != INFINITE)
		&&  (job_ptr->end_time <= now)) {
			last_job_update = now;
			info("Time limit exhausted for JobId=%u",
				job_ptr->job_id);
			_job_timed_out(job_ptr);
			continue;
		}

		/* Give srun command warning message about pending timeout */
		if (job_ptr->end_time <= (now + PERIODIC_TIMEOUT * 2))
			srun_timeout (job_ptr);
	}

	list_iterator_destroy(job_iterator);
}

/* Terminate a job that has exhausted its time limit */
static void _job_timed_out(struct job_record *job_ptr)
{
	xassert(job_ptr);

	if (job_ptr->details) {
		time_t now      = time(NULL);
		job_ptr->end_time           = now;
		job_ptr->time_last_active   = now;
		job_ptr->job_state          = JOB_TIMEOUT | JOB_COMPLETING;
		job_ptr->exit_code = MAX(job_ptr->exit_code, 1);
		deallocate_nodes(job_ptr, true, false);
		job_completion_logger(job_ptr);
	} else
		job_signal(job_ptr->job_id, SIGKILL, 0, 0);
	return;
}

/* _validate_job_desc - validate that a job descriptor for job submit or 
 *	allocate has valid data, set values to defaults as required 
 * IN/OUT job_desc_msg - pointer to job descriptor, modified as needed
 * IN allocate - if clear job to be queued, if set allocate for user now 
 * IN submit_uid - who request originated
 */
static int _validate_job_desc(job_desc_msg_t * job_desc_msg, int allocate, 
			      uid_t submit_uid)
{	
	if ((job_desc_msg->num_procs == NO_VAL) &&
	    (job_desc_msg->min_nodes == NO_VAL) &&
	    (job_desc_msg->req_nodes == NULL)) {
		info("Job specified no num_procs, min_nodes or req_nodes");
		return ESLURM_JOB_MISSING_SIZE_SPECIFICATION;
	}
	if ((allocate == SLURM_CREATE_JOB_FLAG_NO_ALLOCATE_0) &&
	    (job_desc_msg->script == NULL)) {
		info("_validate_job_desc: job failed to specify Script");
		return ESLURM_JOB_SCRIPT_MISSING;
	}
	if (job_desc_msg->user_id == NO_VAL) {
		info("_validate_job_desc: job failed to specify User");
		return ESLURM_USER_ID_MISSING;
	}
	if ( job_desc_msg->group_id == NO_VAL ) {
		debug("_validate_job_desc: job failed to specify group");
		job_desc_msg->group_id = 0;	/* uses user default */
	}
	if ((job_desc_msg->name) &&
	    (strlen(job_desc_msg->name) >= MAX_JOBNAME_LEN)) {
		job_desc_msg->name[MAX_JOBNAME_LEN-1] = '\0';
	}
	if (job_desc_msg->contiguous == (uint16_t) NO_VAL)
		job_desc_msg->contiguous = 0;

        if (job_desc_msg->task_dist == (uint32_t) NO_VAL)
		info("_validate_job_desc: job failed to specify distribution ");
        if (job_desc_msg->plane_size == (uint32_t) NO_VAL)
                job_desc_msg->plane_size = 0;

	if (job_desc_msg->kill_on_node_fail == (uint16_t) NO_VAL)
		job_desc_msg->kill_on_node_fail = 1;

	if (job_desc_msg->job_id != NO_VAL) {
		struct job_record *dup_job_ptr;
		if ((submit_uid != 0) && 
		    (submit_uid != slurmctld_conf.slurm_user_id)) {
			info("attempt by uid %u to set job_id", submit_uid);
			return ESLURM_INVALID_JOB_ID;
		}
		if (job_desc_msg->job_id == 0) {
			info("attempt by uid %u to set zero job_id", submit_uid);
			return ESLURM_INVALID_JOB_ID;
		}
		dup_job_ptr = find_job_record((uint32_t) job_desc_msg->job_id);
		if (dup_job_ptr && 
		    (!(IS_JOB_FINISHED(dup_job_ptr)))) {
			info("attempt re-use active job_id %u", 
			     job_desc_msg->job_id);
			return ESLURM_DUPLICATE_JOB_ID;
		}
		if (dup_job_ptr)	/* Purge the record for re-use */
			_purge_job_record(job_desc_msg->job_id);
	}

	if ((submit_uid != 0) 	/* only root or SlurmUser can set job prio */
	&&  (submit_uid != slurmctld_conf.slurm_user_id)) {
		if (job_desc_msg->priority != 0)
			job_desc_msg->priority = NO_VAL;
		if (job_desc_msg->nice < NICE_OFFSET)
			job_desc_msg->nice = NICE_OFFSET;
	}

	if (job_desc_msg->num_procs == NO_VAL)
		job_desc_msg->num_procs = 1;	/* default cpu count of 1 */
	if (job_desc_msg->min_sockets == NO_VAL)
		job_desc_msg->min_sockets = 1;	/* default socket count of 1 */
	if (job_desc_msg->min_cores == NO_VAL)
		job_desc_msg->min_cores = 1;	/* default core count of 1 */
	if (job_desc_msg->min_threads == NO_VAL)
		job_desc_msg->min_threads = 1;	/* default thread count of 1 */
	if (job_desc_msg->min_nodes == NO_VAL)
		job_desc_msg->min_nodes = 1;	/* default node count of 1 */
	if (job_desc_msg->job_min_procs == NO_VAL)
		job_desc_msg->job_min_procs = 1;   /* default 1 cpu per node */
	if (job_desc_msg->job_min_sockets == NO_VAL)
		job_desc_msg->job_min_sockets = 1; /* default 1 socket per node */
	if (job_desc_msg->job_min_cores == NO_VAL)
		job_desc_msg->job_min_cores = 1;   /* default 1 core per socket */
	if (job_desc_msg->job_min_threads == NO_VAL)
		job_desc_msg->job_min_threads = 1; /* default 1 thread per core */
	if (job_desc_msg->job_min_memory == NO_VAL)
		job_desc_msg->job_min_memory = 1;  /* default 1MB mem per node */
	if (job_desc_msg->job_max_memory == NO_VAL)
		job_desc_msg->job_max_memory = 1;  /* default 1MB mem per node */
	if (job_desc_msg->job_min_tmp_disk == NO_VAL)
		job_desc_msg->job_min_tmp_disk = 1;/* default 1MB disk per node */

	return SLURM_SUCCESS;
}

/* 
 * _list_delete_job - delete a job record and its corresponding job_details,
 *	see common/list.h for documentation
 * IN job_entry - pointer to job_record to delete
 * global: job_list - pointer to global job list
 *	job_count - count of job list entries
 *	job_hash - hash table into job records
 */
static void _list_delete_job(void *job_entry)
{
	struct job_record *job_ptr = (struct job_record *) job_entry;
	struct job_record **job_pptr;

	xassert(job_entry);
	xassert (job_ptr->magic == JOB_MAGIC);

	/* Remove the record from the hash table */
	job_pptr = &job_hash[JOB_HASH_INX(job_ptr->job_id)];
	while ((job_pptr != NULL) && 
	       ((job_ptr = *job_pptr) != (struct job_record *) job_entry)) {
		job_pptr = &job_ptr->job_next;		
	}
	if (job_pptr == NULL)
		fatal("job hash error"); 
	*job_pptr = job_ptr->job_next;

	delete_job_details(job_ptr);
	xfree(job_ptr->alloc_node);
	xfree(job_ptr->nodes);
	FREE_NULL_BITMAP(job_ptr->node_bitmap);
	xfree(job_ptr->cpus_per_node);
	xfree(job_ptr->cpu_count_reps);
	xfree(job_ptr->node_addr);
	xfree(job_ptr->alloc_resp_host);
	xfree(job_ptr->other_host);
	xfree(job_ptr->account);
	xfree(job_ptr->mail_user);
	xfree(job_ptr->network);
	xfree(job_ptr->alloc_lps);
	xfree(job_ptr->comment);
	select_g_free_jobinfo(&job_ptr->select_jobinfo);
	if (job_ptr->step_list) {
		delete_all_step_records(job_ptr);
		list_destroy(job_ptr->step_list);
	}
	job_count--;
	xfree(job_ptr);
}


/*
 * _list_find_job_id - find specific job_id entry in the job list,  
 *	see common/list.h for documentation, key is job_id_ptr 
 * global- job_list - the global partition list
 */
static int _list_find_job_id(void *job_entry, void *key)
{
	uint32_t *job_id_ptr = (uint32_t *) key;

	if (((struct job_record *) job_entry)->job_id == *job_id_ptr)
		return 1;
	else
		return 0;
}


/*
 * _list_find_job_old - find old entries in the job list,  
 *	see common/list.h for documentation, key is ignored 
 * global- job_list - the global partition list
 */
static int _list_find_job_old(void *job_entry, void *key)
{
	time_t now      = time(NULL);
	time_t kill_age = now - (slurmctld_conf.kill_wait + 20);
	time_t min_age  = now - slurmctld_conf.min_job_age;
	struct job_record *job_ptr = (struct job_record *)job_entry;

	if ( (job_ptr->job_state & JOB_COMPLETING) &&
	     (job_ptr->end_time < kill_age) ) {
		re_kill_job(job_ptr);
		return 0;       /* Job still completing */
	}

	if (slurmctld_conf.min_job_age == 0)
		return 0;	/* No job record purging */

	if (job_ptr->end_time > min_age)
		return 0;	/* Too new to purge */

	if (!(IS_JOB_FINISHED(job_ptr))) 
		return 0;	/* Job still active */

	return 1;		/* Purge the job */
}


/* 
 * pack_all_jobs - dump all job information for all jobs in 
 *	machine independent form (for network transmission)
 * OUT buffer_ptr - the pointer is set to the allocated buffer.
 * OUT buffer_size - set to size of the buffer in bytes
 * IN show_flags - job filtering options
 * IN uid - uid of user making request (for partition filtering)
 * global: job_list - global list of job records
 * NOTE: the buffer at *buffer_ptr must be xfreed by the caller
 * NOTE: change _unpack_job_desc_msg() in common/slurm_protocol_pack.c 
 *	whenever the data format changes
 */
extern void pack_all_jobs(char **buffer_ptr, int *buffer_size,
		uint16_t show_flags, uid_t uid)
{
	ListIterator job_iterator;
	struct job_record *job_ptr;
	uint32_t jobs_packed = 0, tmp_offset;
	Buf buffer;
	time_t now = time(NULL);

	buffer_ptr[0] = NULL;
	*buffer_size = 0;

	buffer = init_buf(HUGE_BUF_SIZE);

	/* write message body header : size and time */
	/* put in a place holder job record count of 0 for now */
	pack32((uint32_t) jobs_packed, buffer);
	pack_time(now, buffer);

	/* write individual job records */
	part_filter_set(uid);
	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		xassert (job_ptr->magic == JOB_MAGIC);

		if (((show_flags & SHOW_ALL) == 0) &&
		    (job_ptr->part_ptr) && 
		    (job_ptr->part_ptr->hidden))
			continue;

		pack_job(job_ptr, buffer);
		jobs_packed++;
	}
	part_filter_clear();
	list_iterator_destroy(job_iterator);

	/* put the real record count in the message body header */
	tmp_offset = get_buf_offset(buffer);
	set_buf_offset(buffer, 0);
	pack32((uint32_t) jobs_packed, buffer);
	set_buf_offset(buffer, tmp_offset);

	*buffer_size = get_buf_offset(buffer);
	buffer_ptr[0] = xfer_buf_data(buffer);
}


/* 
 * pack_job - dump all configuration information about a specific job in 
 *	machine independent form (for network transmission)
 * IN dump_job_ptr - pointer to job for which information is requested
 * IN/OUT buffer - buffer in which data is placed, pointers automatically 
 *	updated
 * NOTE: change _unpack_job_info_members() in common/slurm_protocol_pack.c
 *	  whenever the data format changes
 */
void pack_job(struct job_record *dump_job_ptr, Buf buffer)
{
	struct job_details *detail_ptr;
	uint32_t size_tmp;

	pack32((uint32_t)dump_job_ptr->job_id, buffer);
	pack32((uint32_t)dump_job_ptr->user_id, buffer);
	pack32((uint32_t)dump_job_ptr->group_id, buffer);

	pack16((uint16_t) dump_job_ptr->job_state, buffer);
	pack16((uint16_t) dump_job_ptr->batch_flag, buffer);
	pack32((uint32_t)dump_job_ptr->alloc_sid, buffer);
	if ((dump_job_ptr->time_limit == NO_VAL) && dump_job_ptr->part_ptr)
		pack32((uint32_t)dump_job_ptr->part_ptr->max_time, buffer);
	else
		pack32((uint32_t)dump_job_ptr->time_limit, buffer);

	if (dump_job_ptr->details) {
		pack_time(dump_job_ptr->details->submit_time, buffer);
	} else {
		pack_time((time_t) 0, buffer);
	}
	if (IS_JOB_PENDING(dump_job_ptr)) {
		if (dump_job_ptr->details)
			pack_time(dump_job_ptr->details->begin_time,
				buffer);
		else
			pack_time((time_t) 0, buffer);
	} else
		pack_time(dump_job_ptr->start_time, buffer);
	pack_time(dump_job_ptr->end_time, buffer);
	pack_time(dump_job_ptr->suspend_time, buffer);
	pack_time(dump_job_ptr->pre_sus_time, buffer);

	pack32((uint32_t)dump_job_ptr->priority, buffer);

	packstr(dump_job_ptr->nodes, buffer);
	packstr(dump_job_ptr->partition, buffer);
	packstr(dump_job_ptr->account, buffer);
	packstr(dump_job_ptr->network, buffer);
	packstr(dump_job_ptr->comment, buffer);
	pack32((uint32_t)dump_job_ptr->dependency, buffer);
	pack32((uint32_t)dump_job_ptr->exit_code, buffer);

        pack16(dump_job_ptr->num_cpu_groups, buffer);
	size_tmp = dump_job_ptr->num_cpu_groups;
	if (size_tmp < 0) {
	    	size_tmp = 0;
	}
	pack32_array(dump_job_ptr->cpus_per_node, size_tmp, buffer);
	pack32_array(dump_job_ptr->cpu_count_reps, size_tmp, buffer);

	packstr(dump_job_ptr->name, buffer);
	packstr(dump_job_ptr->alloc_node, buffer);
	pack_bit_fmt(dump_job_ptr->node_bitmap, buffer);
	pack32((uint32_t)dump_job_ptr->num_procs, buffer);
	
	select_g_pack_jobinfo(dump_job_ptr->select_jobinfo, buffer);

	detail_ptr = dump_job_ptr->details;
	/* A few details are always dumped here */
	_pack_default_job_details(detail_ptr, buffer);

	/* other job details are only dumped until the job starts 
	 * running (at which time they become meaningless) */
	if (detail_ptr && dump_job_ptr->job_state == JOB_PENDING)
		_pack_pending_job_details(detail_ptr, buffer);
	else
		_pack_pending_job_details(NULL, buffer);
}

/* pack default job details for "get_job_info" RPC */
static void _pack_default_job_details(struct job_details *detail_ptr, Buf buffer)
{
	if (detail_ptr) {
		packstr(detail_ptr->features, buffer);

		pack32((uint32_t) detail_ptr->min_nodes, buffer);
		pack32((uint32_t) detail_ptr->max_nodes, buffer);
		pack32((uint32_t) detail_ptr->min_sockets, buffer);
		pack32((uint32_t) detail_ptr->max_sockets, buffer);
		pack32((uint32_t) detail_ptr->min_cores, buffer);
		pack32((uint32_t) detail_ptr->max_cores, buffer);
		pack32((uint32_t) detail_ptr->min_threads, buffer);
		pack32((uint32_t) detail_ptr->max_threads, buffer);
	} else {
		packnull(buffer);

		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
	}
}

/* pack pending job details for "get_job_info" RPC */
static void _pack_pending_job_details(struct job_details *detail_ptr, Buf buffer)
{
	if (detail_ptr) {		
		pack16((uint16_t) detail_ptr->shared, buffer);
		pack16((uint16_t) detail_ptr->contiguous, buffer);
		pack16((uint16_t) detail_ptr->cpus_per_task, buffer);
		pack16((uint16_t) detail_ptr->ntasks_per_node, buffer);
		pack16((uint16_t) detail_ptr->ntasks_per_socket, buffer);
		pack16((uint16_t) detail_ptr->ntasks_per_core, buffer);

		pack32((uint32_t) detail_ptr->job_min_procs, buffer);
		pack32((uint32_t) detail_ptr->job_min_sockets, buffer);
		pack32((uint32_t) detail_ptr->job_min_cores, buffer);
		pack32((uint32_t) detail_ptr->job_min_threads, buffer);
		pack32((uint32_t) detail_ptr->job_min_memory, buffer);
		pack32((uint32_t) detail_ptr->job_max_memory, buffer);
		pack32((uint32_t) detail_ptr->job_min_tmp_disk, buffer);
		pack16((uint16_t) detail_ptr->wait_reason, buffer);

		packstr(detail_ptr->req_nodes, buffer);
		pack_bit_fmt(detail_ptr->req_node_bitmap, buffer);
		packstr(detail_ptr->exc_nodes, buffer);
		pack_bit_fmt(detail_ptr->exc_node_bitmap, buffer);
	} 

	else {
		pack16((uint16_t) 0, buffer);
		pack16((uint16_t) 0, buffer);
		pack16((uint16_t) 0, buffer);
		pack16((uint16_t) 0, buffer);
		pack16((uint16_t) 0, buffer);
		pack16((uint16_t) 0, buffer);

		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack32((uint32_t) 0, buffer);
		pack16((uint16_t) 0, buffer);

		packnull(buffer);
		packnull(buffer);
		packnull(buffer);
		packnull(buffer);
	}
}

/*
 * purge_old_job - purge old job records. 
 *	the jobs must have completed at least MIN_JOB_AGE minutes ago
 * global: job_list - global job table
 *	last_job_update - time of last job table update
 * NOTE: READ lock_slurmctld config before entry
 */
void purge_old_job(void)
{
	int i;

	i = list_delete_all(job_list, &_list_find_job_old, "");
	if (i) {
		debug2("purge_old_job: purged %d old job records", i);
/*		last_job_update = time(NULL);	don't worry about state save */
	}
}


/*
 * _purge_job_record - purge specific job record
 * IN job_id - job_id of job record to be purged
 * RET int - count of job's purged
 * global: job_list - global job table
 */
static int _purge_job_record(uint32_t job_id)
{
	return list_delete_all(job_list, &_list_find_job_id, (void *) &job_id);
}


/* 
 * reset_job_bitmaps - reestablish bitmaps for existing jobs. 
 *	this should be called after rebuilding node information, 
 *	but before using any job entries.
 * global: last_job_update - time of last job table update
 *	job_list - pointer to global job list
 */
void reset_job_bitmaps(void)
{
	ListIterator job_iterator;
	struct job_record  *job_ptr;
	struct part_record *part_ptr;
	bool job_fail = false;

	xassert(job_list);

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		xassert (job_ptr->magic == JOB_MAGIC);
		job_fail = false;
		part_ptr = list_find_first(part_list, &list_find_part,
					   job_ptr->partition);
		if (part_ptr == NULL) {
			error("Invalid partition (%s) for job_id %u", 
		    	      job_ptr->partition, job_ptr->job_id);
			job_fail = true;
		}
		job_ptr->part_ptr = part_ptr;

		FREE_NULL_BITMAP(job_ptr->node_bitmap);
		if ((job_ptr->nodes) && 
		    (node_name2bitmap(job_ptr->nodes, false, 
				      &job_ptr->node_bitmap))) {
			error("Invalid nodes (%s) for job_id %u", 
		    	      job_ptr->nodes, job_ptr->job_id);
			job_fail = true;
		}
		build_node_details(job_ptr);	/* set: num_cpu_groups, 
						 * cpu_count_reps, node_cnt, 
						 * cpus_per_node, node_addr */
		if (_reset_detail_bitmaps(job_ptr))
			job_fail = true;

		_reset_step_bitmaps(job_ptr);

		if ((job_ptr->kill_on_step_done) &&
		    (list_count(job_ptr->step_list) <= 1))
			job_fail = true;

		if (job_fail) {
			if (job_ptr->job_state == JOB_PENDING) {
				job_ptr->start_time = 
					job_ptr->end_time = time(NULL);
				job_ptr->job_state = JOB_NODE_FAIL;
			} else if (job_ptr->job_state == JOB_RUNNING) {
				job_ptr->end_time = time(NULL);
				job_ptr->job_state = JOB_NODE_FAIL | 
						     JOB_COMPLETING;
			} else if (job_ptr->job_state == JOB_SUSPENDED) {
				job_ptr->end_time = job_ptr->suspend_time;
				job_ptr->job_state = JOB_NODE_FAIL |
						     JOB_COMPLETING;
			}
			job_ptr->exit_code = MAX(job_ptr->exit_code, 1);
			job_completion_logger(job_ptr);
		}
	}

	list_iterator_destroy(job_iterator);
	last_job_update = time(NULL);
}

static int _reset_detail_bitmaps(struct job_record *job_ptr)
{
	if (job_ptr->details == NULL) 
		return SLURM_SUCCESS;

	FREE_NULL_BITMAP(job_ptr->details->req_node_bitmap);
	if ((job_ptr->details->req_nodes) && 
	    (node_name2bitmap(job_ptr->details->req_nodes, false,  
			      &job_ptr->details->req_node_bitmap))) {
		error("Invalid req_nodes (%s) for job_id %u", 
	    	      job_ptr->details->req_nodes, job_ptr->job_id);
		return SLURM_ERROR;
	}

	FREE_NULL_BITMAP(job_ptr->details->exc_node_bitmap);
	if ((job_ptr->details->exc_nodes) && 
	    (node_name2bitmap(job_ptr->details->exc_nodes, true,
			      &job_ptr->details->exc_node_bitmap))) {
		error("Invalid exc_nodes (%s) for job_id %u", 
	    	      job_ptr->details->exc_nodes, job_ptr->job_id);
		return SLURM_ERROR;
	}

	return SLURM_SUCCESS;
}

static void _reset_step_bitmaps(struct job_record *job_ptr)
{
	ListIterator step_iterator;
	struct step_record *step_ptr;

	step_iterator = list_iterator_create (job_ptr->step_list);
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		FREE_NULL_BITMAP(step_ptr->step_node_bitmap);
		if (step_ptr->step_layout &&
		    step_ptr->step_layout->node_list && 		
		    (node_name2bitmap(step_ptr->step_layout->node_list, false, 
				      &step_ptr->step_node_bitmap))) {
			error("Invalid step_node_list (%s) for step_id %u.%u", 
	   	 	      step_ptr->step_layout->node_list, 
			      job_ptr->job_id, step_ptr->step_id);
			delete_step_record (job_ptr, step_ptr->step_id);
		}
	}		

	list_iterator_destroy (step_iterator);
	return;
}

/* update first assigned job id as needed on reconfigure
 * NOTE: READ lock_slurmctld config before entry */
void reset_first_job_id(void)
{
	if (job_id_sequence < slurmctld_conf.first_job_id)
		job_id_sequence = slurmctld_conf.first_job_id;
}

/*
 * get_next_job_id - return the job_id to be used by default for 
 *	the next job
 */
extern uint32_t get_next_job_id(void)
{
	uint32_t next_id;

	if (job_id_sequence == 0)
		job_id_sequence = slurmctld_conf.first_job_id;
	next_id = job_id_sequence + 1;
	if (next_id >= MIN_NOALLOC_JOBID)
		next_id = slurmctld_conf.first_job_id;
	return next_id;
}

/*
 * _set_job_id - set a default job_id, insure that it is unique
 * IN job_ptr - pointer to the job_record
 */
static void _set_job_id(struct job_record *job_ptr)
{
	uint32_t new_id;

	if (job_id_sequence == 0)
		job_id_sequence = slurmctld_conf.first_job_id;

	xassert(job_ptr);
	xassert (job_ptr->magic == JOB_MAGIC);
	if ((job_ptr->partition == NULL)
	    || (strlen(job_ptr->partition) == 0))
		fatal("_set_job_id: partition not set");

	/* Insure no conflict in job id if we roll over 32 bits */
	while (1) {
		if (++job_id_sequence >= MIN_NOALLOC_JOBID)
			job_id_sequence = slurmctld_conf.first_job_id;
		new_id = job_id_sequence;
		if (find_job_record(new_id) == NULL) {
			job_ptr->job_id = new_id;
			break;
		}
	}
}


/*
 * _set_job_prio - set a default job priority
 * IN job_ptr - pointer to the job_record
 * NOTE: this is a simple prototype, we need to re-establish value on restart
 */
static void _set_job_prio(struct job_record *job_ptr)
{
	xassert(job_ptr);
	xassert (job_ptr->magic == JOB_MAGIC);
	job_ptr->priority = slurm_sched_initial_priority(maximum_prio);
	if (job_ptr->priority > 0)
		maximum_prio = MIN(job_ptr->priority, maximum_prio);
}


/* After a node is returned to service, reset the priority of jobs 
 * which may have been held due to that node being unavailable */
void reset_job_priority(void)
{
	ListIterator job_iterator;
	struct job_record *job_ptr;
	int count = 0;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		if (job_ptr->priority == 1) {
			_set_job_prio(job_ptr);
			count++;
		}
	}
	list_iterator_destroy(job_iterator);
	if (count)
		last_job_update = time(NULL);
}

/* 
 * _top_priority - determine if any other job for this partition has a 
 *	higher priority than specified job
 * IN job_ptr - pointer to selected partition
 * RET true if selected job has highest priority
 */
static bool _top_priority(struct job_record *job_ptr)
{
#ifdef HAVE_BG
	/* On BlueGene, all jobs run ASAP. 
	 * Priority only matters within a specific job size. */
	return true;

#else
	struct job_details *detail_ptr = job_ptr->details;
	bool top;

	if (job_ptr->priority == 0)	/* user held */
		top = false;
	else {
		ListIterator job_iterator;
		struct job_record *job_ptr2;

		top = true;	/* assume top priority until found otherwise */
		job_iterator = list_iterator_create(job_list);
		while ((job_ptr2 = (struct job_record *) 
				list_next(job_iterator))) {
			if (job_ptr2 == job_ptr)
				continue;
			if (job_ptr2->job_state != JOB_PENDING)
				continue;
			if (!job_independent(job_ptr2))
				continue;
			if ((job_ptr2->priority >  job_ptr->priority) &&
			    (job_ptr2->part_ptr == job_ptr->part_ptr)) {
				top = false;
				break;
			}
		}
		list_iterator_destroy(job_iterator);
	}

	if ((!top) && detail_ptr) {	/* not top prio */
		if (job_ptr->priority == 0)		/* user/admin hold */
			detail_ptr->wait_reason = WAIT_HELD;
		else if (job_ptr->priority != 1)	/* not system hold */
			detail_ptr->wait_reason = WAIT_PRIORITY;
	}
	return top;
#endif
}


/*
 * update_job - update a job's parameters per the supplied specifications
 * IN job_specs - a job's specification
 * IN uid - uid of user issuing RPC
 * RET returns an error code from slurm_errno.h
 * global: job_list - global list of job entries
 *	last_job_update - time of last job table update
 */
int update_job(job_desc_msg_t * job_specs, uid_t uid)
{
	int error_code = SLURM_SUCCESS;
	int super_user = 0;
	struct job_record *job_ptr;
	struct job_details *detail_ptr;
	struct part_record *tmp_part_ptr;
	bitstr_t *req_bitmap = NULL;
	time_t now = time(NULL);

	job_ptr = find_job_record(job_specs->job_id);
	if (job_ptr == NULL) {
		error("update_job: job_id %u does not exist.",
		      job_specs->job_id);
		return ESLURM_INVALID_JOB_ID;
	}
	if ((uid == 0) || (uid == slurmctld_conf.slurm_user_id))
		super_user = 1;
	if ((job_ptr->user_id != uid) && (super_user == 0)) {
		error("Security violation, JOB_UPDATE RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}

	detail_ptr = job_ptr->details;
	last_job_update = now;

	if ((job_specs->time_limit != NO_VAL) && (!IS_JOB_FINISHED(job_ptr))) {
		if (super_user ||
		    (job_ptr->time_limit > job_specs->time_limit)) {
			time_t old_time =  job_ptr->time_limit;
			job_ptr->time_limit = job_specs->time_limit;
			if (job_ptr->time_limit == INFINITE)	/* one year */
				job_ptr->end_time = now +
						(365 * 24 * 60 * 60);
			else {
				/* Update end_time based upon change
				 * to preserve suspend time info */
				job_ptr->end_time = job_ptr->end_time +
						((job_ptr->time_limit -
						  old_time) * 60);
			}
			if (job_ptr->end_time < now)
				job_ptr->end_time = now;
			if ((job_ptr->job_state == JOB_RUNNING) &&
			    (list_is_empty(job_ptr->step_list) == 0))
				_xmit_new_end_time(job_ptr);
			info("update_job: setting time_limit to %u for "
				"job_id %u", job_specs->time_limit, 
				job_specs->job_id);
		} else {
			error("Attempt to increase time limit for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->priority != NO_VAL) {
		if (super_user ||
		    (job_ptr->priority > job_specs->priority)) {
			job_ptr->priority = job_specs->priority;
			info("update_job: setting priority to %u for "
				"job_id %u", job_ptr->priority, 
				job_specs->job_id);
		} else {
			error("Attempt to increase priority for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->nice != NICE_OFFSET) {
		if (super_user || (job_specs->nice < NICE_OFFSET)) {
			job_ptr->priority -= ((int)job_specs->nice - 
					NICE_OFFSET);
			info("update_job: setting priority to %u for "
				"job_id %u", job_ptr->priority,
				job_specs->job_id);
		} else {
			error("Attempt to increase priority for job %u",
				job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}
 
	if (job_specs->job_min_procs != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->job_min_procs > job_specs->job_min_procs)) {
			detail_ptr->job_min_procs = job_specs->job_min_procs;
			info("update_job: setting job_min_procs to %u for "
				"job_id %u", job_specs->job_min_procs, 
				job_specs->job_id);
		} else {
			error("Attempt to increase job_min_procs for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->job_min_sockets != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->job_min_sockets > job_specs->job_min_sockets)) {
			detail_ptr->job_min_sockets = job_specs->job_min_sockets;
			info("update_job: setting job_min_sockets to %u for "
				"job_id %u", job_specs->job_min_sockets, 
				job_specs->job_id);
		} else {
			error("Attempt to increase job_min_sockets for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->job_min_cores != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->job_min_cores > job_specs->job_min_cores)) {
			detail_ptr->job_min_cores = job_specs->job_min_cores;
			info("update_job: setting job_min_cores to %u for "
				"job_id %u", job_specs->job_min_cores, 
				job_specs->job_id);
		} else {
			error("Attempt to increase job_min_cores for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}


	if (job_specs->job_min_threads != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->job_min_threads > job_specs->job_min_threads)) {
			detail_ptr->job_min_threads = job_specs->job_min_threads;
			info("update_job: setting job_min_threads to %u for "
				"job_id %u", job_specs->job_min_threads, 
				job_specs->job_id);
		} else {
			error("Attempt to increase job_min_threads for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->job_min_memory != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->job_min_memory > job_specs->job_min_memory)) {
			detail_ptr->job_min_memory = job_specs->job_min_memory;
			info("update_job: setting job_min_memory to %u for "
				"job_id %u", job_specs->job_min_memory, 
				job_specs->job_id);
		} else {
			error("Attempt to increase job_min_memory for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->job_min_tmp_disk != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->job_min_tmp_disk > job_specs->job_min_tmp_disk)) {
			detail_ptr->job_min_tmp_disk = job_specs->job_min_tmp_disk;
			info("update_job: setting job_min_tmp_disk to %u for "
				"job_id %u", job_specs->job_min_tmp_disk, 
				job_specs->job_id);
		} else {
			error
			    ("Attempt to increase job_min_tmp_disk for job %u",
			     job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->num_procs != NO_VAL) {
		if (super_user || 
		    (job_ptr->num_procs > job_specs->num_procs)) {
			job_ptr->num_procs = job_specs->num_procs;
			info("update_job: setting num_procs to %u for "
				"job_id %u", job_specs->num_procs, 
				job_specs->job_id);
		} else {
			error("Attempt to increase num_procs for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->min_nodes != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->min_nodes > job_specs->min_nodes)) {
			detail_ptr->min_nodes = job_specs->min_nodes;
			info("update_job: setting min_nodes to %u for "
				"job_id %u", job_specs->min_nodes, 
				job_specs->job_id);
		} else {
			error("Attempt to increase min_nodes for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->min_sockets != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->min_sockets > job_specs->min_sockets)) {
			detail_ptr->min_sockets = job_specs->min_sockets;
			info("update_job: setting min_sockets to %u for "
				"job_id %u", job_specs->min_sockets, 
				job_specs->job_id);
		} else {
			error("Attempt to increase min_sockets for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->min_cores != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->min_cores > job_specs->min_cores)) {
			detail_ptr->min_cores = job_specs->min_cores;
			info("update_job: setting min_cores to %u for "
				"job_id %u", job_specs->min_cores, 
				job_specs->job_id);
		} else {
			error("Attempt to increase min_cores for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->min_threads != NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->min_threads > job_specs->min_threads)) {
			detail_ptr->min_threads = job_specs->min_threads;
			info("update_job: setting min_threads to %u for "
				"job_id %u", job_specs->min_threads, 
				job_specs->job_id);
		} else {
			error("Attempt to increase min_threads for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->shared != (uint16_t) NO_VAL && detail_ptr) {
		if (super_user || (detail_ptr->shared > job_specs->shared)) {
			detail_ptr->shared = job_specs->shared;
			info("update_job: setting shared to %u for job_id %u", 
			     job_specs->shared, job_specs->job_id);
		} else {
			error("Attempt to remove sharing for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->contiguous != (uint16_t) NO_VAL && detail_ptr) {
		if (super_user ||
		    (detail_ptr->contiguous > job_specs->contiguous)) {
			detail_ptr->contiguous = job_specs->contiguous;
			info("update_job: setting contiguous to %u for "
				"job_id %u", job_specs->contiguous, 
				job_specs->job_id);
		} else {
			error("Attempt to add contiguous for job %u",
				job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->kill_on_node_fail != (uint16_t) NO_VAL) {
		job_ptr->kill_on_node_fail = job_specs->kill_on_node_fail;
		info("update_job: setting kill_on_node_fail to %u for "
			"job_id %u", job_specs->kill_on_node_fail, 
			job_specs->job_id);
	}

	if (job_specs->features && detail_ptr) {
		if (super_user) {
			xfree(detail_ptr->features);
			if (job_specs->features[0] != '\0') {
				detail_ptr->features = job_specs->features;
				job_specs->features = NULL;
				info("update_job: setting features to %s for "
					"job_id %u", job_specs->features, 
					job_specs->job_id);
			}
		} else {
			error("Attempt to change features for job %u",
				job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->name) {
		strncpy(job_ptr->name, job_specs->name, MAX_JOBNAME_LEN);
		info("update_job: setting name to %s for job_id %u",
		     job_specs->name, job_specs->job_id);
	}

	if (job_specs->partition) {
		tmp_part_ptr = find_part_record(job_specs->partition);
		if (tmp_part_ptr == NULL)
			error_code = ESLURM_INVALID_PARTITION_NAME;
		if ((super_user && tmp_part_ptr)) {
			strncpy(job_ptr->partition, job_specs->partition,
				MAX_SLURM_NAME);
			job_ptr->part_ptr = tmp_part_ptr;
			info("update_job: setting partition to %s for "
				"job_id %u", job_specs->partition, 
				job_specs->job_id);
		} else {
			error("Attempt to change partition for job %u",
				job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->req_nodes && detail_ptr) {
		if (job_specs->req_nodes[0] == '\0') {
			xfree(detail_ptr->req_nodes);
			FREE_NULL_BITMAP(detail_ptr->req_node_bitmap);
		} else if (super_user) {
			if (node_name2bitmap(job_specs->req_nodes, false, 
						&req_bitmap)) {
				error("Invalid node list for job_update: %s",
					job_specs->req_nodes);
				FREE_NULL_BITMAP(req_bitmap);
				req_bitmap = NULL;
				error_code = ESLURM_INVALID_NODE_NAME;
			}
			if (req_bitmap) {
				xfree(detail_ptr->req_nodes);
				detail_ptr->req_nodes =
				    job_specs->req_nodes;
				FREE_NULL_BITMAP(detail_ptr->req_node_bitmap);
				detail_ptr->req_node_bitmap = req_bitmap;
				info("update_job: setting req_nodes to %s "
					"for job_id %u", job_specs->req_nodes, 
					job_specs->job_id);
				job_specs->req_nodes = NULL;
			}
		} else {
			error("Attempt to change req_nodes for job %u",
			      job_specs->job_id);
			error_code = ESLURM_ACCESS_DENIED;
		}
	}

	if (job_specs->account) {
		xfree(job_ptr->account);
		if (job_specs->account[0] != '\0') {
			job_ptr->account = job_specs->account ;
			info("update_job: setting account to %s for job_id %u",
				job_ptr->account, job_specs->job_id);
			job_specs->account = NULL;
		}
	}

	if (job_specs->dependency != NO_VAL) {
		if (job_specs->dependency == job_ptr->job_id)
			error_code = ESLURM_DEPENDENCY;
		else {
			job_ptr->dependency = job_specs->dependency;
			info("update_job: setting dependency to %u for " 
				"job_id %u",  job_ptr->dependency, 
				job_ptr->job_id);
		}
	}

	if (job_specs->begin_time) {
		if (IS_JOB_PENDING(job_ptr) && detail_ptr)
			detail_ptr->begin_time = job_specs->begin_time;
		else
			error_code = ESLURM_DISABLED;
	}

	return error_code;
}


/*
 * validate_jobs_on_node - validate that any jobs that should be on the node 
 *	are actually running, if not clean up the job records and/or node 
 *	records
 * IN node_name - node which should have jobs running
 * IN/OUT job_count - number of jobs which should be running on specified node
 * IN job_id_ptr - pointer to array of job_ids that should be on this node
 * IN step_id_ptr - pointer to array of job step ids that should be on node
 */
void
validate_jobs_on_node(char *node_name, uint32_t * job_count,
		      uint32_t * job_id_ptr, uint16_t * step_id_ptr)
{
	int i, node_inx, jobs_on_node;
	struct node_record *node_ptr;
	struct job_record *job_ptr;
	time_t now = time(NULL);

	node_ptr = find_node_record(node_name);
	if (node_ptr == NULL) {
		error("slurmd registered on unknown node %s", node_name);
		return;
	}
	node_inx = node_ptr - node_record_table_ptr;

	/* Check that jobs running are really supposed to be there */
	for (i = 0; i < *job_count; i++) {
		if ( (job_id_ptr[i] >= MIN_NOALLOC_JOBID) && 
		     (job_id_ptr[i] <= MAX_NOALLOC_JOBID) ) {
			info("NoAllocate job %u.%u reported on node %s",
				job_id_ptr[i], step_id_ptr[i], node_name);
			continue;
		}

		job_ptr = find_job_record(job_id_ptr[i]);
		if (job_ptr == NULL) {
			error("Orphan job %u.%u reported on node %s",
			      job_id_ptr[i], step_id_ptr[i], node_name);
			kill_job_on_node(job_id_ptr[i], job_ptr, node_ptr);
		}

		else if ((job_ptr->job_state == JOB_RUNNING) ||
				(job_ptr->job_state == JOB_SUSPENDED)) {
			if (bit_test(job_ptr->node_bitmap, node_inx)) {
				debug3("Registered job %u.%u on node %s ",
				       job_id_ptr[i], step_id_ptr[i], 
				       node_name);
				if ((job_ptr->batch_flag) &&
				    (node_inx == bit_ffs(
						job_ptr->node_bitmap))) {
					/* NOTE: Used for purging defunct
					 * batch jobs */
					job_ptr->time_last_active = now;
				}
			} else {
				error
				    ("Registered job %u.%u on wrong node %s ",
				     job_id_ptr[i], step_id_ptr[i], node_name);
				kill_job_on_node(job_id_ptr[i], job_ptr, 
						node_ptr);
			}
		}

		else if (job_ptr->job_state & JOB_COMPLETING) {
			/* Re-send kill request as needed, 
			 * not necessarily an error */
			kill_job_on_node(job_id_ptr[i], job_ptr, node_ptr);
		}


		else if (job_ptr->job_state == JOB_PENDING) {
			error("Registered PENDING job %u.%u on node %s ",
			      job_id_ptr[i], step_id_ptr[i], node_name);
			job_ptr->job_state = JOB_FAILED;
			job_ptr->exit_code = 1;
			last_job_update    = now;
			job_ptr->start_time = job_ptr->end_time  = now;
			kill_job_on_node(job_id_ptr[i], job_ptr, node_ptr);
			job_completion_logger(job_ptr);
			delete_job_details(job_ptr);
		}

		else {		/* else job is supposed to be done */
			error
			    ("Registered job %u.%u in state %s on node %s ",
			     job_id_ptr[i], step_id_ptr[i], 
			     job_state_string(job_ptr->job_state),
			     node_name);
			kill_job_on_node(job_id_ptr[i], job_ptr, node_ptr);
		}
	}

	jobs_on_node = node_ptr->run_job_cnt + node_ptr->comp_job_cnt;
	if (jobs_on_node)
		_purge_lost_batch_jobs(node_inx, now);

	if (jobs_on_node != *job_count) {
		/* slurmd will not know of a job unless the job has
		 * steps active at registration time, so this is not 
		 * an error condition, slurmd is also reporting steps 
		 * rather than jobs */
		debug3("resetting job_count on node %s from %d to %d", 
		     node_name, *job_count, jobs_on_node);
		*job_count = jobs_on_node;
	}

	return;
}

/* Purge any batch job that should have its script running on node 
 * node_inx, but is not (i.e. its time_last_active != now) */
static void _purge_lost_batch_jobs(int node_inx, time_t now)
{
	ListIterator job_iterator;
	struct job_record *job_ptr;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		bool job_active = ((job_ptr->job_state == JOB_RUNNING) ||
				   (job_ptr->job_state == JOB_SUSPENDED));
		if ((!job_active)                       ||
		    (job_ptr->batch_flag == 0)          ||
		    (job_ptr->time_last_active == now)  ||
		    (node_inx != bit_ffs(job_ptr->node_bitmap)))
			continue;

		info("Master node lost JobId=%u, killing it", 
			job_ptr->job_id);
		job_complete(job_ptr->job_id, 0, false, 0);
	}
	list_iterator_destroy(job_iterator);
}

/*
 * kill_job_on_node - Kill the specific job_id on a specific node,
 *	the request is not processed immediately, but queued. 
 *	This is to prevent a flood of pthreads if slurmctld restarts 
 *	without saved state and slurmd daemons register with a 
 *	multitude of running jobs. Slurmctld will not recognize 
 *	these jobs and use this function to kill them - one 
 *	agent request per node as they register.
 * IN job_id - id of the job to be killed
 * IN job_ptr - pointer to terminating job (NULL if unknown, e.g. orphaned)
 * IN node_ptr - pointer to the node on which the job resides
 */
extern void
kill_job_on_node(uint32_t job_id, struct job_record *job_ptr, 
		struct node_record *node_ptr)
{
	agent_arg_t *agent_info;
	kill_job_msg_t *kill_req;

	debug("Killing job %u on node %s", job_id, node_ptr->name);

	kill_req = xmalloc(sizeof(kill_job_msg_t));
	kill_req->job_id	= job_id;
	kill_req->time          = time(NULL);
	kill_req->nodes	        = xstrdup(node_ptr->name);
	if (job_ptr) {  /* NULL if unknown */
		kill_req->select_jobinfo = 
			select_g_copy_jobinfo(job_ptr->select_jobinfo);
	}

	agent_info = xmalloc(sizeof(agent_arg_t));
	agent_info->node_count	= 1;
	agent_info->retry	= 0;
	agent_info->hostlist	= hostlist_create(node_ptr->name);
	agent_info->msg_type	= REQUEST_TERMINATE_JOB;
	agent_info->msg_args	= kill_req;

	agent_queue_request(agent_info);
}


/*
 * job_alloc_info - get details about an existing job allocation
 * IN uid - job issuing the code
 * IN job_id - ID of job for which info is requested
 * OUT job_pptr - set to pointer to job record
 */
extern int
job_alloc_info(uint32_t uid, uint32_t job_id, struct job_record **job_pptr)
{
	struct job_record *job_ptr;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL)
		return ESLURM_INVALID_JOB_ID;
	if ((job_ptr->user_id != uid) && 
	    (uid != 0) && (uid != slurmctld_conf.slurm_user_id))
		return ESLURM_ACCESS_DENIED;
	if (IS_JOB_PENDING(job_ptr))
		return ESLURM_JOB_PENDING;
	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;

	*job_pptr = job_ptr;
	return SLURM_SUCCESS;
}

/*
 * Synchronize the batch job in the system with their files.
 * All pending batch jobs must have script and environment files
 * No other jobs should have such files
 * NOTE: READ lock_slurmctld config before entry
 */
int sync_job_files(void)
{
	List batch_dirs;

	batch_dirs = list_create(_del_batch_list_rec);
	_get_batch_job_dir_ids(batch_dirs);
	_validate_job_files(batch_dirs);
	_remove_defunct_batch_dirs(batch_dirs);
	list_destroy(batch_dirs);
	return SLURM_SUCCESS;
}

/* Append to the batch_dirs list the job_id's associated with 
 *	every batch job directory in existence
 * NOTE: READ lock_slurmctld config before entry
 */
static void _get_batch_job_dir_ids(List batch_dirs)
{
	DIR *f_dir;
	struct dirent *dir_ent;
	long long_job_id;
	uint32_t *job_id_ptr;
	char *endptr;

	xassert(slurmctld_conf.state_save_location);
	f_dir = opendir(slurmctld_conf.state_save_location);
	if (!f_dir) {
		error("opendir(%s): %m", 
		      slurmctld_conf.state_save_location);
		return;
	}

	while ((dir_ent = readdir(f_dir))) {
		if (strncmp("job.#", dir_ent->d_name, 4))
			continue;
		long_job_id = strtol(&dir_ent->d_name[4], &endptr, 10);
		if ((long_job_id == 0) || (endptr[0] != '\0'))
			continue;
		debug3("found batch directory for job_id %ld",long_job_id);
		job_id_ptr = xmalloc(sizeof(uint32_t));
		*job_id_ptr = long_job_id;
		list_append (batch_dirs, job_id_ptr);
	}

	closedir(f_dir);
}

/* All pending batch jobs must have a batch_dir entry, 
 *	otherwise we flag it as FAILED and don't schedule
 * If the batch_dir entry exists for a PENDING or RUNNING batch job, 
 *	remove it the list (of directories to be deleted) */
static void _validate_job_files(List batch_dirs)
{
	ListIterator job_iterator;
	struct job_record *job_ptr;
	int del_cnt;

	job_iterator = list_iterator_create(job_list);
	while ((job_ptr = (struct job_record *) list_next(job_iterator))) {
		if (!job_ptr->batch_flag)
			continue;
		if (IS_JOB_FINISHED(job_ptr))
			continue;
		/* Want to keep this job's files */
		del_cnt = list_delete_all(batch_dirs, _find_batch_dir, 
					  &(job_ptr->job_id));
		if ((del_cnt == 0) && 
		    (job_ptr->job_state == JOB_PENDING)) {
			error("Script for job %u lost, state set to FAILED",
			      job_ptr->job_id);
			job_ptr->job_state = JOB_FAILED;
			job_ptr->exit_code = 1;
			job_ptr->start_time = job_ptr->end_time = time(NULL);
			job_completion_logger(job_ptr);
		}
	}
	list_iterator_destroy(job_iterator);
}

/* List matching function, see common/list.h */
static int _find_batch_dir(void *x, void *key)
{
	uint32_t *key1 = x;
	uint32_t *key2 = key;
	return (int)(*key1 == *key2);
}
/* List entry deletion function, see common/list.h */
static void _del_batch_list_rec(void *x)
{
	xfree(x);
}

/* Remove all batch_dir entries in the list
 * NOTE: READ lock_slurmctld config before entry */
static void _remove_defunct_batch_dirs(List batch_dirs)
{
	ListIterator batch_dir_inx;
	uint32_t *job_id_ptr;

	batch_dir_inx = list_iterator_create(batch_dirs);
	while ((job_id_ptr = list_next(batch_dir_inx))) {
		error("Purging files for defunct batch job %u",
		      *job_id_ptr);
		_delete_job_desc_files(*job_id_ptr);
	}
	list_iterator_destroy(batch_dir_inx);
}

/*
 *  _xmit_new_end_time
 *	Tell all slurmd's associated with a job of its new end time
 * IN job_ptr - pointer to terminating job
 * globals: node_record_count - number of nodes in the system
 *	node_record_table_ptr - pointer to global node table
 */
static void 
_xmit_new_end_time(struct job_record *job_ptr)
{
	job_time_msg_t *job_time_msg_ptr;
	agent_arg_t *agent_args;
	int i;

	agent_args = xmalloc(sizeof(agent_arg_t));
	agent_args->msg_type = REQUEST_UPDATE_JOB_TIME;
	agent_args->retry = 1;
	agent_args->hostlist = hostlist_create("");
	job_time_msg_ptr = xmalloc(sizeof(job_time_msg_t));
	job_time_msg_ptr->job_id          = job_ptr->job_id;
	job_time_msg_ptr->expiration_time = job_ptr->end_time;

	for (i = 0; i < node_record_count; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		hostlist_push(agent_args->hostlist,
			      node_record_table_ptr[i].name);
		agent_args->node_count++;
#ifdef HAVE_FRONT_END		/* operate only on front-end node */
		break;
#endif
	}

	agent_args->msg_args = job_time_msg_ptr;
	agent_queue_request(agent_args);
	return;
}


/*
 * job_epilog_complete - Note the completion of the epilog script for a 
 *	given job
 * IN job_id      - id of the job for which the epilog was executed
 * IN node_name   - name of the node on which the epilog was executed
 * IN return_code - return code from epilog script
 * RET true if job is COMPLETED, otherwise false
 */
extern bool job_epilog_complete(uint32_t job_id, char *node_name, 
		uint32_t return_code)
{
	struct job_record  *job_ptr = find_job_record(job_id);

	if (job_ptr == NULL)
		return true;

	/* There is a potential race condition this handles.
	 * If slurmctld cold-starts while slurmd keeps running, 
	 * slurmd could notify slurmctld of a job epilog completion 
	 * before getting synced up with slurmctld state. If 
	 * a new job arrives and the job_id is reused, we 
	 * could try to note the termination of a job that 
	 * hasn't really started. Very rare obviously. */
	if ((job_ptr->job_state == JOB_PENDING)
	||  (job_ptr->node_bitmap == NULL)) {
		error("Epilog complete request for non-running job %u, "
			"slurmctld and slurmd out of sync", job_id);
		return false;
	}

#ifdef HAVE_FRONT_END		/* operate only on front-end node */
{
	int i;
	struct node_record *node_ptr;

	if (return_code)
		error("Epilog error on %s, setting DOWN", 
			job_ptr->nodes);
	for (i=0; i<node_record_count; i++) {
		if (!bit_test(job_ptr->node_bitmap, i))
			continue;
		node_ptr = &node_record_table_ptr[i];
		if (return_code)
			set_node_down(node_ptr->name, "Epilog error");
		else
			make_node_idle(node_ptr, job_ptr);
	}
}
#else
	if (return_code) {
		error("Epilog error on %s, setting DOWN", node_name);
		set_node_down(node_name, "Epilog error");
	} else {
		struct node_record *node_ptr = find_node_record(node_name);
		if (node_ptr)
			make_node_idle(node_ptr, job_ptr);
	}
#endif

	step_epilog_complete(job_ptr, node_name);
	if (!(job_ptr->job_state & JOB_COMPLETING)) {	/* COMPLETED */
		if ((job_ptr->job_state == JOB_PENDING)
		&&  (job_ptr->batch_flag)) {
			info("requeue batch job %u", job_ptr->job_id);
			if (job_ptr->details) {
				/* the time stamp on the new batch launch 
				 * credential must be larger than the time 
				 * stamp on the revoke request, so delay 
				 * for at least two seconds. */
				job_ptr->details->begin_time = time(NULL) + 2;
			}
		}
		return true;
	} else
		return false;
}

/* job_fini - free all memory associated with job records */
void job_fini (void) 
{
	if (job_list) {
		list_destroy(job_list);
		job_list = NULL;
	}
	xfree(job_hash);
}

/* log the completion of the specified job */
extern void job_completion_logger(struct job_record  *job_ptr)
{
	int base_state;
	xassert(job_ptr);

	base_state = job_ptr->job_state & (~JOB_COMPLETING);
	if ((base_state == JOB_COMPLETE) || (base_state == JOB_CANCELLED)) {
		if (job_ptr->mail_type & MAIL_JOB_END)
			mail_job_info(job_ptr, MAIL_JOB_END);
	} else {	/* JOB_FAILED, JOB_NODE_FAIL, or JOB_TIMEOUT */
		if (job_ptr->mail_type & MAIL_JOB_FAIL)
			mail_job_info(job_ptr, MAIL_JOB_FAIL);
	}

	jobacct_g_job_complete_slurmctld(job_ptr);
	g_slurm_jobcomp_write(job_ptr);
	srun_complete(job_ptr);
}

/*
 * job_independent - determine if this job has a depenendent job pending
 *	or if the job's scheduled begin time is in the future
 * IN job_ptr - pointer to job being tested
 * RET - true if job no longer must be defered for another job
 */
extern bool job_independent(struct job_record *job_ptr)
{
	struct job_record *dep_ptr;
	struct job_details *detail_ptr = job_ptr->details;

	if (detail_ptr && (detail_ptr->begin_time > time(NULL))) {
		detail_ptr->wait_reason = WAIT_TIME;
		return false;	/* not yet time */
	}
		
	if (job_ptr->dependency == 0)
		return true;

	dep_ptr = find_job_record(job_ptr->dependency);
	if (dep_ptr == NULL)
		return true;

	if (((dep_ptr->job_state & JOB_COMPLETING) == 0) &&
	    (dep_ptr->job_state >= JOB_COMPLETE))
		return true;

	if (detail_ptr)
		detail_ptr->wait_reason = WAIT_DEPENDENCY;
	return false;	/* job exists and incomplete */
}
/*
 * determine if job is ready to execute per the node select plugin
 * IN job_id - job to test
 * OUT ready - 1 if job is ready to execute 0 otherwise
 * RET SLURM error code
 */
extern int job_node_ready(uint32_t job_id, int *ready)
{
	int rc;
	struct job_record *job_ptr;
	xassert(ready);

	*ready = 0;
	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL)
		return ESLURM_INVALID_JOB_ID;

	rc = select_g_job_ready(job_ptr);
	if (rc == READY_JOB_FATAL)
		return ESLURM_INVALID_PARTITION_NAME;
	if (rc == READY_JOB_ERROR)
		return EAGAIN;

	if (rc)
		rc = READY_NODE_STATE;
	if (job_ptr->job_state == JOB_RUNNING)
		rc |= READY_JOB_STATE;

	*ready = rc;
	return SLURM_SUCCESS;
}

/* Send specified signal to all steps associated with a job */
static void _signal_job(struct job_record *job_ptr, int signal)
{
	agent_arg_t *agent_args = NULL;
	signal_job_msg_t *signal_job_msg = NULL;
	int i;

	agent_args = xmalloc(sizeof(agent_arg_t));
	agent_args->msg_type = REQUEST_SIGNAL_JOB;
	agent_args->retry = 1;
	agent_args->hostlist = hostlist_create("");
	signal_job_msg = xmalloc(sizeof(kill_tasks_msg_t));
	signal_job_msg->job_id = job_ptr->job_id;
	signal_job_msg->signal = signal;

	for (i = 0; i < node_record_count; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		hostlist_push(agent_args->hostlist,
			      node_record_table_ptr[i].name);
		agent_args->node_count++;
#ifdef HAVE_FRONT_END	/* Operate only on front-end */
		break;
#endif
	}

	if (agent_args->node_count == 0) {
		xfree(signal_job_msg);
		xfree(agent_args);
		return;
	}

	agent_args->msg_args = signal_job_msg;
	agent_queue_request(agent_args);
	return;
}

/* Send suspend request to slumrd of all nodes associated with a job */
static void _suspend_job(struct job_record *job_ptr, uint16_t op)
{
	agent_arg_t *agent_args;
	suspend_msg_t *sus_ptr;
	int i;

	agent_args = xmalloc(sizeof(agent_arg_t));
	agent_args->msg_type = REQUEST_SUSPEND;
	agent_args->retry = 1;
	agent_args->hostlist = hostlist_create("");
	sus_ptr = xmalloc(sizeof(suspend_msg_t));
	sus_ptr->job_id = job_ptr->job_id;
	sus_ptr->op = op;

	for (i = 0; i < node_record_count; i++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		hostlist_push(agent_args->hostlist,
			      node_record_table_ptr[i].name);
		agent_args->node_count++;
#ifdef HAVE_FRONT_END	/* Operate only on front-end */
		break;
#endif
	}

	if (agent_args->node_count == 0) {
		xfree(sus_ptr);
		xfree(agent_args);
		return;
	}

	agent_args->msg_args = sus_ptr;
	agent_queue_request(agent_args);
	return;
}
/* Specified job is being suspended, release allocated nodes */
static int _suspend_job_nodes(struct job_record *job_ptr)
{
	int i, rc;
	struct node_record *node_ptr = node_record_table_ptr;
	uint16_t base_state, node_flags;

	if ((rc = select_g_job_suspend(job_ptr)) != SLURM_SUCCESS)
		return rc;

	for (i=0; i<node_record_count; i++, node_ptr++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;

		if (node_ptr->run_job_cnt)
			(node_ptr->run_job_cnt)--;
		else {
			error("Node %s run_job_cnt underflow", 
				node_ptr->name);
		}
		if (job_ptr->details
		&&  (job_ptr->details->shared == 0)) {
			if (node_ptr->no_share_job_cnt)
				(node_ptr->no_share_job_cnt)--;
			else {
				error("Node %s no_share_job_cnt "
					"underflow", node_ptr->name);
			}
			if (node_ptr->no_share_job_cnt == 0)
				bit_set(share_node_bitmap, i);
		}
		base_state = node_ptr->node_state & NODE_STATE_BASE;
		node_flags = node_ptr->node_state & NODE_STATE_FLAGS;
		if ((node_ptr->run_job_cnt  == 0)
		&&  (node_ptr->comp_job_cnt == 0)) {
			bit_set(idle_node_bitmap, i);
		}
		if (base_state == NODE_STATE_DOWN) {
			debug3("_suspend_job_nodes: Node %s left DOWN",
				node_ptr->name);
		} else if (node_ptr->run_job_cnt) {
			node_ptr->node_state = NODE_STATE_ALLOCATED | 
					node_flags;
		} else {
			node_ptr->node_state = NODE_STATE_IDLE | 
					node_flags;
		}
	}
	last_job_update = last_node_update = time(NULL);
	return rc;
}

/* Specified job is being resumed, re-allocate the nodes */
static int _resume_job_nodes(struct job_record *job_ptr)
{
	int i, rc;
	struct node_record *node_ptr = node_record_table_ptr;
	uint16_t base_state, node_flags;

	if ((rc = select_g_job_resume(job_ptr)) != SLURM_SUCCESS)
		return rc;

	for (i=0; i<node_record_count; i++, node_ptr++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;
		base_state = node_ptr->node_state & NODE_STATE_BASE;
		if (base_state == NODE_STATE_DOWN)
			return SLURM_ERROR;
	}

	node_ptr = node_record_table_ptr;
	for (i=0; i<node_record_count; i++, node_ptr++) {
		if (bit_test(job_ptr->node_bitmap, i) == 0)
			continue;

		node_ptr->run_job_cnt++;
		if (job_ptr->details
		&&  (job_ptr->details->shared == 0)) {
			node_ptr->no_share_job_cnt++;
			if (node_ptr->no_share_job_cnt)
				bit_clear(share_node_bitmap, i);
		}
		bit_clear(idle_node_bitmap, i);
		node_flags = node_ptr->node_state & NODE_STATE_FLAGS;
		node_ptr->node_state = NODE_STATE_ALLOCATED |
				node_flags;
	}
	last_job_update = last_node_update = time(NULL);
	return rc;
}


/*
 * job_suspend - perform some suspend/resume operation
 * IN sus_ptr - suspend/resume request message
 * IN uid - user id of the user issuing the RPC
 * IN conn_fd - file descriptor on which to send reply, 
 *              -1 if none
 * RET 0 on success, otherwise ESLURM error code
 */
extern int job_suspend(suspend_msg_t *sus_ptr, uid_t uid, 
		slurm_fd conn_fd)
{
	int rc = SLURM_SUCCESS;
	time_t now = time(NULL);
	struct job_record *job_ptr = NULL;
	slurm_msg_t resp_msg;
	return_code_msg_t rc_msg;

	/* test if this system configuration
	 * supports job suspend/resume */
	if (strcasecmp(slurmctld_conf.switch_type,
			"switch/federation") == 0) {
		/* Work is needed to support the
		 * release and reuse of switch
		 * windows associated with a job */
		rc = ESLURM_NOT_SUPPORTED;
	}
#ifdef HAVE_BG
	rc = ESLURM_NOT_SUPPORTED;
#endif
	if (rc)
		goto reply;

	/* find the job */
	job_ptr = find_job_record (sus_ptr->job_id);
	if (job_ptr == NULL) {
		rc = ESLURM_INVALID_JOB_ID;
		goto reply;
	}

	/* validate the request */
	if ((uid != 0) && (uid != getuid())) {
		rc = ESLURM_ACCESS_DENIED;
		goto reply;
	}
	if (job_ptr->job_state == JOB_PENDING) {
		rc = ESLURM_JOB_PENDING;
		goto reply;
	}
	if (IS_JOB_FINISHED(job_ptr)) {
		rc = ESLURM_ALREADY_DONE;
		goto reply;
	}

	/* perform the operation */
	if (sus_ptr->op == SUSPEND_JOB) {
		if (job_ptr->job_state != JOB_RUNNING) {
			rc = ESLURM_DISABLED;
			goto reply;
		}
		rc = _suspend_job_nodes(job_ptr);
		if (rc != SLURM_SUCCESS)
			goto reply;
		_suspend_job(job_ptr, sus_ptr->op);
		job_ptr->job_state = JOB_SUSPENDED;
		if (job_ptr->suspend_time) {
			job_ptr->pre_sus_time +=
				difftime(now, 
				job_ptr->suspend_time);
		} else {
			job_ptr->pre_sus_time +=
				difftime(now,
				job_ptr->start_time);
		}
	} else if (sus_ptr->op == RESUME_JOB) {
		if (job_ptr->job_state != JOB_SUSPENDED) {
			rc = ESLURM_DISABLED;
			goto reply;
		}
		rc = _resume_job_nodes(job_ptr);
		if (rc != SLURM_SUCCESS)
			goto reply;
		_suspend_job(job_ptr, sus_ptr->op);
		job_ptr->job_state = JOB_RUNNING;
		if (job_ptr->time_limit != INFINITE) {
			/* adjust effective time_limit */
			job_ptr->end_time = now +
				(job_ptr->time_limit * 60)
				- job_ptr->pre_sus_time;
		}
	}

	job_ptr->time_last_active = now;
	job_ptr->suspend_time = now;
	
    reply:
	if(job_ptr)
		jobacct_g_suspend_slurmctld(job_ptr);

	if (conn_fd >= 0) {
		slurm_msg_t_init(&resp_msg);
		resp_msg.msg_type  = RESPONSE_SLURM_RC;
		rc_msg.return_code = rc;
		resp_msg.data      = &rc_msg;
		slurm_send_node_msg(conn_fd, &resp_msg);
	}
	return rc;
}

/*
 * job_requeue - Requeue a running or pending batch job
 * IN uid - user id of user issuing the RPC
 * IN job_id - id of the job to be requeued
 * IN conn_fd - file descriptor on which to send reply
 * RET 0 on success, otherwise ESLURM error code
 */
extern int job_requeue (uid_t uid, uint32_t job_id, slurm_fd conn_fd)
{
	int rc = SLURM_SUCCESS;
	struct job_record *job_ptr = NULL;
	bool super_user = false, suspended = false;
	slurm_msg_t resp_msg;
	return_code_msg_t rc_msg;
	time_t now = time(NULL);

	/* find the job */
	job_ptr = find_job_record (job_id);
	if (job_ptr == NULL) {
		rc = ESLURM_INVALID_JOB_ID;
		goto reply;
	}

	/* validate the request */
	if ((uid == 0) || (uid == slurmctld_conf.slurm_user_id))
		super_user = 1;
	if ((uid != job_ptr->user_id) && (!super_user)) {
		rc = ESLURM_ACCESS_DENIED;
		goto reply;
	}
	if (IS_JOB_FINISHED(job_ptr)) {
		rc = ESLURM_ALREADY_DONE;
		goto reply;
	}
	if (job_ptr->details && job_ptr->details->no_requeue) {
		rc = ESLURM_DISABLED;
		goto reply;
	}
	if (job_ptr->job_state & JOB_COMPLETING) {
		rc = ESLURM_TRANSITION_STATE_NO_UPDATE;
		goto reply;
	}

	/* reset the priority */
	_set_job_prio(job_ptr);
	last_job_update = now;

	/* nothing else to do if pending */
	if (job_ptr->job_state == JOB_PENDING)
		goto reply;

	if (job_ptr->batch_flag == 0) {
		rc = ESLURM_BATCH_ONLY;
		goto reply;
	}

	if ((job_ptr->job_state != JOB_SUSPENDED)
	&&  (job_ptr->job_state != JOB_RUNNING)) {
		error("job_requeue job %u state is bad %s", job_id,
			job_state_string(job_ptr->job_state));
		rc = EINVAL;
		goto reply;
	}

	if (job_ptr->job_state == JOB_SUSPENDED)
		suspended = true;
	job_ptr->time_last_active  = now;
	job_ptr->job_state         = JOB_PENDING | JOB_COMPLETING;
	if (suspended)
		job_ptr->end_time = job_ptr->suspend_time;
	else
		job_ptr->end_time = now;
	deallocate_nodes(job_ptr, false, suspended);
	job_completion_logger(job_ptr);
//FIXME: Test accounting

    reply:
	if (conn_fd >= 0) {
		slurm_msg_t_init(&resp_msg);
		resp_msg.msg_type  = RESPONSE_SLURM_RC;
		rc_msg.return_code = rc;
		resp_msg.data      = &rc_msg;
		slurm_send_node_msg(conn_fd, &resp_msg);
	}
	return rc;
}

/*
 * job_end_time - Process JOB_END_TIME
 * IN time_req_msg - job end time request
 * OUT timeout_msg - job timeout response to be sent
 * RET SLURM_SUCESS or an error code
 */
extern int job_end_time(job_alloc_info_msg_t *time_req_msg,
			srun_timeout_msg_t *timeout_msg)
{
	struct job_record *job_ptr;
	xassert(timeout_msg);

	job_ptr = find_job_record(time_req_msg->job_id);
	if (!job_ptr)
		return ESLURM_INVALID_JOB_ID;

	timeout_msg->job_id  = time_req_msg->job_id;
	timeout_msg->step_id = NO_VAL;
	timeout_msg->timeout = job_ptr->end_time;
	return SLURM_SUCCESS;
}
