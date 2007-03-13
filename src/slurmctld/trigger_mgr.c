/*****************************************************************************\
 *  trigger_mgr.c - Event trigger management
 *****************************************************************************
 *  Copyright (C) 2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> et. al.
 *  UCRL-CODE-226842.
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
#endif

#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "src/common/bitstring.h"
#include "src/common/list.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/trigger_mgr.h"

#define _DEBUG 1
#define MAX_PROG_TIME 300	/* maximum run time for program */

/* Change TRIGGER_STATE_VERSION value when changing the state save format */
#define TRIGGER_STATE_VERSION      "VER001"

List trigger_list;
uint32_t next_trigger_id = 1;
static pthread_mutex_t trigger_mutex = PTHREAD_MUTEX_INITIALIZER;
bitstr_t *trigger_down_nodes_bitmap = NULL;
bitstr_t *trigger_up_nodes_bitmap   = NULL;
static bool trigger_node_reconfig = false;

typedef struct trig_mgr_info {
	uint32_t trig_id;	/* trigger ID */
	uint8_t  res_type;	/* TRIGGER_RES_TYPE_* */
	char *   res_id;	/* node name or job_id (string) */
	bitstr_t *nodes_bitmap;	/* bitmap of requested nodes (if applicable) */
	uint32_t job_id;	/* job ID (if applicable) */
	struct job_record *job_ptr; /* pointer to job record (if applicable) */
	uint16_t trig_type;	/* TRIGGER_TYPE_* */
	time_t   trig_time;	/* offset (pending) or time stamp (complete) */
	uint32_t user_id;	/* user requesting trigger */
	uint32_t group_id;	/* user's group id (pending) or pid (complete) */
	char *   program;	/* program to execute */
	uint8_t  state;		/* 0=pending, 1=pulled, 2=completed */
} trig_mgr_info_t;

/* Prototype for ListDelF */
void _trig_del(void *x) {
	trig_mgr_info_t * tmp = (trig_mgr_info_t *) x;
	xfree(tmp->res_id);
	xfree(tmp->program);
	FREE_NULL_BITMAP(tmp->nodes_bitmap);
	xfree(tmp);
}

static char *_res_type(uint8_t res_type)
{
	if      (res_type == TRIGGER_RES_TYPE_JOB)
		return "job";
	else if (res_type == TRIGGER_RES_TYPE_NODE)
		return "node";
	else
		return "unknown";
}

static char *_trig_type(uint16_t trig_type)
{
	if      (trig_type == TRIGGER_TYPE_UP)
		return "up";
	else if (trig_type == TRIGGER_TYPE_DOWN)
		return "down";
	else if (trig_type == TRIGGER_TYPE_TIME)
		return "time";
	else if (trig_type == TRIGGER_TYPE_FINI)
		return "fini";
	else if (trig_type == TRIGGER_TYPE_RECONFIG)
		return "reconfig";
	else
		return "unknown";
}

static int _trig_offset(uint16_t offset)
{
	static int rc;
	rc  = offset;
	rc -= 0x8000;
	return rc;
}

static void _dump_trigger_msg(char *header, trigger_info_msg_t *msg)
{
#if _DEBUG
	int i;

	info(header);
	if ((msg == NULL) || (msg->record_count == 0)) {
		info("Trigger has no entries");
		return;
	}

	info("INDEX TRIG_ID RES_TYPE RES_ID TRIG_TYPE OFFSET UID PROGRAM");
	for (i=0; i<msg->record_count; i++) {
		info("trigger[%u] %u %s %s %s %d %u %s", i,
			msg->trigger_array[i].trig_id,
			_res_type(msg->trigger_array[i].res_type),
			msg->trigger_array[i].res_id,
			_trig_type(msg->trigger_array[i].trig_type),
			_trig_offset(msg->trigger_array[i].offset),
			msg->trigger_array[i].user_id,
			msg->trigger_array[i].program);
	}
#endif
}

extern int trigger_clear(uid_t uid, trigger_info_msg_t *msg)
{
	int rc = ESRCH;
	ListIterator trig_iter;
	trigger_info_t *trig_in;
	trig_mgr_info_t *trig_test;
	uint32_t job_id = 0;

	slurm_mutex_lock(&trigger_mutex);
	if (trigger_list == NULL)
		trigger_list = list_create(_trig_del);

	/* validate the request, need a job_id and/or trigger_id */
	_dump_trigger_msg("trigger_clear", msg);
	if (msg->record_count != 1)
		goto fini;
	trig_in = msg->trigger_array;
	if (trig_in->res_type == TRIGGER_RES_TYPE_JOB) {
		job_id = (uint32_t) atol(trig_in->res_id);
		if (job_id == 0) {
			rc = ESLURM_INVALID_JOB_ID;
			goto fini;
		}
	} else if (trig_in->trig_id == 0) {
		rc = EINVAL;
		goto fini;
	}

	/* now look for a valid request, matching uid */
	trig_iter = list_iterator_create(trigger_list);
	while ((trig_test = list_next(trig_iter))) {
		if ((trig_test->user_id != (uint32_t) uid)
		&&  (uid != 0))
			continue;
		if (trig_in->trig_id
		&&  (trig_in->trig_id != trig_test->trig_id))
			continue;
		if (job_id
		&& (job_id != trig_test->job_id))
			continue;
		list_delete(trig_iter);
		rc = SLURM_SUCCESS;
	}
	list_iterator_destroy(trig_iter);

fini:	slurm_mutex_unlock(&trigger_mutex);
	return rc;
}

extern trigger_info_msg_t * trigger_get(uid_t uid, trigger_info_msg_t *msg)
{
	trigger_info_msg_t *resp_data;
	ListIterator trig_iter;
	trigger_info_t *trig_out;
	trig_mgr_info_t *trig_in;
	int recs_written = 0;

	slurm_mutex_lock(&trigger_mutex);
	if (trigger_list == NULL)
		trigger_list = list_create(_trig_del);

	_dump_trigger_msg("trigger_get", NULL);
	resp_data = xmalloc(sizeof(trigger_info_msg_t));
	resp_data->record_count = list_count(trigger_list);
	resp_data->trigger_array = xmalloc(sizeof(trigger_info_t) *
			resp_data->record_count);
	trig_iter = list_iterator_create(trigger_list);
	trig_out = resp_data->trigger_array;
	while ((trig_in = list_next(trig_iter))) {
		/* Note: Filtering currently done by strigger */
		if (trig_in->state > 1)
			continue;	/* no longer pending */
		trig_out->trig_id   = trig_in->trig_id;
		trig_out->res_type  = trig_in->res_type;
		trig_out->res_id    = xstrdup(trig_in->res_id);
		trig_out->trig_type = trig_in->trig_type;
		trig_out->offset    = trig_in->trig_time;
		trig_out->user_id   = trig_in->user_id;
		trig_out->program   = xstrdup(trig_in->program);
		trig_out++;
		recs_written++;
	}
	list_iterator_destroy(trig_iter);
	slurm_mutex_unlock(&trigger_mutex);
	resp_data->record_count = recs_written;

	_dump_trigger_msg("trigger_got", resp_data);
	return resp_data;
}

extern int trigger_set(uid_t uid, gid_t gid, trigger_info_msg_t *msg)
{
	int i;
	int rc = SLURM_SUCCESS;
	uint32_t job_id;
	bitstr_t *bitmap = NULL;
	trig_mgr_info_t * trig_add;
	struct job_record *job_ptr;
	slurmctld_lock_t job_read_lock =
		{ NO_LOCK, READ_LOCK, NO_LOCK, NO_LOCK };

	lock_slurmctld(job_read_lock);
	slurm_mutex_lock(&trigger_mutex);
	if (trigger_list == NULL) {
		trigger_list = list_create(_trig_del);
	} else if ((uid != 0) &&
	           (list_count(trigger_list) >= slurmctld_conf.max_job_cnt)) {
		rc = EAGAIN;
		goto fini;
	}

	_dump_trigger_msg("trigger_set", msg);
	for (i=0; i<msg->record_count; i++) {
		if (msg->trigger_array[i].res_type ==
				TRIGGER_RES_TYPE_JOB) {
			job_id = (uint32_t) atol(
				msg->trigger_array[i].res_id);
			job_ptr = find_job_record(job_id);
			if (job_ptr == NULL) {
				rc = ESLURM_INVALID_JOB_ID;
				continue;
			}
			if (IS_JOB_FINISHED(job_ptr)) {
				rc = ESLURM_ALREADY_DONE;
				continue;
			}
		} else {
			job_id = 0;
			job_ptr = NULL;
			if ((msg->trigger_array[i].res_id != NULL)
			&&  (msg->trigger_array[i].res_id[0] != '*')
			&&  (node_name2bitmap(msg->trigger_array[i].res_id,
						false, &bitmap) != 0)) {
				rc = ESLURM_INVALID_NODE_NAME;
				continue;
			}
		}
		trig_add = xmalloc(sizeof(trig_mgr_info_t));
		msg->trigger_array[i].trig_id = next_trigger_id;
		trig_add->trig_id = next_trigger_id;
		next_trigger_id++;
		trig_add->res_type = msg->trigger_array[i].res_type;
		trig_add->nodes_bitmap = bitmap;
		trig_add->job_id = job_id;
		trig_add->job_ptr = job_ptr;
		/* move don't copy "res_id" */
		trig_add->res_id = msg->trigger_array[i].res_id;
		msg->trigger_array[i].res_id = NULL;
		trig_add->trig_type = msg->trigger_array[i].trig_type;
		trig_add->trig_time = msg->trigger_array[i].offset;
		trig_add->user_id = (uint32_t) uid;
		trig_add->group_id = (uint32_t) gid;
		/* move don't copy "program" */
		trig_add->program = msg->trigger_array[i].program;
		msg->trigger_array[i].program = NULL;
		list_append(trigger_list, trig_add);
	}

fini:	slurm_mutex_unlock(&trigger_mutex);
	unlock_slurmctld(job_read_lock);
	return rc;
}

extern void trigger_node_down(struct node_record *node_ptr)
{
        int inx = node_ptr - node_record_table_ptr;

	slurm_mutex_lock(&trigger_mutex);
	if (trigger_down_nodes_bitmap == NULL)
		trigger_down_nodes_bitmap = bit_alloc(node_record_count);
	bit_set(trigger_down_nodes_bitmap, inx);
	slurm_mutex_unlock(&trigger_mutex);
}

extern void trigger_node_up(struct node_record *node_ptr)
{
        int inx = node_ptr - node_record_table_ptr;

	slurm_mutex_lock(&trigger_mutex);
	if (trigger_up_nodes_bitmap == NULL)
		trigger_up_nodes_bitmap = bit_alloc(node_record_count);
	bit_set(trigger_up_nodes_bitmap, inx);
	slurm_mutex_unlock(&trigger_mutex);
}

extern void trigger_reconfig(void)
{
	slurm_mutex_lock(&trigger_mutex);
	trigger_node_reconfig = true;
	slurm_mutex_unlock(&trigger_mutex);
}

static void _dump_trigger_state(trig_mgr_info_t *trig_ptr, Buf buffer)
{
	pack32   (trig_ptr->trig_id,   buffer);
	pack8    (trig_ptr->res_type,  buffer);
	packstr  (trig_ptr->res_id,    buffer);
	/* rebuild nodes_bitmap as needed from res_id */
	/* rebuild job_id as needed from res_id */
	/* rebuild job_ptr as needed from res_id */
	pack16   (trig_ptr->trig_type, buffer);
	pack_time(trig_ptr->trig_time, buffer);
	pack32   (trig_ptr->user_id,   buffer);
	pack32   (trig_ptr->group_id,  buffer);
	packstr  (trig_ptr->program,   buffer);
	pack8    (trig_ptr->state,     buffer);
}

static int _load_trigger_state(Buf buffer)
{
	trig_mgr_info_t *trig_ptr;
	uint16_t str_len;

	trig_ptr = xmalloc(sizeof(trig_mgr_info_t));
	safe_unpack32   (&trig_ptr->trig_id,   buffer);
	safe_unpack8    (&trig_ptr->res_type,  buffer);
	safe_unpackstr_xmalloc(&trig_ptr->res_id, &str_len, buffer);
	/* rebuild nodes_bitmap as needed from res_id */
	/* rebuild job_id as needed from res_id */
	/* rebuild job_ptr as needed from res_id */
	safe_unpack16   (&trig_ptr->trig_type, buffer);
	safe_unpack_time(&trig_ptr->trig_time, buffer);
	safe_unpack32   (&trig_ptr->user_id,   buffer);
	safe_unpack32   (&trig_ptr->group_id,  buffer);
	safe_unpackstr_xmalloc(&trig_ptr->program, &str_len, buffer);
	safe_unpack8    (&trig_ptr->state,     buffer);
	if ((trig_ptr->res_type < TRIGGER_RES_TYPE_JOB)
	||  (trig_ptr->res_type > TRIGGER_RES_TYPE_NODE)
	||  (trig_ptr->state > 2))
		goto unpack_error;
	if (trig_ptr->res_type == TRIGGER_RES_TYPE_JOB) {
		trig_ptr->job_id = (uint32_t) atol(trig_ptr->res_id);
		trig_ptr->job_ptr = find_job_record(trig_ptr->job_id);
		if ((trig_ptr->job_id == 0)
		||  (trig_ptr->job_ptr == NULL)
		||  (IS_JOB_FINISHED(trig_ptr->job_ptr)))
			goto unpack_error;
	} else {
		trig_ptr->job_id = 0;
		trig_ptr->job_ptr = NULL;
		if ((trig_ptr->res_id != NULL)
		&&  (trig_ptr->res_id[0] != '*')
		&&  (node_name2bitmap(trig_ptr->res_id, false,
				&trig_ptr->nodes_bitmap) != 0))
			goto unpack_error;
	}

	slurm_mutex_lock(&trigger_mutex);
	if (trigger_list == NULL)
		trigger_list = list_create(_trig_del);
	list_append(trigger_list, trig_ptr);
	slurm_mutex_unlock(&trigger_mutex);

	return SLURM_SUCCESS;

unpack_error:
	error("Incomplete trigger record");
	xfree(trig_ptr->res_id);
	xfree(trig_ptr->program);
	FREE_NULL_BITMAP(trig_ptr->nodes_bitmap);
	xfree(trig_ptr);
	return SLURM_FAILURE;
}
extern int trigger_state_save(void)
{
	static int high_buffer_size = (1024 * 1024);
	int error_code = 0, log_fd;
	char *old_file, *new_file, *reg_file;
	Buf buffer = init_buf(high_buffer_size);
	ListIterator trig_iter;
	trig_mgr_info_t *trig_in;
	/* Locks: Read config */
	slurmctld_lock_t config_read_lock =
		{ READ_LOCK, NO_LOCK, NO_LOCK, NO_LOCK };

	/* write header: version, time */
	packstr(TRIGGER_STATE_VERSION, buffer);
	pack_time(time(NULL), buffer);

	/* write individual trigger records */
	slurm_mutex_lock(&trigger_mutex);
	if (trigger_list == NULL)
		trigger_list = list_create(_trig_del);

	trig_iter = list_iterator_create(trigger_list);
	while ((trig_in = list_next(trig_iter)))
		_dump_trigger_state(trig_in, buffer);
	list_iterator_destroy(trig_iter);
	slurm_mutex_unlock(&trigger_mutex);

	/* write the buffer to file */
	lock_slurmctld(config_read_lock);
	old_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(old_file, "/trigger_state.old");
	reg_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(reg_file, "/trigger_state");
	new_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(new_file, "/trigger_state.new");
	unlock_slurmctld(config_read_lock);

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
	return error_code;
}

extern int trigger_state_restore(void)
{
	int data_allocated, data_read = 0, error_code = 0;
	uint32_t data_size = 0;
	int state_fd, trigger_cnt = 0;
	char *data = NULL, *state_file;
	Buf buffer;
	time_t buf_time;
	char *ver_str = NULL;
	uint16_t ver_str_len;

	/* read the file */
	state_file = xstrdup(slurmctld_conf.state_save_location);
	xstrcat(state_file, "/trigger_state");
	lock_state_files();
	state_fd = open(state_file, O_RDONLY);
	if (state_fd < 0) {
		info("No trigger state file (%s) to recover", state_file);
		error_code = ENOENT;
	} else {
		data_allocated = BUF_SIZE;
		data = xmalloc(data_allocated);
		while (1) {
			data_read = read(state_fd, &data[data_size],
					 BUF_SIZE);
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

	buffer = create_buf(data, data_size);
	if (size_buf(buffer) >= sizeof(uint16_t) + 
			strlen(TRIGGER_STATE_VERSION)) {
	        char *ptr = get_buf_data(buffer);

	        if (!memcmp(&ptr[sizeof(uint16_t)], TRIGGER_STATE_VERSION, 3)) {
		        safe_unpackstr_xmalloc(&ver_str, &ver_str_len, buffer);
		        debug3("Version string in trigger_state header is %s",
			       ver_str);
		}
	}
	if (ver_str && (strcmp(ver_str, TRIGGER_STATE_VERSION) != 0)) {
		error("Can't recover trigger state, data version incompatable");
		xfree(ver_str);
		free_buf(buffer);
		return EFAULT;
	}
	xfree(ver_str);

	safe_unpack_time(&buf_time, buffer);

	while (remaining_buf(buffer) > 0) {
		error_code = _load_trigger_state(buffer);
		if (error_code != SLURM_SUCCESS)
			goto unpack_error;
		trigger_cnt++;
	}
	goto fini;

unpack_error:
	error("Incomplete trigger data checkpoint file");
fini:	verbose("State of %d triggers recovered", trigger_cnt);
	free_buf(buffer);
	return SLURM_FAILURE;
}

/* Test if the event has been triggered, change trigger state as needed */
static void _trigger_job_event(trig_mgr_info_t *trig_in, time_t now)
{
	if ((trig_in->job_ptr == NULL)
	||  (trig_in->job_ptr->job_id == trig_in->job_id))
		trig_in->job_ptr = find_job_record(trig_in->job_ptr->job_id);
	if ((trig_in->trig_type & TRIGGER_TYPE_FINI)
	&&  ((trig_in->job_ptr == NULL) ||
	     (IS_JOB_FINISHED(trig_in->job_ptr)))) {
#if _DEBUG
		info("trigger[%u] event for job %u fini",
			trig_in->trig_id, trig_in->job_id);
#endif
		trig_in->state = 1;
		return;
	}
	if (trig_in->job_ptr == NULL) {
#if _DEBUG
		info("trigger[%u] for defunct job %u",
			trig_in->trig_id, trig_in->job_id);
#endif
		trig_in->state = 2;
		trig_in->trig_time = now;
		return;
	}
	if (trig_in->trig_type & TRIGGER_TYPE_TIME) {
		long rem_time = (trig_in->job_ptr->end_time - now);
		if (rem_time <= (0x8000 - trig_in->trig_time)) {
#if _DEBUG
			info("trigger[%u] for job %u time",
				trig_in->trig_id, trig_in->job_id);
#endif
			trig_in->state = 1;
			return;
		}
	}
	if (trig_in->trig_type & TRIGGER_TYPE_DOWN) {
		if (trigger_down_nodes_bitmap
		&&  bit_overlap(trig_in->job_ptr->node_bitmap, 
				trigger_down_nodes_bitmap)) {
#if _DEBUG
			info("trigger[%u] for job %u down",
				trig_in->trig_id, trig_in->job_id);
#endif
			trig_in->state = 1;
			return;
		}
	}
	if (trig_in->trig_type & TRIGGER_TYPE_UP) {
		if (trigger_down_nodes_bitmap
		&&  bit_overlap(trig_in->job_ptr->node_bitmap, 
				trigger_up_nodes_bitmap)) {
#if _DEBUG
			info("trigger[%u] for job %u up",
				trig_in->trig_id, trig_in->job_id);
#endif
			trig_in->state = 1;
			return;
		}
	}
}

static void _trigger_node_event(trig_mgr_info_t *trig_in, time_t now)
{
	if ((trig_in->trig_type & TRIGGER_TYPE_DOWN)
	&&   trigger_down_nodes_bitmap
	&&   (bit_ffs(trigger_down_nodes_bitmap) != -1)) {
		if (trig_in->nodes_bitmap == NULL) {	/* all nodes */
			xfree(trig_in->res_id);
			trig_in->res_id = bitmap2node_name(
					trigger_down_nodes_bitmap);
			trig_in->state = 1;
		} else if (bit_overlap(trig_in->nodes_bitmap, 
					trigger_down_nodes_bitmap)) {
			bit_and(trig_in->nodes_bitmap, 
					trigger_down_nodes_bitmap);
			xfree(trig_in->res_id);
			trig_in->res_id = bitmap2node_name(
					trig_in->nodes_bitmap);
			trig_in->state = 1;
		}
		if (trig_in->state == 1) {
#if _DEBUG
			info("trigger[%u] for node %s down",
				trig_in->trig_id, trig_in->res_id);
#endif
			return;
		}
	}
	if ((trig_in->trig_type & TRIGGER_TYPE_UP)
	&&   trigger_up_nodes_bitmap
	&&   (bit_ffs(trigger_up_nodes_bitmap) != -1)) {
		if (trig_in->nodes_bitmap == NULL) {	/* all nodes */
			xfree(trig_in->res_id);
			trig_in->res_id = bitmap2node_name(
					trigger_up_nodes_bitmap);
			trig_in->state = 1;
		} else if (bit_overlap(trig_in->nodes_bitmap, 
					trigger_up_nodes_bitmap)) {
			bit_and(trig_in->nodes_bitmap, 
					trigger_up_nodes_bitmap);
			xfree(trig_in->res_id);
			trig_in->res_id = bitmap2node_name(
					trig_in->nodes_bitmap);
			trig_in->state = 1;
		}
		if (trig_in->state == 1) {
#if _DEBUG
			info("trigger[%u] for node %s up",
				trig_in->trig_id, trig_in->res_id);
#endif
			return;
		}
	}
	if ((trig_in->trig_type & TRIGGER_TYPE_RECONFIG)
	&&   trigger_node_reconfig) {
		trig_in->state = 1;
		xfree(trig_in->res_id);
		trig_in->res_id = xstrdup("reconfig");
#if _DEBUG
		info("trigger[%u] for reconfig", trig_in->trig_id);
#endif
		return;
	}
}

static void _trigger_run_program(trig_mgr_info_t *trig_in)
{
	char program[1024], arg0[1024], arg1[1024], *pname;
	uid_t uid;
	gid_t gid;
	pid_t child;

	strncpy(program, trig_in->program, sizeof(program));
	pname = strrchr(program, '/');
	if (pname == NULL)
		pname = program;
	else
		pname++;
	strncpy(arg0, pname, sizeof(arg0));
	strncpy(arg1, trig_in->res_id, sizeof(arg1));
	uid = trig_in->user_id;
	gid = trig_in->group_id;
	child = fork();
	if (child > 0) {
		trig_in->group_id = child;
	} else if (child == 0) {
		int i;
		for (i=0; i<128; i++)
			close(i);
		setpgrp();
		setsid();
		setuid(uid);
		setgid(gid);
		execl(program, arg0, arg1, NULL);
		exit(1);
	} else
		error("fork: %m");
}

static void _clear_event_triggers(void)
{
	if (trigger_down_nodes_bitmap)
		bit_nclear(trigger_down_nodes_bitmap, 0, (node_record_count-1));
	if (trigger_up_nodes_bitmap)
		bit_nclear(trigger_up_nodes_bitmap,   0, (node_record_count-1));
	trigger_node_reconfig = false;
}

extern void trigger_process(void)
{
	ListIterator trig_iter;
	trig_mgr_info_t *trig_in;
	time_t now = time(NULL);
	slurmctld_lock_t job_node_read_lock =
		{ NO_LOCK, READ_LOCK, READ_LOCK, NO_LOCK };

	lock_slurmctld(job_node_read_lock);
	slurm_mutex_lock(&trigger_mutex);
	if (trigger_list == NULL)
		trigger_list = list_create(_trig_del);

	trig_iter = list_iterator_create(trigger_list);
	while ((trig_in = list_next(trig_iter))) {
		if (trig_in->state == 0) {
			if (trig_in->res_type == TRIGGER_RES_TYPE_JOB)
				_trigger_job_event(trig_in, now);
			else
				_trigger_node_event(trig_in, now);
		}
		if (trig_in->state == 1) {
#if _DEBUG
			info("launching program for trigger[%u]",
				trig_in->trig_id);
//			info("  program=%s arg=%s", trig_in->program, 
//				trig_in->res_id);
#endif
			trig_in->state = 2;
			trig_in->trig_time = now;
			_trigger_run_program(trig_in);
		} else if ((trig_in->state == 2) && 
			   (difftime(now, trig_in->trig_time) > 
					MAX_PROG_TIME)) {
#if _DEBUG
			info("purging trigger[%u]", trig_in->trig_id);
#endif
			if (trig_in->group_id != 0)
				kill(-trig_in->group_id, SIGKILL);
			list_delete(trig_iter);
		}
	}
	list_iterator_destroy(trig_iter);
	_clear_event_triggers();
	slurm_mutex_unlock(&trigger_mutex);
	unlock_slurmctld(job_node_read_lock);
}
