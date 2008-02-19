/****************************************************************************\
 *  slurmdbd_defs.c - functions for use with Slurm DBD RPCs
 *****************************************************************************
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#if HAVE_CONFIG_H 
#  include "config.h"
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  else
#    if HAVE_STDINT_H
#      include <stdint.h>
#    endif
#  endif			/* HAVE_INTTYPES_H */
#else				/* !HAVE_CONFIG_H */
#  include <inttypes.h>
#endif				/*  HAVE_CONFIG_H */

#include <arpa/inet.h>
#include <pthread.h>
#include <sys/poll.h>

#include "slurm/slurm_errno.h"
#include "src/common/fd.h"
#include "src/common/pack.h"
#include "src/common/slurmdbd_defs.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"

static slurm_fd slurmdbd_fd   = -1;
static pthread_mutex_t slurmdbd_lock = PTHREAD_MUTEX_INITIALIZER;


static void   _close_slurmdbd_fd(void);
static bool   _fd_readable(slurm_fd fd);
static int    _fd_writeable(slurm_fd fd);
static int    _get_return_code(void);
static void   _open_slurmdbd_fd(void);
static Buf    _recv_msg(void);
static void   _reopen_slurmdbd_fd(void);
static int    _send_init_msg(void);
static int    _send_msg(Buf buffer);
static int    _tot_wait (struct timeval *start_time);

/****************************************************************************
 * Socket open/close/read/write functions
 ****************************************************************************/
extern int slurm_open_slurmdbd_conn(void)
{
	slurm_mutex_lock(&slurmdbd_lock);
	if (slurmdbd_fd < 0)
		_open_slurmdbd_fd();
	slurm_mutex_unlock(&slurmdbd_lock);

	return SLURM_SUCCESS;
}

extern int slurm_close_slurmdbd_conn(void)
{
	slurm_mutex_lock(&slurmdbd_lock);
	_close_slurmdbd_fd();
	slurm_mutex_unlock(&slurmdbd_lock);
	return SLURM_SUCCESS;
}

extern int slurm_send_recv_slurmdbd_rc_msg(slurmdbd_msg_t *req, int *resp_code)
{
	int rc;
	Buf buffer;

	xassert(resp_code);
	*resp_code = SLURM_ERROR;

	if (slurmdbd_fd < 0) {
		/* slurm_open_slurmdbd_conn() should have been called first,
		 * but we'll be accomodating and open the connection here
		 * too. However this will slow the RPC down. */
		slurm_open_slurmdbd_conn();
		if (slurmdbd_fd < 0)
			return SLURM_ERROR;
	}

	buffer = init_buf(1024);
	pack16(req->msg_type, buffer);
	switch (req->msg_type) {
		case DBD_INIT:
			slurm_dbd_pack_init_msg(
				(dbd_init_msg_t *) req->data, buffer);
			break;
		case DBD_JOB_COMPLETE:
			slurm_dbd_pack_job_complete_msg(
				(dbd_job_comp_msg_t *) req->data, buffer);
			break;
		case DBD_JOB_START:
			slurm_dbd_pack_job_start_msg(
				(dbd_job_start_msg_t *) req->data, buffer);
			break;
		case DBD_JOB_SUSPEND:
			slurm_dbd_pack_job_suspend_msg(
				(dbd_job_suspend_msg_t *) req->data, buffer);
			break;
		case DBD_STEP_COMPLETE:
			slurm_dbd_pack_step_complete_msg(
				(dbd_step_comp_msg_t *) req->data, buffer);
			break;
		case DBD_STEP_START:
			slurm_dbd_pack_step_start_msg(
				(dbd_step_start_msg_t *) req->data, buffer);
			break;
		default:
			error("slurmdbd: Invalid message type %u",
			      req->msg_type);
			free_buf(buffer);
			return SLURM_ERROR;
	}

	rc = _send_msg(buffer);
	free_buf(buffer);
	if (rc != SLURM_SUCCESS) {
		error("slurmdbd: Sending message type %u", req->msg_type);
		return SLURM_ERROR;
	}

	*resp_code = _get_return_code();
	return SLURM_SUCCESS;
}

/* Open a connection to the Slurm DBD and set slurmdbd_fd */
static void _open_slurmdbd_fd(void)
{
	slurm_addr dbd_addr;
	uint16_t slurmdbd_port;
	char *   slurmdbd_addr;

	if (slurmdbd_fd >= 0) {
		debug("Attempt to re-open slurmdbd socket");
		return;
	}

	slurmdbd_addr = slurm_get_slurmdbd_addr();
	slurmdbd_port = slurm_get_slurmdbd_port();
	if ((slurmdbd_addr == NULL) || (slurmdbd_port == 0)) {
		error("Invalid SlurmDbd address %s:%u",
			slurmdbd_addr, slurmdbd_port);
		xfree(slurmdbd_addr);
		return;
	}

	slurm_set_addr(&dbd_addr, slurmdbd_port, slurmdbd_addr);
	if (dbd_addr.sin_port == 0)
		error("Unable to locate SlurmDBD host %s:%s", 
		      slurmdbd_addr, slurmdbd_addr);
	else {
		slurmdbd_fd = slurm_open_msg_conn(&dbd_addr);
		if (slurmdbd_fd < 0)
			error("slurmdbd: slurm_open_msg_conn: %m");
		else {
			fd_set_nonblocking(slurmdbd_fd);
			if (_send_init_msg() != SLURM_SUCCESS)
				error("slurmdbd: Sending DdbInit msg: %m");
			else
				debug("slurmdbd: Sent DbdInit msg");
		}
	}
	xfree(slurmdbd_addr);
}

static int _send_init_msg(void)
{
	int rc;
	slurmdbd_msg_t msg;
	dbd_init_msg_t req;

	req.version  = SLURM_DBD_VERSION;
	msg.msg_type = DBD_INIT;
	msg.data = &req;

	if (slurm_send_recv_slurmdbd_rc_msg(&msg, &rc) < 0)
		return SLURM_ERROR;

	if (rc)
		slurm_seterrno_ret(rc);

	return SLURM_SUCCESS;
}

/* Close the SlurmDbd connection */
static void _close_slurmdbd_fd(void)
{
	if (slurmdbd_fd >= 0) {
		close(slurmdbd_fd);
		slurmdbd_fd = -1;
	}
}

/* Reopen the Slurm DBD connection due to some error */
static void _reopen_slurmdbd_fd(void)
{
	info("slurmdbd: reopening connection");
	_close_slurmdbd_fd();
	_open_slurmdbd_fd();
}

static int _send_msg(Buf buffer)
{
	uint32_t msg_size, nw_size;
	char *msg;
	ssize_t msg_wrote;
	int rc;

	if (slurmdbd_fd < 0)
		return SLURM_ERROR;

	rc =_fd_writeable(slurmdbd_fd);
	if (rc == -1) {
re_open:	/* SlurmDBD shutdown, try to reopen a connection now */
		_reopen_slurmdbd_fd();
		rc = _fd_writeable(slurmdbd_fd);
	}
	if (rc < 1)
		return SLURM_ERROR;

	msg_size = get_buf_offset(buffer);
	nw_size = htonl(msg_size);
	msg_wrote = write(slurmdbd_fd, &nw_size, sizeof(nw_size));
	if (msg_wrote != sizeof(nw_size))
		return SLURM_ERROR;

	msg = get_buf_data(buffer);
	while (msg_size > 0) {
		rc = _fd_writeable(slurmdbd_fd);
		if (rc == -1)
			goto re_open;
		if (rc < 1)
			return SLURM_ERROR;
		msg_wrote = write(slurmdbd_fd, msg, msg_size);
		if (msg_wrote <= 0)
			return SLURM_ERROR;
		msg += msg_wrote;
		msg_size -= msg_wrote;
	}

	return SLURM_SUCCESS;
}

static int _get_return_code(void)
{
	Buf buffer;
	uint16_t msg_type;
	dbd_rc_msg_t *msg;
	int rc = SLURM_ERROR;

	buffer = _recv_msg();
	if (buffer == NULL)
		return rc;

	safe_unpack16(&msg_type, buffer);
	if (msg_type != DBD_RC)
		error("slurmdbd: bad message type %d != DBD_RC", msg_type);
	else if (slurm_dbd_unpack_rc_msg(&msg, buffer) == SLURM_SUCCESS) {
		rc = msg->return_code;
		slurm_dbd_free_rc_msg(msg);
		if (rc != SLURM_SUCCESS)
			error("slurmdbd: DBD_RC is %d", rc);
	} else
		error("slurmdbd: unpack message error");

 unpack_error:
	free_buf(buffer);
	return rc;
}

static Buf _recv_msg(void)
{
	uint32_t msg_size, nw_size;
	char *msg;
	ssize_t msg_read, offset;
	Buf buffer;

	if (slurmdbd_fd < 0)
		return NULL;

	if (!_fd_readable(slurmdbd_fd))
		return NULL;
	msg_read = read(slurmdbd_fd, &nw_size, sizeof(nw_size));
	if (msg_read != sizeof(nw_size))
		return NULL;
	msg_size = ntohl(nw_size);
	if ((msg_size < 2) || (msg_size > 1000000)) {
		error("slurmdbd: Invalid msg_size (%u)");
		return NULL;
	}

	msg = xmalloc(msg_size);
	offset = 0;
	while (msg_size > offset) {
		if (!_fd_readable(slurmdbd_fd))
			break;		/* problem with this socket */
		msg_read = read(slurmdbd_fd, (msg + offset), 
				(msg_size - offset));
		if (msg_read <= 0) {
			error("slurmdbd: read: %m");
			break;
		}
		offset += msg_read;
	}
	if (msg_size != offset) {
		error("slurmdbd: only read %d of %d bytes", offset, msg_size);
		xfree(msg);
		return NULL;
	}

	buffer = create_buf(msg, msg_size);
	if (buffer == NULL)
		fatal("create_buf: malloc failure");
	return buffer;
}

/* Return time in msec since "start time" */
static int _tot_wait (struct timeval *start_time)
{
	struct timeval end_time;
	int msec_delay;

	gettimeofday(&end_time, NULL);
	msec_delay =   (end_time.tv_sec  - start_time->tv_sec ) * 1000;
	msec_delay += ((end_time.tv_usec - start_time->tv_usec + 500) / 1000);
	return msec_delay;
}

/* Wait until a file is readable, 
 * RET false if can not be read */
static bool _fd_readable(slurm_fd fd)
{
	struct pollfd ufds;
	static int msg_timeout = -1;
	int rc, time_left;
	struct timeval tstart;

	if (msg_timeout == -1)
		msg_timeout = slurm_get_msg_timeout() * 1000;

	ufds.fd     = fd;
	ufds.events = POLLIN;
	gettimeofday(&tstart, NULL);
	while (1) {
		time_left = msg_timeout - _tot_wait(&tstart);
		rc = poll(&ufds, 1, time_left);
		if ((rc == 0) && ((errno == EINTR) || (errno == EAGAIN)))
			continue;
		if (ufds.revents & POLLHUP) {
			debug2("SlurmDBD connection closed");
			return false;
		}
		if (ufds.revents & POLLNVAL) {
			error("SlurmDBD connection is invalid");
			return false;
		}
		if (ufds.revents & POLLERR) {
			error("SlurmDBD connection experienced an error");
			return false;
		}
		if ((ufds.revents & POLLIN) == 0) {
			error("SlurmDBD connection %d events %d", 
				fd, ufds.revents);
			return false;
		}
		break;
	}
	return true;
}

/* Wait until a file is writable, 
 * RET 1 if file can be written now,
 *     0 if can not be written to within 5 seconds
 *     -1 if file has been closed POLLHUP
 */
static int _fd_writeable(slurm_fd fd)
{
	struct pollfd ufds;
	int msg_timeout = 5000;
	int rc, time_left;
	struct timeval tstart;

	ufds.fd     = fd;
	ufds.events = POLLOUT;
	gettimeofday(&tstart, NULL);
	while (1) {
		time_left = msg_timeout - _tot_wait(&tstart);
		rc = poll(&ufds, 1, time_left);
		if ((rc == 0) && ((errno == EINTR) || (errno == EAGAIN)))
			continue;
		if (ufds.revents & POLLHUP) {
			debug2("SlurmDBD connection is closed");
			return -1;
		}
		if (ufds.revents & POLLNVAL) {
			error("SlurmDBD connection is invalid");
			return 0;
		}
		if (ufds.revents & POLLERR) {
			error("SlurmDBD connection experienced an error: %m");
			return 0;
		}
		if ((ufds.revents & POLLOUT) == 0) {
			error("SlurmDBD connection %d events %d", 
				fd, ufds.revents);
			return 0;
		}
		break;
	}
	return 1;
}

/****************************************************************************
 * Free data structures
 ****************************************************************************/
void inline slurm_dbd_free_get_jobs_msg(dbd_get_jobs_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_init_msg(dbd_init_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_job_complete_msg(dbd_job_comp_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_job_start_msg(dbd_job_start_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_job_submit_msg(dbd_job_submit_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_job_suspend_msg(dbd_job_suspend_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_rc_msg(dbd_rc_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_step_complete_msg(dbd_step_comp_msg_t *msg)
{
	xfree(msg);
}

void inline slurm_dbd_free_step_start_msg(dbd_step_start_msg_t *msg)
{
	xfree(msg);
}

/****************************************************************************
 * Pack and unpack data structures
 ****************************************************************************/
void inline 
slurm_dbd_pack_get_jobs_msg(dbd_get_jobs_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
}

int inline 
slurm_dbd_unpack_get_jobs_msg(dbd_get_jobs_msg_t **msg, Buf buffer)
{
	dbd_get_jobs_msg_t *msg_ptr = xmalloc(sizeof(dbd_get_jobs_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_init_msg(dbd_init_msg_t *msg, Buf buffer)
{
	pack16(msg->version, buffer);
}

int inline 
slurm_dbd_unpack_init_msg(dbd_init_msg_t **msg, Buf buffer)
{
	dbd_init_msg_t *msg_ptr = xmalloc(sizeof(dbd_init_msg_t));
	*msg = msg_ptr;
	safe_unpack16(&msg_ptr->version, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_job_complete_msg(dbd_job_comp_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
}

int inline 
slurm_dbd_unpack_job_complete_msg(dbd_job_comp_msg_t **msg, Buf buffer)
{
	dbd_job_comp_msg_t *msg_ptr = xmalloc(sizeof(dbd_job_comp_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_job_start_msg(dbd_job_start_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
}

int inline 
slurm_dbd_unpack_job_start_msg(dbd_job_start_msg_t **msg, Buf buffer)
{
	dbd_job_start_msg_t *msg_ptr = xmalloc(sizeof(dbd_job_start_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_job_submit_msg(dbd_job_submit_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
}

int inline 
slurm_dbd_unpack_job_submit_msg(dbd_job_submit_msg_t **msg, Buf buffer)
{
	dbd_job_submit_msg_t *msg_ptr = xmalloc(sizeof(dbd_job_submit_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_job_suspend_msg(dbd_job_suspend_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
}

int inline 
slurm_dbd_unpack_job_suspend_msg(dbd_job_suspend_msg_t **msg, Buf buffer)
{
	dbd_job_suspend_msg_t *msg_ptr = xmalloc(sizeof(dbd_job_suspend_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_rc_msg(dbd_rc_msg_t *msg, Buf buffer)
{
	pack32(msg->return_code, buffer);
}

int inline 
slurm_dbd_unpack_rc_msg(dbd_rc_msg_t **msg, Buf buffer)
{
	dbd_rc_msg_t *msg_ptr = xmalloc(sizeof(dbd_rc_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->return_code, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_step_complete_msg(dbd_step_comp_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
	pack32(msg->step_id, buffer);
}

int inline 
slurm_dbd_unpack_step_complete_msg(dbd_step_comp_msg_t **msg, Buf buffer)
{
	dbd_step_comp_msg_t *msg_ptr = xmalloc(sizeof(dbd_step_comp_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	safe_unpack32(&msg_ptr->step_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}

void inline 
slurm_dbd_pack_step_start_msg(dbd_step_start_msg_t *msg, Buf buffer)
{
	pack32(msg->job_id, buffer);
	pack32(msg->step_id, buffer);
}

int inline 
slurm_dbd_unpack_step_start_msg(dbd_step_start_msg_t **msg, Buf buffer)
{
	dbd_step_start_msg_t *msg_ptr = xmalloc(sizeof(dbd_step_start_msg_t));
	*msg = msg_ptr;
	safe_unpack32(&msg_ptr->job_id, buffer);
	safe_unpack32(&msg_ptr->step_id, buffer);
	return SLURM_SUCCESS;

unpack_error:
	xfree(msg_ptr);
	*msg = NULL;
	return SLURM_ERROR;
}
