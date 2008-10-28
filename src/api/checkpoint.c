/*****************************************************************************\
 *  checkpoint.c - Process checkpoint related functions.
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> et. al.
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

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <slurm/slurm.h>

#include "src/common/checkpoint.h"
#include "src/common/slurm_protocol_api.h"

#ifdef HAVE_AIX
  char *__progname = "PROGRAM";
#else
  extern char * __progname;
#endif

static int _handle_rc_msg(slurm_msg_t *msg);
static int _checkpoint_op (uint16_t op, uint16_t data,
		uint32_t job_id, uint32_t step_id);
/*
 * _checkpoint_op - perform many checkpoint operation for some job step.
 * IN op      - operation to perform
 * IN data    - operation-specific data
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * RET 0 or a slurm error code
 */
static int _checkpoint_op (uint16_t op, uint16_t data,
		uint32_t job_id, uint32_t step_id)
{
	int rc;
	checkpoint_msg_t ckp_req;
	slurm_msg_t req_msg;

	slurm_msg_t_init(&req_msg);
	ckp_req.op       = op;
	ckp_req.data     = data;
	ckp_req.job_id   = job_id;
	ckp_req.step_id  = step_id;
	req_msg.msg_type = REQUEST_CHECKPOINT;
	req_msg.data     = &ckp_req;

	if (slurm_send_recv_controller_rc_msg(&req_msg, &rc) < 0)
		return SLURM_ERROR;

	slurm_seterrno(rc);
	return rc;
}

/*
 * slurm_checkpoint_able - determine if the specified job step can presently
 *	be checkpointed
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * OUT start_time - time at which checkpoint request was issued
 * RET 0 (can be checkpoined) or a slurm error code
 */
extern int slurm_checkpoint_able (uint32_t job_id, uint32_t step_id, 
		time_t *start_time)
{
	int rc;
	slurm_msg_t req_msg, resp_msg;
	checkpoint_msg_t ckp_req;
	checkpoint_resp_msg_t *resp;

	ckp_req.op       = CHECK_ABLE;
	ckp_req.job_id   = job_id;
	ckp_req.step_id  = step_id;
	slurm_msg_t_init(&req_msg);
	slurm_msg_t_init(&resp_msg);
	req_msg.msg_type = REQUEST_CHECKPOINT;
	req_msg.data     = &ckp_req;

	if (slurm_send_recv_controller_msg(&req_msg, &resp_msg) < 0)
		return SLURM_ERROR;

	switch(resp_msg.msg_type) {
	 case RESPONSE_CHECKPOINT:
		resp = (checkpoint_resp_msg_t *) resp_msg.data;
		*start_time = resp->event_time;
		slurm_free_checkpoint_resp_msg(resp_msg.data);
		rc = SLURM_SUCCESS;
		break;
	 case RESPONSE_SLURM_RC:
		rc = _handle_rc_msg(&resp_msg);
		break;
	 default:
		*start_time = (time_t) NULL;
		rc = SLURM_ERROR;
	}
	return rc;
}

/*
 * slurm_checkpoint_disable - disable checkpoint requests for some job step
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_disable (uint32_t job_id, uint32_t step_id)
{
	return _checkpoint_op (CHECK_DISABLE, 0, job_id, step_id);
}


/*
 * slurm_checkpoint_enable - enable checkpoint requests for some job step
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_enable (uint32_t job_id, uint32_t step_id)
{
	return _checkpoint_op (CHECK_ENABLE, 0, job_id, step_id);
}

/*
 * slurm_checkpoint_create - initiate a checkpoint requests for some job step.
 *	the job will continue execution after the checkpoint operation completes
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * IN max_wait - maximum wait for operation to complete, in seconds
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_create (uint32_t job_id, uint32_t step_id, 
		uint16_t max_wait)
{
	return _checkpoint_op (CHECK_CREATE, max_wait, job_id, step_id);
}

/*
 * slurm_checkpoint_vacate - initiate a checkpoint requests for some job step.
 *	the job will terminate after the checkpoint operation completes
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * IN max_wait - maximum wait for operation to complete, in seconds
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_vacate (uint32_t job_id, uint32_t step_id, 
		uint16_t max_wait)
{
	return _checkpoint_op (CHECK_VACATE, max_wait, job_id, step_id);
}

/*
 * slurm_checkpoint_restart - restart execution of a checkpointed job step.
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_restart (uint32_t job_id, uint32_t step_id)
{
	return _checkpoint_op (CHECK_RESTART, 0, job_id, step_id);
}

/*
 * slurm_checkpoint_complete - note the completion of a job step's checkpoint
 *	operation.
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * IN begin_time - time at which checkpoint began
 * IN error_code - error code, highest value for all complete calls is preserved
 * IN error_msg - error message, preserved for highest error_code
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_complete (uint32_t job_id, uint32_t step_id,
		time_t begin_time, uint32_t error_code, char *error_msg)
{
	int rc;
	slurm_msg_t msg;
	checkpoint_comp_msg_t req;

	slurm_msg_t_init(&msg);
	req.job_id       = job_id;
	req.step_id      = step_id;
	req.begin_time   = begin_time;
	req.error_code   = error_code;
	req.error_msg    = error_msg;
	msg.msg_type     = REQUEST_CHECKPOINT_COMP;
	msg.data         = &req;

	if (slurm_send_recv_controller_rc_msg(&msg, &rc) < 0)
		return SLURM_ERROR;
	if (rc)
		slurm_seterrno_ret(rc);
	return SLURM_SUCCESS;
}

/*
 * slurm_checkpoint_error - gather error information for the last checkpoint 
 *	operation for some job step
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * OUT error_code - error number associated with the last checkpoint operation,
 *	this value is dependent upon the checkpoint plugin used and may be
 *	completely unrelated to slurm error codes, the highest value for all
 *	complete calls is preserved
 * OUT error_msg - error message, preserved for highest error_code, value 
 *	must be freed by the caller to prevent memory leak
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_error ( uint32_t job_id, uint32_t step_id, 
		uint32_t *error_code, char **error_msg)
{
	int rc;
	slurm_msg_t msg;
	checkpoint_msg_t req;
	slurm_msg_t resp_msg;
	checkpoint_resp_msg_t *ckpt_resp;

	if ((error_code == NULL) || (error_msg == NULL))
		return EINVAL;

	/*
	 * Request message:
	 */
	req.op       = CHECK_ERROR;
	req.job_id   = job_id;
	req.step_id  = step_id;
	slurm_msg_t_init(&msg);
	slurm_msg_t_init(&resp_msg);
	msg.msg_type = REQUEST_CHECKPOINT;
	msg.data     = &req;

	rc = slurm_send_recv_controller_msg(&msg, &resp_msg);

	if (rc == SLURM_SOCKET_ERROR) 
		return rc;

	switch (resp_msg.msg_type) {
	 case RESPONSE_SLURM_RC:
		*error_code = 0;
		*error_msg = strdup("");
		rc = _handle_rc_msg(&resp_msg);
		break;
	 case RESPONSE_CHECKPOINT:
		ckpt_resp = (checkpoint_resp_msg_t *) resp_msg.data;
		*error_code = ckpt_resp->error_code;
		if (ckpt_resp->error_msg)
			*error_msg = strdup(ckpt_resp->error_msg);
		else
			*error_msg = strdup("");
		slurm_free_checkpoint_resp_msg(ckpt_resp);
		rc = SLURM_SUCCESS;
		break;
	 default:
		rc = SLURM_UNEXPECTED_MSG_ERROR;
	}

	return rc;
}

/*
 *  Handle a return code message type. 
 *    Sets errno to return code and returns it
 */
static int
_handle_rc_msg(slurm_msg_t *msg)
{
	int rc = ((return_code_msg_t *) msg->data)->return_code;
	slurm_free_return_code_msg(msg->data);
	slurm_seterrno(rc);
	return rc;
}

/*
 * slurm_checkpoint_task_complete - note the completion of a task's checkpoint
 *	operation.
 * IN job_id  - job on which to perform operation
 * IN step_id - job step on which to perform operation
 * IN task_id - task which completed the operation
 * IN begin_time - time at which checkpoint began
 * IN error_code - error code, highest value for all complete calls is preserved
 * IN error_msg - error message, preserved for highest error_code
 * RET 0 or a slurm error code
 */
extern int slurm_checkpoint_task_complete (uint32_t job_id, uint32_t step_id,
		uint32_t task_id, time_t begin_time, uint32_t error_code, char *error_msg)
{
	int rc;
	slurm_msg_t msg;
	checkpoint_task_comp_msg_t req;

	slurm_msg_t_init(&msg);
	req.job_id       = job_id;
	req.step_id      = step_id;
	req.task_id      = task_id;
	req.begin_time   = begin_time;
	req.error_code   = error_code;
	req.error_msg    = error_msg;
	msg.msg_type     = REQUEST_CHECKPOINT_TASK_COMP;
	msg.data         = &req;

	if (slurm_send_recv_controller_rc_msg(&msg, &rc) < 0)
		return SLURM_ERROR;
	if (rc)
		slurm_seterrno_ret(rc);
	return SLURM_SUCCESS;
}

/*
 * slurm_get_checkpoint_file_path - return the checkpoint file
 *      path of this process, creating the directory if needed.
 * IN len: length of the file path buffer
 * OUT buf: buffer to store the checkpoint file path
 * RET: 0 on success, -1 on failure with errno set
 */
extern int
slurm_get_checkpoint_file_path(size_t len, char *buf)
{
       char *ckpt_path, *job_id, *step_id, *proc_id;
       struct stat mystat;
       int idx;

       len --;                 /* for a terminating 0 */

       ckpt_path = getenv("SLURM_CHECKPOINT_PATH");
       if (ckpt_path == NULL) { /* this should not happen since the program may chdir */
               ckpt_path = getcwd(buf, len);
               if (ckpt_path == NULL)  /* ERANGE: len is too short */
                       return -1;
       } else {
               if (snprintf(buf, len, "%s", ckpt_path) >= len) { /* glibc >= 2.1 */
                       errno = ERANGE;
                       return -1;
               }
               ckpt_path = buf;
       }
       idx = strlen(ckpt_path) - 1;
       while (idx > 0 && ckpt_path[idx] == '/')
               ckpt_path[idx --] = 0;

       if (stat(ckpt_path, &mystat) < 0)
               return -1;
       if (! S_ISDIR(mystat.st_mode)) {
               errno = ENOTDIR;
               return -1;
       }

       job_id = getenv("SLURM_JOBID");
       step_id = getenv("SLURM_STEPID");
       proc_id = getenv("SLURM_PROCID");
       if (job_id == NULL || step_id == NULL || proc_id == NULL) {
               errno = ENODATA;
               return -1;
       }
       idx = strlen(buf);
       if (snprintf(buf + idx, len - idx, "/%s.%s", job_id, step_id) >= len - idx) {
               errno = ERANGE;
               return -1;
       }

       if (stat(buf, &mystat) < 0) {
               if (errno == ENOENT) { /* dir does not exists */
                       if (mkdir(buf, 0750) < 0 && errno != EEXIST)
                               return -1;
                       if (stat(buf, &mystat) < 0)
                               return -1;
               }
               else
                       return -1;
       }
       if (! S_ISDIR(mystat.st_mode)) {
               errno = ENOTDIR;
               return -1;
       }

       idx = strlen(buf);
       if (snprintf(buf + idx, len - idx, "/%s.%s.ckpt", __progname, proc_id) >= len - idx) {
               errno = ERANGE;
               return -1;
       }

       return 0;
}
