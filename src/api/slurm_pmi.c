/****************************************************************************\
 *  slurm_pmi.c - PMI support functions internal to SLURM
 *****************************************************************************
 *  Copyright (C) 2005-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>.
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

#include <stdlib.h>
#include <slurm/slurm.h>
#include <slurm/slurm_errno.h>

#include "src/api/slurm_pmi.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_protocol_defs.h"
#include "src/common/forward.h"
#include "src/common/read_config.h"
#include "src/common/xmalloc.h"
#include "src/common/fd.h"
#include "src/common/slurm_auth.h"

#define MAX_RETRIES 5
#define PMI_TIME    500	/* spacing between RPCs, usec */

int pmi_fd = -1;
uint16_t srun_port = 0;
slurm_addr srun_addr;

static int _forward_comm_set(struct kvs_comm_set *kvs_set_ptr);
static int _get_addr(void);

static int _get_addr(void)
{
	char *env_host, *env_port;

	if (srun_port)
		return SLURM_SUCCESS;

	env_host = getenv("SLURM_SRUN_COMM_HOST");
	env_port = getenv("SLURM_SRUN_COMM_PORT");
	if (!env_host || !env_port)
		return SLURM_ERROR;

	srun_port = (uint16_t) atol(env_port);
	slurm_set_addr(&srun_addr, srun_port, env_host);
	return SLURM_SUCCESS;
}

/* Transmit PMI Keyval space data */
int slurm_send_kvs_comm_set(struct kvs_comm_set *kvs_set_ptr, 
		int pmi_rank, int pmi_size)
{
	slurm_msg_t msg_send;
	int rc, retries = 0, timeout = 0;

	if (kvs_set_ptr == NULL)
		return EINVAL;

	if ((rc = _get_addr()) != SLURM_SUCCESS)
		return rc; 

	slurm_msg_t_init(&msg_send);
	msg_send.address = srun_addr;
	msg_send.msg_type = PMI_KVS_PUT_REQ;
	msg_send.data = (void *) kvs_set_ptr;
	
	/* Send the RPC to the local srun communcation manager.
	 * Since the srun can be sent thousands of messages at 
	 * the same time and refuse some connections, retry as 
	 * needed. Spread out messages by task's rank. Also 
	 * increase the timeout if many tasks since the srun 
	 * command is very overloaded.
	 * We also increase the timeout (default timeout is
	 * 10 secs). */
	usleep(pmi_rank * PMI_TIME);
	if      (pmi_size > 1000)	/* 100 secs */
		timeout = slurm_get_msg_timeout() * 10000;
	else if (pmi_size > 100)	/* 50 secs */
		timeout = slurm_get_msg_timeout() * 5000;
	else if (pmi_size > 10)		/* 20 secs */
		timeout = slurm_get_msg_timeout() * 2000;

	while (slurm_send_recv_rc_msg_only_one(&msg_send, &rc, timeout) < 0) {
		if (retries++ > MAX_RETRIES) {
			error("slurm_send_kvs_comm_set: %m");
			return SLURM_ERROR;
		}
		usleep(pmi_rank * PMI_TIME);
	}

	return rc;
}

/* Wait for barrier and get full PMI Keyval space data */
int  slurm_get_kvs_comm_set(struct kvs_comm_set **kvs_set_ptr, 
		int pmi_rank, int pmi_size)
{
	int rc, srun_fd, retries = 0, timeout = 0;
	slurm_msg_t msg_send, msg_rcv;
	slurm_addr slurm_addr, srun_reply_addr;
	char hostname[64];
	uint16_t port;
	kvs_get_msg_t data;
	char *env_pmi_ifhn;
	
	if (kvs_set_ptr == NULL)
		return EINVAL;

	if ((rc = _get_addr()) != SLURM_SUCCESS) {
		error("_get_addr: %m");
		return rc;
	}
	if (pmi_fd < 0) {
		if ((pmi_fd = slurm_init_msg_engine_port(0)) < 0) {
			error("slurm_init_msg_engine_port: %m");
			return SLURM_ERROR;
		}
		fd_set_blocking(pmi_fd);
	}
	if (slurm_get_stream_addr(pmi_fd, &slurm_addr) < 0) {
		error("slurm_get_stream_addr: %m");
		return SLURM_ERROR;
	}
	/* hostname is not set here, so slurm_get_addr fails
	slurm_get_addr(&slurm_addr, &port, hostname, sizeof(hostname)); */
	port = ntohs(slurm_addr.sin_port); 
	if ((env_pmi_ifhn = getenv("SLURM_PMI_RESP_IFHN"))) {
		strncpy(hostname, env_pmi_ifhn, sizeof(hostname));
		hostname[sizeof(hostname)-1] = 0;
	} else
		gethostname_short(hostname, sizeof(hostname));

	data.task_id = pmi_rank;
	data.size = pmi_size;
	data.port = port;
	data.hostname = hostname;
	slurm_msg_t_init(&msg_send);
	slurm_msg_t_init(&msg_rcv);
	msg_send.address = srun_addr;
	msg_send.msg_type = PMI_KVS_GET_REQ;
	msg_send.data = &data;

	/* Send the RPC to the local srun communcation manager.
	 * Since the srun can be sent thousands of messages at 
	 * the same time and refuse some connections, retry as 
	 * needed. Spread out messages by task's rank. Also
	 * increase the timeout if many tasks since the srun
	 * command is very overloaded.
	 * We also increase the timeout (default timeout is
	 * 10 secs). */
	usleep(pmi_rank * PMI_TIME);
	if      (pmi_size > 1000)	/* 100 secs */
		timeout = slurm_get_msg_timeout() * 10000;
	else if (pmi_size > 100)	/* 50 secs */
		timeout = slurm_get_msg_timeout() * 5000;
	else if (pmi_size > 10)		/* 20 secs */
		timeout = slurm_get_msg_timeout() * 2000;

	while (slurm_send_recv_rc_msg_only_one(&msg_send, &rc, timeout) < 0) {
		if (retries++ > MAX_RETRIES) {
			error("slurm_get_kvs_comm_set: %m");
			return SLURM_ERROR;
		}
		usleep(pmi_rank * PMI_TIME);
	}
	if (rc != SLURM_SUCCESS) {
		error("slurm_get_kvs_comm_set error_code=%d", rc);
		return rc;
	}

	/* get the message after all tasks reach the barrier */
	srun_fd = slurm_accept_msg_conn(pmi_fd, &srun_reply_addr);
	if (srun_fd < 0) {
		error("slurm_accept_msg_conn: %m");
		return errno;
	}

	while ((rc = slurm_receive_msg(srun_fd, &msg_rcv, 0)) != 0) {
		if (errno == EINTR)
			continue;
		error("slurm_receive_msg: %m");
		slurm_close_accepted_conn(srun_fd);
		return errno;
	}
	if(msg_rcv.auth_cred)
		(void)g_slurm_auth_destroy(msg_rcv.auth_cred);
	
	if (msg_rcv.msg_type != PMI_KVS_GET_RESP) {
		error("slurm_get_kvs_comm_set msg_type=%d", msg_rcv.msg_type);
		slurm_close_accepted_conn(srun_fd);
		return SLURM_UNEXPECTED_MSG_ERROR;
	}
	if (slurm_send_rc_msg(&msg_rcv, SLURM_SUCCESS) < 0)
		error("slurm_send_rc_msg: %m");
	
	slurm_close_accepted_conn(srun_fd);
	*kvs_set_ptr = msg_rcv.data;

	rc = _forward_comm_set(*kvs_set_ptr);
	return rc;
}

/* Forward keypair info to other tasks as required.
* Clear message forward structure upon completion. */
static int _forward_comm_set(struct kvs_comm_set *kvs_set_ptr)
{
	int i, rc = SLURM_SUCCESS;
	int tmp_host_cnt = kvs_set_ptr->host_cnt;
	slurm_msg_t msg_send;
	int msg_rc;

	kvs_set_ptr->host_cnt = 0;
	for (i=0; i<tmp_host_cnt; i++) {
		slurm_msg_t_init(&msg_send);
		msg_send.msg_type = PMI_KVS_GET_RESP;
		msg_send.data = (void *) kvs_set_ptr;
		slurm_set_addr(&msg_send.address,
			kvs_set_ptr->kvs_host_ptr[i].port,
			kvs_set_ptr->kvs_host_ptr[i].hostname);
		if (slurm_send_recv_rc_msg_only_one(&msg_send,
				&msg_rc, 0) < 0) {
			error("Could not forward msg to %s",
				kvs_set_ptr->kvs_host_ptr[i].hostname);
			msg_rc = 1;
		}
		rc = MAX(rc, msg_rc);
		xfree(kvs_set_ptr->kvs_host_ptr[i].hostname);
	}
	xfree(kvs_set_ptr->kvs_host_ptr);
	return rc;
}

static void _free_kvs_comm(struct kvs_comm *kvs_comm_ptr)
{
	int i;

	if (kvs_comm_ptr == NULL)
		return;

	for (i=0; i<kvs_comm_ptr->kvs_cnt; i++) {
		xfree(kvs_comm_ptr->kvs_keys[i]);
		xfree(kvs_comm_ptr->kvs_values[i]);
	}
	xfree(kvs_comm_ptr->kvs_name);
	xfree(kvs_comm_ptr->kvs_keys);
	xfree(kvs_comm_ptr->kvs_values);
	xfree(kvs_comm_ptr);
}

/* Free kvs_comm_set returned by slurm_get_kvs_comm_set() */
void slurm_free_kvs_comm_set(struct kvs_comm_set *kvs_set_ptr)
{
	int i;

	if (kvs_set_ptr == NULL)
		return;

	for (i=0; i<kvs_set_ptr->host_cnt; i++)
		xfree(kvs_set_ptr->kvs_host_ptr[i].hostname);
	xfree(kvs_set_ptr->kvs_host_ptr);

	for (i=0; i<kvs_set_ptr->kvs_comm_recs; i++)
		_free_kvs_comm(kvs_set_ptr->kvs_comm_ptr[i]);
	xfree(kvs_set_ptr->kvs_comm_ptr);
	xfree(kvs_set_ptr);
}

/* Finalization processing */
void slurm_pmi_finalize(void)
{
	if (pmi_fd >= 0) {
		slurm_shutdown_msg_engine(pmi_fd);
		pmi_fd = -1;
	}
	srun_port = 0;
}
