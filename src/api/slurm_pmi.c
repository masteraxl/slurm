/****************************************************************************\
 *  slurm_pmi.c - PMI support functions internal to SLURM
 *****************************************************************************
 *  Copyright (C) 2005-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>.
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
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#include <stdlib.h>
#include <slurm/slurm.h>
#include <slurm/slurm_errno.h>

#include "src/api/slurm_pmi.h"
#include "src/common/slurm_protocol_defs.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/forward.h"
#include "src/common/read_config.h"
#include "src/common/xmalloc.h"
#include "src/common/fd.h"

#define MAX_RETRIES 5

int pmi_fd = -1;
uint16_t srun_port = 0;
slurm_addr srun_addr;

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
		int pmi_rank)
{
	slurm_msg_t msg_send;
	int rc, retries = 0;

	if (kvs_set_ptr == NULL)
		return EINVAL;

	if ((rc = _get_addr()) != SLURM_SUCCESS)
		return rc; 

	msg_send.address = srun_addr;
	msg_send.msg_type = PMI_KVS_PUT_REQ;
	msg_send.data = (void *) kvs_set_ptr;
	forward_init(&msg_send.forward, NULL);
	msg_send.ret_list = NULL;
	msg_send.forward_struct_init = 0;
	
	/* Send the RPC to the local srun communcation manager.
	 * Since the srun can be sent thousands of messages at 
	 * the same time and refuse some connections, retry as 
	 * needed. Spread out messages by task's rank.*/
	while (slurm_send_recv_rc_msg_only_one(&msg_send, &rc, 0) < 0) {
		if (retries++ > MAX_RETRIES) {
			error("slurm_get_kvs_comm_set: %m");
			return SLURM_ERROR;
		}
		usleep(pmi_rank * 1000);
	}

	return rc;
}

/* Wait for barrier and get full PMI Keyval space data */
int  slurm_get_kvs_comm_set(struct kvs_comm_set **kvs_set_ptr, 
		int pmi_rank, int pmi_size)
{
	int rc, srun_fd, retries = 0;
	slurm_msg_t msg_send, msg_rcv;
	slurm_addr slurm_addr, srun_reply_addr;
	char hostname[64];
	uint16_t port;
	kvs_get_msg_t data;
	char *env_pmi_ifhn;
	List ret_list = NULL;

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
		getnodename(hostname, sizeof(hostname));

	data.task_id = pmi_rank;
	data.size = pmi_size;
	data.port = port;
	data.hostname = hostname;
	msg_send.address = srun_addr;
	msg_send.msg_type = PMI_KVS_GET_REQ;
	msg_send.data = &data;

	/* Send the RPC to the srun communcation manager */
	forward_init(&msg_send.forward, NULL);
	msg_send.ret_list = NULL;
	msg_send.forward_struct_init = 0;

	/* Send the RPC to the local srun communcation manager.
	 * Since the srun can be sent thousands of messages at 
	 * the same time and refuse some connections, retry as 
	 * needed. Spread out messages by task's rank. */
	while (slurm_send_recv_rc_msg_only_one(&msg_send, &rc, 0) < 0) {
		if (retries++ > MAX_RETRIES) {
			error("slurm_get_kvs_comm_set: %m");
			return SLURM_ERROR;
		}
		usleep(pmi_rank * 1000);
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

	while ((ret_list = slurm_receive_msg(srun_fd, &msg_rcv, 0)) == NULL) {
		if (errno == EINTR)
			continue;
		error("slurm_receive_msg: %m");
		slurm_close_accepted_conn(srun_fd);
		return errno;
	}
	if(ret_list) {
		if(list_count(ret_list)>0) {
			error("We didn't do things correctly "
			      "got %d responses didn't expect any",
			      list_count(ret_list));
		}
		list_destroy(ret_list);
	}
	msg_rcv.conn_fd = srun_fd;
	if (msg_rcv.msg_type != PMI_KVS_GET_RESP) {
		error("slurm_get_kvs_comm_set msg_type=%d", msg_rcv.msg_type);
		slurm_close_accepted_conn(srun_fd);
		return SLURM_UNEXPECTED_MSG_ERROR;
	}
	if (slurm_send_rc_msg(&msg_rcv, SLURM_SUCCESS) < 0)
		error("slurm_send_rc_msg: %m");
	slurm_close_accepted_conn(srun_fd);
	*kvs_set_ptr = msg_rcv.data;
	return SLURM_SUCCESS;
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

	for (i=0; i<kvs_set_ptr->kvs_comm_recs; i++)
		_free_kvs_comm(kvs_set_ptr->kvs_comm_ptr[i]);
	xfree(kvs_set_ptr->kvs_comm_ptr);
	xfree(kvs_set_ptr);
}

