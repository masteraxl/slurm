/*****************************************************************************\
 *  node_select_info.c - get the node select plugin state information of slurm
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <syslog.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

#include <slurm/slurm.h>

#include "src/api/node_select_info.h"
#include "src/common/slurm_protocol_api.h"

/*
 * slurm_load_node_select - issue RPC to get slurm all node select plugin 
 *	information if changed since update_time 
 * IN update_time - time of current configuration data
 * IN node_select_info_msg_pptr - place to store a node select configuration 
 *	pointer
 * RET 0 or a slurm error code
 * NOTE: free the response using slurm_free_node_select_info_msg
 */
extern int slurm_load_node_select (time_t update_time, 
		node_select_info_msg_t **node_select_info_msg_pptr)
{
        int rc;
        slurm_msg_t req_msg;
        slurm_msg_t resp_msg;
	node_info_select_request_msg_t req;

        req.last_update  = update_time;
        req_msg.msg_type = REQUEST_NODE_SELECT_INFO;
        req_msg.data     = &req;

	if (slurm_send_recv_controller_msg(&req_msg, &resp_msg) < 0)
		return SLURM_ERROR;

	switch (resp_msg.msg_type) {
	case RESPONSE_NODE_SELECT_INFO:
		*node_select_info_msg_pptr = (node_select_info_msg_t *) 
			resp_msg.data;
		break;
	case RESPONSE_SLURM_RC:
		rc = ((return_code_msg_t *) resp_msg.data)->return_code;
		slurm_free_return_code_msg(resp_msg.data);	
		if (rc) 
			slurm_seterrno_ret(rc);
		*node_select_info_msg_pptr = NULL;
		break;
	default:
		*node_select_info_msg_pptr = NULL;
		slurm_seterrno_ret(SLURM_UNEXPECTED_MSG_ERROR);
		break;
	}

        return SLURM_SUCCESS;
}

/*
 * slurm_free_node_select_info_msg - free buffer returned by
 *      slurm_load_node_select
 * IN node_select_info_msg_pptr - data is freed and pointer is set to NULL
 * RET 0 or a slurm error code
 */
extern int slurm_free_node_select_info_msg (node_select_info_msg_t **
                node_select_info_msg_pptr)
{
	if (node_select_info_msg_pptr == NULL)
		return EINVAL;

	//free it
	*node_select_info_msg_pptr = NULL;
	return SLURM_SUCCESS;
}
