/****************************************************************************\
 *  update_config.c - request that slurmctld update its configuration
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> and Kevin Tew <tew1@llnl.gov>.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <slurm/slurm.h>

#include "src/common/slurm_protocol_api.h"

static int _slurm_update (void * data, slurm_msg_type_t msg_type);

/*
 * slurm_update_job - issue RPC to a job's configuration per request, 
 *	only usable by user root or (for some parameters) the job's owner
 * IN job_msg - description of job updates
 * RET 0 on success or slurm error code
 */
int 
slurm_update_job ( job_desc_msg_t * job_msg ) 
{
	return _slurm_update ((void *) job_msg, REQUEST_UPDATE_JOB);
}

/*
 * slurm_update_node - issue RPC to a node's configuration per request, 
 *	only usable by user root
 * IN node_msg - description of node updates
 * RET 0 on success or slurm error code
 */
int 
slurm_update_node ( update_node_msg_t * node_msg ) 
{
	return _slurm_update ((void *) node_msg, REQUEST_UPDATE_NODE);
}

/*
 * slurm_update_partition - issue RPC to a partition's configuration per  
 *	request, only usable by user root
 * IN part_msg - description of partition updates
 * RET 0 on success or slurm error code
 */
int 
slurm_update_partition ( update_part_msg_t * part_msg ) 
{
	return _slurm_update ((void *) part_msg, REQUEST_UPDATE_PARTITION);
}

/*
 * slurm_delete_partition - issue RPC to delete a partition, only usable 
 *	by user root
 * IN part_msg - description of partition updates
 * RET 0 on success or slurm error code
 */
int 
slurm_delete_partition ( delete_part_msg_t * part_msg ) 
{
	return _slurm_update ((void *) part_msg, REQUEST_DELETE_PARTITION);
}

/* _slurm_update - issue RPC for all update requests */
static int 
_slurm_update (void *data, slurm_msg_type_t msg_type)
{
	int rc;
	slurm_msg_t req_msg;

	req_msg.msg_type = msg_type;
	req_msg.data     = data; 

	if (slurm_send_recv_controller_rc_msg(&req_msg, &rc) < 0)
		return SLURM_ERROR;

	if (rc != SLURM_SUCCESS)
		slurm_seterrno_ret(rc);

        return SLURM_PROTOCOL_SUCCESS;
}
