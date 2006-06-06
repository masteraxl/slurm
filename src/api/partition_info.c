/*****************************************************************************\
 *  partition_info.c - get/print the partition state information of slurm
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> et. al.
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

#include "src/api/job_info.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"

#define HUGE_BUF 4096

/*
 * slurm_print_partition_info_msg - output information about all Slurm 
 *	partitions based upon message as loaded using slurm_load_partitions
 * IN out - file to write to
 * IN part_info_ptr - partitions information message pointer
 * IN one_liner - print as a single line if true
 */
void slurm_print_partition_info_msg ( FILE* out, 
		partition_info_msg_t * part_info_ptr, int one_liner )
{
	int i ;
	partition_info_t * part_ptr = part_info_ptr->partition_array ;
	char time_str[32];

	slurm_make_time_str ((time_t *)&part_info_ptr->last_update, time_str, 
		sizeof(time_str));
	fprintf( out, "Partition data as of %s, record count %d\n",
		time_str, part_info_ptr->record_count);

	for (i = 0; i < part_info_ptr->record_count; i++) {
		slurm_print_partition_info ( out, & part_ptr[i], one_liner ) ;
	}

}

/*
 * slurm_print_partition_info - output information about a specific Slurm 
 *	partition based upon message as loaded using slurm_load_partitions
 * IN out - file to write to
 * IN part_ptr - an individual partition information record pointer
 * IN one_liner - print as a single line if true
 */
void slurm_print_partition_info ( FILE* out, partition_info_t * part_ptr, 
				  int one_liner )
{
	char *print_this = slurm_sprint_partition_info(part_ptr, one_liner);
	fprintf ( out, "%s", print_this);
	xfree(print_this);
}


/*
 * slurm_sprint_partition_info - output information about a specific Slurm 
 *	partition based upon message as loaded using slurm_load_partitions
 * IN part_ptr - an individual partition information record pointer
 * IN one_liner - print as a single line if true
 * RET out - char * containing formatted output (must be freed after call)
 *           NULL is returned on failure.
 */
char *slurm_sprint_partition_info ( partition_info_t * part_ptr, 
				    int one_liner )
{
	int j;
	char tmp1[7];
	char tmp2[100];
	int len = 0;
	int lentmp2 = 0;
	int tmplen = 0;
	char *out = xmalloc(HUGE_BUF);
	/****** Line 1 ******/
	sprintf ( tmp2, "PartitionName=%s ", part_ptr->name);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out, "%s", tmp2);
	len += lentmp2;

	convert_to_kilo(part_ptr->total_nodes, tmp1);
	sprintf ( tmp2, "TotalNodes=%s ", tmp1);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	convert_to_kilo(part_ptr->total_cpus, tmp1);
	sprintf ( tmp2, "TotalCPUs=%s ", tmp1);	
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (part_ptr->root_only)
		sprintf ( tmp2, "RootOnly=YES");
	else
		sprintf ( tmp2, "RootOnly=NO");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (one_liner)
		sprintf ( tmp2, " ");
	else
		sprintf ( tmp2, "\n   ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	/****** Line 2 ******/
	if (part_ptr->default_part)
		sprintf ( tmp2, "Default=YES ");
	else
		sprintf ( tmp2, "Default=NO ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (part_ptr->shared == SHARED_NO)
		sprintf ( tmp2, "Shared=NO ");
	else if (part_ptr->shared == SHARED_YES)
		sprintf ( tmp2, "Shared=YES ");
	else
		sprintf ( tmp2, "Shared=FORCE ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (part_ptr->state_up)
		sprintf ( tmp2, "State=UP ");
	else
		sprintf ( tmp2, "State=DOWN ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (part_ptr->max_time == INFINITE)
		sprintf ( tmp2, "MaxTime=UNLIMITED ");
	else
		sprintf ( tmp2, "MaxTime=%u ", part_ptr->max_time);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (part_ptr->hidden)
		sprintf ( tmp2, "Hidden=YES");
	else
		sprintf ( tmp2, "Hidden=NO");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (one_liner)
		sprintf ( tmp2, " ");
	else
		sprintf ( tmp2, "\n   ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	/****** Line 3 ******/
	convert_to_kilo(part_ptr->min_nodes, tmp1);
	sprintf ( tmp2, "MinNodes=%s ", tmp1);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (part_ptr->max_nodes == INFINITE)
		sprintf ( tmp2, "MaxNodes=UNLIMITED ");
	else {
		convert_to_kilo(part_ptr->max_nodes, tmp1);
		sprintf ( tmp2, "MaxNodes=%s ", tmp1);
	}
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if ((part_ptr->allow_groups == NULL) || 
	    (part_ptr->allow_groups[0] == '\0'))
		sprintf ( tmp2, "AllowGroups=ALL");
	else
		sprintf ( tmp2, "AllowGroups=%s", part_ptr->allow_groups);
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	if (one_liner)
		sprintf ( tmp2, " ");
	else
		sprintf ( tmp2, "\n   ");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	
	/****** Line 4 ******/
#ifdef HAVE_BG
	sprintf ( tmp2, "BasePartitions=%s BPIndices=", part_ptr->nodes);
#else
	sprintf ( tmp2, "Nodes=%s NodeIndices=", part_ptr->nodes);
#endif
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	len += lentmp2;
	
	for (j = 0; part_ptr->node_inx; j++) {
		if (j > 0)
			sprintf( tmp2, ",%d", part_ptr->node_inx[j]);
		else
			sprintf( tmp2, "%d", part_ptr->node_inx[j]);
		lentmp2 = strlen(tmp2);
		tmplen += lentmp2;
		if(tmplen>HUGE_BUF) {
			j = len + HUGE_BUF;
			xrealloc(out, j);
		}
		sprintf ( out+len, "%s", tmp2);
		len += lentmp2;	
		if (part_ptr->node_inx[j] == -1)
			break;
	}
	sprintf( tmp2, "\n\n");
	lentmp2 = strlen(tmp2);
	tmplen += lentmp2;
	if(tmplen>HUGE_BUF) {
		j = len + HUGE_BUF;
		xrealloc(out, j);
	}
	sprintf ( out+len, "%s", tmp2);
	
	return out;
}



/*
 * slurm_load_partitions - issue RPC to get slurm all partition configuration  
 *	information if changed since update_time 
 * IN update_time - time of current configuration data
 * IN partition_info_msg_pptr - place to store a partition configuration 
 *	pointer
 * IN show_flags - partition filtering options
 * RET 0 or a slurm error code
 * NOTE: free the response using slurm_free_partition_info_msg
 */
extern int slurm_load_partitions (time_t update_time, 
		partition_info_msg_t **resp, uint16_t show_flags)
{
        int rc;
        slurm_msg_t req_msg;
        slurm_msg_t resp_msg;
        part_info_request_msg_t req;

        req.last_update  = update_time;
	req.show_flags   = show_flags;
        req_msg.msg_type = REQUEST_PARTITION_INFO;
        req_msg.data     = &req;
	
	if (slurm_send_recv_controller_msg(&req_msg, &resp_msg) < 0)
		return SLURM_ERROR;
	
	switch (resp_msg.msg_type) {
	case RESPONSE_PARTITION_INFO:
		*resp = (partition_info_msg_t *) resp_msg.data;
		break;
	case RESPONSE_SLURM_RC:
		rc = ((return_code_msg_t *) resp_msg.data)->return_code;
		slurm_free_return_code_msg(resp_msg.data);	
		if (rc) 
			slurm_seterrno_ret(rc);
		*resp = NULL;
		break;
	default:
		slurm_seterrno_ret(SLURM_UNEXPECTED_MSG_ERROR);
		break;
	}

	return SLURM_PROTOCOL_SUCCESS;
}
