/*****************************************************************************\
 *  select_job_res.c - functions to manage data structure identifying specific
 *	CPUs allocated to a job, step or partition
 *****************************************************************************
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Written by Morris Jette <jette1@llnl.gov>.
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
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <slurm/slurm_errno.h>

#include "src/common/hostlist.h"
#include "src/common/log.h"
#include "src/common/select_job_res.h"
#include "src/common/xmalloc.h"
#include "src/common/xassert.h"
#include "src/slurmctld/slurmctld.h"


/* Create an empty select_job_res data structure */
extern select_job_res_t create_select_job_res(void)
{
	select_job_res_t select_job_res;

	select_job_res = xmalloc(sizeof(struct select_job_res));
	return select_job_res;
}


extern int build_select_job_res(select_job_res_t select_job_res,
				 char *hosts, uint16_t fast_schedule,
				 void *node_finder)
{
	hostset_t hs;
	char *host_name;
	int core_cnt = 0, host_inx = 0, sock_inx = -1;
	struct node_record * (*node_finder_local) (char *host_name);
	struct node_record *node_ptr;

	xassert(hosts);
	hs = hostset_create(hosts);
	if (!hs) {
		error("build_select_job_res: Invalid hostlist: %s", hosts);
		return SLURM_ERROR;
	}
	select_job_res->nhosts = hostset_count(hs);
	xfree(select_job_res->sockets_per_node);
	xfree(select_job_res->cores_per_socket);
	xfree(select_job_res->sock_core_rep_count);
	select_job_res->sockets_per_node = xmalloc(sizeof(uint32_t) * 
						   select_job_res->nhosts);
	select_job_res->cores_per_socket = xmalloc(sizeof(uint32_t) * 
						   select_job_res->nhosts);
	select_job_res->sock_core_rep_count = xmalloc(sizeof(uint32_t) * 
						      select_job_res->nhosts);

	node_finder_local = (struct node_record * (*) (char *host_name)) 
			    node_finder;
	while ((host_name = hostset_shift(hs))) {
		node_ptr = node_finder_local(host_name);
		if (++host_inx > select_job_res->nhosts) {
			error("build_select_job_res: "
			      "hostlist parsing problem: %s",
			      hosts);
			free(host_name);
			return SLURM_ERROR;
		} else if (node_ptr) {
			uint32_t cores, socks;
			if (fast_schedule) {
				socks = node_ptr->config_ptr->sockets;
				cores = node_ptr->config_ptr->cores;
			} else {
				socks = node_ptr->sockets;
				cores = node_ptr->cores;
			}
			if ((sock_inx < 0) ||
			    (socks != select_job_res->
				      sockets_per_node[sock_inx]) ||
			    (cores != select_job_res->
				      cores_per_socket[sock_inx])){
				sock_inx++;
				select_job_res->sockets_per_node[sock_inx] = 
						socks;
				select_job_res->cores_per_socket[sock_inx] = 
						cores;
			}
			select_job_res->sock_core_rep_count[sock_inx]++;
			core_cnt += (cores * socks);
		} else {
			error("build_select_job_res: Invalid host: %s", 
			      host_name);
			free(host_name);
			return SLURM_ERROR;
		}
		free(host_name);
	}
	hostset_destroy(hs);
	select_job_res->alloc_core_bitmap = bit_alloc(core_cnt);
	return SLURM_SUCCESS;
}

extern select_job_res_t copy_select_job_res(select_job_res_t
					    select_job_res_ptr)
{
	int i, mem_inx = 0, sock_inx = 0;
	select_job_res_t new_layout = xmalloc(sizeof(struct select_job_res));

	xassert(select_job_res_ptr);
	new_layout->nhosts = select_job_res_ptr->nhosts;
	new_layout->nprocs = select_job_res_ptr->nprocs;
	new_layout->node_req = select_job_res_ptr->node_req;
	new_layout->alloc_core_bitmap = bit_copy(select_job_res_ptr->
						 alloc_core_bitmap);

	/* Copy memory_allocated and memory_rep_count */
	new_layout->memory_allocated = xmalloc(sizeof(uint32_t) * 
					       new_layout->nhosts);
	new_layout->memory_rep_count = xmalloc(sizeof(uint32_t) * 
					       new_layout->nhosts);
	for (i=0; i<new_layout->nhosts; i++) {
		if (select_job_res_ptr->memory_rep_count[i] ==  0) {
			error("copy_select_job_res: memory_rep_count=0");
			break;
		}
		mem_inx += select_job_res_ptr->memory_rep_count[i];
		if (mem_inx >= select_job_res_ptr->nhosts) {
			i++;
			break;
		}
	}
	memcpy(new_layout->memory_allocated, 
	       select_job_res_ptr->memory_allocated, (sizeof(uint32_t) * i));
	memcpy(new_layout->memory_rep_count, 
	       select_job_res_ptr->memory_rep_count, (sizeof(uint32_t) * i));

	/* Copy sockets_per_node, cores_per_socket and core_sock_rep_count */
	new_layout->sockets_per_node = xmalloc(sizeof(uint32_t) * 
					       new_layout->nhosts);	
	new_layout->cores_per_socket = xmalloc(sizeof(uint32_t) * 
					       new_layout->nhosts);	
	new_layout->sock_core_rep_count = xmalloc(sizeof(uint32_t) * 
						  new_layout->nhosts);	
	for (i=0; i<new_layout->nhosts; i++) {
		if (select_job_res_ptr->sock_core_rep_count[i] ==  0) {
			error("copy_select_job_res: sock_core_rep_count=0");
			break;
		}
		sock_inx += select_job_res_ptr->sock_core_rep_count[i];
		if (sock_inx >= select_job_res_ptr->nhosts) {
			i++;
			break;
		}
	}
	memcpy(new_layout->sockets_per_node, 
	       select_job_res_ptr->sockets_per_node, (sizeof(uint32_t) * i));
	memcpy(new_layout->cores_per_socket, 
	       select_job_res_ptr->cores_per_socket, (sizeof(uint32_t) * i));
	memcpy(new_layout->sock_core_rep_count, 
	       select_job_res_ptr->sock_core_rep_count, 
	       (sizeof(uint32_t) * i));

	return new_layout;
}

extern void free_select_job_res(select_job_res_t *select_job_res_pptr)
{
	if (select_job_res_pptr) {
		select_job_res_t select_job_res_ptr = *select_job_res_pptr;
		xfree(select_job_res_ptr->memory_allocated);
		xfree(select_job_res_ptr->memory_rep_count);
		xfree(select_job_res_ptr->sockets_per_node);
		xfree(select_job_res_ptr->cores_per_socket);
		xfree(select_job_res_ptr->sock_core_rep_count);
		if (select_job_res_ptr->alloc_core_bitmap)
			bit_free(select_job_res_ptr->alloc_core_bitmap);
		xfree(select_job_res_ptr);
		*select_job_res_pptr = NULL;
	}
}

/* Log the contents of a select_job_res data structure using info() */
extern void log_select_job_res(select_job_res_t select_job_res_ptr)
{
	int bit_inx = 0, bit_reps, i;
	int mem_inx = 0, mem_reps = 0;
	int node_inx;
	int sock_inx = 0, sock_reps = 0;

	xassert(select_job_res_ptr);
	info("====================");
	info("nhosts:%u nprocs:%u node_req:%u", 
	     select_job_res_ptr->nhosts, select_job_res_ptr->nprocs,
	     select_job_res_ptr->node_req);
	for (node_inx=0; node_inx<select_job_res_ptr->nhosts; node_inx++) {
		info("Node[%d]:", node_inx);

		if (mem_reps >= 
		    select_job_res_ptr->memory_rep_count[mem_inx]) {
			mem_inx++;
			mem_reps = 0;
		}
		mem_reps++;

		if (sock_reps >= 
		    select_job_res_ptr->sock_core_rep_count[sock_inx]) {
			sock_inx++;
			sock_reps = 0;
		}
		sock_reps++;

		info("  Mem(MB):%u  Sockets:%u  Cores:%u", 
		     select_job_res_ptr->memory_allocated[mem_inx],
		     select_job_res_ptr->sockets_per_node[sock_inx],
		     select_job_res_ptr->cores_per_socket[sock_inx]);

		bit_reps = select_job_res_ptr->sockets_per_node[sock_inx] *
			   select_job_res_ptr->cores_per_socket[sock_inx];
		for (i=0; i<bit_reps; i++) {
			if (bit_test(select_job_res_ptr->alloc_core_bitmap,
				     bit_inx)) {
				info("  Socket[%d] Core[%d] in use",
				     (i / select_job_res_ptr->
				          cores_per_socket[sock_inx]),
				     (i % select_job_res_ptr->
					  cores_per_socket[sock_inx]));
			}
			bit_inx++;
		}
	}
	info("====================");
}

extern void pack_select_job_res(select_job_res_t select_job_res_ptr, 
				Buf buffer)
{
	int i;
	uint32_t core_cnt = 0, mem_recs = 0, sock_recs = 0;

	xassert(select_job_res_ptr);
	pack32(select_job_res_ptr->nhosts, buffer);
	pack32(select_job_res_ptr->nprocs, buffer);
	pack8(select_job_res_ptr->node_req, buffer);
	for (i=0; i<select_job_res_ptr->nhosts; i++) {
		mem_recs += select_job_res_ptr->memory_rep_count[i];
		if (mem_recs >= select_job_res_ptr->nhosts)
			break;
	}
	i++;
	pack32_array(select_job_res_ptr->memory_allocated,  
		     (uint32_t) i, buffer);
	pack32_array(select_job_res_ptr->memory_rep_count, 
		     (uint32_t) i, buffer);

	for (i=0; i<select_job_res_ptr->nhosts; i++) {
		core_cnt += select_job_res_ptr->sockets_per_node[i] *
			    select_job_res_ptr->cores_per_socket[i] *
			    select_job_res_ptr->sock_core_rep_count[i];
		sock_recs += select_job_res_ptr->sock_core_rep_count[i];
		if (sock_recs >= select_job_res_ptr->nhosts)
			break;
	}
	i++;
	pack32_array(select_job_res_ptr->sockets_per_node,
		     (uint32_t) i, buffer);
	pack32_array(select_job_res_ptr->cores_per_socket,
		     (uint32_t) i, buffer);
	pack32_array(select_job_res_ptr->sock_core_rep_count, 
		     (uint32_t) i, buffer);
	pack32(core_cnt, buffer);
	xassert (core_cnt == bit_size(select_job_res_ptr->alloc_core_bitmap));
	pack_bit_fmt(select_job_res_ptr->alloc_core_bitmap, buffer);
}

extern int  unpack_select_job_res(select_job_res_t *select_job_res_pptr, Buf buffer)
{
	char *bit_fmt = NULL;
	uint32_t core_cnt, tmp32;
	select_job_res_t select_job_res;

	xassert(select_job_res_pptr);
	select_job_res = xmalloc(sizeof(struct select_job_res));
	safe_unpack32(&select_job_res->nhosts, buffer);
	safe_unpack32(&select_job_res->nprocs, buffer);
	safe_unpack8(&select_job_res->node_req, buffer);
	safe_unpack32_array(&select_job_res->memory_allocated,
			    &tmp32, buffer);
	safe_unpack32_array(&select_job_res->memory_rep_count,
			    &tmp32, buffer);
	safe_unpack32_array(&select_job_res->sockets_per_node,
			    &tmp32, buffer);
	safe_unpack32_array(&select_job_res->cores_per_socket,
			    &tmp32, buffer);
	safe_unpack32_array(&select_job_res->sock_core_rep_count,
			    &tmp32, buffer);
	safe_unpack32(&core_cnt, buffer);    /* NOTE: Not part of struct */
	safe_unpackstr_xmalloc(&bit_fmt, &tmp32, buffer);
	select_job_res->alloc_core_bitmap = bit_alloc((bitoff_t) core_cnt);
	if (bit_unfmt(select_job_res->alloc_core_bitmap, bit_fmt))
		goto unpack_error;
	xfree(bit_fmt);
	*select_job_res_pptr = select_job_res;
	return SLURM_SUCCESS;

  unpack_error:
	xfree(select_job_res);
	xfree(bit_fmt);
	*select_job_res_pptr = NULL;
	return SLURM_ERROR;
}

extern int get_select_job_res_bit(select_job_res_t select_job_res_ptr, 
				  uint32_t node_id, uint32_t socket_id, 
				  uint32_t core_id)
{
	int i, bit_inx = 0;

	xassert(select_job_res_ptr);

	for (i=0; i<select_job_res_ptr->nhosts; i++) {
		if (select_job_res_ptr->sock_core_rep_count[i] <= node_id) {
			bit_inx += select_job_res_ptr->sockets_per_node[i] *
				   select_job_res_ptr->cores_per_socket[i] *
				   select_job_res_ptr->sock_core_rep_count[i];
			node_id -= select_job_res_ptr->sock_core_rep_count[i];
		} else {
			bit_inx += select_job_res_ptr->sockets_per_node[i] *
				   select_job_res_ptr->cores_per_socket[i] *
				   node_id;
			bit_inx += select_job_res_ptr->cores_per_socket[i] *
				   socket_id;
			bit_inx += core_id;
			break;
		}
	}
	i = bit_size(select_job_res_ptr->alloc_core_bitmap);
	if (bit_inx >= i) {
		error("get_select_job_res_bit: offset >= bitmap size "
		      "(%d >= %d)", bit_inx, i);
		return 0;
	}

	return bit_test(select_job_res_ptr->alloc_core_bitmap, bit_inx);
}

extern int set_select_job_res_bit(select_job_res_t select_job_res_ptr, 
				  uint32_t node_id, uint32_t socket_id, 
				  uint32_t core_id)
{
	int i, bit_inx = 0;

	xassert(select_job_res_ptr);

	for (i=0; i<select_job_res_ptr->nhosts; i++) {
		if (select_job_res_ptr->sock_core_rep_count[i] <= node_id) {
			bit_inx += select_job_res_ptr->sockets_per_node[i] *
				   select_job_res_ptr->cores_per_socket[i] *
				   select_job_res_ptr->sock_core_rep_count[i];
			node_id -= select_job_res_ptr->sock_core_rep_count[i];
		} else {
			bit_inx += select_job_res_ptr->sockets_per_node[i] *
				   select_job_res_ptr->cores_per_socket[i] *
				   node_id;
			bit_inx += select_job_res_ptr->cores_per_socket[i] *
				   socket_id;
			bit_inx += core_id;
			break;
		}
	}
	i = bit_size(select_job_res_ptr->alloc_core_bitmap);
	if (bit_inx >= i) {
		error("set_select_job_res_bit: offset >= bitmap size "
		      "(%d >= %d)", bit_inx, i);
		return SLURM_ERROR;
	}

	bit_set(select_job_res_ptr->alloc_core_bitmap, bit_inx);
	return SLURM_SUCCESS;
}
