/*****************************************************************************\
 *  select_job_res.h - functions to manage data structure identifying specific
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

#ifndef _SELECT_JOB_RES_H
#define _SELECT_JOB_RES_H

#if HAVE_CONFIG_H
#  include "config.h"
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  else
#    if HAVE_STDINT_H
#      include <stdint.h>
#    endif
#  endif			/* HAVE_INTTYPES_H */
#endif

#include "src/common/bitstring.h"
#include "src/common/pack.h"
#include "src/slurmctld/slurmctld.h"

/* struct select_job_res defines exactly which resources are allocated
 *	to a job, step, partition, etc.
 *
 * node_cnt		- Number of nodes in the allocation
 * memory_reserved	- MB per node reserved
 * memory_rep_count	- How many consecutive nodes that memory_reserved
 *			  applies to
 * sockets_per_node	- Count of sockets on this node
 * cores_per_socket	- Count of cores per socket on this node
 * sock_core_rep_count	- How many consecutive nodes that sockets_per_node
 *			  and cores_per_socket apply to
 * allocated_cores	- bitmap of selected cores for all nodes and sockets
 *
 * Sample layout:
 *   |               Node_0              |               Node_1              |
 *   |      Sock_0     |      Sock_1     |      Sock_0     |      Sock_1     |
 *   | Core_0 | Core_1 | Core_0 | Core_1 | Core_0 | Core_1 | Core_0 | Core_1 |
 *   | Bit_0  | Bit_1  | Bit_2  | Bit_3  | Bit_4  | Bit_5  | Bit_6  | Bit_7  |
 */
struct select_job_res {
	uint32_t	node_cnt;
	uint32_t *	memory_reserved;
	uint32_t *	memory_rep_count;
	uint32_t *	sockets_per_node;
	uint32_t *	cores_per_socket;
	uint32_t *	sock_core_rep_count;
	bitstr_t *	allocated_cores;
};

/* Create an empty select_job_res data structure */
extern select_job_res_t create_select_job_res(void);

/* Set the socket and core counts associated with a set of selected
 * nodes of a select_job_res data structure based upon slurmctld state.
 * Call this ONLY from slurmctld. We pass a pointer to slurmctld's
 * find_node_record function so this module can be loaded in libslurm
 * and the other functions used from slurmd. Example of use:
 * select_job_res_t select_job_res_ptr = create_select_job_res();
 * rc = build_select_job_res(select_job_res_ptr,
 *			     "tux[2,5,10-12,16]", 
 *			     slurmctld_conf.fast_schedule,
 *			     find_node_record);
 */
extern int build_select_job_res(select_job_res_t select_job_res_ptr,
				char *hosts, uint16_t fast_schedule,
				void *node_finder);

/* Make a copy of a select_job_res data structure, free using free_select_job_res() */
extern select_job_res_t copy_select_job_res(select_job_res_t 
					    select_job_res_ptr);

/* Free select_job_res data structure created using copy_select_job_res() or
 *	unpack_select_job_res() */
extern void free_select_job_res(select_job_res_t *select_job_res_pptr);

/* Log the contents of a select_job_res data structure using info() */
extern void log_select_job_res(select_job_res_t select_job_res_ptr);

/* Un/pack full select_job_res data structure */
extern void pack_select_job_res(select_job_res_t select_job_res_ptr, 
				Buf buffer);
extern int unpack_select_job_res(select_job_res_t *select_job_res_pptr, 
				 Buf buffer);

/* Get/set bit value at specified location.
 *	node_id, socket_id and core_id are all zero origin */
extern int get_select_job_res_bit(select_job_res_t select_job_res_ptr, 
				  uint32_t node_id,
				  uint32_t socket_id, uint32_t core_id);
extern int set_select_job_res_bit(select_job_res_t select_job_res_ptr, 
				  uint32_t node_id,
				  uint32_t socket_id, uint32_t core_id);
#endif /* !_SELECT_JOB_RES_H */
