/*****************************************************************************\
 *  bgl_job_place.c - blue gene job placement (e.g. base partition selection)
 *  functions. 
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Dan Phung <phung4@llnl.gov> and Morris Jette <jette1@llnl.gov>
 *  and Danny Auble <da@llnl.gov>
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

#include "src/common/node_select.h"
#include "bluegene.h"

#define BUFSIZE 4096
#define BITSIZE 128

#define _DEBUG 0

#define SWAP(a,b,t)	\
_STMT_START {		\
	(t) = (a);	\
	(a) = (b);	\
	(b) = (t);	\
} _STMT_END

static int  _find_best_partition_match(struct job_record* job_ptr,
				bitstr_t* slurm_part_bitmap,
				int min_nodes, int max_nodes,
				int spec, bgl_record_t** found_bgl_record);
static void _rotate_geo(uint16_t *req_geometry, int rot_cnt);

/* Rotate a 3-D geometry array through its six permutations */
static void _rotate_geo(uint16_t *req_geometry, int rot_cnt)
{
	uint16_t tmp;

	switch (rot_cnt) {
		case 0:		/* ABC -> ACB */
			SWAP(req_geometry[Y], req_geometry[Z], tmp);
			break;
		case 1:		/* ACB -> CAB */
			SWAP(req_geometry[X], req_geometry[Y], tmp);
			break;
		case 2:		/* CAB -> CBA */
			SWAP(req_geometry[Y], req_geometry[Z], tmp);
			break;
		case 3:		/* CBA -> BCA */
			SWAP(req_geometry[X], req_geometry[Y], tmp);
			break;
		case 4:		/* BCA -> BAC */
			SWAP(req_geometry[Y], req_geometry[Z], tmp);
			break;
		case 5:		/* BAC -> ABC */
			SWAP(req_geometry[X], req_geometry[Y], tmp);
			break;
	}
}

/*
 * finds the best match for a given job request 
 * 
 * IN - int spec right now holds the place for some type of
 * specification as to the importance of certain job params, for
 * instance, geometry, type, size, etc.
 * 
 * OUT - part_id of matched partition, NULL otherwise
 * returns 1 for error (no match)
 * 
 */
static int _find_best_partition_match(struct job_record* job_ptr, 
		bitstr_t* slurm_part_bitmap, int min_nodes, int max_nodes,
		int spec, bgl_record_t** found_bgl_record)
{
	ListIterator itr;
	bgl_record_t* record;
	int i;
	uint16_t req_geometry[SYSTEM_DIMENSIONS];
	uint16_t conn_type, node_use, rotate, target_size = 1;

	sort_bgl_record_inc_size(bgl_list);

	select_g_get_jobinfo(job_ptr->select_jobinfo,
		SELECT_DATA_CONN_TYPE, &conn_type);
	select_g_get_jobinfo(job_ptr->select_jobinfo,
		SELECT_DATA_GEOMETRY, req_geometry);
	select_g_get_jobinfo(job_ptr->select_jobinfo,
		SELECT_DATA_NODE_USE, &node_use);
	select_g_get_jobinfo(job_ptr->select_jobinfo,
		SELECT_DATA_ROTATE, &rotate);
	for (i=0; i<SYSTEM_DIMENSIONS; i++)
		target_size *= req_geometry[i];
	if (target_size == 0)	/* no geometry specified */
		target_size = min_nodes;

	/* this is where we should have the control flow depending on
	 * the spec arguement */
	itr = list_iterator_create(bgl_list);
	*found_bgl_record = NULL;

	/*
	 * FIXME: NEED TO PUT THIS LOGIC IN: 
	 * if RM_NAV, then the partition with both the TORUS and the
	 * dims should be favored over the MESH and the dims, but
	 * foremost is the correct num of dims. 
	 */
	debug("number of partitions to check: %d", list_count(bgl_list));
	while ((record = (bgl_record_t*) list_next(itr))) {
		/*
		 * check that the number of nodes is suitable
		 */
 		if ((record->bp_count < min_nodes)
		||  (max_nodes != 0 && record->bp_count > max_nodes)
		||  (record->bp_count < target_size)) {
			debug("partition %s node count not suitable",
				record->bgl_part_id);
			continue;
		}
		
		/*
		 * Next we check that this partition's bitmap is within 
		 * the set of nodes which the job can use. 
		 * Nodes not available for the job could be down,
		 * drained, allocated to some other job, or in some 
		 * SLURM partition not available to this job.
		 */
		char bitstring[BITSIZE];
		char bitstring2[BITSIZE];
		bit_fmt(bitstring, BITSIZE, record->bitmap);
		bit_fmt(bitstring2, BITSIZE, slurm_part_bitmap);
		printf("These are the bitmaps, is %s inside of %s?\n",bitstring,bitstring2);
		if (!bit_super_set(record->bitmap, slurm_part_bitmap)) {
			debug("bgl partition %s has nodes not usable by this "
				"job", record->bgl_part_id);
			continue;
		}

		/*
		 * Insure that any required nodes are in this BGL partition
		 */
		if (job_ptr->details->req_node_bitmap
		&& (!bit_super_set(job_ptr->details->req_node_bitmap,
				record->bitmap))) {
			info("bgl partition %s lacks required nodes",
				record->bgl_part_id);
			continue;
		}

		/***********************************************/
		/* check the connection type specified matches */
		/***********************************************/
		if ((conn_type != record->conn_type)
		&&  (conn_type != SELECT_NAV)) {
			debug("bgl partition %s conn-type not usable", 
				record->bgl_part_id);
			continue;
		} 

		/***********************************************/
		/* check the node_use specified matches        */
		/***********************************************/
		if ((node_use != record->node_use) 
		&&  (node_use != SELECT_NAV)) {
			debug("bgl partition %s node-use not usable", 
					record->bgl_part_id);
			continue;
		}

		/*****************************************/
		/* match up geometry as "best" possible  */
		/*****************************************/
		if (req_geometry[0] == 0)
			;	/* Geometry not specified */
		else {	/* match requested geometry */
			bool match = false;
			int rot_cnt;	/* attempt six rotations  */

			for (rot_cnt=0; rot_cnt<6; rot_cnt++) {
				if ((record->coord[X] >= req_geometry[X])
				&&  (record->coord[Y] >= req_geometry[Y])
				&&  (record->coord[Z] >= req_geometry[Z])) {
					match = true;
					break;
				}
				if (!rotate)
					break;
				_rotate_geo(req_geometry, rot_cnt);
			}

			if (!match) 
				continue;	/* Not usable */
		}

		if ((*found_bgl_record == NULL)
		||  (record->bp_count < (*found_bgl_record)->bp_count)) {
			*found_bgl_record = record;
			if (record->bp_count == target_size)
				break;
		}
	}
	list_iterator_destroy(itr);	
	
	/* set the bitmap and do other allocation activities */
	if (*found_bgl_record) {
		debug("_find_best_partition_match %s <%s>", 
			(*found_bgl_record)->bgl_part_id, 
			(*found_bgl_record)->nodes);
		bit_and(slurm_part_bitmap, (*found_bgl_record)->bitmap);
		return SLURM_SUCCESS;
	}
	
	debug("_find_best_partition_match none found");
	return SLURM_ERROR;
}

/*
 * Try to find resources for a given job request
 * IN job_ptr - pointer to job record in slurmctld
 * IN/OUT bitmap - nodes availble for assignment to job, clear those not to
 *	be used
 * IN min_nodes, max_nodes  - minimum and maximum number of nodes to allocate
 *	to this job (considers slurm partition limits)
 * RET - SLURM_SUCCESS if job runnable now, error code otherwise
 */
extern int submit_job(struct job_record *job_ptr, bitstr_t *slurm_part_bitmap,
		      int min_nodes, int max_nodes)
{
	int spec = 1; /* this will be like, keep TYPE a priority, etc,  */
	bgl_record_t* record;
	char buf[100];

	select_g_sprint_jobinfo(job_ptr->select_jobinfo, buf, sizeof(buf), 
		SELECT_PRINT_MIXED);
	debug("bluegene:submit_job: %s nodes=%d-%d", buf, min_nodes, max_nodes);
	
	if ((_find_best_partition_match(job_ptr, slurm_part_bitmap, min_nodes, 
					max_nodes, spec, &record)) != SLURM_SUCCESS) {
		return SLURM_ERROR;
	} else {
		/* now we place the part_id into the env of the script to run */
		char bgl_part_id[BITSIZE];
#ifdef HAVE_BGL_FILES
		snprintf(bgl_part_id, BITSIZE, "%s", record->bgl_part_id);
#else
		snprintf(bgl_part_id, BITSIZE, "UNDEFINED");
#endif
		select_g_set_jobinfo(job_ptr->select_jobinfo,
			SELECT_DATA_PART_ID, bgl_part_id);
	}

	return SLURM_SUCCESS;
}
