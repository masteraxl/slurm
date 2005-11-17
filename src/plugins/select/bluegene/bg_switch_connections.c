/*****************************************************************************\
 *  bg_switch_connections.c - Blue Gene switch management functions, 
 *  establish switch connections
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Dan Phung <phung4@llnl.gov> and Danny Auble <da@llnl.gov>
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

#include "bluegene.h"

#ifdef HAVE_BG_FILES

List bg_bp_list;

static int _get_bp_by_location(rm_BGL_t* my_bg, 
			       int* curr_coord, 
			       rm_BP_t** bp);
//static int _set_switch(rm_switch_t* curr_switch, pa_connection_t *int_wire);
static int _add_switch_conns(rm_switch_t* curr_switch, 
			     bg_switch_t *bg_switch);
static int _lookat_path(bg_bp_t *bg_bp, 
			pa_switch_t *curr_switch, 
			int source, 
			int target, 
			int dim);
static int _destroy_bg_bp_list(List bg_bp);
/** 
 * this is just stupid.  there are some implicit rules for where
 * "NextBP" goes to, but we don't know, so we have to do this.
 */
static int _get_bp_by_location(rm_BGL_t* my_bg, int* curr_coord, rm_BP_t** bp)
{
	int i, bp_num, rc;
	rm_location_t loc;

	if ((rc = rm_get_data(my_bg, RM_BPNum, &bp_num)) != STATUS_OK) {
		fatal("rm_get_data: RM_BPNum: %s", bg_err_str(rc));
		return SLURM_ERROR;
	}

	for (i=0; i<bp_num; i++){
		if(i) {
			if ((rc = rm_get_data(my_bg, RM_NextBP, bp)) 
			    != STATUS_OK) {
				fatal("rm_get_data: RM_NextBP: %s", 
				      bg_err_str(rc));
				return SLURM_ERROR;
			}	
		} else {
			if ((rc = rm_get_data(my_bg, RM_FirstBP, bp)) 
			    != STATUS_OK) {
				fatal("rm_get_data: RM_FirstBP: %s", 
				      bg_err_str(rc));
				return SLURM_ERROR;
			}
		}	
		if ((rc = rm_get_data(*bp, RM_BPLoc, &loc)) != STATUS_OK) {
			fatal("rm_get_data: RM_BPLoc: %s", bg_err_str(rc));
			return SLURM_ERROR;
		}

		if ((loc.X == curr_coord[X])
		&&  (loc.Y == curr_coord[Y])
		&&  (loc.Z == curr_coord[Z])) {
			return SLURM_SUCCESS;
		}
	}

	// error("_get_bp_by_location: could not find specified bp.");
	return SLURM_ERROR;
}

static int _add_switch_conns(rm_switch_t* curr_switch, 
			     bg_switch_t *bg_switch)
{
	ListIterator itr;
	bg_conn_t *bg_conn;
	
	int firstconnect=1;
	rm_connection_t conn;
	int j, rc;
	int conn_num=0;
	int port = 0;
	
	itr = list_iterator_create(bg_switch->conn_list);
	while((bg_conn = list_next(itr)) != NULL) {
		if(bg_conn->source == bg_conn->target)
			continue;
		
		for(j=0;j<2;j++) {
			switch(j) {
			case 0:
				port = bg_conn->source;
				break;
			case 1:
				port = bg_conn->target;
				break;
			}
			switch(port) {
			case 0:
				conn.p2 = RM_PORT_S0; 
				break;
			case 1:
				conn.p1 = RM_PORT_S1;
				break;
			case 2:
				conn.p1 = RM_PORT_S2;
				break;
			case 3:
				conn.p2 = RM_PORT_S3; 
				break;
			case 4:
				conn.p1 = RM_PORT_S4;
				break;
			case 5:
				conn.p2 = RM_PORT_S5; 
				break;
			}
		}
		conn.part_state = RM_PARTITION_READY;
		
		if(firstconnect) {
			if ((rc = rm_set_data(curr_switch, 
					RM_SwitchFirstConnection, &conn)) 
					!= STATUS_OK) {
				fatal("rm_set_data(RM_SwitchFirstConnection):"
					" %s", bg_err_str(rc));
				list_iterator_destroy(itr);
				return SLURM_ERROR;
			}
			firstconnect=0;
		} else {
			if ((rc = rm_set_data(curr_switch, 
					RM_SwitchNextConnection, &conn)) 
					!= STATUS_OK) {
				fatal("rm_set_data(RM_SwitchNextConnection):"
					" %s", bg_err_str(rc));
				list_iterator_destroy(itr);
				return SLURM_ERROR;
			}
		} 
		conn_num++;
		debug2("adding %d -> %d",bg_conn->source, bg_conn->target);
	}
	list_iterator_destroy(itr);
	if ((rc = rm_set_data(curr_switch, RM_SwitchConnNum, &conn_num)) 
	    != STATUS_OK) {
		fatal("rm_set_data: RM_SwitchConnNum: %s", bg_err_str(rc));
		return SLURM_ERROR;
	} 	
	return SLURM_SUCCESS;
}

static int _lookat_path(bg_bp_t *bg_bp, pa_switch_t *curr_switch, 
			int source, int target, int dim) 
{
	ListIterator bg_itr, switch_itr, conn_itr;
	bg_switch_t *bg_switch;
	bg_conn_t *bg_conn;
	int *node_tar;
	int port_tar;
	int port_tar1;
	int *node_src;
	pa_switch_t *next_switch; 
	
	switch_itr = list_iterator_create(bg_bp->switch_list);
	while((bg_switch = list_next(switch_itr)) != NULL) {
		if(bg_switch->dim == dim)
			break;
	}
	list_iterator_destroy(switch_itr);
	
	if(bg_switch == NULL) {
		bg_switch = xmalloc(sizeof(bg_switch_t));
		bg_switch->dim=dim;
		bg_switch->conn_list = list_create(NULL);
		list_append(bg_bp->switch_list, bg_switch);
	}
		
	port_tar = curr_switch->int_wire[source].port_tar;
	
	conn_itr = list_iterator_create(bg_switch->conn_list);
	while((bg_conn = list_next(conn_itr)) != NULL) {
		if(port_tar == curr_switch->ext_wire[port_tar].port_tar) {
			//list_delete(conn_itr);
			//continue;
			debug3("I found these %d %d",port_tar, 
			      curr_switch->ext_wire[port_tar].port_tar);
		}
		if(((bg_conn->source == port_tar)
		    && (bg_conn->target == source))
		   || ((bg_conn->source == source)
		       && (bg_conn->target == port_tar)))
			break;
	}
	list_iterator_destroy(conn_itr);
	
	if(bg_conn == NULL) {
		bg_conn = xmalloc(sizeof(bg_conn_t));
		bg_conn->source = source;
		bg_conn->target = port_tar;
		
		list_append(bg_switch->conn_list, bg_conn);
	} else {		
		return SLURM_SUCCESS;	
	}
	if(port_tar==target) {
		return SLURM_SUCCESS;
	}
	/* keep this around to tell us where we are coming from */
	port_tar1 = port_tar;
	/* set port target to to where the external wire is 
	   going to on the next node */
	port_tar = curr_switch->ext_wire[port_tar1].port_tar;
	/* set node target to where the external wire is going to */
	node_tar = curr_switch->ext_wire[port_tar1].node_tar;
	/* set source to the node you are on */
	node_src = curr_switch->ext_wire[0].node_tar;

	debug("dim %d trying from %d%d%d %d -> %d%d%d %d",
	      dim,
	      node_src[X], 
	      node_src[Y], 
	      node_src[Z],
	      port_tar1,
	      node_tar[X], 
	      node_tar[Y], 
	      node_tar[Z],
	      port_tar);


	bg_itr = list_iterator_create(bg_bp_list);
	while((bg_bp = list_next(bg_itr)) != NULL) {
		if((bg_bp->coord[X] == node_tar[X]) 
		   && (bg_bp->coord[Y] == node_tar[Y]) 
		   && (bg_bp->coord[Z] == node_tar[Z]))
			break;
	}
	list_iterator_destroy(bg_itr);
	/* It appears this is a past through node */
	if(bg_bp == NULL) {
		bg_bp = xmalloc(sizeof(bg_bp_t));
		bg_bp->coord = node_tar;
		bg_bp->switch_list = list_create(NULL);
		list_append(bg_bp_list, bg_bp);
		bg_bp->used = 0;
	}
	
	next_switch = &pa_system_ptr->
		grid[node_tar[X]][node_tar[Y]][node_tar[Z]].axis_switch[dim];
	
	if(_lookat_path(bg_bp, next_switch, port_tar, target, dim) 
			== SLURM_ERROR)
		return SLURM_ERROR;
	
	return SLURM_SUCCESS;
}

static int _destroy_bg_bp_list(List bg_bp_list)
{
	bg_switch_t *bg_switch;
	bg_conn_t *bg_conn;
	bg_bp_t *bg_bp;
	
	if(bg_bp_list) {
		while((bg_bp = list_pop(bg_bp_list)) != NULL) {
			while((bg_switch = list_pop(bg_bp->switch_list)) 
			      != NULL) {
				while((bg_conn = list_pop(
					       bg_switch->conn_list)) 
				      != NULL) {
					if(bg_conn)
						xfree(bg_conn);
				}
				list_destroy(bg_switch->conn_list);
				if(bg_switch)
					xfree(bg_switch);
			}
			list_destroy(bg_bp->switch_list);
			if(bg_bp)
				xfree(bg_bp);
		}
		list_destroy(bg_bp_list);
	}
	return SLURM_SUCCESS;
}

extern int configure_small_partition(bg_record_t *bg_record)
{
	bool small = true;
	ListIterator itr;
	pa_node_t* pa_node = NULL;
	int rc = SLURM_SUCCESS;
	rm_BP_t *curr_bp;
	rm_bp_id_t bp_id = NULL;
	int num_ncards = 4;
	rm_nodecard_t *ncard;
	rm_nodecard_list_t *ncard_list;
	rm_quarter_t quarter;
	int num, i;

	if(bg_record->bp_count != 1) {
		error("Requesting small partition with %d bps, needs to be 1.",
		      bg_record->bp_count);
		return SLURM_ERROR;
	}
	
	/* set that we are doing a small partition */
	if ((rc = rm_set_data(bg_record->bg_part, RM_PartitionSmall, 
			      &small)) != STATUS_OK) {
		fatal("rm_set_data(RM_PartitionPsetsPerBP)", bg_err_str(rc));
	}

	if ((rc = rm_set_data(bg_record->bg_part,
			      RM_PartitionNodeCardNum,
			      &num_ncards))
	    != STATUS_OK) {
		fatal("rm_set_data: RM_PartitionBPNum: %s", bg_err_str(rc));
	}

	itr = list_iterator_create(bg_record->bg_part_list);
	pa_node = list_next(itr);
	list_iterator_destroy(itr);

	if (_get_bp_by_location(bg, pa_node->coord, &curr_bp) 
	    == SLURM_ERROR) {
		fatal("_get_bp_by_location()");
	}
	
	/* Set the one BP */
	if ((rc = rm_set_data(bg_record->bg_part,
			      RM_PartitionBPNum,
			      &bg_record->bp_count)) 
	    != STATUS_OK) {
		fatal("rm_set_data: RM_PartitionBPNum: %s", bg_err_str(rc));
		return SLURM_ERROR;
	}	
	if ((rc = rm_set_data(bg_record->bg_part,
			      RM_PartitionFirstBP, 
			      curr_bp)) 
	    != STATUS_OK) {
		fatal("rm_set_data("
		      "RM_PartitionFirstBP): %s", 
		      bg_err_str(rc));
		return SLURM_ERROR;
	}
	
	/* find the bp_id of the bp to get the nodecards */
	if ((rc = rm_get_data(curr_bp, RM_BPID, &bp_id))
	    != STATUS_OK) {
		error("rm_get_data(): %d", rc);
		return SLURM_ERROR;
	}

	if ((rc = rm_get_nodecards(bp_id, &ncard_list))
	    != STATUS_OK) {
		error("rm_get_nodecards(%s): %d",
		       bp_id, rc);
		return SLURM_ERROR;
	}
	
	if((rc = rm_get_data(ncard_list, RM_NodeCardListSize, &num))
	   != STATUS_OK) {
		error("rm_get_data(RM_NodeCardListSize): %s", bg_err_str(rc));
		return SLURM_ERROR;
	}
	num_ncards = 0;
	for(i=0; i<num; i++) {
		if (i) {
			if ((rc = rm_get_data(ncard_list, 
					      RM_NodeCardListNext, 
					      &ncard)) != STATUS_OK) {
				error("rm_get_data(RM_NodeCardListNext): %s",
				      rc);
				rc = SLURM_ERROR;
				goto cleanup;
			}
		} else {
			if ((rc = rm_get_data(ncard_list, 
					      RM_NodeCardListFirst, 
					      &ncard)) != STATUS_OK) {
				error("rm_get_data(RM_NodeCardListFirst: %s",
				      rc);
				rc = SLURM_ERROR;
				goto cleanup;
			}
		}
		
		if ((rc = rm_get_data(ncard, 
				      RM_NodeCardQuarter, 
				      &quarter)) != STATUS_OK) {
			error("rm_get_data(PartitionID): %d",rc);
			rc = SLURM_ERROR;
			goto cleanup;
		}
		if(bg_record->quarter != quarter)
			continue;
		if (num_ncards) {
			if ((rc = rm_set_data(bg_record->bg_part,
					      RM_PartitionNextNodeCard, 
					      ncard)) 
			    != STATUS_OK) {
				fatal("rm_set_data("
				      "RM_PartitionNextNodeCard): %s", 
				      bg_err_str(rc));
				
			}
		} else {
			if ((rc = rm_set_data(bg_record->bg_part,
					      RM_PartitionFirstNodeCard, 
					      ncard)) 
			    != STATUS_OK) {
				fatal("rm_set_data("
				      "RM_PartitionFirstNodeCard): %s", 
				      bg_err_str(rc));
			}
		}
		num_ncards++;
		if(num_ncards == 4)
			break;
	}
cleanup:
	if ((rc = rm_free_nodecard_list(ncard_list)) != STATUS_OK) {
		error("rm_free_nodecard_list(): %s", bg_err_str(rc));
		return SLURM_ERROR;
	}
	
	debug("making the small partition");
	return rc;
}

/**
 * connect the given switch up with the given connections
 */
extern int configure_partition_switches(bg_record_t * bg_record)
{
	int i, rc = SLURM_SUCCESS;
	ListIterator itr, switch_itr, bg_itr;
	pa_node_t* pa_node;
	char *name2;
	rm_BP_t *curr_bp;
	rm_switch_t *coord_switch[PA_SYSTEM_DIMENSIONS];
	rm_switch_t *curr_switch;
	pa_switch_t *pa_switch;
	char *bpid, *curr_bpid;
	int found_bpid = 0;
	int switch_count;
	bg_bp_t *bg_bp;
	bg_switch_t *bg_switch;
	int first_bp=1;
	int first_switch=1;
	
	bg_bp_list = list_create(NULL);
	bg_record->switch_count = 0;
	bg_record->bp_count = 0;
		
	itr = list_iterator_create(bg_record->bg_part_list);
	while ((pa_node = (pa_node_t *) list_next(itr)) != NULL) {
		debug("node %d%d%d",
		      pa_node->coord[X], 
		      pa_node->coord[Y], 
		      pa_node->coord[Z]);
		bg_itr = list_iterator_create(bg_bp_list);
		while((bg_bp = list_next(bg_itr)) != NULL) {
			if((bg_bp->coord[X] == pa_node->coord[X])
			   && (bg_bp->coord[Y] == pa_node->coord[Y])
			   && (bg_bp->coord[Z] == pa_node->coord[Z]))
				break;
		}
		list_iterator_destroy(bg_itr);
		
		if(bg_bp == NULL) {
			bg_bp = xmalloc(sizeof(bg_bp_t));
			bg_bp->coord = pa_node->coord;
			bg_bp->switch_list = list_create(NULL);
			list_append(bg_bp_list, bg_bp);
		}
		bg_record->bp_count++;
		bg_bp->used = 1;
		for(i=0;i<PA_SYSTEM_DIMENSIONS;i++) {
			
			pa_switch = &pa_node->axis_switch[i];
			if(pa_switch->int_wire[0].used) {
				_lookat_path(bg_bp, pa_switch, 0, 1, i);
			}
			
		/* 	if(pa_switch->int_wire[1].used) { */
/* 				_lookat_path(bg_bp, pa_switch, 1, 0, i); */
/* 			} */
		}
	}
	list_iterator_destroy(itr);
	
	bg_itr = list_iterator_create(bg_bp_list);
	while((bg_bp = list_next(bg_itr)) != NULL) {
		debug3("node %d%d%d",
		      bg_bp->coord[X], 
		      bg_bp->coord[Y], 
		      bg_bp->coord[Z]);
		itr = list_iterator_create(bg_bp->switch_list);
		while((bg_switch = list_next(itr)) != NULL) {
			bg_record->switch_count++;
		}
		list_iterator_destroy(itr);
	}
	list_iterator_destroy(bg_itr);

	if ((rc = rm_set_data(bg_record->bg_part,
			      RM_PartitionBPNum,
			      &bg_record->bp_count)) 
	    != STATUS_OK) {
		fatal("rm_set_data: RM_PartitionBPNum: %s", bg_err_str(rc));
		rc = SLURM_ERROR;
		goto cleanup;
	}
	debug3("BP count %d",bg_record->bp_count);
	if ((rc = rm_set_data(bg_record->bg_part,
			      RM_PartitionSwitchNum,
			      &bg_record->switch_count)) 
	    != STATUS_OK) {
		fatal("rm_set_data: RM_PartitionSwitchNum: %s", 
		      bg_err_str(rc));
		rc = SLURM_ERROR;
		goto cleanup;
	}
	debug3("switch count %d",bg_record->switch_count);
		
	first_bp = 1;
	first_switch = 1;
	
	if ((rc = rm_get_data(bg, RM_SwitchNum, &switch_count)) 
	    != STATUS_OK) {
		fatal("rm_get_data: RM_SwitchNum: %s", bg_err_str(rc));
		rc = SLURM_ERROR;
		goto cleanup;
	}
	
	bg_itr = list_iterator_create(bg_bp_list);
	while((bg_bp = list_next(bg_itr)) != NULL) {
			
		if (_get_bp_by_location(bg, bg_bp->coord, &curr_bp) 
		    == SLURM_ERROR) {
			list_iterator_destroy(bg_itr);
			rc = SLURM_ERROR;
			goto cleanup;
		}
		
		if(bg_bp->used) {
			if (first_bp){
				if ((rc = rm_set_data(bg_record->bg_part,
						      RM_PartitionFirstBP, 
						      curr_bp)) 
				    != STATUS_OK) {
					list_iterator_destroy(bg_itr);
					fatal("rm_set_data("
					      "RM_PartitionFirstBP): %s", 
					      bg_err_str(rc));
				}
				first_bp = 0;
			} else {
				if ((rc = rm_set_data(bg_record->bg_part,
						      RM_PartitionNextBP, 
						      curr_bp)) 
				    != STATUS_OK) {
					list_iterator_destroy(bg_itr);
					fatal("rm_set_data(RM_PartitionNextBP)"
					      ": %s", bg_err_str(rc));
				}
			}
		}

		if ((rc = rm_get_data(curr_bp,  RM_BPID, &bpid)) 
		    != STATUS_OK) {
			list_iterator_destroy(bg_itr);
			fatal("rm_get_data: RM_BPID: %s", bg_err_str(rc));
		}		

		if(!bpid) {
			error("No BP ID was returned from database");
			continue;
		}

		found_bpid = 0;
		for (i=0; i<switch_count; i++) {
			if(i) {
				if ((rc = rm_get_data(bg, RM_NextSwitch, 
						      &curr_switch)) 
				    != STATUS_OK) {
					list_iterator_destroy(bg_itr);
					fatal("rm_get_data: RM_NextSwitch: %s",
					      bg_err_str(rc));
				}
			} else {
				if ((rc = rm_get_data(bg, RM_FirstSwitch, 
						      &curr_switch)) 
				    != STATUS_OK) {
					list_iterator_destroy(bg_itr);
					fatal("rm_get_data: "
					      "RM_FirstSwitch: %s",
					      bg_err_str(rc));
				}
			}
			if ((rc = rm_get_data(curr_switch, RM_SwitchBPID, 
					      &curr_bpid)) != STATUS_OK) {
				list_iterator_destroy(bg_itr);
				fatal("rm_get_data: RM_SwitchBPID: %s", 
				      bg_err_str(rc));
			}

			if(!curr_bpid) {
				error("No BP ID was returned from database");
				continue;
			}

			if (!strcasecmp((char *)bpid, (char *)curr_bpid)) {
				coord_switch[found_bpid] = curr_switch;
				found_bpid++;
				if(found_bpid==PA_SYSTEM_DIMENSIONS) {
					free(curr_bpid);
					break;
				}
			}
			free(curr_bpid);
		}

		free(bpid);

		if(found_bpid==PA_SYSTEM_DIMENSIONS) {
						
			debug2("adding midplane %d%d%d",
			       bg_bp->coord[X],
			       bg_bp->coord[Y],
			       bg_bp->coord[Z]);
			switch_itr = list_iterator_create(bg_bp->switch_list);
			while((bg_switch = list_next(switch_itr)) != NULL) {
				
				debug2("adding switch dim %d",
				       bg_switch->dim);
				     
				if (_add_switch_conns(coord_switch
						      [bg_switch->dim],
						      bg_switch) 
				    == SLURM_ERROR) {
					list_iterator_destroy(switch_itr);
					list_iterator_destroy(bg_itr);
					rc = SLURM_ERROR;
					goto cleanup;
				}
				
				if (first_switch){
					if ((rc = rm_set_data(
						     bg_record->bg_part,
						     RM_PartitionFirstSwitch,
						     coord_switch
						     [bg_switch->dim])) 
					    != STATUS_OK) {
						list_iterator_destroy(
							switch_itr);
						list_iterator_destroy(bg_itr);
						fatal("rm_set_data("
						      "RM_PartitionFirst"
						      "Switch): %s", 
						      bg_err_str(rc));
					}
					
					first_switch = 0;
				} else {
					if ((rc = rm_set_data(
						     bg_record->bg_part,
						     RM_PartitionNextSwitch,
						     coord_switch
						     [bg_switch->dim])) 
					    != STATUS_OK) {
						list_iterator_destroy(
							switch_itr);
						list_iterator_destroy(bg_itr);
						fatal("rm_set_data("
						      "RM_PartitionNext"
						      "Switch:) %s", 
						      bg_err_str(rc));
					}
				}
			}
			list_iterator_destroy(switch_itr);
		}
	}
	rc = SLURM_SUCCESS;
cleanup:
	if (_destroy_bg_bp_list(bg_bp_list) == SLURM_ERROR)
		return SLURM_ERROR;	
	
	return rc;	
}


#endif
