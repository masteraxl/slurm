/*****************************************************************************\
 *  block_allocator.c - Assorted functions for layout of bglblocks, 
 *	 wiring, mapping for smap, etc.
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Dan Phung <phung4@llnl.gov>, Danny Auble <da@llnl.gov>
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "block_allocator.h"
#include "src/common/uid.h"

#define DEBUG_PA
#define BEST_COUNT_INIT 20

/* Global */
bool _initialized = false;
bool _wires_initialized = false;
bool _bp_map_initialized = false;

/* _ba_system is the "current" system that the structures will work
 *  on */
ba_system_t *ba_system_ptr = NULL;
List path = NULL;
List best_path = NULL;
int best_count;
int color_count = 0;
bool *passthrough = NULL;

/* extern Global */
List bp_map_list = NULL;
char letters[62];
char colors[6];
#ifdef HAVE_BG
int DIM_SIZE[BA_SYSTEM_DIMENSIONS] = {0,0,0};
#else
int DIM_SIZE[BA_SYSTEM_DIMENSIONS] = {0};
#endif

s_p_options_t bg_conf_file_options[] = {
	{"BlrtsImage", S_P_STRING}, 
	{"LinuxImage", S_P_STRING},
	{"MloaderImage", S_P_STRING},
	{"RamDiskImage", S_P_STRING},
	{"BridgeAPILogFile", S_P_STRING},
	{"RamDiskImage", S_P_STRING},
	{"LayoutMode", S_P_STRING},
	{"BridgeAPIVerbose", S_P_UINT16},
	{"BasePartitionNodeCnt", S_P_UINT16},
	{"NodeCardNodeCnt", S_P_UINT16},
	{"Numpsets", S_P_UINT16},
	{"BPs", S_P_ARRAY, parse_blockreq, destroy_blockreq},
	/* these are just going to be put into a list that will be
	   freed later don't free them after reading them */
	{"AltBlrtsImage", S_P_ARRAY, parse_image, NULL}, 
	{"AltLinuxImage", S_P_ARRAY, parse_image, NULL},
	{"AltMloaderImage", S_P_ARRAY, parse_image, NULL},
	{"AltRamDiskImage", S_P_ARRAY, parse_image, NULL},
	{NULL}
};

#ifdef HAVE_BG
/** internal helper functions */
#ifdef HAVE_BG_FILES
/** */
static void _bp_map_list_del(void *object);
static int _port_enum(int port);
#endif
/* */
static int _check_for_options(ba_request_t* ba_request); 
/* */
static int _append_geo(int *geo, List geos, int rotate);
/* */
static int _fill_in_coords(List results, List start_list, 
			   int *geometry, int conn_type);
/* */
static int _copy_the_path(List nodes, ba_switch_t *curr_switch,
			  ba_switch_t *mark_switch, 
			  int source, int dim);
/* */
static int _find_yz_path(ba_node_t *ba_node, int *first, 
			 int *geometry, int conn_type);
#endif

#ifndef HAVE_BG_FILES
#ifdef HAVE_BG
/* */
static int _create_config_even(ba_node_t ***grid);
#else
/* */
static int _create_config_even(ba_node_t *grid);
#endif
#endif

/** */
static void _new_ba_node(ba_node_t *ba_node, int *coord,
			 bool track_down_nodes);
/** */
static int _reset_the_path(ba_switch_t *curr_switch, int source, 
			   int target, int dim);
/** */
static void _create_ba_system(void);
/* */
static void _delete_ba_system(void);
/* */
static void _delete_path_list(void *object);

/* find the first block match in the system */
static int _find_match(ba_request_t* ba_request, List results);

static bool _node_used(ba_node_t* ba_node, int *geometry);

/* */
static void _switch_config(ba_node_t* source, ba_node_t* target, int dim, 
			   int port_src, int port_tar);
/* */
static int _set_external_wires(int dim, int count, ba_node_t* source, 
				ba_node_t* target);
/* */
static char *_set_internal_wires(List nodes, int size, int conn_type);
/* */
static int _find_x_path(List results, ba_node_t *ba_node, 
			int *start, int *first, 
			int *geometry, int found, int conn_type);
/* */
static int _find_x_path2(List results, ba_node_t *ba_node, 
			 int *start, int *first, 
			 int *geometry, int found, int conn_type);
/* */
static int _remove_node(List results, int *node_tar);
/* */
static int _find_next_free_using_port_2(ba_switch_t *curr_switch, 
					int source_port, 
					List nodes, int dim, 
					int count);
/* */
static int _find_passthrough(ba_switch_t *curr_switch, int source_port, 
			     List nodes, int dim, 
			     int count, int highest_phys_x); 
/* */
static int _finish_torus(ba_switch_t *curr_switch, int source_port, 
			   List nodes, int dim, int count, int *start);
/* */
static int *_set_best_path();

/* */
static int _set_one_dim(int *start, int *end, int *coord);

/* */
static void _destroy_geo(void *object);

extern char *bg_block_state_string(rm_partition_state_t state)
{
	static char tmp[16];

#ifdef HAVE_BG
	switch (state) {
		case RM_PARTITION_BUSY: 
			return "BUSY";
		case RM_PARTITION_CONFIGURING:
			return "CONFIG";
		case RM_PARTITION_DEALLOCATING:
			return "DEALLOC";
		case RM_PARTITION_ERROR:
			return "ERROR";
		case RM_PARTITION_FREE:
			return "FREE";
		case RM_PARTITION_NAV:
			return "NAV";
		case RM_PARTITION_READY:
			return "READY";
	}
#endif

	snprintf(tmp, sizeof(tmp), "%d", state);
	return tmp;
}

extern int parse_blockreq(void **dest, slurm_parser_enum_t type,
			  const char *key, const char *value, 
			  const char *line, char **leftover)
{
	s_p_options_t block_options[] = {
		{"Type", S_P_STRING},
		{"Nodecards", S_P_UINT16},
		{"Quarters", S_P_UINT16},
		{"BlrtsImage", S_P_STRING},
		{"LinuxImage", S_P_STRING},
		{"MloaderImage", S_P_STRING},
		{"RamDiskImage", S_P_STRING},
		{NULL}
	};
	s_p_hashtbl_t *tbl;
	char *tmp = NULL;
	blockreq_t *n = NULL;
	hostlist_t hl = NULL;
	char temp[BUFSIZE];
	tbl = s_p_hashtbl_create(block_options);
	s_p_parse_line(tbl, *leftover, leftover);
	if(!value) {
		return 0;
	}
	n = xmalloc(sizeof(blockreq_t));
	hl = hostlist_create(value);
	hostlist_ranged_string(hl, BUFSIZE, temp);
	hostlist_destroy(hl);

	n->block = xstrdup(temp);
	s_p_get_string(&n->blrtsimage, "BlrtsImage", tbl);
	s_p_get_string(&n->linuximage, "LinuxImage", tbl);
	s_p_get_string(&n->mloaderimage, "MloaderImage", tbl);
	s_p_get_string(&n->ramdiskimage, "RamDiskImage", tbl);
	
	s_p_get_string(&tmp, "Type", tbl);
	if (!tmp || !strcasecmp(tmp,"TORUS"))
		n->conn_type = SELECT_TORUS;
	else if(!strcasecmp(tmp,"MESH"))
		n->conn_type = SELECT_MESH;
	else
		n->conn_type = SELECT_SMALL;
	xfree(tmp);
	
	if (!s_p_get_uint16(&n->nodecards, "Nodecards", tbl))
		n->nodecards = 0;
	if (!s_p_get_uint16(&n->quarters, "Quarters", tbl))
		n->quarters = 0;

	s_p_hashtbl_destroy(tbl);

	*dest = (void *)n;
	return 1;
}

extern void destroy_blockreq(void *ptr)
{
	blockreq_t *n = (blockreq_t *)ptr;
	if(n) {
		xfree(n->block);
		xfree(n->blrtsimage);
		xfree(n->linuximage);
		xfree(n->mloaderimage);
		xfree(n->ramdiskimage);
		xfree(n);
	}
}

extern int parse_image(void **dest, slurm_parser_enum_t type,
		       const char *key, const char *value, 
		       const char *line, char **leftover)
{
	s_p_options_t image_options[] = {
		{"GROUPS", S_P_STRING},
		{NULL}
	};
	s_p_hashtbl_t *tbl = NULL;
	char *tmp = NULL;
	image_t *n = NULL;
	image_group_t *image_group = NULL;
	int i = 0, j = 0;

	tbl = s_p_hashtbl_create(image_options);
	s_p_parse_line(tbl, *leftover, leftover);
	
	n = xmalloc(sizeof(image_t));
	n->name = xstrdup(value);
	n->def = false;
	debug3("image %s", n->name);
	n->groups = list_create(destroy_image_group_list);
	s_p_get_string(&tmp, "Groups", tbl);
	if(tmp) {
		for(i=0; i<strlen(tmp); i++) {
			if(tmp[i] == ':') {
				image_group = xmalloc(sizeof(image_group));
				image_group->name = xmalloc(i-j+2);
				snprintf(image_group->name,
					 (i-j)+1, "%s", tmp+j);
				image_group->gid =
					gid_from_string(image_group->name);
				debug3("adding group %s %d", image_group->name,
				       image_group->gid);
				list_append(n->groups, image_group);
				j=i;
				j++;
			} 		
		}
		if(j != i) {
			image_group = xmalloc(sizeof(image_group));
			image_group->name = xmalloc(i-j+2);
			snprintf(image_group->name, (i-j)+1, "%s", tmp+j);
			image_group->gid = gid_from_string(image_group->name);
			debug3("adding group %s %d", image_group->name,
			       image_group->gid);
			list_append(n->groups, image_group);
		}
		xfree(tmp);
	}
	s_p_hashtbl_destroy(tbl);

	*dest = (void *)n;
	return 1;
}

extern void destroy_image_group_list(void *ptr)
{
	image_group_t *image_group = (image_group_t *)ptr;
	if(image_group) {
		xfree(image_group->name);
		xfree(image_group);
	}
}

extern void destroy_image(void *ptr)
{
	image_t *n = (image_t *)ptr;
	if(n) {
		xfree(n->name);
		if(n->groups) {
			list_destroy(n->groups);
			n->groups = NULL;
		}
		xfree(n);
	}
}

extern void destroy_ba_node(void *ptr)
{
	ba_node_t *ba_node = (ba_node_t *)ptr;
	if(ba_node) {
		xfree(ba_node);
	}
}

/**
 * create a block request.  Note that if the geometry is given,
 * then size is ignored.  
 * 
 * IN/OUT - ba_request: structure to allocate and fill in.  
 * return SUCCESS of operation.
 */
extern int new_ba_request(ba_request_t* ba_request)
{
	int i=0;
#ifdef HAVE_BG
	float sz=1;
	int geo[BA_SYSTEM_DIMENSIONS] = {0,0,0};
	int i2, i3, picked, total_sz=1 , size2, size3;
	ListIterator itr;
	int checked[8];
	int *geo_ptr;
	int messed_with = 0;
	
	ba_request->save_name= NULL;
	ba_request->rotate_count= 0;
	ba_request->elongate_count = 0;
	ba_request->elongate_geos = list_create(_destroy_geo);
	geo[X] = ba_request->geometry[X];
	geo[Y] = ba_request->geometry[Y];
	geo[Z] = ba_request->geometry[Z];
	passthrough = &ba_request->passthrough;

	if(geo[X] != (uint16_t)NO_VAL) { 
		for (i=0; i<BA_SYSTEM_DIMENSIONS; i++){
			if ((geo[i] < 1) 
			    ||  (geo[i] > DIM_SIZE[i])){
				error("new_ba_request Error, "
				      "request geometry is invalid %d "
				      "DIMS are %d%d%d", 
				      geo[i],
				      DIM_SIZE[X],
				      DIM_SIZE[Y],
				      DIM_SIZE[Z]); 
				return 0;
			}
		}
		_append_geo(geo, ba_request->elongate_geos, 0);
		sz=1;
		for (i=0; i<BA_SYSTEM_DIMENSIONS; i++)
			sz *= ba_request->geometry[i];
		ba_request->size = sz;
		sz=0;
	}
	
	if(ba_request->elongate || sz) {
		sz=1;
		/* decompose the size into a cubic geometry */
		ba_request->rotate= 1;
		ba_request->elongate = 1;		
		
		for (i=0; i<BA_SYSTEM_DIMENSIONS; i++) { 
			total_sz *= DIM_SIZE[i];
			geo[i] = 1;
		}
		
		if(ba_request->size==1) {
			_append_geo(geo, 
				    ba_request->elongate_geos,
				    ba_request->rotate);
			goto endit;
		}
		
		if(ba_request->size<=DIM_SIZE[Y]) {
			geo[X] = 1;
			geo[Y] = ba_request->size;
			geo[Z] = 1;
			sz=ba_request->size;
			_append_geo(geo, 
				    ba_request->elongate_geos, 
				    ba_request->rotate);
		}
		
		i = ba_request->size/4;
		if(!(ba_request->size%2)
		   && i <= DIM_SIZE[Y]
		   && i <= DIM_SIZE[Z]
		   && i*i == ba_request->size) {
			geo[X] = 1;
			geo[Y] = i;
			geo[Z] = i;
			sz=ba_request->size;
			_append_geo(geo,
				    ba_request->elongate_geos,
				    ba_request->rotate);
		}
	
		if(ba_request->size>total_sz || ba_request->size<1) {
			return 0;			
		}
		sz = ba_request->size % (DIM_SIZE[Y] * DIM_SIZE[Z]);
		if(!sz) {
		      i = ba_request->size / (DIM_SIZE[Y] * DIM_SIZE[Z]);
		      geo[X] = i;
		      geo[Y] = DIM_SIZE[Y];
		      geo[Z] = DIM_SIZE[Z];
		      sz=ba_request->size;
		      _append_geo(geo,
				  ba_request->elongate_geos,
				  ba_request->rotate);
		}	
	startagain:		
		picked=0;
		for(i=0;i<8;i++)
			checked[i]=0;
		
		size3=ba_request->size;
		
		for (i=0; i<BA_SYSTEM_DIMENSIONS; i++) {
			total_sz *= DIM_SIZE[i];
			geo[i] = 1;
		}
	
		sz = 1;
		size3=ba_request->size;
		picked=0;
	tryagain:	
		if(size3!=ba_request->size)
			size2=size3;
		else
			size2=ba_request->size;
		//messedup:

		for (i=picked; i<BA_SYSTEM_DIMENSIONS; i++) { 
			if(size2<=1) 
				break;
			sz = size2%DIM_SIZE[i];
			if(!sz) {
				geo[i] = DIM_SIZE[i];	
				size2 /= DIM_SIZE[i];
			} else if (size2 > DIM_SIZE[i]){
				for(i2=(DIM_SIZE[i]-1);i2>1;i2--) {
					/* go through each number to see if 
					   the size is divisable by a smaller 
					   number that is 
					   good in the other dims. */
					if (!(size2%i2) && !checked[i2]) {
						size2 /= i2;
									
						if(i==0)
							checked[i2]=1;
							
						if(i2<DIM_SIZE[i]) 
							geo[i] = i2;
						else {
							goto tryagain;
						}
						if((i2-1)!=1 && 
						   i!=(BA_SYSTEM_DIMENSIONS-1))
							break;
					}		
				}				
				if(i2==1) {
					ba_request->size +=1;
					goto startagain;
				}
						
			} else {
				geo[i] = sz;	
				break;
			}					
		}
		
		if((geo[X]*geo[Y]) <= DIM_SIZE[Y]) {
			ba_request->geometry[X] = 1;
			ba_request->geometry[Y] = geo[X] * geo[Y];
			ba_request->geometry[Z] = geo[Z];
			_append_geo(ba_request->geometry, 
				    ba_request->elongate_geos, 
				    ba_request->rotate);

		}
		if((geo[X]*geo[Z]) <= DIM_SIZE[Y]) {
			ba_request->geometry[X] = 1;
			ba_request->geometry[Y] = geo[Y];
			ba_request->geometry[Z] = geo[X] * geo[Z];	
			_append_geo(ba_request->geometry, 
				    ba_request->elongate_geos, 
				    ba_request->rotate);		
	
		}
		if((geo[X]/2) <= DIM_SIZE[Y]) {
			if(geo[Y] == 1) {
				ba_request->geometry[Y] = geo[X]/2;
				messed_with = 1;
			} else  
				ba_request->geometry[Y] = geo[Y];
			if(!messed_with && geo[Z] == 1) {
				messed_with = 1;
				ba_request->geometry[Z] = geo[X]/2;
			} else
				ba_request->geometry[Z] = geo[Z];
			if(messed_with) {
				messed_with = 0;
				ba_request->geometry[X] = 2;
				_append_geo(ba_request->geometry, 
					    ba_request->elongate_geos, 
					    ba_request->rotate);
			}
		}
		if(geo[X] == DIM_SIZE[X]
		   && (geo[Y] < DIM_SIZE[Y] 
		       || geo[Z] < DIM_SIZE[Z])) {
			if(DIM_SIZE[Y]<DIM_SIZE[Z]) {
				i = DIM_SIZE[Y];
				DIM_SIZE[Y] = DIM_SIZE[Z];
				DIM_SIZE[Z] = i;
			}
			ba_request->geometry[X] = geo[X];
			ba_request->geometry[Y] = geo[Y];
			ba_request->geometry[Z] = geo[Z];
			if(ba_request->geometry[Y] < DIM_SIZE[Y]) {
				i = (DIM_SIZE[Y] - ba_request->geometry[Y]);
				ba_request->geometry[Y] +=i;
			}
			if(ba_request->geometry[Z] < DIM_SIZE[Z]) {
				i = (DIM_SIZE[Z] - ba_request->geometry[Z]);
				ba_request->geometry[Z] +=i;
			}
			for(i = DIM_SIZE[X]; i>0; i--) {
				ba_request->geometry[X]--;
				i2 = (ba_request->geometry[X]
				      * ba_request->geometry[Y]
				      * ba_request->geometry[Z]);
				if(i2 < ba_request->size) {
					ba_request->geometry[X]++;
					messed_with = 1;
					break;
				}					
			}			
			if(messed_with) {
				messed_with = 0;
				_append_geo(ba_request->geometry, 
					    ba_request->elongate_geos, 
					    ba_request->rotate);
			}
		}
		
		_append_geo(geo, 
			    ba_request->elongate_geos, 
			    ba_request->rotate);
	
		/* see if We can find a cube or square root of the 
		   size to make an easy cube */
		for(i=0;i<BA_SYSTEM_DIMENSIONS-1;i++) {
			sz = powf((float)ba_request->size,
				  (float)1/(BA_SYSTEM_DIMENSIONS-i));
			if(pow(sz,(BA_SYSTEM_DIMENSIONS-i))==ba_request->size)
				break;
		}
	
		if(i<BA_SYSTEM_DIMENSIONS-1) {
			/* we found something that looks like a cube! */
			i3=i;
			for (i=0; i<i3; i++) 
				geo[i] = 1;			
			
			for (i=i3; i<BA_SYSTEM_DIMENSIONS; i++)  
				if(sz<=DIM_SIZE[i]) 
					geo[i] = sz;	
				else
					goto endit;
				
			_append_geo(geo, 
				    ba_request->elongate_geos, 
				    ba_request->rotate);
		} 
	}
	
endit:
	itr = list_iterator_create(ba_request->elongate_geos);
	geo_ptr = list_next(itr);
	list_iterator_destroy(itr);
	
	if(geo_ptr == NULL)
		return 0;

	ba_request->elongate_count++;
	ba_request->geometry[X] = geo_ptr[X];
	ba_request->geometry[Y] = geo_ptr[Y];
	ba_request->geometry[Z] = geo_ptr[Z];
	sz=1;
	for (i=0; i<BA_SYSTEM_DIMENSIONS; i++)
		sz *= ba_request->geometry[i];
	ba_request->size = sz;
	
#else
	int geo[BA_SYSTEM_DIMENSIONS] = {0};
	
	ba_request->rotate_count= 0;
	ba_request->elongate_count = 0;
	ba_request->elongate_geos = list_create(_destroy_geo);
	geo[X] = ba_request->geometry[X];
		
	if(geo[X] != NO_VAL) { 
		for (i=0; i<BA_SYSTEM_DIMENSIONS; i++){
			if ((geo[i] < 1) 
			    ||  (geo[i] > DIM_SIZE[i])){
				error("new_ba_request Error, "
				      "request geometry is invalid %d", 
				      geo[i]); 
				return 0;
			}
		}
		
		ba_request->size = ba_request->geometry[X];

	} else if (ba_request->size) {
		ba_request->geometry[X] = ba_request->size;
	} else
		return 0;
			
#endif
	return 1;
}

/**
 * delete a block request 
 */
extern void delete_ba_request(void *arg)
{
	ba_request_t *ba_request = (ba_request_t *)arg;
	if(ba_request) {
		xfree(ba_request->save_name);
		if(ba_request->elongate_geos)
			list_destroy(ba_request->elongate_geos);
		xfree(ba_request->blrtsimage);
		xfree(ba_request->linuximage);
		xfree(ba_request->mloaderimage);
		xfree(ba_request->ramdiskimage);
		
		xfree(ba_request);
	}
}

/**
 * print a block request 
 */
extern void print_ba_request(ba_request_t* ba_request)
{
	int i;

	if (ba_request == NULL){
		error("print_ba_request Error, request is NULL");
		return;
	}
	debug("  ba_request:");
	debug("    geometry:\t");
	for (i=0; i<BA_SYSTEM_DIMENSIONS; i++){
		debug("%d", ba_request->geometry[i]);
	}
	debug("        size:\t%d", ba_request->size);
	debug("   conn_type:\t%d", ba_request->conn_type);
	debug("      rotate:\t%d", ba_request->rotate);
	debug("    elongate:\t%d", ba_request->elongate);
}

/**
 * empty a list that we don't want to destroy the memory of the
 * elements always returns 1
*/
extern int empty_null_destroy_list(void *arg, void *key)
{
	return 1;
}

/**
 * Initialize internal structures by either reading previous block
 * configurations from a file or by running the graph solver.
 * 
 * IN: node_info_msg_t * can be null, 
 *     should be from slurm_load_node().
 * 
 * return: void.
 */
extern void ba_init(node_info_msg_t *node_info_ptr)
{
	int x,y,z;

#ifdef HAVE_BG
	node_info_t *node_ptr = NULL;
	int start, temp;
	char *numeric = NULL;
	int i, j=0;
	slurm_conf_node_t *node = NULL, **ptr_array;
	int count, number;
	int end[BA_SYSTEM_DIMENSIONS];
	
#ifdef HAVE_BG_FILES
	rm_BGL_t *bg = NULL;
	rm_size3D_t bp_size;
	int rc = 0;
#endif /* HAVE_BG_FILES */
	
#endif /* HAVE_BG */

	/* We only need to initialize once, so return if already done so. */
	if (_initialized){
		return;
	}
	
#ifdef HAVE_BG_FILES
	bridge_init();
#endif	
	y = 65;
	for (x = 0; x < 62; x++) {
		if (y == 91)
			y = 97;
		else if(y == 123)
			y = 48;
		else if(y == 58)
			y = 65;
		letters[x] = y;
		y++;
	}

	z=1;
	for (x = 0; x < 6; x++) {
		if(z == 4)
			z++;
		colors[x] = z;				
		z++;
	}
		
	best_count=BEST_COUNT_INIT;

	if(ba_system_ptr)
		_delete_ba_system();
	
	ba_system_ptr = (ba_system_t *) xmalloc(sizeof(ba_system_t));
	
	ba_system_ptr->xcord = 1;
	ba_system_ptr->ycord = 1;
	ba_system_ptr->num_of_proc = 0;
	ba_system_ptr->resize_screen = 0;
	
	if(node_info_ptr!=NULL) {
#ifdef HAVE_BG
		for (i = 0; i < node_info_ptr->record_count; i++) {
			node_ptr = &node_info_ptr->node_array[i];
			start = 0;
			
			if(!node_ptr->name) {
				DIM_SIZE[X] = 0;
				DIM_SIZE[Y] = 0;
				DIM_SIZE[Z] = 0;
				goto node_info_error;
			}

			numeric = node_ptr->name;
			while (numeric) {
				if ((numeric[0] < '0')
				||  (numeric[0] > '9')) {
					numeric++;
					continue;
				}
				start = atoi(numeric);
				break;
			}
			
			temp = start / 100;
			if (DIM_SIZE[X] < temp)
				DIM_SIZE[X] = temp;
			temp = (start / 10) % 10;
			if (DIM_SIZE[Y] < temp)
				DIM_SIZE[Y] = temp;
			temp = start % 10;
			if (DIM_SIZE[Z] < temp)
				DIM_SIZE[Z] = temp;
		}
		DIM_SIZE[X]++;
		DIM_SIZE[Y]++;
		DIM_SIZE[Z]++;
#else
		DIM_SIZE[X] = node_info_ptr->record_count;
#endif
		ba_system_ptr->num_of_proc = node_info_ptr->record_count;
	} 
#ifdef HAVE_BG
node_info_error:

#ifdef HAVE_BG_FILES
	if (have_db2
	    && ((DIM_SIZE[X]==0) || (DIM_SIZE[Y]==0) || (DIM_SIZE[Z]==0))) {
		if ((rc = bridge_get_bg(&bg)) != STATUS_OK) {
			error("bridge_get_BGL(): %d", rc);
			return;
		}
		
		if ((bg != NULL)
		&&  ((rc = bridge_get_data(bg, RM_Msize, &bp_size)) 
		     == STATUS_OK)) {
			DIM_SIZE[X]=bp_size.X;
			DIM_SIZE[Y]=bp_size.Y;
			DIM_SIZE[Z]=bp_size.Z;
		} else {
			error("bridge_get_data(RM_Msize): %d", rc);	
		}
		if ((rc = bridge_free_bg(bg)) != STATUS_OK)
			error("bridge_free_BGL(): %d", rc);
	}
#endif

	if ((DIM_SIZE[X]==0) || (DIM_SIZE[Y]==0) || (DIM_SIZE[Z]==0)) {
		debug("Setting dimensions from slurm.conf file");
		count = slurm_conf_nodename_array(&ptr_array);
		if (count == 0)
			fatal("No NodeName information available!");
		
		for (i = 0; i < count; i++) {
			node = ptr_array[i];
			j = 0;
			while (node->nodenames[j] != '\0') {
				if ((node->nodenames[j] == '['
				     || node->nodenames[j] == ',')
				    && (node->nodenames[j+8] == ']' 
					|| node->nodenames[j+8] == ',')
				    && (node->nodenames[j+4] == 'x'
					|| node->nodenames[j+4] == '-')) {
					j+=5;
				} else if((node->nodenames[j] < 58 
					   && node->nodenames[j] > 47)) {
				} else {
					j++;
					continue;
				}
				number = strtol(node->nodenames + j,
						NULL, BG_BASE);
				
				end[X] = number / (BG_BASE * BG_BASE);
				end[Y] = (number % (BG_BASE * BG_BASE))
					/ BG_BASE;
				end[Z] = (number % BG_BASE);
				DIM_SIZE[X] = MAX(DIM_SIZE[X], end[X]);
				DIM_SIZE[Y] = MAX(DIM_SIZE[Y], end[Y]);
				DIM_SIZE[Z] = MAX(DIM_SIZE[Z], end[Z]);
				break;
			}
				
		}
		if ((DIM_SIZE[X]==0) && (DIM_SIZE[Y]==0) && (DIM_SIZE[Z]==0)) 
			info("are you sure you only have 1 midplane? %s",
			      node->nodenames);
		DIM_SIZE[X]++;
		DIM_SIZE[Y]++;
		DIM_SIZE[Z]++;
	}
	debug("DIM_SIZE = %dx%dx%d", DIM_SIZE[X], DIM_SIZE[Y], DIM_SIZE[Z]);
	
#else 
	if (DIM_SIZE[X]==0) {
		debug("Setting default system dimensions");
		DIM_SIZE[X]=100;
	}	
#endif
	if(!ba_system_ptr->num_of_proc)
		ba_system_ptr->num_of_proc = 
			DIM_SIZE[X] 
#ifdef HAVE_BG
			* DIM_SIZE[Y] 
			* DIM_SIZE[Z]
#endif 
			;

	_create_ba_system();
	
#ifndef HAVE_BG_FILES
	_create_config_even(ba_system_ptr->grid);
#endif
	path = list_create(_delete_path_list);
	best_path = list_create(_delete_path_list);

	_initialized = true;
	init_grid(node_info_ptr);
}

extern void init_wires()
{
	int x, y, z, i;
	ba_node_t *source = NULL;
	if(_wires_initialized)
		return;

	for(x=0;x<DIM_SIZE[X];x++) {
		for(y=0;y<DIM_SIZE[Y];y++) {
			for(z=0;z<DIM_SIZE[Z];z++) {
#ifdef HAVE_BG
				source = &ba_system_ptr->grid[x][y][z];
#else
				source = &ba_system_ptr->grid[x];
#endif
				for(i=0; i<6; i++) {
					_switch_config(source, source, 
						       X, i, i);
					_switch_config(source, source, 
						       Y, i, i);
					_switch_config(source, source, 
						       Z, i, i);
				}
			}
		}
	}
#ifdef HAVE_BG_FILES	
	_set_external_wires(0,0,NULL,NULL);
	if(!bp_map_list) {
		if(set_bp_map() == -1) {
			return;
		}
	}
#endif
	
	_wires_initialized = true;
	return;
}


/** 
 * destroy all the internal (global) data structs.
 */
extern void ba_fini()
{
	if (!_initialized){
		return;
	}

	if (path) {
		list_destroy(path);
		path = NULL;
	}
	if (best_path) {
		list_destroy(best_path);
		best_path = NULL;
	}
#ifdef HAVE_BG_FILES
	if (bp_map_list) {
		list_destroy(bp_map_list);
		bp_map_list = NULL;
		_bp_map_initialized = false;
	}
	bridge_fini();
#endif
	_delete_ba_system();
//	debug2("pa system destroyed");
}


/** 
 * set the node in the internal configuration as unusable
 * 
 * IN ba_node: ba_node_t to put down
 */
extern void ba_update_node_state(ba_node_t *ba_node, uint16_t state)
{
	uint16_t node_base_state = state & NODE_STATE_BASE;

	if (!_initialized){
		error("Error, configuration not initialized, "
		      "calling ba_init(NULL)");
		ba_init(NULL);
	}

#ifdef HAVE_BG
	debug("ba_update_node_state: new state of node[%d%d%d] is %s", 
	      ba_node->coord[X], ba_node->coord[Y], ba_node->coord[Z],
	      node_state_string(state)); 
#else
	debug("ba_update_node_state: new state of node[%d] is %s", 
	      ba_node->coord[X],
	      node_state_string(state)); 
#endif

	/* basically set the node as used */
	if((node_base_state == NODE_STATE_DOWN)
	   || (ba_node->state & NODE_STATE_DRAIN)) 
		ba_node->used = true;
	else
		ba_node->used = false;
	ba_node->state = state;
}
/** 
 * copy info from a ba_node
 * 
 * IN ba_node: node to be copied
 * OUT ba_node_t *: copied info must be freed with destroy_ba_node
 */
extern ba_node_t *ba_copy_node(ba_node_t *ba_node)
{
	ba_node_t *new_ba_node = xmalloc(sizeof(ba_node_t));
	
	memcpy(new_ba_node, ba_node, sizeof(ba_node_t));
	return new_ba_node;
}
/** 
 * Try to allocate a block.
 * 
 * IN - ba_request: allocation request
 * OUT - results: List of results of the allocation request.  Each
 * list entry will be a coordinate.  allocate_block will create the
 * list, but the caller must destroy it.
 * 
 * return: success or error of request
 */
extern int allocate_block(ba_request_t* ba_request, List results)
{

	if (!_initialized){
		error("Error, configuration not initialized, "
		      "calling ba_init(NULL)");
	}

	if (!ba_request){
		error("allocate_block Error, request not initialized");
		return 0;
	}
	
	// _backup_ba_system();
	if (_find_match(ba_request, results)){
		return 1;
	} else {
		return 0;
	}
}


/** 
 * Doh!  Admin made a boo boo.  
 *
 * returns SLURM_SUCCESS if undo was successful.
 */
extern int remove_block(List nodes, int new_count)
{
	int dim;
	ba_node_t* ba_node = NULL;
	ba_switch_t *curr_switch = NULL; 
	ListIterator itr;
	
	itr = list_iterator_create(nodes);
	while((ba_node = (ba_node_t*) list_next(itr)) != NULL) {
		ba_node->used = false;
		ba_node->color = 7;
		ba_node->letter = '.';
		for(dim=0;dim<BA_SYSTEM_DIMENSIONS;dim++) {		
			curr_switch = &ba_node->axis_switch[dim];
			if(curr_switch->int_wire[0].used) {
				_reset_the_path(curr_switch, 0, 1, dim);
			}
		}
	}
	list_iterator_destroy(itr);
	if(new_count == -1)
		color_count--;
	else
		color_count=new_count;			
	if(color_count < 0)
		color_count = 0;
	return 1;
}

/** 
 * Doh!  Admin made a boo boo.  Note: Undo only has one history
 * element, so two consecutive undo's will fail.
 *
 * returns SLURM_SUCCESS if undo was successful.
 */
extern int alter_block(List nodes, int conn_type)
{
	/* int dim; */
/* 	ba_node_t* ba_node = NULL; */
/* 	ba_switch_t *curr_switch = NULL;  */
/* 	int size=0; */
/* 	char *name = NULL; */
/* 	ListIterator results_i;	 */	

	return SLURM_ERROR;
	/* results_i = list_iterator_create(nodes); */
/* 	while ((ba_node = list_next(results_i)) != NULL) { */
/* 		ba_node->used = false; */
		
/* 		for(dim=0;dim<BA_SYSTEM_DIMENSIONS;dim++) { */
/* 			curr_switch = &ba_node->axis_switch[dim]; */
/* 			if(curr_switch->int_wire[0].used) { */
/* 				_reset_the_path(curr_switch, 0, 1, dim); */
/* 			} */
/* 		} */
/* 		size++; */
/* 	} */
/* 	list_iterator_destroy(results_i); */
/* 	if((name = _set_internal_wires(nodes, size, conn_type)) == NULL) */
/* 		return SLURM_ERROR; */
/* 	else { */
/* 		xfree(name); */
/* 		return SLURM_SUCCESS; */
/* 	} */
}

/** 
 * After a block is deleted or altered following allocations must
 * be redone to make sure correct path will be used in the real system
 *
 */
extern int redo_block(List nodes, int *geo, int conn_type, int new_count)
{
       	ba_node_t* ba_node;
	char *name = NULL;

	ba_node = list_peek(nodes);
	if(!ba_node)
		return SLURM_ERROR;

	remove_block(nodes, new_count);
	list_delete_all(nodes, &empty_null_destroy_list, "");
		
	name = set_bg_block(nodes, ba_node->coord, geo, conn_type);
	if(!name)
		return SLURM_ERROR;
	else {
		xfree(name);
		return SLURM_SUCCESS;
	}
}

extern int copy_node_path(List nodes, List dest_nodes)
{
	int rc = SLURM_ERROR;
	
#ifdef HAVE_BG
	ListIterator itr = NULL;
	ListIterator itr2 = NULL;
	ba_node_t *ba_node = NULL, *new_ba_node = NULL;
	int dim;
	ba_switch_t *curr_switch = NULL, *new_switch = NULL; 
	
	if(!nodes)
		return SLURM_ERROR;
	if(!dest_nodes)
		dest_nodes = list_create(destroy_ba_node);

	itr = list_iterator_create(nodes);
	while((ba_node = list_next(itr))) {
		itr2 = list_iterator_create(dest_nodes);
		while((new_ba_node = list_next(itr2))) {
			if (ba_node->coord[X] == new_ba_node->coord[X] &&
			    ba_node->coord[Y] == new_ba_node->coord[Y] &&
			    ba_node->coord[Z] == new_ba_node->coord[Z]) 
				break;	/* we found it */
		}
		list_iterator_destroy(itr2);
	
		if(!new_ba_node) {
			debug2("adding %d%d%d as a new node",
			       ba_node->coord[X], 
			       ba_node->coord[Y],
			       ba_node->coord[Z]);
			new_ba_node = ba_copy_node(ba_node);
			_new_ba_node(new_ba_node, ba_node->coord, false);
			list_push(dest_nodes, new_ba_node);
			
		}
		new_ba_node->used = true;
		for(dim=0;dim<BA_SYSTEM_DIMENSIONS;dim++) {		
			curr_switch = &ba_node->axis_switch[dim];
			new_switch = &new_ba_node->axis_switch[dim];
			if(curr_switch->int_wire[0].used) {
				_copy_the_path(dest_nodes, 
					       curr_switch, new_switch,
					       0, dim);
			}
		}
		
	}
	list_iterator_destroy(itr);
	rc = SLURM_SUCCESS;
#endif	
	return rc;
}
extern int check_and_set_node_list(List nodes)
{
	int rc = SLURM_ERROR;

#ifdef HAVE_BG
	int i, j;
	ba_switch_t *ba_switch = NULL, *curr_ba_switch = NULL; 
	ba_node_t *ba_node = NULL, *curr_ba_node = NULL;
	ListIterator itr = NULL;

	if(!nodes)
		return rc;

	itr = list_iterator_create(nodes);
	while((ba_node = list_next(itr))) { 
		/* info("checking %d%d%d", */
/* 		     ba_node->coord[X],  */
/* 		     ba_node->coord[Y], */
/* 		     ba_node->coord[Z]); */
					      
		curr_ba_node = &ba_system_ptr->
			grid[ba_node->coord[X]]
			[ba_node->coord[Y]]
			[ba_node->coord[Z]];
		if(ba_node->used && curr_ba_node->used) {
			debug3("I have already been to "
			       "this node %d%d%d",
			       ba_node->coord[X], 
			       ba_node->coord[Y],
			       ba_node->coord[Z]);
			rc = SLURM_ERROR;
			goto end_it;
		}
		
		if(ba_node->used) 
			curr_ba_node->used = true;		
		for(i=0; i<BA_SYSTEM_DIMENSIONS; i++) {
			ba_switch = &ba_node->axis_switch[i];
			curr_ba_switch = &curr_ba_node->axis_switch[i];
			//info("checking dim %d", i);
		
			for(j=0; j<BA_SYSTEM_DIMENSIONS; j++) {
				//info("checking port %d", j);
		
				if(ba_switch->int_wire[j].used 
				   && curr_ba_switch->int_wire[j].used) {
					debug3("%d%d%d dim %d port %d "
					       "is already in use",
					       ba_node->coord[X], 
					       ba_node->coord[Y],
					       ba_node->coord[Z], 
					       i,
					       j);
					rc = SLURM_ERROR;
					goto end_it;
				}
				if(!ba_switch->int_wire[j].used)
					continue;
				curr_ba_switch->int_wire[j].used = 1;
				curr_ba_switch->int_wire[j].port_tar 
					= ba_switch->int_wire[j].port_tar;
			}
		}
	}
	rc = SLURM_SUCCESS;
end_it:
	list_iterator_destroy(itr);
#endif
	return rc;
}

extern char *set_bg_block(List results, int *start, 
			  int *geometry, int conn_type)
{
	char *name = NULL;
	ba_node_t* ba_node = NULL;
	int size = 0;
	int send_results = 0;
	int found = 0;


#ifdef HAVE_BG
	if(start[X]>=DIM_SIZE[X] 
	   || start[Y]>=DIM_SIZE[Y]
	   || start[Z]>=DIM_SIZE[Z])
		return NULL;
	if(geometry[X]<=0 
	   || geometry[Y]<=0
	   || geometry[Z]<=0) {
		error("problem with geometry %d%d%d, needs to be at least 111",
		      geometry[X],
		      geometry[Y],
		      geometry[Z]);		      
		return NULL;
	}
	size = geometry[X] * geometry[Y] * geometry[Z];
	ba_node = &ba_system_ptr->
		grid[start[X]][start[Y]][start[Z]];
#else
	if(start[X]>=DIM_SIZE[X])
		return NULL;
	size = geometry[X];
	ba_node = &ba_system_ptr->
			grid[start[X]];	
#endif
	

	if(!ba_node)
		return NULL;

	if(!results)
		results = list_create(NULL);
	else
		send_results = 1;
		
	list_append(results, ba_node);
	if(conn_type == SELECT_SMALL) {
		/* adding the ba_node and ending */
		ba_node->used = true;
		name = xmalloc(4);
		snprintf(name, 4, "%d%d%d",
			 ba_node->coord[X],
			 ba_node->coord[Y],
			 ba_node->coord[Z]);
		goto end_it; 
	}
	found = _find_x_path(results, ba_node,
			     ba_node->coord, 
			     ba_node->coord, 
			     geometry, 
			     1,
			     conn_type);

	if(!found) {
		debug2("trying less efficient code");
		remove_block(results, color_count);
		list_delete_all(results, &empty_null_destroy_list, "");
		list_append(results, ba_node);
		found = _find_x_path2(results, ba_node,
				      ba_node->coord,
				      ba_node->coord,
				      geometry,
				      1,
				      conn_type);
	}
	if(found) {
#ifdef HAVE_BG
		List start_list = NULL;
		ListIterator itr;

		start_list = list_create(NULL);
		itr = list_iterator_create(results);
		while((ba_node = (ba_node_t*) list_next(itr))) {
			list_append(start_list, ba_node);
		}
		list_iterator_destroy(itr);
		
		if(!_fill_in_coords(results, 
				    start_list, 
				    geometry, 
				    conn_type)) {
			list_destroy(start_list);
			goto end_it;
		}
		list_destroy(start_list);			
#endif		
	} else {
		goto end_it;
	}

	name = _set_internal_wires(results,
				   size,
				   conn_type);
end_it:
	if(!send_results && results) {
		list_destroy(results);
	}
	if(name!=NULL) {
		debug2("name = %s", name);
	} else {
		debug2("can't allocate");
		xfree(name);
	}

	return name;	
}

extern int reset_ba_system(bool track_down_nodes)
{
	int x;
#ifdef HAVE_BG
	int y, z;
#endif
	int coord[BA_SYSTEM_DIMENSIONS];

	for (x = 0; x < DIM_SIZE[X]; x++) {
#ifdef HAVE_BG
		for (y = 0; y < DIM_SIZE[Y]; y++)
			for (z = 0; z < DIM_SIZE[Z]; z++) {
				coord[X] = x;
				coord[Y] = y;
				coord[Z] = z;
				_new_ba_node(&ba_system_ptr->grid[x][y][z], 
					     coord, track_down_nodes);
			}
#else
		coord[X] = x;
		_new_ba_node(&ba_system_ptr->grid[x], coord, track_down_nodes);

#endif
	}
				
	return 1;
}
/* init_grid - set values of every grid point */
extern void init_grid(node_info_msg_t * node_info_ptr)
{
	node_info_t *node_ptr = NULL;
	int x, i = 0;
	uint16_t node_base_state;
	/* For systems with more than 62 active jobs or BG blocks, 
	 * we just repeat letters */

#ifdef HAVE_BG
	int y,z;
	for (x = 0; x < DIM_SIZE[X]; x++)
		for (y = 0; y < DIM_SIZE[Y]; y++)
			for (z = 0; z < DIM_SIZE[Z]; z++) {
				if(node_info_ptr!=NULL) {
					node_ptr = 
						&node_info_ptr->node_array[i];
					node_base_state = node_ptr->node_state 
						& NODE_STATE_BASE;
					ba_system_ptr->grid[x][y][z].color = 7;
					if ((node_base_state 
					     == NODE_STATE_DOWN) || 
					    (node_ptr->node_state &
					     NODE_STATE_DRAIN)) {
						ba_system_ptr->
							grid[x][y][z].color 
							= 0;
						ba_system_ptr->
							grid[x][y][z].letter 
							= '#';
						if(_initialized) {
							ba_update_node_state(
							&ba_system_ptr->
							grid[x][y][z],
							node_ptr->node_state);
						}
					} else {
						ba_system_ptr->grid[x][y][z].
							color = 7;
						ba_system_ptr->grid[x][y][z].
							letter = '.';
					}
					ba_system_ptr->grid[x][y][z].state 
						= node_ptr->node_state;
				} else {
					ba_system_ptr->grid[x][y][z].color = 7;
					ba_system_ptr->grid[x][y][z].letter 
						= '.';
					ba_system_ptr->grid[x][y][z].state = 
						NODE_STATE_IDLE;
				}
				ba_system_ptr->grid[x][y][z].index = i++;
			}
#else
	for (x = 0; x < DIM_SIZE[X]; x++) {
		if(node_info_ptr!=NULL) {
			node_ptr = &node_info_ptr->node_array[i];
			node_base_state = node_ptr->node_state 
				& NODE_STATE_BASE;
			ba_system_ptr->grid[x].color = 7;
			if ((node_base_state == NODE_STATE_DOWN) || 
			    (node_ptr->node_state & NODE_STATE_DRAIN)) {
				ba_system_ptr->grid[x].color = 0;
				ba_system_ptr->grid[x].letter = '#';
				if(_initialized) {
					ba_update_node_state(
						&ba_system_ptr->grid[x],
						node_ptr->node_state);
				}
			} else {
				ba_system_ptr->grid[x].color = 7;
				ba_system_ptr->grid[x].letter = '.';
			}
			ba_system_ptr->grid[x].state = node_ptr->node_state;
		} else {
			ba_system_ptr->grid[x].color = 7;
			ba_system_ptr->grid[x].letter = '.';
			ba_system_ptr->grid[x].state = NODE_STATE_IDLE;
		}
		ba_system_ptr->grid[x].index = i++;
	}
#endif
	return;
}

extern int *find_bp_loc(char* bp_id)
{
#ifdef HAVE_BG_FILES
	ba_bp_map_t *bp_map = NULL;
	ListIterator itr;
	
	if(!bp_map_list) {
		if(set_bp_map() == -1)
			return NULL;
	}
	itr = list_iterator_create(bp_map_list);
	while ((bp_map = list_next(itr)) != NULL)
		if (!strcasecmp(bp_map->bp_id, bp_id)) 
			break;	/* we found it */
	
	list_iterator_destroy(itr);
	if(bp_map != NULL)
		return bp_map->coord;
	else
		return NULL;

#else
	return NULL;
#endif
}

extern char *find_bp_rack_mid(char* xyz)
{
#ifdef HAVE_BG_FILES
	ba_bp_map_t *bp_map = NULL;
	ListIterator itr;
	int number;
	int coord[BA_SYSTEM_DIMENSIONS];
	int len = strlen(xyz);
	len -= 3;
	if(len<0)
		return NULL;
	number = atoi(&xyz[X]+len);
	coord[X] = number / 100;
	coord[Y] = (number % 100) / 10;
	coord[Z] = (number % 10);
	if(!bp_map_list) {
		if(set_bp_map() == -1)
			return NULL;
	}
	
	itr = list_iterator_create(bp_map_list);
	while ((bp_map = list_next(itr)) != NULL)
		if (bp_map->coord[X] == coord[X] &&
		    bp_map->coord[Y] == coord[Y] &&
		    bp_map->coord[Z] == coord[Z]) 
			break;	/* we found it */
	
	list_iterator_destroy(itr);
	if(bp_map != NULL)
		return bp_map->bp_id;
	else
		return NULL;

#else
	return NULL;
#endif
}

extern int load_block_wiring(char *bg_block_id)
{
#ifdef HAVE_BG_FILES
	int rc, i, j;
	rm_partition_t *block_ptr = NULL;
	int cnt = 0;
	int switch_cnt = 0;
	rm_switch_t *curr_switch = NULL;
	rm_BP_t *curr_bp = NULL;
	char *switchid = NULL;
	rm_connection_t curr_conn;
	int dim;
	ba_switch_t *ba_switch = NULL; 
	int *geo = NULL;
	
	debug2("getting info for block %s\n", bg_block_id);
	
	if ((rc = bridge_get_block(bg_block_id,  &block_ptr)) != STATUS_OK) {
		error("bridge_get_block(%s): %s", 
		      bg_block_id, 
		      bg_err_str(rc));
		return SLURM_ERROR;
	}	
	
	if ((rc = bridge_get_data(block_ptr, RM_PartitionSwitchNum,
				  &switch_cnt)) != STATUS_OK) {
		error("bridge_get_data(RM_PartitionSwitchNum): %s",
		      bg_err_str(rc));
		return SLURM_ERROR;
	} 
	if(!switch_cnt) {
		debug3("no switch_cnt");
		if ((rc = bridge_get_data(block_ptr, 
					  RM_PartitionFirstBP, 
					  &curr_bp)) 
		    != STATUS_OK) {
			error("bridge_get_data: "
			      "RM_PartitionFirstBP: %s",
			      bg_err_str(rc));
			return SLURM_ERROR;
		}
		if ((rc = bridge_get_data(curr_bp, RM_BPID, &switchid))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchBPID: %s",
			      bg_err_str(rc));
			return SLURM_ERROR;
		} 

		geo = find_bp_loc(switchid);	
		if(!geo) {
			error("find_bp_loc: bpid %s not known", switchid);
			return SLURM_ERROR;
		}
		ba_system_ptr->grid[geo[X]][geo[Y]][geo[Z]].used = true;
		return SLURM_SUCCESS;
	}
	for (i=0; i<switch_cnt; i++) {
		if(i) {
			if ((rc = bridge_get_data(block_ptr, 
						  RM_PartitionNextSwitch, 
						  &curr_switch)) 
			    != STATUS_OK) {
				error("bridge_get_data: "
				      "RM_PartitionNextSwitch: %s",
				      bg_err_str(rc));
				return SLURM_ERROR;
			}
		} else {
			if ((rc = bridge_get_data(block_ptr, 
						  RM_PartitionFirstSwitch, 
						  &curr_switch)) 
			    != STATUS_OK) {
				error("bridge_get_data: "
				      "RM_PartitionFirstSwitch: %s",
				      bg_err_str(rc));
				return SLURM_ERROR;
			}
		}
		if ((rc = bridge_get_data(curr_switch, RM_SwitchDim, &dim))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchDim: %s",
			      bg_err_str(rc));
			return SLURM_ERROR;
		} 
		if ((rc = bridge_get_data(curr_switch, RM_SwitchBPID, 
					  &switchid))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchBPID: %s",
			      bg_err_str(rc));
			return SLURM_ERROR;
		} 

		geo = find_bp_loc(switchid);
		if(!geo) {
			error("find_bp_loc: bpid %s not known", switchid);
			return SLURM_ERROR;
		}
		
		if ((rc = bridge_get_data(curr_switch, RM_SwitchConnNum, &cnt))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchBPID: %s",
			      bg_err_str(rc));
			return SLURM_ERROR;
		}
		debug2("switch id = %s dim %d conns = %d", 
		       switchid, dim, cnt);
		ba_switch = &ba_system_ptr->
			grid[geo[X]][geo[Y]][geo[Z]].axis_switch[dim];
		for (j=0; j<cnt; j++) {
			if(j) {
				if ((rc = bridge_get_data(
					     curr_switch, 
					     RM_SwitchNextConnection, 
					     &curr_conn)) 
				    != STATUS_OK) {
					error("bridge_get_data: "
					      "RM_SwitchNextConnection: %s",
					       bg_err_str(rc));
					return SLURM_ERROR;
				}
			} else {
				if ((rc = bridge_get_data(
					     curr_switch, 
					     RM_SwitchFirstConnection,
					     &curr_conn)) 
				    != STATUS_OK) {
					error("bridge_get_data: "
					      "RM_SwitchFirstConnection: %s",
					      bg_err_str(rc));
					return SLURM_ERROR;
				}
			}
			switch(curr_conn.p1) {
			case RM_PORT_S1:
				curr_conn.p1 = 1;
				break;
			case RM_PORT_S2:
				curr_conn.p1 = 2;
				break;
			case RM_PORT_S4:
				curr_conn.p1 = 4;
				break;
			default:
				error("1 unknown port %d", 
				      _port_enum(curr_conn.p1));
				return SLURM_ERROR;
			}
			
			switch(curr_conn.p2) {
			case RM_PORT_S0:
				curr_conn.p2 = 0;
				break;
			case RM_PORT_S3:
				curr_conn.p2 = 3;
				break;
			case RM_PORT_S5:
				curr_conn.p2 = 5;
				break;
			default:
				error("2 unknown port %d", 
				      _port_enum(curr_conn.p2));
				return SLURM_ERROR;
			}

			if(curr_conn.p1 == 1 && dim == X) {
				if(ba_system_ptr->
				   grid[geo[X]][geo[Y]][geo[Z]].used) {
					debug("I have already been to "
					      "this node %d%d%d",
					      geo[X], geo[Y], geo[Z]);
					return SLURM_ERROR;
				}
				ba_system_ptr->grid[geo[X]][geo[Y]][geo[Z]].
					used = true;		
			}
			debug3("connection going from %d -> %d",
			      curr_conn.p1, curr_conn.p2);
			
			if(ba_switch->int_wire[curr_conn.p1].used) {
				debug("%d%d%d dim %d port %d "
				      "is already in use",
				      geo[X],
				      geo[Y],
				      geo[Z],
				      dim,
				      curr_conn.p1);
				return SLURM_ERROR;
			}
			ba_switch->int_wire[curr_conn.p1].used = 1;
			ba_switch->int_wire[curr_conn.p1].port_tar 
				= curr_conn.p2;
		
			if(ba_switch->int_wire[curr_conn.p2].used) {
				debug("%d%d%d dim %d port %d "
				      "is already in use",
				      geo[X],
				      geo[Y],
				      geo[Z],
				      dim,
				      curr_conn.p2);
				return SLURM_ERROR;
			}
			ba_switch->int_wire[curr_conn.p2].used = 1;
			ba_switch->int_wire[curr_conn.p2].port_tar 
				= curr_conn.p1;
		}
	}
	return SLURM_SUCCESS;

#else
	return SLURM_ERROR;
#endif
	
}

extern List get_and_set_block_wiring(char *bg_block_id)
{
#ifdef HAVE_BG_FILES
	int rc, i, j;
	rm_partition_t *block_ptr = NULL;
	int cnt = 0;
	int switch_cnt = 0;
	rm_switch_t *curr_switch = NULL;
	rm_BP_t *curr_bp = NULL;
	char *switchid = NULL;
	rm_connection_t curr_conn;
	int dim;
	ba_node_t *ba_node = NULL; 
	ba_switch_t *ba_switch = NULL; 
	int *geo = NULL;
	List results = list_create(destroy_ba_node);
	ListIterator itr = NULL;
	
	debug2("getting info for block %s\n", bg_block_id);
	
	if ((rc = bridge_get_block(bg_block_id,  &block_ptr)) != STATUS_OK) {
		error("bridge_get_block(%s): %s", 
		      bg_block_id, 
		      bg_err_str(rc));
		goto end_it;
	}	
	
	if ((rc = bridge_get_data(block_ptr, RM_PartitionSwitchNum,
				  &switch_cnt)) != STATUS_OK) {
		error("bridge_get_data(RM_PartitionSwitchNum): %s",
		      bg_err_str(rc));
		goto end_it;
	} 
	if(!switch_cnt) {
		debug3("no switch_cnt");
		if ((rc = bridge_get_data(block_ptr, 
					  RM_PartitionFirstBP, 
					  &curr_bp)) 
		    != STATUS_OK) {
			error("bridge_get_data: "
			      "RM_PartitionFirstBP: %s",
			      bg_err_str(rc));
			goto end_it;
		}
		if ((rc = bridge_get_data(curr_bp, RM_BPID, &switchid))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchBPID: %s",
			      bg_err_str(rc));
			goto end_it;
		} 

		geo = find_bp_loc(switchid);	
		if(!geo) {
			error("find_bp_loc: bpid %s not known", switchid);
			goto end_it;
		}
		ba_node = xmalloc(sizeof(ba_node_t));
		list_push(results, ba_node);
		ba_node->coord[X] = geo[X];
		ba_node->coord[Y] = geo[Y];
		ba_node->coord[Z] = geo[Z];
		
		ba_node->used = TRUE;
		return results;
	}
	for (i=0; i<switch_cnt; i++) {
		if(i) {
			if ((rc = bridge_get_data(block_ptr, 
						  RM_PartitionNextSwitch, 
						  &curr_switch)) 
			    != STATUS_OK) {
				error("bridge_get_data: "
				      "RM_PartitionNextSwitch: %s",
				      bg_err_str(rc));
				goto end_it;
			}
		} else {
			if ((rc = bridge_get_data(block_ptr, 
						  RM_PartitionFirstSwitch, 
						  &curr_switch)) 
			    != STATUS_OK) {
				error("bridge_get_data: "
				      "RM_PartitionFirstSwitch: %s",
				      bg_err_str(rc));
				goto end_it;
			}
		}
		if ((rc = bridge_get_data(curr_switch, RM_SwitchDim, &dim))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchDim: %s",
			      bg_err_str(rc));
			goto end_it;
		} 
		if ((rc = bridge_get_data(curr_switch, RM_SwitchBPID, 
					  &switchid))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchBPID: %s",
			      bg_err_str(rc));
			goto end_it;
		} 

		geo = find_bp_loc(switchid);
		if(!geo) {
			error("find_bp_loc: bpid %s not known", switchid);
			goto end_it;
		}
		
		if ((rc = bridge_get_data(curr_switch, RM_SwitchConnNum, &cnt))
		    != STATUS_OK) { 
			error("bridge_get_data: RM_SwitchBPID: %s",
			      bg_err_str(rc));
			goto end_it;
		}
		debug2("switch id = %s dim %d conns = %d", 
		       switchid, dim, cnt);
		
		itr = list_iterator_create(results);
		while((ba_node = list_next(itr))) {
			if (ba_node->coord[X] == geo[X] &&
			    ba_node->coord[Y] == geo[Y] &&
			    ba_node->coord[Z] == geo[Z]) 
				break;	/* we found it */
		}
		list_iterator_destroy(itr);
		if(!ba_node) {
			ba_node = xmalloc(sizeof(ba_node_t));
			
			list_push(results, ba_node);
			ba_node->coord[X] = geo[X];
			ba_node->coord[Y] = geo[Y];
			ba_node->coord[Z] = geo[Z];
		}
		ba_switch = &ba_node->axis_switch[dim];
		for (j=0; j<cnt; j++) {
			if(j) {
				if ((rc = bridge_get_data(
					     curr_switch, 
					     RM_SwitchNextConnection, 
					     &curr_conn)) 
				    != STATUS_OK) {
					error("bridge_get_data: "
					      "RM_SwitchNextConnection: %s",
					       bg_err_str(rc));
					goto end_it;
				}
			} else {
				if ((rc = bridge_get_data(
					     curr_switch, 
					     RM_SwitchFirstConnection,
					     &curr_conn)) 
				    != STATUS_OK) {
					error("bridge_get_data: "
					      "RM_SwitchFirstConnection: %s",
					      bg_err_str(rc));
					goto end_it;
				}
			}
			switch(curr_conn.p1) {
			case RM_PORT_S1:
				curr_conn.p1 = 1;
				break;
			case RM_PORT_S2:
				curr_conn.p1 = 2;
				break;
			case RM_PORT_S4:
				curr_conn.p1 = 4;
				break;
			default:
				error("1 unknown port %d", 
				      _port_enum(curr_conn.p1));
				goto end_it;
			}
			
			switch(curr_conn.p2) {
			case RM_PORT_S0:
				curr_conn.p2 = 0;
				break;
			case RM_PORT_S3:
				curr_conn.p2 = 3;
				break;
			case RM_PORT_S5:
				curr_conn.p2 = 5;
				break;
			default:
				error("2 unknown port %d", 
				      _port_enum(curr_conn.p2));
				goto end_it;
			}

			if(curr_conn.p1 == 1 && dim == X) {
				if(ba_node->used) {
					debug("I have already been to "
					      "this node %d%d%d",
					      geo[X], geo[Y], geo[Z]);
					goto end_it;
				}
				ba_node->used = true;		
			}
			debug3("connection going from %d -> %d",
			      curr_conn.p1, curr_conn.p2);
			
			if(ba_switch->int_wire[curr_conn.p1].used) {
				debug("%d%d%d dim %d port %d "
				      "is already in use",
				      geo[X],
				      geo[Y],
				      geo[Z],
				      dim,
				      curr_conn.p1);
				goto end_it;
			}
			ba_switch->int_wire[curr_conn.p1].used = 1;
			ba_switch->int_wire[curr_conn.p1].port_tar 
				= curr_conn.p2;
		
			if(ba_switch->int_wire[curr_conn.p2].used) {
				debug("%d%d%d dim %d port %d "
				      "is already in use",
				      geo[X],
				      geo[Y],
				      geo[Z],
				      dim,
				      curr_conn.p2);
				goto end_it;
			}
			ba_switch->int_wire[curr_conn.p2].used = 1;
			ba_switch->int_wire[curr_conn.p2].port_tar 
				= curr_conn.p1;
		}
	}
	return results;
end_it:
	list_destroy(results);
	return NULL;
#else
	return NULL;
#endif
	
}

/********************* Local Functions *********************/

#ifdef HAVE_BG

#ifdef HAVE_BG_FILES
static void _bp_map_list_del(void *object)
{
	ba_bp_map_t *bp_map = (ba_bp_map_t *)object;
	
	if (bp_map) {
		xfree(bp_map->bp_id);
		xfree(bp_map);		
	}
}

static int _port_enum(int port)
{
	switch(port) {
	case RM_PORT_S0:
		return 0;
		break;
	case RM_PORT_S1:
		return 1;
		break;
	case RM_PORT_S2:
		return 2;
		break;
	case RM_PORT_S3:
		return 3;
		break;
	case RM_PORT_S4:
		return 4;
		break;
	case RM_PORT_S5:
		return 5;
		break;
	default:
		return -1;
	}
}

#endif

static int _check_for_options(ba_request_t* ba_request) 
{
	int temp;
	int set=0;
	int *geo = NULL;
	ListIterator itr;

	if(ba_request->rotate) {
	rotate_again:
		debug2("Rotating! %d",ba_request->rotate_count);
		
		if (ba_request->rotate_count==(BA_SYSTEM_DIMENSIONS-1)) {
			temp=ba_request->geometry[X];
			ba_request->geometry[X]=ba_request->geometry[Z];
			ba_request->geometry[Z]=temp;
			ba_request->rotate_count++;
			set=1;
		
		} else if(ba_request->rotate_count<(BA_SYSTEM_DIMENSIONS*2)) {
			temp=ba_request->geometry[X];
			ba_request->geometry[X]=ba_request->geometry[Y];
			ba_request->geometry[Y]=ba_request->geometry[Z];
			ba_request->geometry[Z]=temp;
			ba_request->rotate_count++;
			set=1;
		} else 
			ba_request->rotate = false;
		if(set) {
			if(ba_request->geometry[X]<=DIM_SIZE[X] 
			   && ba_request->geometry[Y]<=DIM_SIZE[Y]
			   && ba_request->geometry[Z]<=DIM_SIZE[Z])
				return 1;
			else {
				set = 0;
				goto rotate_again;
			}
		}
	}
	if(ba_request->elongate) {
	elongate_again:
		debug2("Elongating! %d",ba_request->elongate_count);
		ba_request->rotate_count=0;
		ba_request->rotate = true;
		
		set = 0;
		itr = list_iterator_create(ba_request->elongate_geos);
		for(set=0; set<=ba_request->elongate_count; set++)
			geo = list_next(itr);
		list_iterator_destroy(itr);
		if(geo == NULL)
			return 0;
		ba_request->elongate_count++;
		ba_request->geometry[X] = geo[X];
		ba_request->geometry[Y] = geo[Y];
		ba_request->geometry[Z] = geo[Z];
		if(ba_request->geometry[X]<=DIM_SIZE[X] 
		   && ba_request->geometry[Y]<=DIM_SIZE[Y]
		   && ba_request->geometry[Z]<=DIM_SIZE[Z]) {
			return 1;
		} else 				
			goto elongate_again;
		
	}
	return 0;
}

static int _append_geo(int *geometry, List geos, int rotate) 
{
	ListIterator itr;
	int *geo_ptr = NULL;
	int *geo = NULL;
	int temp_geo;
	int i, j;
	
	if(rotate) {
		for (i = (BA_SYSTEM_DIMENSIONS - 1); i >= 0; i--) {
			for (j = 1; j <= i; j++) {
				if (geometry[j-1] > geometry[j]) {
					temp_geo = geometry[j-1];
					geometry[j-1] = geometry[j];
					geometry[j] = temp_geo;
				}
			}
		}
	}
	itr = list_iterator_create(geos);
	while ((geo_ptr = list_next(itr)) != NULL) {
		if(geometry[X] == geo_ptr[X]
		   && geometry[Y] == geo_ptr[Y]
		   && geometry[Z] == geo_ptr[Z])
			break;
		
	}
	list_iterator_destroy(itr);
	
	if(geo_ptr == NULL) { 
		geo = xmalloc(sizeof(int)*BA_SYSTEM_DIMENSIONS);
		geo[X] = geometry[X];
		geo[Y] = geometry[Y];
		geo[Z] = geometry[Z];
		debug3("adding geo %d%d%d",geo[X],geo[Y],geo[Z]);
		list_append(geos, geo);
	}
	return 1;
}

static int _fill_in_coords(List results, List start_list,
			    int *geometry, int conn_type)
{
	ba_node_t *ba_node = NULL;
	ba_node_t *check_node = NULL;
	int rc = 1;
	ListIterator itr = NULL;
	int y=0, z=0;
	ba_switch_t *curr_switch = NULL; 
	ba_switch_t *next_switch = NULL; 
	
	if(!start_list)
		return 0;
	itr = list_iterator_create(start_list);
	while((check_node = (ba_node_t*) list_next(itr))) {		
		curr_switch = &check_node->axis_switch[X];
	
		for(y=0; y<geometry[Y]; y++) {
			if((check_node->coord[Y]+y) 
			   >= DIM_SIZE[Y]) {
				rc = 0;
				goto failed;
			}
			for(z=0; z<geometry[Z]; z++) {
				if((check_node->coord[Z]+z) 
				   >= DIM_SIZE[Z]) {
					rc = 0;
					goto failed;
				}
				ba_node = &ba_system_ptr->grid
					[check_node->coord[X]]
					[check_node->coord[Y]+y]
					[check_node->coord[Z]+z];
				if(ba_node->coord[Y] 
				   == check_node->coord[Y]
				   && ba_node->coord[Z] 
				   == check_node->coord[Z])
					continue;
				if (!_node_used(ba_node,geometry)) {
					debug3("here Adding %d%d%d",
					       ba_node->coord[X],
					       ba_node->coord[Y],
					       ba_node->coord[Z]);
					list_append(results, ba_node);
					next_switch = &ba_node->axis_switch[X];
					_copy_the_path(NULL, curr_switch, 
						       next_switch, 
						       0, X);
				} else {
					rc = 0;
					goto failed;
				}
			}
		}
		
	}
	list_iterator_destroy(itr);
	itr = list_iterator_create(start_list);
	check_node = (ba_node_t*) list_next(itr);
	list_iterator_destroy(itr);
	
	itr = list_iterator_create(results);
	while((ba_node = (ba_node_t*) list_next(itr))) {	
		if(!_find_yz_path(ba_node, 
				  check_node->coord, 
				  geometry, 
				  conn_type)){
			rc = 0;
			goto failed;
		}
	}
failed:
	list_iterator_destroy(itr);				
				
	return rc;
}

static int _copy_the_path(List nodes, ba_switch_t *curr_switch, 
			  ba_switch_t *mark_switch, 
			  int source, int dim)
{
	int *node_tar;
	int *mark_node_tar;
	int *node_curr;
	int port_tar, port_tar1;
	ba_switch_t *next_switch = NULL; 
	ba_switch_t *next_mark_switch = NULL; 
	/*set the switch to not be used */
	mark_switch->int_wire[source].used = 
		curr_switch->int_wire[source].used;
	mark_switch->int_wire[source].port_tar = 
		curr_switch->int_wire[source].port_tar;

	port_tar = curr_switch->int_wire[source].port_tar;
	
	mark_switch->int_wire[port_tar].used = 
		curr_switch->int_wire[port_tar].used;
	mark_switch->int_wire[port_tar].port_tar = 
		curr_switch->int_wire[port_tar].port_tar;
	port_tar1 = port_tar;
	
	/* follow the path */
	node_curr = curr_switch->ext_wire[0].node_tar;
	node_tar = curr_switch->ext_wire[port_tar].node_tar;
	if(mark_switch->int_wire[source].used)
		debug2("setting dim %d %d%d%d %d-> %d%d%d %d",
		       dim,
		       node_curr[X],
		       node_curr[Y],
		       node_curr[Z],
		       source, 
		       node_tar[X],
		       node_tar[Y],
		       node_tar[Z],
		       port_tar);	
	
	if(port_tar == 1) {
		mark_switch->int_wire[1].used = 
			curr_switch->int_wire[1].used;
		mark_switch->int_wire[1].port_tar = 
			curr_switch->int_wire[1].port_tar;
		return 1;
	}

	mark_node_tar = mark_switch->ext_wire[port_tar].node_tar;
	port_tar = curr_switch->ext_wire[port_tar].port_tar;
	
	if(node_curr[X] == node_tar[X]
	   && node_curr[Y] == node_tar[Y]
	   && node_curr[Z] == node_tar[Z]) {
		debug4("something bad happened!!");
		return 0;
	}
	next_switch = &ba_system_ptr->
		grid[node_tar[X]][node_tar[Y]][node_tar[Z]].axis_switch[dim];
	if(!nodes) {
		next_mark_switch = &ba_system_ptr->
			grid[mark_node_tar[X]]
			[mark_node_tar[Y]]
			[mark_node_tar[Z]]
			.axis_switch[dim];
	} else {
		ba_node_t *ba_node = NULL;
		ListIterator itr = list_iterator_create(nodes);
		while((ba_node = list_next(itr))) {
			if (ba_node->coord[X] == mark_node_tar[X] &&
			    ba_node->coord[Y] == mark_node_tar[Y] &&
			    ba_node->coord[Z] == mark_node_tar[Z]) 
				break;	/* we found it */
		}
		list_iterator_destroy(itr);
		if(!ba_node) {
			ba_node = ba_copy_node(&ba_system_ptr->
					       grid[mark_node_tar[X]]
					       [mark_node_tar[Y]]
					       [mark_node_tar[Z]]);
			_new_ba_node(ba_node, mark_node_tar, false);
			list_push(nodes, ba_node);
			debug3("adding %d%d%d as a pass through",
			       ba_node->coord[X], 
			       ba_node->coord[Y],
			       ba_node->coord[Z]);
		}
		next_mark_switch = &ba_node->axis_switch[dim];
			
	}
	_copy_the_path(nodes, next_switch, next_mark_switch,
		       port_tar, dim);
	return 1;
}

static int _find_yz_path(ba_node_t *ba_node, int *first, 
			 int *geometry, int conn_type)
{
	ba_node_t *next_node = NULL;
	int *node_tar = NULL;
	ba_switch_t *dim_curr_switch = NULL; 
	ba_switch_t *dim_next_switch = NULL; 
	int i2;
	int count = 0;

	for(i2=1;i2<=2;i2++) {
		if(geometry[i2] > 1) {
			debug3("%d node %d%d%d"
			       " port 2 -> ",
			       i2,
			       ba_node->coord[X],
			       ba_node->coord[Y],
			       ba_node->coord[Z]);
							       
			dim_curr_switch = 
				&ba_node->
				axis_switch[i2];
			if(dim_curr_switch->int_wire[2].used) {
				debug4("returning here");
				return 0;
			}
							
			node_tar = dim_curr_switch->
				ext_wire[2].node_tar;
							
			next_node = &ba_system_ptr->
				grid[node_tar[X]][node_tar[Y]][node_tar[Z]];
			dim_next_switch = &next_node->axis_switch[i2];
			debug3("%d%d%d port 5",
			       next_node->coord[X],
			       next_node->coord[Y],
			       next_node->coord[Z]);
							  
			if(dim_next_switch->int_wire[5].used) {
				debug2("returning here 2");
				return 0;
			}
			debug4("%d %d %d %d",i2, node_tar[i2],
			       first[i2], geometry[i2]);
			if(node_tar[i2] < first[i2])
				count = DIM_SIZE[i2]-first[i2]+node_tar[i2];
			else
				count = node_tar[i2]+first[i2];
			if((count) == (geometry[i2])) {
				debug4("found end of me %d%d%d",
				       node_tar[X],
				       node_tar[Y],
				       node_tar[Z]);
				if(conn_type == SELECT_TORUS) {
					dim_curr_switch->
						int_wire[0].used = 1;
					dim_curr_switch->
						int_wire[0].port_tar
						= 2;
					dim_curr_switch->
						int_wire[2].used
						= 1;
					dim_curr_switch->
						int_wire[2].
						port_tar = 0;
					dim_curr_switch = dim_next_switch;
									
					while(node_tar[i2] != first[i2]) {
						debug3("on dim %d at %d "
						       "looking for %d",
						       i2,
						       node_tar[i2],
						       first[i2]);
						
						if(dim_curr_switch->
						   int_wire[2].used) {
							debug3("returning "
							       "here 3");
							return 0;
						} 
						dim_curr_switch->
							int_wire[2].used = 1;
						dim_curr_switch->
							int_wire[2].port_tar
							= 5;
						dim_curr_switch->
							int_wire[5].used
							= 1;
						dim_curr_switch->
							int_wire[5].
							port_tar = 2;
						
						
						node_tar = dim_curr_switch->
							ext_wire[2].node_tar;
						next_node = &ba_system_ptr->
							grid
							[node_tar[X]]
							[node_tar[Y]]
							[node_tar[Z]];
						dim_curr_switch = 
							&next_node->
							axis_switch[i2];
					}
									
					debug3("back to first on dim %d "
					       "at %d looking for %d",
					       i2,
					       node_tar[i2],
					       first[i2]);
									
					dim_curr_switch->
						int_wire[5].used = 1;
					dim_curr_switch->
						int_wire[5].port_tar
						= 1;
					dim_curr_switch->
						int_wire[1].used
						= 1;
					dim_curr_switch->
						int_wire[1].
						port_tar = 5;
				}
								
			} else {
				if(conn_type == SELECT_TORUS || 
				   (conn_type == SELECT_MESH && 
				    (node_tar[i2] != first[i2]))) {
					dim_curr_switch->
						int_wire[0].used = 1;
					dim_curr_switch->
						int_wire[0].port_tar
						= 2;
					dim_curr_switch->
						int_wire[2].used
						= 1;
					dim_curr_switch->
						int_wire[2].
						port_tar = 0;
								
					dim_next_switch->int_wire[5].used
						= 1;
					dim_next_switch->
						int_wire[5].port_tar
						= 1;
					dim_next_switch->
						int_wire[1].used = 1;
					dim_next_switch->
						int_wire[1].port_tar
						= 5;
				}
			}
		}
	}
	return 1;
}

#endif

#ifndef HAVE_BG_FILES
/** */
#ifdef HAVE_BG
static int _create_config_even(ba_node_t ***grid)
#else
static int _create_config_even(ba_node_t *grid)
#endif
{
	int x;
	ba_node_t *source = NULL, *target = NULL;

#ifdef HAVE_BG
	int y,z;
	init_wires();
	
	for(x=0;x<DIM_SIZE[X];x++) {
		for(y=0;y<DIM_SIZE[Y];y++) {
			for(z=0;z<DIM_SIZE[Z];z++) {
				source = &grid[x][y][z];
				
				if(x<(DIM_SIZE[X]-1)) {
					target = &grid[x+1][y][z];
				} else
					target = &grid[0][y][z];

				_set_external_wires(X, x, source, 
						    target);
				
				if(y<(DIM_SIZE[Y]-1)) 
					target = &grid[x][y+1][z];
				else 
					target = &grid[x][0][z];
				
				_set_external_wires(Y, y, source, 
						    target);
				if(z<(DIM_SIZE[Z]-1)) 
					target = &grid[x][y][z+1];
				else 
					target = &grid[x][y][0];
				
				_set_external_wires(Z, z, source, 
						    target);
			}
		}
	}
#else
	for(x=0;x<DIM_SIZE[X];x++) {
		source = &grid[x];
		target = &grid[x+1];
		_set_external_wires(X, x, source, 
				    target);
	}
#endif
	return 1;
}
#endif


static int _reset_the_path(ba_switch_t *curr_switch, int source, 
			   int target, int dim)
{
	int *node_tar;
	int *node_curr;
	int port_tar, port_tar1;
	ba_switch_t *next_switch = NULL; 

	if(source < 0 || source > NUM_PORTS_PER_NODE) {
		fatal("source port was %d can only be 0->%d",
		      source, NUM_PORTS_PER_NODE);
	}
	if(target < 0 || target > NUM_PORTS_PER_NODE) {
		fatal("target port was %d can only be 0->%d",
		      target, NUM_PORTS_PER_NODE);
	}
	/*set the switch to not be used */
	if(!curr_switch->int_wire[source].used) {
		debug("I reached the end, the source isn't used");
		return 1;
	}
	curr_switch->int_wire[source].used = 0;
	port_tar = curr_switch->int_wire[source].port_tar;
	if(port_tar < 0 || port_tar > NUM_PORTS_PER_NODE) {
		fatal("port_tar port was %d can only be 0->%d",
		      source, NUM_PORTS_PER_NODE);
	}
	
	port_tar1 = port_tar;
	curr_switch->int_wire[source].port_tar = source;
	curr_switch->int_wire[port_tar].used = 0;
	curr_switch->int_wire[port_tar].port_tar = port_tar;
	if(port_tar==target) {
		return 1;
	}
	/* follow the path */
	node_curr = curr_switch->ext_wire[0].node_tar;
	node_tar = curr_switch->ext_wire[port_tar].node_tar;
	port_tar = curr_switch->ext_wire[port_tar].port_tar;
	if(source == port_tar1) {
		debug("got this bad one %d%d%d %d %d -> %d%d%d %d",
		      node_curr[X],
		      node_curr[Y],
		      node_curr[Z],
		      source,
		      port_tar1,
		      node_tar[X],
		      node_tar[Y],
		      node_tar[Z],
		      port_tar);
		return 0;
	}
	debug4("from %d%d%d %d %d -> %d%d%d %d",
	       node_curr[X],
	       node_curr[Y],
	       node_curr[Z],
	       source,
	       port_tar1,
	       node_tar[X],
	       node_tar[Y],
	       node_tar[Z],
	       port_tar);
	if(node_curr[X] == node_tar[X]
	   && node_curr[Y] == node_tar[Y]
	   && node_curr[Z] == node_tar[Z]) {
		debug4("%d something bad happened!!", dim);
		return 0;
	}
	next_switch = &ba_system_ptr->
		grid[node_tar[X]]
#ifdef HAVE_BG
		[node_tar[Y]]
		[node_tar[Z]]
#endif
		.axis_switch[dim];

	_reset_the_path(next_switch, port_tar, target, dim);
	return 1;
}

/*
 * Convert a BG API error code to a string
 * IN inx - error code from any of the BG Bridge APIs
 * RET - string describing the error condition
 */
extern char *bg_err_str(status_t inx)
{
#ifdef HAVE_BG_FILES
	switch (inx) {
	case STATUS_OK:
		return "Status OK";
	case PARTITION_NOT_FOUND:
		return "Partition not found";
	case JOB_NOT_FOUND:
		return "Job not found";
	case BP_NOT_FOUND:
		return "Base partition not found";
	case SWITCH_NOT_FOUND:
		return "Switch not found";
	case JOB_ALREADY_DEFINED:
		return "Job already defined";
	case CONNECTION_ERROR:
		return "Connection error";
	case INTERNAL_ERROR:
		return "Internal error";
	case INVALID_INPUT:
		return "Invalid input";
	case INCOMPATIBLE_STATE:
		return "Incompatible state";
	case INCONSISTENT_DATA:
		return "Inconsistent data";
	}
#endif

	return "?";
}

/** */
extern int set_bp_map(void)
{
#ifdef HAVE_BG_FILES
	static rm_BGL_t *bg = NULL;
	int rc;
	rm_BP_t *my_bp = NULL;
	ba_bp_map_t *bp_map = NULL;
	int bp_num, i;
	char *bp_id = NULL;
	rm_location_t bp_loc;
	int number = 0;

	if(_bp_map_initialized)
		return 1;

	bp_map_list = list_create(_bp_map_list_del);

	if (!have_db2) {
		fatal("Can't access DB2 library, run from service node");
		return -1;
	}

	if (!getenv("DB2INSTANCE") || !getenv("VWSPATH")) {
		fatal("Missing DB2INSTANCE or VWSPATH env var."
			"Execute 'db2profile'");
		return -1;
	}
	
	if ((rc = bridge_get_bg(&bg)) != STATUS_OK) {
		error("bridge_get_BGL(): %d", rc);
		return -1;
	}
	
	if ((rc = bridge_get_data(bg, RM_BPNum, &bp_num)) != STATUS_OK) {
		error("bridge_get_data(RM_BPNum): %d", rc);
		bp_num = 0;
	}

	for (i=0; i<bp_num; i++) {

		if (i) {
			if ((rc = bridge_get_data(bg, RM_NextBP, &my_bp))
			    != STATUS_OK) {
				error("bridge_get_data(RM_NextBP): %d", rc);
				break;
			}
		} else {
			if ((rc = bridge_get_data(bg, RM_FirstBP, &my_bp))
			    != STATUS_OK) {
				error("bridge_get_data(RM_FirstBP): %d", rc);
				break;
			}
		}
		
		bp_map = (ba_bp_map_t *) xmalloc(sizeof(ba_bp_map_t));
		
		if ((rc = bridge_get_data(my_bp, RM_BPID, &bp_id))
		    != STATUS_OK) {
			xfree(bp_map);
			error("bridge_get_data(RM_BPID): %d", rc);
			continue;
		}

		if(!bp_id) {
			error("No BP ID was returned from database");
			continue;
		}
			
		if ((rc = bridge_get_data(my_bp, RM_BPLoc, &bp_loc))
		    != STATUS_OK) {
			xfree(bp_map);
			error("bridge_get_data(RM_BPLoc): %d", rc);
			continue;
		}
		
		bp_map->bp_id = xstrdup(bp_id);
		bp_map->coord[X] = bp_loc.X;
		bp_map->coord[Y] = bp_loc.Y;
		bp_map->coord[Z] = bp_loc.Z;
		number = atoi(bp_id+1);		
		if(DIM_SIZE[X] > bp_loc.X
		   && DIM_SIZE[Y] > bp_loc.Y
		   && DIM_SIZE[Z] > bp_loc.Z)
			ba_system_ptr->grid
				[bp_loc.X]
				[bp_loc.Y]
				[bp_loc.Z].phys_x = number / 100;
		
		list_push(bp_map_list, bp_map);
		
		free(bp_id);		
	}

	if ((rc = bridge_free_bg(bg)) != STATUS_OK)
		error("bridge_free_BGL(): %s", rc);	
	
#endif
	_bp_map_initialized = true;
	return 1;
	
}

static void _new_ba_node(ba_node_t *ba_node, int *coord, bool track_down_nodes)
{
	int i,j;
	uint16_t node_base_state = ba_node->state & NODE_STATE_BASE;
	
	if(((node_base_state != NODE_STATE_DOWN)
	   && !(ba_node->state & NODE_STATE_DRAIN)) || !track_down_nodes) 
		ba_node->used = false;

	for (i=0; i<BA_SYSTEM_DIMENSIONS; i++){
		ba_node->coord[i] = coord[i];
		
		for(j=0;j<NUM_PORTS_PER_NODE;j++) {
			ba_node->axis_switch[i].int_wire[j].used = 0;	
			if(i!=X) {
				if(j==3 || j==4) 
					ba_node->axis_switch[i].int_wire[j].
						used = 1;	
			}
			ba_node->axis_switch[i].int_wire[j].port_tar = j;
		}
	}
}

static void _create_ba_system(void)
{
	int x;
	int coord[BA_SYSTEM_DIMENSIONS];
				
#ifdef HAVE_BG
	int y,z;
	ba_system_ptr->grid = (ba_node_t***) 
		xmalloc(sizeof(ba_node_t**) * DIM_SIZE[X]);
#else
	ba_system_ptr->grid = (ba_node_t*) 
		xmalloc(sizeof(ba_node_t) * DIM_SIZE[X]);
#endif
	for (x=0; x<DIM_SIZE[X]; x++) {
#ifdef HAVE_BG
		ba_system_ptr->grid[x] = (ba_node_t**) 
			xmalloc(sizeof(ba_node_t*) * DIM_SIZE[Y]);
		for (y=0; y<DIM_SIZE[Y]; y++) {
			ba_system_ptr->grid[x][y] = (ba_node_t*) 
				xmalloc(sizeof(ba_node_t) * DIM_SIZE[Z]);
			for (z=0; z<DIM_SIZE[Z]; z++){
				coord[X] = x;
				coord[Y] = y;
				coord[Z] = z;
				_new_ba_node(&ba_system_ptr->grid[x][y][z], 
					     coord, true);
			}
		}
#else
		coord[X] = x;
		_new_ba_node(&ba_system_ptr->grid[x], coord, true);
#endif
	}
}

/** */
static void _delete_ba_system(void)
{
#ifdef HAVE_BG
	int x=0;
	int y;
#endif
	if (!ba_system_ptr){
		return;
	}
	
	if(ba_system_ptr->grid) {
#ifdef HAVE_BG
		for (x=0; x<DIM_SIZE[X]; x++) {
			for (y=0; y<DIM_SIZE[Y]; y++)
				xfree(ba_system_ptr->grid[x][y]);
			
			xfree(ba_system_ptr->grid[x]);
		}
#endif
		
		
		xfree(ba_system_ptr->grid);
	}
	xfree(ba_system_ptr);
}

static void _delete_path_list(void *object)
{
	ba_path_switch_t *path_switch = (ba_path_switch_t *)object;

	if (path_switch) {
		xfree(path_switch);
	}
	return;
}

/** 
 * algorithm for finding match
 */
static int _find_match(ba_request_t *ba_request, List results)
{
	int x=0;
#ifdef HAVE_BG
	int start[BA_SYSTEM_DIMENSIONS] = {0,0,0};
#else
	int start[BA_SYSTEM_DIMENSIONS] = {0};
#endif
	ba_node_t *ba_node = NULL;
	char *name=NULL;
	int startx = (start[X]-1);
	
	if(startx == -1)
		startx = DIM_SIZE[X]-1;
	if(ba_request->start_req) {
		if(ba_request->start[X]>DIM_SIZE[X] 
#ifdef HAVE_BG
		   || ba_request->start[Y]>DIM_SIZE[Y]
		   || ba_request->start[Z]>DIM_SIZE[Z]
#endif
			)
			return 0;
		for(x=0;x<BA_SYSTEM_DIMENSIONS;x++) {
			start[x] = ba_request->start[x];
		}
	}
	x=0;
	
	if(ba_request->geometry[X]>DIM_SIZE[X] 
#ifdef HAVE_BG
	   || ba_request->geometry[Y]>DIM_SIZE[Y]
	   || ba_request->geometry[Z]>DIM_SIZE[Z]
#endif
		)
#ifdef HAVE_BG
		if(!_check_for_options(ba_request))
#endif
			return 0;

#ifdef HAVE_BG
start_again:
#endif
	x=0;
	if(x == startx)
		x = startx-1;
	while(x!=startx) {
		x++;
		debug3("finding %d%d%d try %d",
		       ba_request->geometry[X],
#ifdef HAVE_BG
		       ba_request->geometry[Y],
		       ba_request->geometry[Z],
#endif
		       x);
#ifdef HAVE_BG
	new_node:
#endif
		debug2("starting at %d%d%d",
		       start[X]
#ifdef HAVE_BG
		       , start[Y],
		       start[Z]
#endif
			);
		
		ba_node = &ba_system_ptr->
			grid[start[X]]
#ifdef HAVE_BG
			[start[Y]]
			[start[Z]]
#endif
			;

		if (!_node_used(ba_node, ba_request->geometry)) {
			debug3("trying this node %d%d%d %d%d%d %d",
			       start[X], start[Y], start[Z],
			       ba_request->geometry[X],
			       ba_request->geometry[Y],
			       ba_request->geometry[Z], 
			       ba_request->conn_type);
			name = set_bg_block(results,
					    start, 
					    ba_request->geometry, 
					    ba_request->conn_type);
			if(name) {
				ba_request->save_name = xstrdup(name);
				xfree(name);
				return 1;
			}
			
			if(results) {
				remove_block(results, color_count);
				list_delete_all(results,
						&empty_null_destroy_list, "");
			}
			if(ba_request->start_req) 
				goto requested_end;
			//exit(0);
			debug2("trying something else");
			
		}
		
#ifdef HAVE_BG
		
		if((DIM_SIZE[Z]-start[Z]-1)
		   >= ba_request->geometry[Z])
			start[Z]++;
		else {
			start[Z] = 0;
			if((DIM_SIZE[Y]-start[Y]-1)
			   >= ba_request->geometry[Y])
				start[Y]++;
			else {
				start[Y] = 0;
				if ((DIM_SIZE[X]-start[X]-1)
				    >= ba_request->geometry[X])
					start[X]++;
				else {
					if(ba_request->size == 1)
						goto requested_end;
					if(!_check_for_options(ba_request))
						return 0;
					else {
						start[X]=0;
						start[Y]=0;
						start[Z]=0;
						goto start_again;
					}
				}
			}
		}
		goto new_node;
#endif
	}							
requested_end:
	debug("can't allocate");
	
	return 0;
}

/* bool _node_used(ba_node_t* ba_node, int geometry,  */
static bool _node_used(ba_node_t* ba_node, int *geometry)
{
	int i=0;
	ba_switch_t* ba_switch = NULL;
	
	/* if we've used this node in another block already */
	if (!ba_node || ba_node->used) {
		debug3("node used");
		return true;
	}
	/* if we've used this nodes switches completely in another 
	   block already */
	for(i=0;i<1;i++) {
		if(geometry[i]>1) {
			ba_switch = &ba_node->axis_switch[i];
			
			if(ba_switch->int_wire[3].used 
			   && ba_switch->int_wire[5].used) {
				debug3("switch in use dim %d!",i);
				return true;
			}
		}
	}
		
	return false;

}


static void _switch_config(ba_node_t* source, ba_node_t* target, int dim, 
			   int port_src, int port_tar)
{
	ba_switch_t* config = NULL, *config_tar = NULL;
	int i;

	if (!source || !target)
		return;
	
	config = &source->axis_switch[dim];
	config_tar = &target->axis_switch[dim];
	for(i=0;i<BA_SYSTEM_DIMENSIONS;i++) {
		/* Set the coord of the source target node to the target */
		config->ext_wire[port_src].node_tar[i] = target->coord[i];
	
		/* Set the coord of the target back to the source */
		config_tar->ext_wire[port_tar].node_tar[i] = source->coord[i];
	}

	/* Set the port of the source target node to the target */
	config->ext_wire[port_src].port_tar = port_tar;
	
	/* Set the port of the target back to the source */
	config_tar->ext_wire[port_tar].port_tar = port_src;
}

static int _set_external_wires(int dim, int count, ba_node_t* source, 
			       ba_node_t* target)
{
#ifdef HAVE_BG_FILES
	rm_BGL_t *bg = NULL;
	int rc;
	int i;
	rm_wire_t *my_wire = NULL;
	rm_port_t *my_port = NULL;
	char *wire_id = NULL;
	char from_node[5];
	char to_node[5];
	int from_port, to_port;
	int wire_num;
	int *coord;

	if (!have_db2) {
		error("Can't access DB2 library, run from service node");
		return -1;
	}
	
	if ((rc = bridge_get_bg(&bg)) != STATUS_OK) {
		error("bridge_get_BGL(): %d", rc);
		return -1;
	}
		
	if (bg == NULL) 
		return -1;
	
	if ((rc = bridge_get_data(bg, RM_WireNum, &wire_num)) != STATUS_OK) {
		error("bridge_get_data(RM_BPNum): %d", rc);
		wire_num = 0;
	}
	/* find out system wires on each bp */
	
	for (i=0; i<wire_num; i++) {
		
		if (i) {
			if ((rc = bridge_get_data(bg, RM_NextWire, &my_wire))
			    != STATUS_OK) {
				error("bridge_get_data(RM_NextWire): %d", rc);
				break;
			}
		} else {
			if ((rc = bridge_get_data(bg, RM_FirstWire, &my_wire))
			    != STATUS_OK) {
				error("bridge_get_data(RM_FirstWire): %d", rc);
				break;
			}
		}
		if ((rc = bridge_get_data(my_wire, RM_WireID, &wire_id))
		    != STATUS_OK) {
			error("bridge_get_data(RM_FirstWire): %d", rc);
			break;
		}
		
		if(!wire_id) {
			error("No Wire ID was returned from database");
			continue;
		}

		if(wire_id[7] != '_') 
			continue;
		switch(wire_id[0]) {
		case 'X':
			dim = X;
			break;
		case 'Y':
			dim = Y;
			break;
		case 'Z':
			dim = Z;
			break;
		}
		if(strlen(wire_id)<12) {
			error("Wire_id isn't correct %s",wire_id);
			continue;
		}
		strncpy(from_node, wire_id+2, 4);
		strncpy(to_node, wire_id+8, 4);
		
		free(wire_id);
		
		from_node[4] = '\0';
		to_node[4] = '\0';
		if ((rc = bridge_get_data(my_wire, RM_WireFromPort, &my_port))
		    != STATUS_OK) {
			error("bridge_get_data(RM_FirstWire): %d", rc);
			break;
		}
		if ((rc = bridge_get_data(my_port, RM_PortID, &from_port))
		    != STATUS_OK) {
			error("bridge_get_data(RM_PortID): %d", rc);
			break;
		}
		if ((rc = bridge_get_data(my_wire, RM_WireToPort, &my_port))
		    != STATUS_OK) {
			error("bridge_get_data(RM_WireToPort): %d", rc);
			break;
		}
		if ((rc = bridge_get_data(my_port, RM_PortID, &to_port))
		    != STATUS_OK) {
			error("bridge_get_data(RM_PortID): %d", rc);
			break;
		}

		coord = find_bp_loc(from_node);
		if(!coord) {
			error("1 find_bp_loc: bpid %s not known", from_node);
			continue;
		}
		
		if(coord[X]>=DIM_SIZE[X] 
		   || coord[Y]>=DIM_SIZE[Y]
		   || coord[Z]>=DIM_SIZE[Z]) {
			error("got coord %d%d%d greater than system dims "
			      "%d%d%d",
			      coord[X],
			      coord[Y],
			      coord[Z],
			      DIM_SIZE[X],
			      DIM_SIZE[Y],
			      DIM_SIZE[Z]);
			continue;
		}
		source = &ba_system_ptr->
			grid[coord[X]][coord[Y]][coord[Z]];
		coord = find_bp_loc(to_node);
		if(!coord) {
			error("2 find_bp_loc: bpid %s not known", to_node);
			continue;
		}
		if(coord[X]>=DIM_SIZE[X] 
		   || coord[Y]>=DIM_SIZE[Y]
		   || coord[Z]>=DIM_SIZE[Z]) {
			error("got coord %d%d%d greater than system dims "
			      "%d%d%d",
			      coord[X],
			      coord[Y],
			      coord[Z],
			      DIM_SIZE[X],
			      DIM_SIZE[Y],
			      DIM_SIZE[Z]);
			continue;
		}
		target = &ba_system_ptr->
			grid[coord[X]][coord[Y]][coord[Z]];
		_switch_config(source, 
			       target, 
			       dim, 
			       _port_enum(from_port),
			       _port_enum(to_port));	
		
		debug2("dim %d from %d%d%d %d -> %d%d%d %d",
		       dim,
		       source->coord[X], source->coord[Y], source->coord[Z],
		       _port_enum(from_port),
		       target->coord[X], target->coord[Y], target->coord[Z],
		       _port_enum(to_port));
	}
	if ((rc = bridge_free_bg(bg)) != STATUS_OK)
		error("bridge_free_BGL(): %s", rc);
	
#else

	_switch_config(source, source, dim, 0, 0);
	_switch_config(source, source, dim, 1, 1);
	if(dim!=X) {
		_switch_config(source, target, dim, 2, 5);
		_switch_config(source, source, dim, 3, 3);
		_switch_config(source, source, dim, 4, 4);
		return 1;
	}
	/* always 2->5 of next. If it is the last
		   it will go to the first.*/
		
	
#ifdef HAVE_BG
	_switch_config(source, target, dim, 2, 5);
	if(count == 0 || count==4) {
		/* 0 and 4th Node */
		/* 3->4 of next */
		_switch_config(source, target, dim, 3, 4);
		/* 4 is not in use */
		//_switch_config(source, source, dim, 4, 4);
	} else if( count == 1 || count == 5) {
		/* 1st and 5th Node */
		/* 3 is not in use */
		//_switch_config(source, source, dim, 3, 3);
	} else if(count == 2) {
		/* 2nd Node */
		/* make sure target is the last node */
		target = &ba_system_ptr->grid[DIM_SIZE[X]-1]
			[source->coord[Y]]
			[source->coord[Z]];
		/* 3->4 of last */
		_switch_config(source, target, dim, 3, 4);
		/* 4->3 of last */
		_switch_config(source, target, dim, 4, 3);
	} else if(count == 3) {
		/* 3rd Node */
		/* make sure target is the next to last node */
		target = &ba_system_ptr->grid[DIM_SIZE[X]-2]
			[source->coord[Y]]
			[source->coord[Z]];
		/* 3->4 of next to last */
		_switch_config(source, target, dim, 3, 4);
		/* 4->3 of next to last */
		_switch_config(source, target, dim, 4, 3);
	}

	if(DIM_SIZE[X] <= 4) {
		/* 4 X dim fixes for wires */
		
		if(count == 2) {
			/* 2 not in use */
			_switch_config(source, source, dim, 2, 2);
		} else if(count == 3) {
			/* 5 not in use */
			_switch_config(source, source, dim, 5, 5);
		}
	} else if(DIM_SIZE[X] != 8) {
		fatal("Do don't have a config to do this BG system.");
	}
#else
	if(count == 0)
		_switch_config(source, source, dim, 5, 5);
	else if(count < DIM_SIZE[X]-1)
		_switch_config(source, target, dim, 2, 5);
	else
		_switch_config(source, source, dim, 2, 2);
	_switch_config(source, source, dim, 3, 3);
	_switch_config(source, source, dim, 4, 4);
#endif /* HAVE_BG */
#endif /* HAVE_BG_FILES */
	return 1;
}
				
static char *_set_internal_wires(List nodes, int size, int conn_type)
{
	ba_node_t* ba_node[size+1];
	int count=0, i, set=0;
	int *start = NULL;
	int *end = NULL;
	char *name;
	ListIterator itr;
	hostlist_t hostlist;
	char temp_name[4];

	if(!nodes)
		return NULL;

	name = xmalloc(BUFSIZE);
	hostlist = hostlist_create(NULL);
	itr = list_iterator_create(nodes);
	while((ba_node[count] = (ba_node_t*) list_next(itr))) {
		snprintf(temp_name, sizeof(temp_name), "%d%d%d", 
			 ba_node[count]->coord[X],
			 ba_node[count]->coord[Y],
			 ba_node[count]->coord[Z]);
		debug3("name = %s", temp_name);
		count++;
		hostlist_push(hostlist, temp_name);
	}
	list_iterator_destroy(itr);
		
	start = ba_node[0]->coord;
	end = ba_node[count-1]->coord;	
	hostlist_ranged_string(hostlist, BUFSIZE, name);
	hostlist_destroy(hostlist);

	for(i=0;i<count;i++) {
		if(!ba_node[i]->used) {
			ba_node[i]->used=1;
			if(ba_node[i]->letter == '.') {
				ba_node[i]->letter = letters[color_count%62];
				ba_node[i]->color = colors[color_count%6];
				debug3("count %d setting letter = %c "
				       "color = %d",
				       color_count,
				       ba_node[i]->letter,
				       ba_node[i]->color);
				set=1;
			}
		} else {
			debug("No network connection to create "
			      "bgblock containing %s", name);
			debug("Use smap to define bgblocks in "
			      "bluegene.conf");
			xfree(name);
			return NULL;
		}
	}

	if(conn_type == SELECT_TORUS)
		for(i=0;i<count;i++) {
			_set_one_dim(start, end, ba_node[i]->coord);
		}

	if(set)
		color_count++;		

	return name;
}				

static int _find_x_path(List results, ba_node_t *ba_node, 
	int *start, int *first, int *geometry, 
	int found, int conn_type) 
{
	ba_switch_t *curr_switch = NULL; 
	ba_switch_t *next_switch = NULL; 
	
	int port_tar = 0;
	int source_port=0;
	int target_port=0;
	int broke = 0, not_first = 0;
	int ports_to_try[2] = {3,5};
	int *node_tar = NULL;
	int i = 0;
	ba_node_t *next_node = NULL;
	ba_node_t *check_node = NULL;
	int highest_phys_x = geometry[X] - start[X];
	
	ListIterator itr;

	if(!ba_node)
		return 0;

	if(!source_port) {
		target_port=1;
		ports_to_try[0] = 4;
		ports_to_try[1] = 2;
			
	}
	curr_switch = &ba_node->axis_switch[X];
	if(geometry[X] == 1) {
		goto found_one;
	}
	debug3("found - %d",found);
	for(i=0;i<2;i++) {
		/* check to make sure it isn't used */
		if(!curr_switch->int_wire[ports_to_try[i]].used) {
			/* looking at the next node on the switch 
			   and it's port we are going to */
			node_tar = curr_switch->
				ext_wire[ports_to_try[i]].node_tar;
			port_tar = curr_switch->
				ext_wire[ports_to_try[i]].port_tar;

			/* check to see if we are back at the start of the
			   block */
			if((node_tar[X] == 
			    start[X] && 
			    node_tar[Y] == 
			    start[Y] && 
			    node_tar[Z] == 
			    start[Z])) {
				broke = 1;
				goto broke_it;
			}
			/* check to see if the port points to itself */
			if((node_tar[X] == 
			    ba_node->coord[X] && 
			    node_tar[Y] == 
			    ba_node->coord[Y] && 
			    node_tar[Z] == 
			    ba_node->coord[Z])) {
				continue;
			}
			/* check to see if I am going to a place I have
			   already been before */
			itr = list_iterator_create(results);
			while((next_node = (ba_node_t*) list_next(itr))) {
				debug3("looking at %d%d%d and %d%d%d",
				       next_node->coord[X],
				       next_node->coord[Y],
				       next_node->coord[Z],
				       node_tar[X],
				       node_tar[Y],
				       node_tar[Z]);
				if((node_tar[X] == next_node->coord[X] && 
				    node_tar[Y] == next_node->coord[Y] && 
				    node_tar[Z] == next_node->coord[Z])) {
					not_first = 1;
					break;
				}				
			}
			list_iterator_destroy(itr);
			if(not_first && found<DIM_SIZE[X]) {
				debug2("already been there before");
				not_first = 0;
				continue;
			} 
			not_first = 0;
				
		broke_it:
			next_node = &ba_system_ptr->
				grid[node_tar[X]]
#ifdef HAVE_BG
				[node_tar[Y]]
				[node_tar[Z]]
#endif
				;
			next_switch = &next_node->axis_switch[X];

 			if((conn_type == SELECT_MESH) 
			   && (found == (geometry[X]))) {
				debug2("we found the end of the mesh");
				return 1;
			}
			debug3("Broke = %d Found = %d geometry[X] = %d",
			       broke, found, geometry[X]);
			debug3("Next Phys X %d Highest X %d",
			       next_node->phys_x, highest_phys_x);
			if(next_node->phys_x >= highest_phys_x) {
				debug3("looking for a passthrough");
				if(best_path)
					list_destroy(best_path);
				best_path = list_create(_delete_path_list);
				if(path)
					list_destroy(path);
				path = list_create(_delete_path_list);
	
				_find_passthrough(curr_switch,
						  0,
						  results,
						  X,
						  0,
						  highest_phys_x);
				if(best_count < BEST_COUNT_INIT) {
					debug2("yes found next free %d", 
					       best_count);
					node_tar = _set_best_path();
					next_node = &ba_system_ptr->
						grid[node_tar[X]]
#ifdef HAVE_BG
						[node_tar[Y]]
						[node_tar[Z]]
#endif
						;
					next_switch = 
						&next_node->axis_switch[X];
					
#ifdef HAVE_BG
					debug2("found %d looking at "
					       "%d%d%d going to %d%d%d %d",
					       found,
					       ba_node->coord[X],
					       ba_node->coord[Y],
					       ba_node->coord[Z],
					       node_tar[X],
					       node_tar[Y],
					       node_tar[Z],
					       port_tar);
#endif		
					list_append(results, next_node);
					found++;
					if(_find_x_path(results, next_node, 
							start, first, geometry,
							found, conn_type)) {
						return 1;
					} else {
						found--;
						_reset_the_path(curr_switch, 0,
								1, X);
						_remove_node(results, 
							     next_node->coord);
						return 0;
					}
				}
			}			

			if(broke && (found == geometry[X])) {
				goto found_path;
			} else if(found == geometry[X]) {
				debug2("finishing the torus!");
				if(best_path)
					list_destroy(best_path);
				best_path = list_create(_delete_path_list);
				if(path)
					list_destroy(path);
				path = list_create(_delete_path_list);
				_finish_torus(curr_switch, 
					      0, 
					      results, 
					      X, 
					      0, 
					      start);
				if(best_count < BEST_COUNT_INIT) {
					debug2("Found a best path with %d "
					       "steps.", best_count);
					_set_best_path();
					return 1;
				} else {
					return 0;
				}
			} else if(broke) {
				broke = 0;
				continue;
			}

			if (!_node_used(next_node, geometry)) {
#ifdef HAVE_BG
				debug2("found %d looking at %d%d%d "
				       "%d going to %d%d%d %d",
				       found,
				       ba_node->coord[X],
				       ba_node->coord[Y],
				       ba_node->coord[Z],
				       ports_to_try[i],
				       node_tar[X],
				       node_tar[Y],
				       node_tar[Z],
				       port_tar);
#endif
				itr = list_iterator_create(results);
				while((check_node = 
				       (ba_node_t*) list_next(itr))) {
					if((node_tar[X] == 
					    check_node->coord[X] && 
					    node_tar[Y] == 
					    check_node->coord[Y] && 
					    node_tar[Z] == 
					    check_node->coord[Z])) {
						break;
					}
				}
				list_iterator_destroy(itr);
				if(!check_node) {
#ifdef HAVE_BG
					debug2("add %d%d%d",
					       next_node->coord[X],
					       next_node->coord[Y],
					       next_node->coord[Z]);
#endif					       
					list_append(results, next_node);
				} else {
#ifdef HAVE_BG
					debug2("Hey this is already added "
					       "%d%d%d",
					       node_tar[X],
					       node_tar[Y],
					       node_tar[Z]);
#endif
					continue;
				}
				found++;
				
				if(!_find_x_path(results, next_node, 
						 start, first, geometry, 
						 found, conn_type)) {
					_remove_node(results,
						     next_node->coord);
					found--;
					continue;
				} else {
				found_path:
#ifdef HAVE_BG
					debug2("added node %d%d%d %d %d -> "
					       "%d%d%d %d %d",
					       ba_node->coord[X],
					       ba_node->coord[Y],
					       ba_node->coord[Z],
					       source_port,
					       ports_to_try[i],
					       node_tar[X],
					       node_tar[Y],
					       node_tar[Z],
					       port_tar,
					       target_port);
#endif					
				found_one:			
					if(geometry[X] != 1) {
						curr_switch->
							int_wire
							[source_port].used = 1;
						curr_switch->
							int_wire
							[source_port].port_tar
							= ports_to_try[i];
						curr_switch->
							int_wire
							[ports_to_try[i]].used
							= 1;
						curr_switch->
							int_wire
							[ports_to_try[i]].
							port_tar = source_port;
					
						next_switch->
							int_wire[port_tar].used
							= 1;
						next_switch->
							int_wire
							[port_tar].port_tar
							= target_port;
						next_switch->
							int_wire
							[target_port].used = 1;
						next_switch->
							int_wire
							[target_port].port_tar
							= port_tar;
					}
					return 1;

				}
			} 			
		}
	}

	debug2("couldn't find path");
	return 0;
}

static int _find_x_path2(List results, ba_node_t *ba_node, 
			 int *start, int *first, int *geometry, 
			 int found, int conn_type) 
{
	ba_switch_t *curr_switch = NULL; 
	ba_switch_t *next_switch = NULL; 
	
	int port_tar = 0;
	int source_port=0;
	int target_port=0;
	int broke = 0, not_first = 0;
	int ports_to_try[2] = {3,5};
	int *node_tar = NULL;
	int i = 0;
	ba_node_t *next_node = NULL;
	ba_node_t *check_node = NULL;
	
	ListIterator itr;
	
	if(!ba_node)
		return 0;

	if(!source_port) {
		target_port=1;
		ports_to_try[0] = 2;
		ports_to_try[1] = 4;
			
	}
	curr_switch = &ba_node->axis_switch[X];
	if(geometry[X] == 1) {
		goto found_one;
	}
	debug2("found - %d",found);
	for(i=0;i<2;i++) {
		/* check to make sure it isn't used */
		if(!curr_switch->int_wire[ports_to_try[i]].used) {
			node_tar = curr_switch->
				ext_wire[ports_to_try[i]].node_tar;
			port_tar = curr_switch->
				ext_wire[ports_to_try[i]].port_tar;
			if((node_tar[X] == 
			    start[X] && 
			    node_tar[Y] == 
			    start[Y] && 
			    node_tar[Z] == 
			    start[Z])) {
				broke = 1;
				goto broke_it;
			}
			if((node_tar[X] == 
			    ba_node->coord[X] && 
			    node_tar[Y] == 
			    ba_node->coord[Y] && 
			    node_tar[Z] == 
			    ba_node->coord[Z])) {
				continue;
			}
			itr = list_iterator_create(results);
			while((next_node = (ba_node_t*) list_next(itr))) {
				if((node_tar[X] == 
				    next_node->coord[X] && 
				    node_tar[Y] == 
				    next_node->coord[Y] && 
				    node_tar[Z] == 
				    next_node->coord[Z])) {
					not_first = 1;
					break;
				}
				
			}
			list_iterator_destroy(itr);
			if(not_first && found<DIM_SIZE[X]) {
				not_first = 0;
				continue;
			} 
			not_first = 0;
				
		broke_it:
			next_node = &ba_system_ptr->
				grid[node_tar[X]]
#ifdef HAVE_BG
				[node_tar[Y]]
				[node_tar[Z]]
#endif
				;

			next_switch = &next_node->axis_switch[X];
		
			
 			if((conn_type == SELECT_MESH) 
			   && (found == (geometry[X]))) {
				debug2("we found the end of the mesh");
				return 1;
			}
			debug3("Broke = %d Found = %d geometry[X] = %d",
			       broke, found, geometry[X]);
			if(broke && (found == geometry[X])) {
				goto found_path;
			} else if(found == geometry[X]) {
				debug2("finishing the torus!");
				if(best_path)
					list_destroy(best_path);
				best_path = list_create(_delete_path_list);
				if(path)
					list_destroy(path);
				path = list_create(_delete_path_list);
				_finish_torus(curr_switch, 
					      0, 
					      results, 
					      X, 
					      0, 
					      start);
				if(best_count < BEST_COUNT_INIT) {
					debug2("Found a best path with %d "
					       "steps.", best_count);
					_set_best_path();
					return 1;
				} else {
					return 0;
				}
			} else if(broke) {
				broke = 0;
				continue;
			}

			if (!_node_used(next_node, geometry)) {
#ifdef HAVE_BG
				debug2("found %d looking at %d%d%d "
				       "%d going to %d%d%d %d",
				       found,
				       ba_node->coord[X],
				       ba_node->coord[Y],
				       ba_node->coord[Z],
				       ports_to_try[i],
				       node_tar[X],
				       node_tar[Y],
				       node_tar[Z],
				       port_tar);
#endif
				itr = list_iterator_create(results);
				while((check_node = 
				       (ba_node_t*) list_next(itr))) {
					if((node_tar[X] == 
					    check_node->coord[X] && 
					    node_tar[Y] == 
					    check_node->coord[Y] && 
					    node_tar[Z] == 
					    check_node->coord[Z])) {
						break;
					}
				}
				list_iterator_destroy(itr);
				if(!check_node) {
#ifdef HAVE_BG
					debug2("add %d%d%d",
					       next_node->coord[X],
					       next_node->coord[Y],
					       next_node->coord[Z]);
#endif					       
					list_append(results, next_node);
				} else {
#ifdef HAVE_BG
					debug2("Hey this is already added "
					       "%d%d%d",
					       node_tar[X],
					       node_tar[Y],
					       node_tar[Z]);
#endif
					continue;
				}
				found++;
				
				if(!_find_x_path2(results, next_node, 
						 start, first, geometry, 
						 found, conn_type)) {
					_remove_node(results,
						     next_node->coord);
					found--;
					continue;
				} else {
				found_path:
#ifdef HAVE_BG
					debug2("added node %d%d%d %d %d -> "
					       "%d%d%d %d %d",
					       ba_node->coord[X],
					       ba_node->coord[Y],
					       ba_node->coord[Z],
					       source_port,
					       ports_to_try[i],
					       node_tar[X],
					       node_tar[Y],
					       node_tar[Z],
					       port_tar,
					       target_port);
#endif					
				found_one:			
					if(geometry[X] != 1) {
						curr_switch->
							int_wire
							[source_port].used = 1;
						curr_switch->
							int_wire
							[source_port].port_tar
							= ports_to_try[i];
						curr_switch->
							int_wire
							[ports_to_try[i]].used
							= 1;
						curr_switch->
							int_wire
							[ports_to_try[i]].
							port_tar = source_port;
					
						next_switch->
							int_wire[port_tar].used
							= 1;
						next_switch->
							int_wire
							[port_tar].port_tar
							= target_port;
						next_switch->
							int_wire
							[target_port].used = 1;
						next_switch->
							int_wire
							[target_port].port_tar
							= port_tar;
					}
					return 1;
				}
			} 			
		}
	}
#ifdef HAVE_BG
	debug2("looking for the next free node starting at %d%d%d",
	       ba_node->coord[X],
	       ba_node->coord[Y],
	       ba_node->coord[Z]);
#endif

	if(best_path)
		list_destroy(best_path);
	best_path = list_create(_delete_path_list);
	if(path)
		list_destroy(path);
	path = list_create(_delete_path_list);
	
	_find_next_free_using_port_2(curr_switch, 
				     0, 
				     results, 
				     X, 
				     0);
	if(best_count < BEST_COUNT_INIT) {
		debug2("yes found next free %d", best_count);
		node_tar = _set_best_path();

		next_node = &ba_system_ptr->
			grid[node_tar[X]]
#ifdef HAVE_BG
			[node_tar[Y]]
			[node_tar[Z]]
#endif
			;

		next_switch = &next_node->axis_switch[X];
		
#ifdef HAVE_BG
		debug2("found %d looking at %d%d%d going to %d%d%d %d",
		       found,
		       ba_node->coord[X],
		       ba_node->coord[Y],
		       ba_node->coord[Z],
		       node_tar[X],
		       node_tar[Y],
		       node_tar[Z],
		       port_tar);
#endif		
		list_append(results, next_node);
		found++;
		if(_find_x_path2(results, next_node, 
				start, first, geometry, found, conn_type)) {
			return 1;
		} else {
			found--;
			_reset_the_path(curr_switch, 0, 1, X);
			_remove_node(results, next_node->coord);
			debug2("couldn't finish the path off this one");
		}
	} 
	
	debug2("couldn't find path 2");
	return 0;
}

static int _remove_node(List results, int *node_tar)
{
	ListIterator itr;
	ba_node_t *ba_node = NULL;
	
	itr = list_iterator_create(results);
	while((ba_node = (ba_node_t*) list_next(itr))) {
		
#ifdef HAVE_BG
		if(node_tar[X] == ba_node->coord[X] 
		   && node_tar[Y] == ba_node->coord[Y] 
		   && node_tar[Z] == ba_node->coord[Z]) {
			debug2("removing %d%d%d from list",
			       node_tar[X],
			       node_tar[Y],
			       node_tar[Z]);
			list_remove (itr);
			break;
		}
#else
		if(node_tar[X] == ba_node->coord[X]) {
			debug2("removing %d from list",
			       node_tar[X]);
			list_remove (itr);
			break;
		}
#endif
	}
	list_iterator_destroy(itr);
	return 1;
}

static int _find_next_free_using_port_2(ba_switch_t *curr_switch, 
					int source_port, 
					List nodes, 
					int dim, 
					int count) 
{
	ba_switch_t *next_switch = NULL; 
	ba_path_switch_t *path_add = 
		(ba_path_switch_t *) xmalloc(sizeof(ba_path_switch_t));
	ba_path_switch_t *path_switch = NULL;
	ba_path_switch_t *temp_switch = NULL;
	int port_tar;
	int target_port = 0;
	int port_to_try = 2;
	int *node_tar= curr_switch->ext_wire[0].node_tar;
	int *node_src = curr_switch->ext_wire[0].node_tar;
	int used = 0;
	int broke = 0;
	ba_node_t *ba_node = NULL;
	
	ListIterator itr;
	static bool found = false;

	path_add->geometry[X] = node_src[X];
#ifdef HAVE_BG
	path_add->geometry[Y] = node_src[Y];
	path_add->geometry[Z] = node_src[Z];
#endif
	path_add->dim = dim;
	path_add->in = source_port;

	if(count>=best_count)
		goto return_0;
	
	itr = list_iterator_create(nodes);
	while((ba_node = (ba_node_t*) list_next(itr))) {
		
		if(node_tar[X] == ba_node->coord[X] 
#ifdef HAVE_BG
		   && node_tar[Y] == ba_node->coord[Y] 
		   && node_tar[Z] == ba_node->coord[Z] 
#endif
			)
		{
			broke = 1;
			break;
		}
	}
	list_iterator_destroy(itr);
	
	if(!broke && count>0 &&
	   !ba_system_ptr->grid[node_tar[X]]
#ifdef HAVE_BG
	   [node_tar[Y]]
	   [node_tar[Z]]
#endif
	   .used) {
		
#ifdef HAVE_BG
		debug2("this one not found %d%d%d",
		       node_tar[X],
		       node_tar[Y],
		       node_tar[Z]);
#endif		
		broke = 0;
				
		if((source_port%2))
			target_port=1;
		
		list_destroy(best_path);
		best_path = list_create(_delete_path_list);
		found = true;
		path_add->out = target_port;
		list_push(path, path_add);
		
		itr = list_iterator_create(path);
		while((path_switch = (ba_path_switch_t*) list_next(itr))){
		
			temp_switch = (ba_path_switch_t *) 
				xmalloc(sizeof(ba_path_switch_t));
			 
			temp_switch->geometry[X] = path_switch->geometry[X];
#ifdef HAVE_BG
			temp_switch->geometry[Y] = path_switch->geometry[Y];
			temp_switch->geometry[Z] = path_switch->geometry[Z];
#endif
			temp_switch->dim = path_switch->dim;
			temp_switch->in = path_switch->in;
			temp_switch->out = path_switch->out;
			list_append(best_path,temp_switch);
		}
		list_iterator_destroy(itr);
		best_count = count;
		return 1;
	} 

	used=0;
	if(!curr_switch->int_wire[port_to_try].used) {
		itr = list_iterator_create(path);
		while((path_switch = 
		       (ba_path_switch_t*) list_next(itr))){
				
			if(((path_switch->geometry[X] == node_src[X]) 
#ifdef HAVE_BG
			    && (path_switch->geometry[Y] 
				== node_src[Y])
			    && (path_switch->geometry[Z] 
				== node_tar[Z])
#endif
				   )) {
					
				if( path_switch->out
				    == port_to_try) {
					used = 1;
					break;
				}
			}
		}
		list_iterator_destroy(itr);
			
		if(curr_switch->
		   ext_wire[port_to_try].node_tar[X]
		   == curr_switch->ext_wire[0].node_tar[X]  
#ifdef HAVE_BG
		   && curr_switch->
		   ext_wire[port_to_try].node_tar[Y] 
		   == curr_switch->ext_wire[0].node_tar[Y] 
		   && curr_switch->
		   ext_wire[port_to_try].node_tar[Z] 
		   == curr_switch->ext_wire[0].node_tar[Z]
#endif
			) {
			used = 1;
		}
						
		if(!used) {
			port_tar = curr_switch->
				ext_wire[port_to_try].port_tar;
			node_tar = curr_switch->
				ext_wire[port_to_try].node_tar;
				
			next_switch = &ba_system_ptr->
				grid[node_tar[X]]
#ifdef HAVE_BG
				[node_tar[Y]]
				[node_tar[Z]]
#endif
				.axis_switch[X];
				
			count++;
			path_add->out = port_to_try;
			list_push(path, path_add);
			_find_next_free_using_port_2(next_switch, 
					port_tar, nodes,
					dim, count);
			while((temp_switch = list_pop(path)) != path_add){
				xfree(temp_switch);
				debug3("something here 1");
			}
		}
	}
return_0:
	xfree(path_add);
	return 0;
}

static int _find_passthrough(ba_switch_t *curr_switch, int source_port, 
			     List nodes, int dim, int count, int highest_phys_x) 
{
	ba_switch_t *next_switch = NULL; 
	ba_path_switch_t *path_add = 
		(ba_path_switch_t *) xmalloc(sizeof(ba_path_switch_t));
	ba_path_switch_t *path_switch = NULL;
	ba_path_switch_t *temp_switch = NULL;
	int port_tar;
	int target_port = 0;
	int ports_to_try[2] = {3,5};
	int *node_tar= curr_switch->ext_wire[0].node_tar;
	int *node_src = curr_switch->ext_wire[0].node_tar;
	int i;
	int used=0;
	int broke = 0;
	ba_node_t *ba_node = NULL;
	
	ListIterator itr;
	static bool found = false;

	path_add->geometry[X] = node_src[X];
#ifdef HAVE_BG
	path_add->geometry[Y] = node_src[Y];
	path_add->geometry[Z] = node_src[Z];
#endif
	path_add->dim = dim;
	path_add->in = source_port;
	
	if(count>=best_count) {
		xfree(path_add);
		return 0;
	}

	itr = list_iterator_create(nodes);
	while((ba_node = (ba_node_t*) list_next(itr))) {
		
#ifdef HAVE_BG
		if(node_tar[X] == ba_node->coord[X] 
		   && node_tar[Y] == ba_node->coord[Y] 
		   && node_tar[Z] == ba_node->coord[Z]) {
			broke = 1;
			break;
		}
#else
		if(node_tar[X] == ba_node->coord[X]) {
			broke = 1;
			break;
		}
#endif
		
	}
	list_iterator_destroy(itr);
	ba_node = &ba_system_ptr->
		grid[node_tar[X]]
#ifdef HAVE_BG
		[node_tar[Y]]
		[node_tar[Z]]
#endif
		;
	if(!broke && count>0
	   && !ba_node->used 
	   && (ba_node->phys_x < highest_phys_x)) {
		
		debug3("this one not found %d%d%d",
		       node_tar[X],
		       node_tar[Y],
		       node_tar[Z]);
		
		broke = 0;
				
		if((source_port%2))
			target_port=1;
		
		list_destroy(best_path);
		best_path = list_create(_delete_path_list);
		found = true;
		path_add->out = target_port;
		list_push(path, path_add);
		
		itr = list_iterator_create(path);
		while((path_switch = (ba_path_switch_t*) list_next(itr))){
		
			temp_switch = (ba_path_switch_t *) 
				xmalloc(sizeof(ba_path_switch_t));
			 
			temp_switch->geometry[X] = path_switch->geometry[X];
#ifdef HAVE_BG
			temp_switch->geometry[Y] = path_switch->geometry[Y];
			temp_switch->geometry[Z] = path_switch->geometry[Z];
#endif
			temp_switch->dim = path_switch->dim;
			temp_switch->in = path_switch->in;
			temp_switch->out = path_switch->out;
			list_append(best_path,temp_switch);
		}
		list_iterator_destroy(itr);
		best_count = count;
		return 1;
	} 

	if(source_port==0 || source_port==3 || source_port==5) {
		if(count==0) {
			ports_to_try[0] = 2;
			ports_to_try[1] = 4;	
		} else {
			ports_to_try[0] = 4;
			ports_to_try[1] = 2;	
		}
	}
			
	for(i=0;i<2;i++) {
		used=0;
		if(!curr_switch->int_wire[ports_to_try[i]].used) {
			itr = list_iterator_create(path);
			while((path_switch = 
			       (ba_path_switch_t*) list_next(itr))){
				
				if(((path_switch->geometry[X] == node_src[X]) 
#ifdef HAVE_BG
				    && (path_switch->geometry[Y] 
					== node_src[Y])
				    && (path_switch->geometry[Z] 
					== node_tar[Z])
#endif
					   )) {
					
					if( path_switch->out
					    == ports_to_try[i]) {
						used = 1;
						break;
					}
				}
			}
			list_iterator_destroy(itr);
			
			if(curr_switch->
			   ext_wire[ports_to_try[i]].node_tar[X]
			   == curr_switch->ext_wire[0].node_tar[X]  
#ifdef HAVE_BG
			   && curr_switch->
			   ext_wire[ports_to_try[i]].node_tar[Y] 
			   == curr_switch->ext_wire[0].node_tar[Y] 
			   && curr_switch->
			   ext_wire[ports_to_try[i]].node_tar[Z] 
			   == curr_switch->ext_wire[0].node_tar[Z]
#endif
				) {
				continue;
			}
						
			if(!used) {
				port_tar = curr_switch->
					ext_wire[ports_to_try[i]].port_tar;
				node_tar = curr_switch->
					ext_wire[ports_to_try[i]].node_tar;
				
				next_switch = &ba_system_ptr->
					grid[node_tar[X]]
#ifdef HAVE_BG
					[node_tar[Y]]
					[node_tar[Z]]
#endif
					.axis_switch[X];
				
				count++;
				path_add->out = ports_to_try[i];
				list_push(path, path_add);
				debug3("looking at this one "
				       "%d%d%d %d -> %d%d%d %d",
				       ba_node->coord[X],
				       ba_node->coord[Y],
				       ba_node->coord[Z],
				       ports_to_try[i],
				       node_tar[X],
				       node_tar[Y],
				       node_tar[Z],
				       port_tar);
		
				_find_passthrough(next_switch, port_tar, nodes,
						dim, count, highest_phys_x);
				while((temp_switch = list_pop(path)) 
				      != path_add){
					xfree(temp_switch);
					debug3("something here 2");
				}
			}
		}
	}
	xfree(path_add);
	return 0;
}

static int _finish_torus(ba_switch_t *curr_switch, int source_port, 
			 List nodes, int dim, int count, int *start) 
{
	ba_switch_t *next_switch = NULL; 
	ba_path_switch_t *path_add = 
		(ba_path_switch_t *) xmalloc(sizeof(ba_path_switch_t));
	ba_path_switch_t *path_switch = NULL;
	ba_path_switch_t *temp_switch = NULL;
	int port_tar;
	int target_port=0;
	int ports_to_try[2] = {3,5};
	int *node_tar= curr_switch->ext_wire[0].node_tar;
	int *node_src = curr_switch->ext_wire[0].node_tar;
	int i;
	int used=0;
	ListIterator itr;
	static bool found = false;

	path_add->geometry[X] = node_src[X];
#ifdef HAVE_BG
	path_add->geometry[Y] = node_src[Y];
	path_add->geometry[Z] = node_src[Z];
#endif
	path_add->dim = dim;
	path_add->in = source_port;

	if(count>=best_count) {
		xfree(path_add);
		return 0;
	}
	if(node_tar[X] == start[X] 
#ifdef HAVE_BG
	    && node_tar[Y] == start[Y] 
	    && node_tar[Z] == start[Z]
#endif
		) {
		
		if((source_port%2))
			target_port=1;
		if(!curr_switch->int_wire[target_port].used) {
			
			list_destroy(best_path);
			best_path = list_create(_delete_path_list);
			found = true;
			path_add->out = target_port;
			list_push(path, path_add);
			
			itr = list_iterator_create(path);
			while((path_switch = 
			       (ba_path_switch_t*) list_next(itr))){
				
				temp_switch = (ba_path_switch_t *) 
					xmalloc(sizeof(ba_path_switch_t));
				
				temp_switch->geometry[X] = 
					path_switch->geometry[X];
#ifdef HAVE_BG
				temp_switch->geometry[Y] = 
					path_switch->geometry[Y];
				temp_switch->geometry[Z] = 
					path_switch->geometry[Z];
#endif
				temp_switch->dim = path_switch->dim;
				temp_switch->in = path_switch->in;
				temp_switch->out = path_switch->out;
				list_append(best_path,temp_switch);
			}
			list_iterator_destroy(itr);
			best_count = count;
			return 1;
		} 
	}
	
	if(source_port==0 || source_port==3 || source_port==5) {
		ports_to_try[0] = 4;
		ports_to_try[1] = 2;		
	}
	
	for(i=0;i<2;i++) {
		used=0;
		if(!curr_switch->int_wire[ports_to_try[i]].used) {
			itr = list_iterator_create(path);
			while((path_switch = 
			       (ba_path_switch_t*) list_next(itr))){
				
				if(((path_switch->geometry[X] == node_src[X]) 
#ifdef HAVE_BG
				    && (path_switch->geometry[Y] 
					== node_src[Y])
				    && (path_switch->geometry[Z] 
					== node_tar[Z])
#endif
					)) {
					if( path_switch->out
					    == ports_to_try[i]) {
						used = 1;
						break;
					}
				}
			}
			list_iterator_destroy(itr);
			if((curr_switch->
			    ext_wire[ports_to_try[i]].node_tar[X] == 
			    curr_switch->ext_wire[0].node_tar[X] && 
			    curr_switch->
			    ext_wire[ports_to_try[i]].node_tar[Y] == 
			    curr_switch->ext_wire[0].node_tar[Y] && 
			    curr_switch->
			    ext_wire[ports_to_try[i]].node_tar[Z] == 
			    curr_switch->ext_wire[0].node_tar[Z])) {
				continue;
			}
			if(!used) {
				port_tar = curr_switch->
					ext_wire[ports_to_try[i]].port_tar;
				node_tar = curr_switch->
					ext_wire[ports_to_try[i]].node_tar;
				
				next_switch = &ba_system_ptr->
					grid[node_tar[X]]
#ifdef HAVE_BG
					[node_tar[Y]]
					[node_tar[Z]]
#endif
					.axis_switch[dim];
			
				
				count++;
				path_add->out = ports_to_try[i];
				list_push(path, path_add);
				_finish_torus(next_switch, port_tar, nodes,
						dim, count, start);
				while((temp_switch = list_pop(path)) 
				      != path_add){
					xfree(temp_switch);
					debug3("something here 3");
				} 
			}
		}
       }
       xfree(path_add);
       return 0;
}

static int *_set_best_path()
{
	ListIterator itr;
	ba_path_switch_t *path_switch = NULL;
	ba_switch_t *curr_switch = NULL; 
	int *geo = NULL;
	if(!best_path)
		return NULL;
	itr = list_iterator_create(best_path);
	while((path_switch = (ba_path_switch_t*) list_next(itr))) {
		if(passthrough)
			*passthrough = true;
#ifdef HAVE_BG
		debug3("mapping %d%d%d",path_switch->geometry[X],
		       path_switch->geometry[Y],
		       path_switch->geometry[Z]);
		if(!geo)
			geo = path_switch->geometry;
		curr_switch = &ba_system_ptr->
			grid
			[path_switch->geometry[X]]
			[path_switch->geometry[Y]]
			[path_switch->geometry[Z]].  
			axis_switch[path_switch->dim];
#else
		curr_switch = &ba_system_ptr->
			grid[path_switch->geometry[X]].
			axis_switch[path_switch->dim];
#endif
	
		curr_switch->int_wire[path_switch->in].used = 1;
		curr_switch->int_wire[path_switch->in].port_tar = 
			path_switch->out;
		curr_switch->int_wire[path_switch->out].used = 1;
		curr_switch->int_wire[path_switch->out].port_tar = 
			path_switch->in;
	}
	list_iterator_destroy(itr);

	best_count=BEST_COUNT_INIT;
	return geo;
}

static int _set_one_dim(int *start, int *end, int *coord)
{
	int dim;
	ba_switch_t *curr_switch = NULL; 
	
	for(dim=0;dim<BA_SYSTEM_DIMENSIONS;dim++) {
		if(start[dim]==end[dim]) {
			curr_switch = &ba_system_ptr->grid[coord[X]]
#ifdef HAVE_BG
				[coord[Y]]
				[coord[Z]]
#endif
				.axis_switch[dim];

			if(!curr_switch->int_wire[0].used 
			   && !curr_switch->int_wire[1].used) {
				curr_switch->int_wire[0].used = 1;
				curr_switch->int_wire[0].port_tar = 1;
				curr_switch->int_wire[1].used = 1;
				curr_switch->int_wire[1].port_tar = 0;
			}
		}
	}
	return 1;
}

static void _destroy_geo(void *object) {
	int *geo_ptr = (int *)object;
	xfree(geo_ptr);
}

//#define BUILD_EXE
#ifdef BUILD_EXE
/** */
int main(int argc, char** argv)
{
	ba_request_t *request = (ba_request_t*) xmalloc(sizeof(ba_request_t)); 
	log_options_t log_opts = LOG_OPTS_INITIALIZER;
	int debug_level = 6;

	List results;
//	List results2;
//	int i,j;
	log_opts.stderr_level  = debug_level;
	log_opts.logfile_level = debug_level;
	log_opts.syslog_level  = debug_level;
	
	log_alter(log_opts, LOG_DAEMON, 
		  "/dev/null");
	
	DIM_SIZE[X]=0;
	DIM_SIZE[Y]=0;
	DIM_SIZE[Z]=0;
	ba_init(NULL);
	init_wires(NULL);
						
	/* results = list_create(NULL); */
/* 	request->geometry[0] = 1; */
/* 	request->geometry[1] = 4; */
/* 	request->geometry[2] = 4; */
/* 	request->start[0] = 5; */
/* 	request->start[1] = 0; */
/* 	request->start[2] = 0; */
/* 	request->start_req = 1; */
/* 	request->size = 32; */
/* 	request->rotate = 0; */
/* 	request->elongate = 0; */
/* 	request->conn_type = SELECT_TORUS; */
/* 	new_ba_request(request); */
/* 	print_ba_request(request); */
/* 	if(!allocate_block(request, results)) { */
/*        		debug("couldn't allocate %d%d%d", */
/* 		       request->geometry[0], */
/* 		       request->geometry[1], */
/* 		       request->geometry[2]);	 */
/* 	} */
/* 	list_destroy(results); */

	results = list_create(NULL);
	request->geometry[0] = 1;
	request->geometry[1] = 1;
	request->geometry[2] = 1;
	request->start[0] = 0;
	request->start[1] = 0;
	request->start[2] = 0;
	request->start_req = 1;
	request->size = 1;
	request->rotate = 0;
	request->elongate = 0;
	request->conn_type = SELECT_TORUS;
	new_ba_request(request);
	print_ba_request(request);
	if(!allocate_block(request, results)) {
       		debug("couldn't allocate %d%d%d",
		       request->geometry[0],
		       request->geometry[1],
		       request->geometry[2]);	
	}
	list_destroy(results);

	results = list_create(NULL);
	request->geometry[0] = 1;
	request->geometry[1] = 1;
	request->geometry[2] = 1;
	request->start_req = 0;
	request->size = 1;
	request->conn_type = SELECT_TORUS;
	new_ba_request(request);
	print_ba_request(request);
	if(!allocate_block(request, results)) {
       		debug("couldn't allocate %d%d%d",
		       request->geometry[0],
		       request->geometry[1],
		       request->geometry[2]);
	}
	list_destroy(results);
	
/* 	results = list_create(NULL); */
/* 	request->geometry[0] = 4; */
/* 	request->geometry[1] = 4; */
/* 	request->geometry[2] = 4; */
/* 	//request->size = 2; */
/* 	request->conn_type = SELECT_TORUS; */
/* 	new_ba_request(request); */
/* 	print_ba_request(request); */
/* 	if(!allocate_block(request, results)) { */
/*        		printf("couldn't allocate %d%d%d\n", */
/* 		       request->geometry[0], */
/* 		       request->geometry[1], */
/* 		       request->geometry[2]); */
/* 	} */

/* 	results = list_create(NULL); */
/* 	request->geometry[0] = 1; */
/* 	request->geometry[1] = 4; */
/* 	request->geometry[2] = 4; */
/* 	//request->size = 2; */
/* 	request->conn_type = SELECT_TORUS; */
/* 	new_ba_request(request); */
/* 	print_ba_request(request); */
/* 	if(!allocate_block(request, results)) { */
/*        		printf("couldn't allocate %d%d%d\n", */
/* 		       request->geometry[0], */
/* 		       request->geometry[1], */
/* 		       request->geometry[2]); */
/* 	} */

	
	int dim,j;
	int x,y,z;
	int startx=0;
	int starty=0;
	int startz=0;
	int endx=DIM_SIZE[X];
	int endy=1;//DIM_SIZE[Y];
	int endz=1;//DIM_SIZE[Z];
	for(x=startx;x<endx;x++) {
		for(y=starty;y<endy;y++) {
			for(z=startz;z<endz;z++) {
				info("Node %d%d%d Used = %d Letter = %c",
				       x,y,z,ba_system_ptr->grid[x][y][z].used,
				       ba_system_ptr->grid[x][y][z].letter);
				for(dim=0;dim<1;dim++) {
					info("Dim %d",dim);
					ba_switch_t *wire =
						&ba_system_ptr->
						grid[x][y][z].axis_switch[dim];
					for(j=0;j<6;j++)
						info("\t%d -> %d -> %d%d%d %d "
						     "Used = %d",
						     j, wire->int_wire[j].
						     port_tar,
						     wire->ext_wire[
							     wire->int_wire[j].
							     port_tar].
						     node_tar[X],
						     wire->ext_wire[
							     wire->int_wire[j].
							     port_tar].
						     node_tar[Y],
						     wire->ext_wire[
							     wire->int_wire[j].
							     port_tar].
						     node_tar[Z],
						     wire->ext_wire[
							     wire->int_wire[j].
							     port_tar].
						     port_tar,
						     wire->int_wire[j].used);
				}
			}
		}
	}
	/* list_destroy(results); */

/* 	ba_fini(); */

/* 	delete_ba_request(request); */
	
	return 0;
}


#endif
