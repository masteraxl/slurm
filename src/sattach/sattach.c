/*****************************************************************************\
 *  sattach.c - Attach to a running job step.
 *****************************************************************************
 *  Copyright (C) 2006-2007 The Regents of the University of California.
 *  Copyright (C) 2008-2010 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher J. Morrone <morrone2@llnl.gov>
 *  CODE-OCEC-09-009. All rights reserved.
 *
 *  This file is part of SLURM, a resource management program.
 *  For details, see <https://computing.llnl.gov/linux/slurm/>.
 *  Please also read the included file: DISCLAIMER.
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <netinet/in.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

#include "src/common/bitstring.h"
#include "src/common/list.h"
#include "src/common/timers.h"
#include "src/common/xassert.h"
#include "src/common/xstring.h"
#include "src/common/xmalloc.h"

#define _DEBUG 1
#define BGQ_CNODE_DIM_CNT 5	/* Number of dimensions to manage */
#define BGQ_CNODE_CNT     512	/* Number of c-nodes in a midplane */
#define BGQ_GEO_TABLE_LEN 60	/* Number of possible geometries */

/* Number of elements in each dimension */
static int bgq_cnode_dim_size[BGQ_CNODE_DIM_CNT] = {4, 4, 4, 4, 2};


typedef struct geo_table {
	int size;
	int geo[BGQ_CNODE_DIM_CNT];
} geo_table_t;
geo_table_t geo_table[BGQ_GEO_TABLE_LEN];
int geo_table_cnt = 0;

/* Translate a 5-D coordinate into a 1-D offset in the cnode bitmap */
static void _bgl_cnode_xlate_5d_1d(int *offset_1d, int *offset_5d)
{
	int i, map_offset;

	xassert(offset_1d);
	xassert(offset_5d);
	map_offset = offset_5d[0];
	for (i = 1; i < BGQ_CNODE_DIM_CNT; i++) {
		map_offset *= bgq_cnode_dim_size[i];
		map_offset += offset_5d[i];
	}
	*offset_1d = map_offset;
}

/* Translate a 1-D offset in the cnode bitmap to a 5-D coordinate */
static void _bgl_cnode_xlate_1d_5d(int offset_1d, int *offset_5d)
{
	int i, map_offset;

	xassert(offset_5d);

	map_offset = offset_1d;
	for (i = BGQ_CNODE_DIM_CNT - 1; i >= 0; i--) {
		offset_5d[i] = map_offset % bgq_cnode_dim_size[i];
		map_offset /= bgq_cnode_dim_size[i];
	}
}

/* Allocate a bgq cnode map. Use bgq_cnode_map_free() to free */
extern bitstr_t *bgq_cnode_map_alloc(void)
{
	bitstr_t *cnode_map = bit_alloc( BGQ_CNODE_CNT);
	if (cnode_map == NULL)
		fatal("bit_alloc: malloc failure");
	return cnode_map;
}

/* Free a bgq cnode map created by bgq_cnode_map_alloc() */
extern void bgq_cnode_map_free(bitstr_t *cnode_bitmap)
{
	xassert(bit_size(cnode_bitmap) == BGQ_CNODE_CNT);
	bit_free(cnode_bitmap);
}

extern void bgq_cnode_map_set(bitstr_t *cnode_bitmap, int *offset_5d)
{
	int offset_1d;

	_bgl_cnode_xlate_5d_1d(&offset_1d, offset_5d);
	bit_set(cnode_bitmap, offset_1d);
}

/*
 * Add a new allocation's cnode bitmap to that of the midplane's currently
 *	allocated bitmap
 */
extern void bgq_cnode_map_add(bitstr_t *cnode_bitmap, bitstr_t *alloc_bitmap)
{
	xassert(bit_size(cnode_bitmap) == BGQ_CNODE_CNT);
	xassert(bit_size(alloc_bitmap) == BGQ_CNODE_CNT);
	bit_or(cnode_bitmap, alloc_bitmap);
}

/*
 * Remove a terminating allocation's cnode bitmap from that of the midplane's
 *	currently allocated bitmap
 */
extern void bgq_cnode_map_rm(bitstr_t *cnode_bitmap, bitstr_t *alloc_bitmap)
{
	xassert(bit_size(cnode_bitmap) == BGQ_CNODE_CNT);
	xassert(bit_size(alloc_bitmap) == BGQ_CNODE_CNT);
	bit_not(alloc_bitmap);
	bit_and(cnode_bitmap, alloc_bitmap);
	bit_not(alloc_bitmap);
}

/*
 * Print the contents of a bgq cnode map created by bgq_cnode_map_alloc() or
 *	bgq_cnode_map_test().
 */
extern void bgq_cnode_map_print(bitstr_t *cnode_bitmap)
{
	int i, j, offset[BGQ_CNODE_DIM_CNT];

	xassert(cnode_bitmap);
	xassert(bit_size(cnode_bitmap) == BGQ_CNODE_CNT);

	for (i = 0; i < BGQ_CNODE_CNT; i++) {
		if (bit_test(cnode_bitmap, i)) {
			char dim_buf[16], full_buf[64];
			full_buf[0] = '\0';
			_bgl_cnode_xlate_1d_5d(i, offset);
			for (j = 0; j < BGQ_CNODE_DIM_CNT; j++) {
				snprintf(dim_buf, sizeof(dim_buf),  "%d ",
				 	 offset[j]);
				strcat(full_buf, dim_buf);
			}
			info("%s", full_buf);
		}
	}
}

/*
 * Build bitmaps for all possible configurations of a specific c-node count
 */
extern int bgq_cnode_build_all(int *alloc_dims)
{
	int alloc_size = 1;	/* Total cnodes to be allocated */
	int map_cnt = 0;	/* Number of bitmaps built */
	int i, j;
	int start_offset[BGQ_CNODE_DIM_CNT], next_offset[BGQ_CNODE_DIM_CNT];
	int tmp_offset[BGQ_CNODE_DIM_CNT];
	bitstr_t *cnode_bitmap;

	xassert(alloc_dims);
	for (j = 0; j < BGQ_CNODE_DIM_CNT; j++)
		alloc_size *= alloc_dims[j];
	xassert(alloc_size >= 1);
	/* Start at location 00000 and move through all starting locations */
	for (j = 0; j < BGQ_CNODE_DIM_CNT; j++)
		start_offset[j] = 0;
	for (i = 0; i < BGQ_CNODE_CNT; i++) {
#if 0
		int axis;
		/* Cycle through rotations in each dimension. */
		/* NOTE: This logic is only working for allocations
		 * that have non-zero width in only one dimension. */
		for (axis = 0; axis < BGQ_CNODE_DIM_CNT; axis++) {
			int cnode_cnt;
			if (bgq_cnode_dim_size[axis] < alloc_size)
				continue;	/* too small */
			if ((bgq_cnode_dim_size[axis] == alloc_size) &&
			    (start_offset[axis] != 0))
				continue;	/* already filled, no shift */
			map_cnt++;
			cnode_bitmap = bgq_cnode_map_alloc();
			bgq_cnode_map_set(cnode_bitmap, start_offset);
			for (cnode_cnt = 1; cnode_cnt < alloc_size;
			     cnode_cnt++) {
				for (j = 0; j < BGQ_CNODE_DIM_CNT; j++)
					next_offset[j] = start_offset[j];
				next_offset[axis] += cnode_cnt;
				if (next_offset[axis] >=
				    bgq_cnode_dim_size[axis])
					next_offset[axis] = 0;
				bgq_cnode_map_set(cnode_bitmap, next_offset);
			}
			//if (map_cnt < 10) { bgq_cnode_map_print(cnode_bitmap); info(" "); }
			if (alloc_size == 1)	/* No need to rotate */
				break;
		}
#endif
		map_cnt++;
		cnode_bitmap = bgq_cnode_map_alloc();
		for (j = 0; j < BGQ_CNODE_DIM_CNT; j++)
			tmp_offset[j] = 0;
		while (1) {
			/* Compute location of next entry on the grid */
			for (j = 0; j < BGQ_CNODE_DIM_CNT; j++) {
				next_offset[j] = start_offset[j] +
						 tmp_offset[j];
				if (next_offset[j] >= bgq_cnode_dim_size[j]) {
					next_offset[j] -=
						bgq_cnode_dim_size[j];
				}
			}

			/* Add a point on the grid */
			bgq_cnode_map_set(cnode_bitmap, next_offset);

			/* Increment tmp_offset */
			for (j = 0; j < BGQ_CNODE_DIM_CNT; j++) {
				tmp_offset[j]++;
				if (tmp_offset[j] < alloc_dims[j])
					break;
				tmp_offset[j] = 0;
			}
			if (j >= BGQ_CNODE_DIM_CNT)
				break;

		}
		if (map_cnt < 4) { bgq_cnode_map_print(cnode_bitmap); info(" "); }

		/* Move to next starting location */
		for (j = 0; j < BGQ_CNODE_DIM_CNT; j++) {
			if (alloc_dims[j] == bgq_cnode_dim_size[j])
				continue;
			if (++start_offset[j] < bgq_cnode_dim_size[j])
				break;
			start_offset[j] = 0;
		}
		if (j >= BGQ_CNODE_DIM_CNT)
			break;
	}
	return map_cnt;
}

/*
 * Increment a geometry array, return false after reaching the last entry
 * NOTE: This is 1-origin and only generates unique number by insuring
 * that each number in the array is no larger than the number in the
 * prior dimension.
 */
static bool _incr_geo(int *geo)
{
	int dim, i;
	for (dim = BGQ_CNODE_DIM_CNT - 1; dim >= 1; dim--) {
		if ((geo[dim] < geo[dim-1]) &&
		    (geo[dim] < bgq_cnode_dim_size[dim])) {
			geo[dim]++;
			for (i = dim + 1; i < BGQ_CNODE_DIM_CNT; i++)
				geo[i] = 1;
			return true;
		}
	}

	if (++geo[0] <= bgq_cnode_dim_size[0]) {
		for (dim = 1; dim < BGQ_CNODE_DIM_CNT; dim++) {
			geo[dim] = 1;
		}
		return true;
	}
	
	return false;
}

/*
 * Swap two records in the geometry table
 */
static void _swap_geo_table(int inx1, int inx2)
{
	int i, tmp;

	tmp = geo_table[inx1].size;
	geo_table[inx1].size = geo_table[inx2].size;
	geo_table[inx2].size = tmp;
	for (i = 0; i < BGQ_CNODE_DIM_CNT; i++) {
		tmp = geo_table[inx1].geo[i];
		geo_table[inx1].geo[i] = geo_table[inx2].geo[i];
		geo_table[inx2].geo[i] = tmp;
	}
}

/*
 * Sort a geo_table in order of increasing node counts
 * Do bubble sort to preserve order
 */
static void _sort_geo_table(void)
{
	int i;
	bool sorted = false;

	do {
		sorted = true;
		for (i = 1; i < geo_table_cnt; i++) {
			if (geo_table[i-1].size > geo_table[i].size) {
				_swap_geo_table(i-1, i);
				sorted = false;
			}
		}
	} while (!sorted);
}

/*
 * Print the contents of geo_table
 */
static void _print_geo_table(void)
{
	int i;

	for (i = 0; i < geo_table_cnt; i++) {
		info("%d %d %d %d %d : %d", geo_table[i].geo[0],
		    geo_table[i].geo[1], geo_table[i].geo[2],
		    geo_table[i].geo[3], geo_table[i].geo[4],
		    geo_table[i].size);
	}
}

/*
 * Build a geo_table of possible unique geometries
 */
static void _build_geo_table(void)
{
	int dim, inx[BGQ_CNODE_DIM_CNT], product;

	for (dim = 0; dim < BGQ_CNODE_DIM_CNT; dim++)
		inx[dim] = 1;

	do {
		/* Store new value */
		product = 1;
		for (dim = 0; dim < BGQ_CNODE_DIM_CNT; dim++) {
			geo_table[geo_table_cnt].geo[dim] = inx[dim];
			product *= inx[dim];
		}
		geo_table[geo_table_cnt++].size = product;
		if (geo_table_cnt >= BGQ_GEO_TABLE_LEN)
			fatal("BGQ geometry table overflow");
	} while (_incr_geo(inx));	/* Generate next geometry */
}

static void _geo_list_destroy(void *x)
{
	geo_table_t *my_geo = (geo_table_t *) x;
	xfree(my_geo);
}

static int _geo_list_print(void *x, void *arg)
{
	geo_table_t *my_geo = (geo_table_t *) x;
	info("%d %d %d %d %d : %d", my_geo->geo[0], my_geo->geo[1],
	     my_geo->geo[2], my_geo->geo[3], my_geo->geo[4], my_geo->size);
	return 0;
}

/*
 * Find the possible geometries of a given size.
 * Search from the bottom of the list to prefer geometries that more
 * fully use a dimension (e.g. prefer "4,1,1,1,1" over "2,2,1,1,1").
 */
static List _find_dims(int node_cnt)
{
	int dim, i;
	geo_table_t *my_geo;
	List my_list = NULL;

	for (i = geo_table_cnt-1; i >= 0; i--) {
		if (geo_table[i].size < node_cnt)
			break;
		if (geo_table[i].size > node_cnt)
			continue;
		if (my_list == NULL) {
			my_list = list_create(_geo_list_destroy);
			if (my_list == NULL)
				fatal("list_create: malloc failure");
		}
		my_geo = xmalloc(sizeof(geo_table_t));
		my_geo->size = geo_table[i].size;
		for (dim = 0; dim < BGQ_CNODE_DIM_CNT; dim++) {
			my_geo->geo[dim] = geo_table[i].geo[dim];
		}
		list_append(my_list, my_geo);
	}
	return my_list;
}

///////////////////////////////////////////////////////////////////////////////
// Time to build various maps, with rotation
// Size	Count	USec
//    1   512     88
//    2  2304    220
//    3  2048    225
//    4   512    103
///////////////////////////////////////////////////////////////////////////////
// Time to build various maps, without rotation
// Size	Count	USec
//    1   512     85
//    2   512    114
//    3   512    128
//    4   128     64
//    6   128     65
//   16    32     56  (4x4)
//   16   128    109  (4x2x2)
//   64     8     56  (4x4x4)
///////////////////////////////////////////////////////////////////////////////

int sattach(int argc, char *argv[])
{
	DEF_TIMERS;
	char in_buf[32];
	int map_cnt = 0, node_cnt;
	int alloc_dims[BGQ_CNODE_DIM_CNT] = {1, 1, 1, 1, 1};
	List my_geo_list;

	START_TIMER;
	_build_geo_table();
	_sort_geo_table();
	_print_geo_table();
	END_TIMER;
	info("built table: %s", TIME_STR);

	while (1) {
		printf("node_count: ");
		if (gets(in_buf) == NULL)
			break;
		node_cnt = atoi(in_buf);
		if (node_cnt == 0)
			break;
		my_geo_list = _find_dims(node_cnt);
		if (my_geo_list) {
			list_for_each(my_geo_list, _geo_list_print, NULL);
			list_destroy(my_geo_list);
		} else {
			error("no geometry of size %d", node_cnt);
		}
	}
	exit(0);


	START_TIMER;
	map_cnt += bgq_cnode_build_all(alloc_dims);
	END_TIMER;
	info("built %d maps in %s", map_cnt, TIME_STR);

	exit(0);
}
