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
#include "src/common/timers.h"
#include "src/common/xassert.h"
#include "src/common/xstring.h"
#include "src/common/xmalloc.h"

#define BGQ_CNODE_DIM_CNT 5	/* Number of dimensions to manage */
#define BGQ_CNODE_CNT     512	/* Number of c-nodes in a midplane */

/* Number of elements in each dimension */
static int bgq_cnode_dim_size[BGQ_CNODE_DIM_CNT] = {4, 4, 4, 4, 2};

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
	int i, map_cnt = 0;
	int alloc_dims[BGQ_CNODE_DIM_CNT] = {1, 1, 1, 1, 1};

	if (argc < 2) {
		info("Usage: elements required in each dimension: X[Y[Z[A[B]]]]");
		exit(0);
	}

	for (i = 1; i < argc; i++) {
		alloc_dims[i-1] = atoi(argv[i]);
		if (alloc_dims[i-1] <= 0) {
			info("Invalid input: %s", argv[i]);
			exit(1);
		}
	}

	START_TIMER;
	map_cnt += bgq_cnode_build_all(alloc_dims);
	END_TIMER;
	info("built %d maps in %s", map_cnt, TIME_STR);

	exit(0);
}
