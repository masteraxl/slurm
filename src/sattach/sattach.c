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

#include <slurm/slurm_errno.h>

#include "src/common/bitstring.h"
#include "src/common/list.h"
#include "src/common/timers.h"
#include "src/common/xassert.h"
#include "src/common/xstring.h"
#include "src/common/xmalloc.h"

#define _DEBUG            1	/* Print debugging information */
#define DISPLAY_1D        1	/* Print allocation information using 1-D */
#define DISPLAY_5D        1	/* Print allocation information using 5-D */
#define MAX_ATTEMPT_CNT   1000	/* Maximum number of attempts to place a job.
				 * There are over 500,000 possible placements
				 * for some allocation sizes, which could be
				 * too slow to attempt. Stop after reaching
				 * this number of possible placements. */

typedef struct geo_table {
	int size;	/* Total size */
	int *geometry;	/* Size in each dimension */
	struct geo_table *next_ptr;
} geo_table_t;

typedef struct system_geo {
	int dim_count;	/* Number of system dimensions */
	int *dim_size;	/* System size in each dimension */
	int total_size;	/* Total number of elements in the system */
} system_geo_t;

geo_table_t **geo_table_size_ptr;
int geo_table_size = 0;

/* Local functions */
static void _build_geo_table(system_geo_t *my_system_geo);
static void _free_geo_table(system_geo_t *my_system_geo);
static int  _geo_list_print(geo_table_t *geo_ptr, char *header,
			    system_geo_t *my_system_geo);
static bool _incr_geo(int *geo, system_geo_t *my_system_geo);

/* Translate a 5-D coordinate into a 1-D offset in the cnode bitmap */
static void _bgl_cnode_xlate_5d_1d(int *offset_1d, int *offset_5d,
				   system_geo_t *my_system_geo)
{
	int i, map_offset;

	xassert(offset_1d);
	xassert(offset_5d);
	map_offset = offset_5d[0];
	for (i = 1; i < my_system_geo->dim_count; i++) {
		map_offset *= my_system_geo->dim_size[i];
		map_offset += offset_5d[i];
	}
	*offset_1d = map_offset;
}

#if DISPLAY_5D
/* Translate a 1-D offset in the cnode bitmap to a 5-D coordinate */
static void _bgl_cnode_xlate_1d_5d(int offset_1d, int *offset_5d,
				   system_geo_t *my_system_geo)
{
	int i, map_offset;

	xassert(offset_5d);

	map_offset = offset_1d;
	for (i = my_system_geo->dim_count - 1; i >= 0; i--) {
		offset_5d[i] = map_offset % my_system_geo->dim_size[i];
		map_offset /= my_system_geo->dim_size[i];
	}
}
#endif

/* Allocate a bgq cnode map. Use bgq_cnode_map_free() to free */
extern bitstr_t *bgq_cnode_map_alloc(system_geo_t *my_system_geo)
{
	bitstr_t *cnode_map = bit_alloc(my_system_geo->total_size);
	if (cnode_map == NULL)
		fatal("bit_alloc: malloc failure");
	return cnode_map;
}

/* Free a bgq cnode map created by bgq_cnode_map_alloc() */
extern void bgq_cnode_map_free(bitstr_t *cnode_bitmap,
			       system_geo_t *my_system_geo)
{
	xassert(bit_size(cnode_bitmap) == my_system_geo->total_size);
	bit_free(cnode_bitmap);
}

/*
 * Set the contents of the specified position in the bitmap
 */
extern void bgq_cnode_map_set(bitstr_t *cnode_bitmap, int *offset_5d,
			      system_geo_t *my_system_geo)
{
	int offset_1d;

	_bgl_cnode_xlate_5d_1d(&offset_1d, offset_5d, my_system_geo);
	bit_set(cnode_bitmap, offset_1d);
}

/*
 * Return the contents of the specified position in the bitmap
 */
extern int  bgq_cnode_map_test(bitstr_t *cnode_bitmap, int *offset_5d,
			       system_geo_t *my_system_geo)
{
	int offset_1d;

	_bgl_cnode_xlate_5d_1d(&offset_1d, offset_5d, my_system_geo);
	return bit_test(cnode_bitmap, offset_1d);
}

/*
 * Add a new allocation's cnode bitmap to that of the midplane's currently
 *	allocated bitmap
 */
extern void bgq_cnode_map_add(bitstr_t *cnode_bitmap, bitstr_t *alloc_bitmap,
			      system_geo_t *my_system_geo)
{
	xassert(bit_size(cnode_bitmap) == my_system_geo->total_size);
	xassert(bit_size(alloc_bitmap) == my_system_geo->total_size);
	bit_or(cnode_bitmap, alloc_bitmap);
}

/*
 * Remove a terminating allocation's cnode bitmap from that of the midplane's
 *	currently allocated bitmap
 */
extern void bgq_cnode_map_rm(bitstr_t *cnode_bitmap, bitstr_t *alloc_bitmap,
			     system_geo_t *my_system_geo)
{
	xassert(bit_size(cnode_bitmap) == my_system_geo->total_size);
	xassert(bit_size(alloc_bitmap) == my_system_geo->total_size);
	bit_not(alloc_bitmap);
	bit_and(cnode_bitmap, alloc_bitmap);
	bit_not(alloc_bitmap);
}

/*
 * Print the contents of a bgq cnode map created by bgq_cnode_map_alloc() or
 *	bgq_cnode_map_test().
 */
extern void bgq_cnode_map_print(bitstr_t *cnode_bitmap,
				system_geo_t *my_system_geo)
{
#if DISPLAY_1D
{
	char out_buf[256];
	bit_fmt(out_buf, sizeof(out_buf), cnode_bitmap);
	info("%s", out_buf);
}
#endif
#if DISPLAY_5D
{
	int i, j, offset[my_system_geo->dim_count];

	xassert(cnode_bitmap);
	xassert(bit_size(cnode_bitmap) == my_system_geo->total_size);

	for (i = 0; i < my_system_geo->total_size; i++) {
		if (bit_test(cnode_bitmap, i)) {
			char dim_buf[16], full_buf[64];
			full_buf[0] = '\0';
			_bgl_cnode_xlate_1d_5d(i, offset, my_system_geo);
			for (j = 0; j < my_system_geo->dim_count; j++) {
				snprintf(dim_buf, sizeof(dim_buf),  "%d ",
				 	 offset[j]);
				strcat(full_buf, dim_buf);
			}
			info("%s", full_buf);
		}
	}
}
#endif
}

/*
 * Attempt to place a new allocation into an existing c-node state.
 * Do not rotate or change the requested geometry, but do attempt to place
 * it using all possible starting locations.
 *
 * cnode_bitmap IN - bitmap representing current system state, bits are set
 *                   for currently allocated c-nodes
 * alloc_cnode_bitmap OUT - bitmap representing where to place the allocation
 *                          set only if RET == SLURM_SUCCESS
 * geo_req IN - geometry required for the new allocation
 * attempt_cnt OUT - number of job placements attempted
 * my_system_geo IN - system geometry specification
 * RET - SLURM_SUCCESS if allocation can be made, otherwise SLURM_ERROR
 */
extern int  bgq_cnode_test_all(bitstr_t *cnode_bitmap,
			       bitstr_t **alloc_cnode_bitmap,
			       geo_table_t *geo_req, int *attempt_cnt,
			       system_geo_t *my_system_geo)
{
	int rc = SLURM_ERROR;
	int i, j;
	int start_offset[my_system_geo->dim_count];
	int next_offset[my_system_geo->dim_count];
	int tmp_offset[my_system_geo->dim_count];
	bitstr_t *new_bitmap;

	xassert(cnode_bitmap);
	xassert(alloc_cnode_bitmap);
	xassert(geo_req);
	xassert(attempt_cnt);

	*attempt_cnt = 0;
	/* Start at location 00000 and move through all starting locations */
	for (j = 0; j < my_system_geo->dim_count; j++)
		start_offset[j] = 0;
	for (i = 0; i < my_system_geo->total_size; i++) {
		(*attempt_cnt)++;
		for (j = 0; j < my_system_geo->dim_count; j++)
			tmp_offset[j] = 0;
		while (1) {
			/* Compute location of next entry on the grid */
			for (j = 0; j < my_system_geo->dim_count; j++) {
				next_offset[j] = start_offset[j] +
						 tmp_offset[j];
				next_offset[j] %= my_system_geo->dim_size[j];
			}

			/* Test that point on the grid */
			if (bgq_cnode_map_test(cnode_bitmap, next_offset,
					       my_system_geo))
				break;

			/* Increment tmp_offset */
			for (j = 0; j < my_system_geo->dim_count; j++) {
				tmp_offset[j]++;
				if (tmp_offset[j] < geo_req->geometry[j])
					break;
				tmp_offset[j] = 0;
			}
			if (j >= my_system_geo->dim_count) {
				rc = SLURM_SUCCESS;
				break;
			}
		}
		if (rc == SLURM_SUCCESS)
			break;

		/* Move to next starting location */
		for (j = 0; j < my_system_geo->dim_count; j++) {
			if (geo_req->geometry[j] == my_system_geo->dim_size[j])
				continue;	/* full axis used */
			if (++start_offset[j] < my_system_geo->dim_size[j])
				break;		/* sucess */
			start_offset[j] = 0;	/* move to next dimension */
		}
		if (j >= my_system_geo->dim_count)
			return rc;		/* end of starting locations */
	}

	new_bitmap = bgq_cnode_map_alloc(my_system_geo);
	for (j = 0; j < my_system_geo->dim_count; j++)
		tmp_offset[j] = 0;
	while (1) {
		/* Compute location of next entry on the grid */
		for (j = 0; j < my_system_geo->dim_count; j++) {
			next_offset[j] = start_offset[j] + tmp_offset[j];
			if (next_offset[j] >= my_system_geo->dim_size[j])
				next_offset[j] -= my_system_geo->dim_size[j];
		}

		bgq_cnode_map_set(new_bitmap, next_offset, my_system_geo);

		/* Increment tmp_offset */
		for (j = 0; j < my_system_geo->dim_count; j++) {
			tmp_offset[j]++;
			if (tmp_offset[j] < geo_req->geometry[j])
				break;
			tmp_offset[j] = 0;
		}
		if (j >= my_system_geo->dim_count) {
			rc = SLURM_SUCCESS;
			break;
		}
	}
	*alloc_cnode_bitmap = new_bitmap;

	return rc;
}

/*
 * Increment a geometry index array, return false after reaching the last entry
 */
static bool _incr_geo(int *geo, system_geo_t *my_system_geo)
{
	int dim, i;

	for (dim = my_system_geo->dim_count - 1; dim >= 0; dim--) {
		if (geo[dim] < my_system_geo->dim_size[dim]) {
			geo[dim]++;
			for (i = dim + 1; i < my_system_geo->dim_count; i++)
				geo[i] = 1;
			return true;
		}
	}
	
	return false;
}

#if _DEBUG
/*
 * Print the contents of geo_table
 */
static void _print_geo_table(system_geo_t *my_system_geo)
{
	int i;
	geo_table_t *geo_ptr;

	xassert(geo_table_size_ptr);
	for (i = 1; i <= my_system_geo->total_size; i++) {
		geo_ptr = geo_table_size_ptr[i];
		while (geo_ptr) {
			_geo_list_print(geo_ptr, "", my_system_geo);
			geo_ptr = geo_ptr->next_ptr;
		}
	}
}
#endif

/*
 * Build a geo_table of possible unique geometries
 * Release memory using _free_geo_table()
 */
static void _build_geo_table(system_geo_t *my_system_geo)
{
	geo_table_t *geo_ptr;
	int dim, inx[my_system_geo->dim_count], product;

	if (geo_table_size_ptr)
		fatal("geo_table_size_ptr is already set");

	geo_table_size_ptr = xmalloc(sizeof(geo_table_t *) *
				     (my_system_geo->total_size + 1));
	for (dim = 0; dim < my_system_geo->dim_count; dim++)
		inx[dim] = 1;

	do {
		/* Store new value */
		geo_ptr = xmalloc(sizeof(geo_table_t));
		geo_ptr->geometry = xmalloc(sizeof(int) *
					    my_system_geo->dim_count);
		product = 1;
		for (dim = 0; dim < my_system_geo->dim_count; dim++) {
			geo_ptr->geometry[dim] = inx[dim];
			product *= inx[dim];
		}
		geo_ptr->size = product;
		xassert(product <= my_system_geo->total_size);
		geo_ptr->next_ptr = geo_table_size_ptr[product];
		geo_table_size_ptr[product] = geo_ptr;
		geo_table_size++;
	} while (_incr_geo(inx, my_system_geo));   /* Generate next geometry */
}

/*
 * Free memory allocated by _build_geo_table()
 */
static void _free_geo_table(system_geo_t *my_system_geo)
{
	geo_table_t *geo_ptr, *next_ptr;
	int i;

	for (i = 0; i <= my_system_geo->total_size; i++) {
		geo_ptr = geo_table_size_ptr[i];
		geo_table_size_ptr[i] = NULL;
		while (geo_ptr) {
			next_ptr = geo_ptr->next_ptr;
			xfree(geo_ptr->geometry);
			xfree(geo_ptr);
			geo_ptr = next_ptr;
		}
	}
	geo_table_size = 0;
	xfree(geo_table_size_ptr);
}

/*
 * Print a geo_table_t entry from a list.
 * IN geo_ptr - geo_table entry to print
 * IN header - message header
 * IN my_system_geo - system geometry specification
 */
static int  _geo_list_print(geo_table_t *geo_ptr, char *header,
			    system_geo_t *my_system_geo)
{
	int i;
	char dim_buf[16], full_buf[128];

	full_buf[0] = '\0';
	for (i = 0; i < my_system_geo->dim_count; i++) {
		snprintf(dim_buf, sizeof(dim_buf), "%d ", geo_ptr->geometry[i]);
		strcat(full_buf, dim_buf);
	}
	snprintf(dim_buf, sizeof(dim_buf), ": %d", geo_ptr->size);
	strcat(full_buf, dim_buf);
	info("%s%s", header, full_buf);

	return 0;
}

int sattach(int argc, char *argv[])
{
	DEF_TIMERS;
	bitstr_t *cnode_bitmap, *alloc_cnode_bitmap = NULL;
	char in_buf[32];
	int attempt_cnt = 0, total_attempt_cnt;
	int free_node_cnt, node_cnt, rc;
	geo_table_t *my_geo;
	system_geo_t *my_system_geo;

	/* Initialize system configuration */
	my_system_geo = xmalloc(sizeof(system_geo_t));
	my_system_geo->dim_count = 5;
	my_system_geo->dim_size = xmalloc(sizeof(int) * 5);
	my_system_geo->dim_size[0] = 4;
	my_system_geo->dim_size[1] = 4;
	my_system_geo->dim_size[2] = 4;
	my_system_geo->dim_size[3] = 4;
	my_system_geo->dim_size[4] = 2;
	my_system_geo->total_size = 4 * 4 * 4 * 4 * 2;

	START_TIMER;
	_build_geo_table(my_system_geo);
#if _DEBUG
	_print_geo_table(my_system_geo);
#endif
	END_TIMER;
	info("built table of size %d in time %s", geo_table_size, TIME_STR);

	cnode_bitmap = bgq_cnode_map_alloc(my_system_geo);
	while (1) {
		printf("node_count: ");
		if (gets(in_buf) == NULL)
			break;
		node_cnt = atoi(in_buf);
		if (node_cnt == 0)
			break;

		START_TIMER;
		if (node_cnt > my_system_geo->total_size) {
			info("more nodes requested than exist");
			continue;
		}
		free_node_cnt = bit_clear_count(cnode_bitmap);
		if (node_cnt > free_node_cnt) {
			END_TIMER;
			info("only %d free nodes remain, time %s",
			     free_node_cnt, TIME_STR);
			continue;
		}

		total_attempt_cnt = 0;
		rc = SLURM_ERROR;
		if ((node_cnt >= 1) &&
		    (node_cnt <= my_system_geo->total_size))
			my_geo = geo_table_size_ptr[node_cnt];
		else
			my_geo = NULL;
		while (my_geo) {
			_geo_list_print(my_geo, "testing to allocate: ",
					my_system_geo);
			rc = bgq_cnode_test_all(cnode_bitmap,
						&alloc_cnode_bitmap,
						my_geo, &attempt_cnt,
						my_system_geo);
			if (rc == SLURM_SUCCESS) {
				END_TIMER;
				info("allocation successful at:");
				bgq_cnode_map_print(alloc_cnode_bitmap,
						    my_system_geo);
				bgq_cnode_map_add(cnode_bitmap,
						  alloc_cnode_bitmap,
						  my_system_geo);
				bgq_cnode_map_free(alloc_cnode_bitmap,
						   my_system_geo);
				break;
			}
			total_attempt_cnt += attempt_cnt;
			if (total_attempt_cnt >= MAX_ATTEMPT_CNT)
				break;	/* Abandon effort */
			if (rc == SLURM_SUCCESS)
				break;
			my_geo = my_geo->next_ptr;
		}
		if (rc != SLURM_SUCCESS) {
			END_TIMER;
			info("allocation unsuccessful after %d attempts",
			     total_attempt_cnt);
		}
		info("full system allocation:");
		bgq_cnode_map_print(cnode_bitmap, my_system_geo);
		info("allocation time %s", TIME_STR);
	}

	/* Free memory */
	bgq_cnode_map_free(cnode_bitmap, my_system_geo);
	_free_geo_table(my_system_geo);
	xfree(my_system_geo->dim_size);
	xfree(my_system_geo);

	exit(0);
}
