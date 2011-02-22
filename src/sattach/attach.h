/*****************************************************************************\
 *  geo_bitmaps.h - Functions used to manage multiple dimension bitmaps
 *                  especially for packing resources on a BlueGene system
 *****************************************************************************
 *  Copyright (C) 2011 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#ifndef _GEO_BITMAPS_H
#define _GEO_BITMAPS_H


typedef struct geo_table {
	uint16_t size;			/* Total object count */
	uint16_t *geometry;		/* Size in each dimension */
	struct geo_table *next_ptr;	/* Next geometry of this size */
} geo_table_t;

typedef struct system_geo {
	uint16_t dim_count;		/* Number of system dimensions */
	uint16_t *dim_size;		/* System size in each dimension */
	uint16_t total_size;		/* Total number of nodes in system */

	geo_table_t **geo_table_ptr;	/* Pointers to possible geometries.
					 * Index is request size */
	uint16_t geo_table_size;	/* Number of geo_table_t records */
} system_geo_t;

/*
 * Build a geo_table of possible unique geometries
 * IN/OUT my_system_geo - system geometry specification.
 *		Set dim_count and dim_size. Other fields should be NULL.
 *		This function will set total_size, geo_table_ptr, and
 *		geo_table_size.
 * Release memory using bg_free_geo_table().
 */
extern void bg_build_geo_table(system_geo_t *my_system_geo);

/*
 * Free memory allocated by bg_build_geo_table().
 * IN my_system_geo - System geometry specification.
 */
extern void bg_free_geo_table(system_geo_t *my_system_geo);

/*
 * Print the contents of all geo_table_t entries.
 */
extern void bg_print_geo_table(system_geo_t *my_system_geo);

/*
 * Print a linked list of geo_table_t entries.
 * IN geo_ptr - first geo_table entry to print
 * IN header - message header
 * IN my_system_geo - system geometry specification
 */
extern int  bg_geo_list_print(geo_table_t *geo_ptr, char *header,
			      system_geo_t *my_system_geo);

/*
 * Attempt to place a new allocation into an existing node state.
 * Do not rotate or change the requested geometry, but do attempt to place
 * it using all possible starting locations.
 *
 * IN node_bitmap - bitmap representing current system state, bits are set
 *                  for currently allocated nodes
 * OUT alloc_node_bitmap - bitmap representing where to place the allocation
 *                         set only if RET == SLURM_SUCCESS
 * IN geo_req - geometry required for the new allocation
 * OUT attempt_cnt - number of job placements attempted
 * IN my_system_geo - system geometry specification
 * RET - SLURM_SUCCESS if allocation can be made, otherwise SLURM_ERROR
 */
extern int  bg_geo_test_all(bitstr_t *node_bitmap,
			    bitstr_t **alloc_node_bitmap,
			    geo_table_t *geo_req, int *attempt_cnt,
			    system_geo_t *my_system_geo);

/*
 * Print the contents of a node map created by bg_node_map_alloc() or
 *	bg_geo_test_all(). Output may be in one-dimension or more depending
 *	upon configuration.
 * IN node_bitmap - bitmap representing current system state, bits are set
 *                  for currently allocated nodes
 * IN my_system_geo - system geometry specification
 */
extern void bg_node_map_print(bitstr_t *node_bitmap,
			      system_geo_t *my_system_geo);

/*
 * Allocate a multi-dimensional node bitmap. Use bg_node_map_free() to free
 * IN my_system_geo - system geometry specification
 */
extern bitstr_t *bg_node_map_alloc(system_geo_t *my_system_geo);

/*
 * Free a node map created by bg_node_map_alloc()
 * IN node_bitmap - bitmap of currently allocated nodes
 * IN my_system_geo - system geometry specification
 */
extern void bg_node_map_free(bitstr_t *node_bitmap,
			     system_geo_t *my_system_geo);

/*
 * Return the contents of the specified position in the bitmap
 * IN node_bitmap - bitmap of currently allocated nodes
 * IN full_offset - N-dimension zero-origin offset to test
 * IN my_system_geo - system geometry specification
 */
extern int bg_node_map_test(bitstr_t *node_bitmap, int *full_offset,
			    system_geo_t *my_system_geo);

/*
 * Set the contents of the specified position in the bitmap
 * IN node_bitmap - bitmap of currently allocated nodes
 * IN full_offset - N-dimension zero-origin offset to test
 * IN my_system_geo - system geometry specification
 */
extern void bg_node_map_set(bitstr_t *node_bitmap, int *full_offset,
			    system_geo_t *my_system_geo);

/*
 * Add a new allocation's node bitmap to that of the currently
 *	allocated bitmap
 * IN/OUT node_bitmap - bitmap of currently allocated nodes
 * IN alloc_bitmap - bitmap of nodes to be added fromtonode_bitmap
 * IN my_system_geo - system geometry specification
 */
extern void bg_node_map_add(bitstr_t *node_bitmap, bitstr_t *alloc_bitmap,
			    system_geo_t *my_system_geo);

/*
 * Remove a terminating allocation's node bitmap from that of the currently
 *	allocated bitmap
 * IN/OUT node_bitmap - bitmap of currently allocated nodes
 * IN alloc_bitmap - bitmap of nodes to be removed from node_bitmap
 * IN my_system_geo - system geometry specification
 */
extern void bg_node_map_rm(bitstr_t *node_bitmap, bitstr_t *alloc_bitmap,
			   system_geo_t *my_system_geo);

#endif	/* _GEO_BITMAPS_H */
