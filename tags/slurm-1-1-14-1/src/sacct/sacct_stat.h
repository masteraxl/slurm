/*****************************************************************************\
 *  sacct_stat.h - header file for sacct
 *
 *  $Id: sacct.h 7541 2006-03-18 01:44:58Z da $
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>.
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
#ifndef _SACCT_STAT_H
#define _SACCT_STAT_H

#include "src/common/slurm_protocol_api.h"

typedef struct {
	uint16_t taskid; /* contains which task number it was on */
	uint32_t nodeid; /* contains which node number it was on */	
} jobacct_id_t;

typedef struct sacct_struct {
	uint32_t max_vsize; 
	jobacct_id_t max_vsize_id;
	float ave_vsize;
	uint32_t max_rss;
	jobacct_id_t max_rss_id;
	float ave_rss;
	uint32_t max_pages;
	jobacct_id_t max_pages_id;
	float ave_pages;
	float min_cpu;
	jobacct_id_t min_cpu_id;
	float ave_cpu;	
} sacct_t;

extern int sacct_stat(uint32_t jobid, uint32_t stepid);

#endif