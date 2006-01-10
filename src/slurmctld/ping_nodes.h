/*****************************************************************************\
 *  ping_nodes.h - header to manager node ping
 *****************************************************************************
 *  Copyright (C) 2003 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by moe jette <jette1@llnl.gov> et. al.
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
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#ifndef _HAVE_PING_NODES_H
#define _HAVE_PING_NODES_H

#include "src/common/macros.h"

/*
 * is_ping_done - test if the last node ping cycle has completed.
 *	Use this to avoid starting a new set of ping requests before the 
 *	previous one completes
 * RET true if ping process is done, false otherwise
 */
extern bool is_ping_done (void);

/*
 * ping_begin - record that a ping cycle has begin. This can be called more 
 *	than once (for REQUEST_PING and simultaneous REQUEST_NODE_REGISTRATION 
 *	for selected nodes). Matching ping_end calls must be made for each 
 *	before is_ping_done returns true.
 */
extern void ping_begin (void);

/*
 * ping_end - record that a ping cycle has ended. This can be called more 
 *	than once (for REQUEST_PING and simultaneous REQUEST_NODE_REGISTRATION 
 *	for selected nodes). Matching ping_end calls must be made for each 
 *	before is_ping_done returns true.
 */
extern void ping_end (void);

/*
 * ping_nodes - check that all nodes and daemons are alive,  
 *	get nodes in UNKNOWN state to register
 */
extern void ping_nodes (void);

#endif /* !_HAVE_PING_NODES_H */
