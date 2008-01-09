/*****************************************************************************\
 *  node_acct.h - header to node accounting
 *****************************************************************************
 *  Copyright (C) 2007-2008 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
 *  UCRL-CODE-226842.
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

#ifndef _HAVE_NODE_ACCT_H
#define _HAVE_NODE_ACCT_H

#include <src/slurmctld/slurmctld.h>

/* Note that all nodes entered a DOWN state.
 * Typically done after a cold-start of SLURM */
extern void node_acct_all_down(char *reason);

/* Note that a node has entered a DOWN or DRAINED state */
extern void node_acct_down(struct node_record *node_ptr);

/* Note that a node has exited from a DOWN or DRAINED state */
extern void node_acct_up(struct node_record *node_ptr);

/* Note the total processor count in a cluster */
extern void node_acct_procs(char *cluster_name, uint32_t procs);

/* Note that the cluster is up and ready for work.
 * Generates a record of the cluster's processor count.
 * This should be executed whenever the cluser's processor count changes. */
extern void node_acct_ready(void);

#endif /* !_HAVE_NODE_ACCT_H */
