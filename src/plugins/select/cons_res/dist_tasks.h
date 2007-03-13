/*****************************************************************************\
 *  dist_tasks - Assign task count for each resources
 *
 *  $Id: dist_tasks.h,v 1.2 2006/10/31 19:31:31 palermo Exp $
 *****************************************************************************
 *  Copyright (C) 2006 Hewlett-Packard Development Company, L.P.
 *  Written by Susanne M. Balle, <susanne.balle@hp.com>
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

#ifndef _CONS_RES_DIST_TASKS_H
#define _CONS_RES_DIST_TASKS_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#if HAVE_STRING_H
#  include <string.h>
#endif

#include "select_cons_res.h"

int cr_exclusive_dist(struct select_cr_job *job,
		      const select_type_plugin_info_t cr_type);

int cr_dist(struct select_cr_job *job, int cyclic,
	    const select_type_plugin_info_t cr_type,
	    const uint16_t fast_schedule);

int cr_plane_dist(struct select_cr_job *job, 
		  const uint16_t plane_size,
		  const select_type_plugin_info_t cr_type);

int compute_c_b_task_dist(struct select_cr_job *job, 	    
			  const select_type_plugin_info_t cr_type,
			  const uint16_t fast_schedule);

#endif /* !_CONS_RES_DIST_TASKS_H */
