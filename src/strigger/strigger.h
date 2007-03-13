/****************************************************************************\
 *  strigger.h - definitions used for strigger functions
 *****************************************************************************
 *  Copyright (C) 2007 The Regents of the University of California.
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
\****************************************************************************/

#ifndef _STRIGGER_H
#define _STRIGGER_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <slurm/slurm.h>
#include <src/common/macros.h>
#include <src/common/slurm_protocol_defs.h>

struct strigger_parameters {
	bool     job_fini;
	uint32_t job_id;
	bool     mode_set;
	bool     mode_get;
	bool     mode_clear;
	bool     node_down;
	char *   node_id;
	bool     node_up;
	int      offset;
	char *   program;
	bool     quiet;
	bool     reconfig;
	bool     time_limit;
	uint32_t trigger_id;
	uint32_t user_id;
	int      verbose;
};

extern struct strigger_parameters params;

extern void parse_command_line(int argc, char *argv[]);

#endif
