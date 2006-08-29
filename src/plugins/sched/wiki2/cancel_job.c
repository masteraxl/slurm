/*****************************************************************************\
 *  cancel_job.c - Process Wiki cancel job request
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#include "./msg.h"

#define TYPE_ADMIN   0
#define TYPE_TIMEOUT 1

/* RET 0 on success, -1 on failure */
extern int	cancel_job(char *cmd_ptr, int *err_code, char **err_msg)
{
	char *arg_ptr, *tmp_char;
	int cancel_type;
	uint32_t jobid;
	static char reply_msg[128];

	arg_ptr = strstr(cmd_ptr, "ARG=");
	if (arg_ptr == NULL) {
		*err_code = 300;
		*err_msg = "CANCELJOB lacks ARG";
		error("wiki: CANCELJOB lacks ARG");
		return -1;
	}
	jobid = strtol(arg_ptr+4, &tmp_char, 10);
	if (!isspace(tmp_char[0])) {
		*err_code = 300;
		*err_msg = "Invalid ARG value";
		error("wiki: CANCELJOB has invalid jobid");
		return -1;
	}

	if      (strstr(cmd_ptr, "TYPE=TIMEOUT") != 0)
		cancel_type = TYPE_TIMEOUT;
	else if (strstr(cmd_ptr, "TYPE=WALLCLOCK") != 0)
		cancel_type = TYPE_TIMEOUT;
	else if (strstr(cmd_ptr, "TYPE=ADMIN") != 0)
		cancel_type = TYPE_ADMIN;
	else if (strstr(cmd_ptr, "TYPE=") != 0) {
		*err_code = 300;
		*err_msg = "Invalid TYPE value";
		error("wiki: CANCELJOB has invalid TYPE");
		return -1;
	}
	
	if (sched_cancel_job(jobid) != SLURM_SUCCESS) {
		*err_code = 700;
		*err_msg = "failed to cancel job";
		error("wiki: failed to cancel job %u", jobid);
		return -1;
	}

	snprintf(reply_msg, sizeof(reply_msg), 
		"job %u cancelled successfully", jobid);
	*err_msg = reply_msg;
	return 0;
}
