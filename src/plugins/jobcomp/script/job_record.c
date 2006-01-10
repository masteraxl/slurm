/******************************************************************************
 *  job_record.c - Functions to alter script plugin job record structure.
 ******************************************************************************
 *  Produced at Center for High Performance Computing, North Dakota State
 *  University
 *  Written by Nathan Huff <nhuff@acm.org>
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

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#if HAVE_STDINT_H
#	include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#	include <inttypes.h>
#endif


#include <unistd.h>
#include <sys/types.h>
#include <slurm/slurm.h>
#include "job_record.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"

/* 
 * Create a new job_record containing the job completion information
 */
job_record job_record_create(uint32_t job_id, uint32_t user_id, char *job_name,
    char *job_state, char *partition, uint32_t limit, time_t start, time_t end,
    time_t submit, uint16_t batch_flag, char *node_list)
{
	job_record ret;

	ret = xmalloc(sizeof(struct job_record_));
	ret->job_id = job_id;
	ret->user_id = user_id;
	ret->job_name = xstrdup(job_name);
	ret->job_state = xstrdup(job_state);
	ret->partition = xstrdup(partition);
	ret->limit = limit;
	ret->start = start;
	ret->submit = submit;
	ret->batch_flag = batch_flag;
	ret->end = end;
	ret->node_list = xstrdup(node_list);
	return ret;
}

/*
 * Free the memory from a job_record structure
 */
void job_record_destroy(void *job) {
	job_record j = job;
	if (j == NULL) return;
	xfree(j->job_name);
	xfree(j->partition);
	xfree(j->node_list);
	xfree(j->job_state);
	xfree(j);
}
