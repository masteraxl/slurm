/*****************************************************************************\
 *  flatfile_jobcomp_process.c - functions the processing of
 *                               information from the flatfile jobcomp
 *                               database.
 *****************************************************************************
 *
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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
 *
 *  This file is patterned after jobcomp_linux.c, written by Morris Jette and
 *  Copyright (C) 2002 The Regents of the University of California.
\*****************************************************************************/
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>

#include "src/common/xstring.h"
#include "src/common/xmalloc.h"
#include "src/common/slurm_jobacct.h"
/* Map field names to positions */

/* slurmd uses "(uint32_t) -2" to track data for batch allocations
 * which have no logical jobsteps. */
#define BUFFER_SIZE 1024

typedef struct {
	char *name;
	char *val;
} flatfile_jobcomp_info_t;


static void _destroy_flatfile_jobcomp_info(void *object)
{
	flatfile_jobcomp_info_t *jobcomp_info =
		(flatfile_jobcomp_info_t *)object;
	if(jobcomp_info) {
		xfree(jobcomp_info);
	}
}


/* _open_log_file() -- find the current or specified log file, and open it
 *
 * IN:		Nothing
 * RETURNS:	Nothing
 *
 * Side effects:
 * 	- Sets opt_filein to the current system accounting log unless
 * 	  the user specified another file.
 */

static FILE *_open_log_file(char *logfile)
{
	FILE *fd = fopen(logfile, "r");
	if (fd == NULL) {
		perror(logfile);
		exit(1);
	}
	return fd;
}

static void _do_fdump(List job_info_list, int lc)
{
	flatfile_jobcomp_info_t *jobcomp_info = NULL;
	ListIterator itr = list_iterator_create(job_info_list);

	printf("\n------- Line %d -------\n", lc);
	while((jobcomp_info = list_next(itr))) {
		printf("%12s: %s\n", jobcomp_info->name, jobcomp_info->val);
	}
}

static jobcomp_job_rec_t *_parse_line(List job_info_list)
{
	ListIterator itr = NULL;
	flatfile_jobcomp_info_t *jobcomp_info = NULL;
	jobcomp_job_rec_t *job = xmalloc(sizeof(jobcomp_job_rec_t));
	char *temp = NULL;
	char *temp2 = NULL;

	itr = list_iterator_create(job_info_list);
	while((jobcomp_info = list_next(itr))) {
		if(!strcasecmp("JobID", jobcomp_info->name)) {
			job->jobid = atoi(jobcomp_info->val);
		} else if(!strcasecmp("Partition", jobcomp_info->name)) {
			job->partition = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("StartTime", jobcomp_info->name)) {
			job->start_time = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("EndTime", jobcomp_info->name)) {
			job->end_time = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("Userid", jobcomp_info->name)) {
			temp = strstr(jobcomp_info->val, "(");
			if(!temp) 
				job->uid = atoi(jobcomp_info->val);
			*temp++ = 0; 
			temp2 = temp;
			temp = strstr(temp, ")");
			if(!temp) {
				error("problem getting correct uid from %s",
				      jobcomp_info->val);
			} else {
				*temp = 0;
				job->uid = atoi(temp2);
				job->uid_name = xstrdup(jobcomp_info->val);
			}
		} else if(!strcasecmp("GroupId", jobcomp_info->name)) {
			temp = strstr(jobcomp_info->val, "(");
			if(!temp) 
				job->gid = atoi(jobcomp_info->val);
			*temp++ = 0; 
			temp2 = temp;
			temp = strstr(temp, ")");
			if(!temp) {
				error("problem getting correct gid from %s",
				      jobcomp_info->val);
			} else {
				*temp = 0;
				job->gid = atoi(temp2);
				job->gid_name = xstrdup(jobcomp_info->val);
			}
		} else if(!strcasecmp("Block_Id", jobcomp_info->name)) {
			job->blockid = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("Name", jobcomp_info->name)) {
			job->jobname = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("NodeList", jobcomp_info->name)) {
			job->nodelist = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("NodeCnt", jobcomp_info->name)) {
			job->node_cnt = atoi(jobcomp_info->val);
		} else if(!strcasecmp("MaxProcs", jobcomp_info->name)) {
			job->max_procs = atoi(jobcomp_info->val);
		} else if(!strcasecmp("JobState", jobcomp_info->name)) {
			job->state = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("Timelimit", jobcomp_info->name)) {
			job->timelimit = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("Connection", jobcomp_info->name)) {
			job->connection = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("reboot", jobcomp_info->name)) {
			job->reboot = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("rotate", jobcomp_info->name)) {
			job->rotate = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("geometry", jobcomp_info->name)) {
			job->geo = xstrdup(jobcomp_info->val);
		} else if(!strcasecmp("start", jobcomp_info->name)) {
			job->bg_start_point = xstrdup(jobcomp_info->val);
		} else {
			error("Unknown type %s: %s", jobcomp_info->name,
			      jobcomp_info->val);
		}
	}
	list_iterator_destroy(itr);
	
	return job;
}

extern void flatfile_jobcomp_process_get_jobs(List job_list, 
					      List selected_steps,
					      List selected_parts,
					      sacct_parameters_t *params)
{
	char line[BUFFER_SIZE];
	char *fptr = NULL;
	char *jobid = NULL;
	char *partition = NULL;
	FILE *fd = NULL;
	int lc = 0;
	jobcomp_job_rec_t *job = NULL;
	jobacct_selected_step_t *selected_step = NULL;
	char *selected_part = NULL;
	ListIterator itr = NULL;
	List job_info_list = NULL;
	flatfile_jobcomp_info_t *jobcomp_info = NULL;

	fd = _open_log_file(params->opt_filein);
	
	while (fgets(line, BUFFER_SIZE, fd)) {
		lc++;
		fptr = line;	/* break the record into NULL-
				   terminated strings */
		if(job_info_list) 
			list_destroy(job_info_list);
		jobid = NULL;
		partition = NULL;
		job_info_list = list_create(_destroy_flatfile_jobcomp_info);
		while(fptr) {
			jobcomp_info =
				xmalloc(sizeof(flatfile_jobcomp_info_t));
			list_append(job_info_list, jobcomp_info);
			jobcomp_info->name = fptr;
			fptr = strstr(fptr, "=");
			*fptr++ = 0;
			jobcomp_info->val = fptr;
			fptr = strstr(fptr, " ");
			if(!strcasecmp("JobId", jobcomp_info->name)) 
				jobid = jobcomp_info->val;
			else if(!strcasecmp("Partition",
					    jobcomp_info->name)) 
				partition = jobcomp_info->val;
			
			
			if(!fptr) {
				fptr = strstr(jobcomp_info->val, "\n");
				if (fptr)
					*fptr = 0;
				break; 
			} else {
				*fptr++ = 0;
				if(*fptr == '\n') {
					*fptr = 0;
					break;
				}
			}
		}
				
		if (list_count(selected_steps)) {
			if(!jobid) 
				continue;
			itr = list_iterator_create(selected_steps);
			while((selected_step = list_next(itr))) {
				if (strcmp(selected_step->job, jobid))
					continue;
				/* job matches */
				list_iterator_destroy(itr);
				goto foundjob;
			}
			list_iterator_destroy(itr);
			continue;	/* no match */
		}
	foundjob:
		
		if (list_count(selected_parts)) {
			if(!partition) 
				continue;
			itr = list_iterator_create(selected_parts);
			while((selected_part = list_next(itr))) 
				if (!strcasecmp(selected_part, partition)) {
					list_iterator_destroy(itr);
					goto foundp;
				}
			list_iterator_destroy(itr);
			continue;	/* no match */
		}
	foundp:
		
		if (params->opt_fdump) {
			_do_fdump(job_info_list, lc);
			continue;
		}
		
		
		job = _parse_line(job_info_list);
		
		if(job)
			list_append(job_list, job);
	}
	if(job_info_list) 
		list_destroy(job_info_list);
	
	if (ferror(fd)) {
		perror(params->opt_filein);
		exit(1);
	} 
	fclose(fd);

	return;
}

extern void flatfile_jobcomp_process_archive(List selected_parts,
					     sacct_parameters_t *params)
{
	info("No code to archive jobcomp.");
}
