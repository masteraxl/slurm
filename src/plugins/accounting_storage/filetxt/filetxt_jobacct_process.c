/*****************************************************************************\
 *  filetxt_jobacct_process.c - functions the processing of
 *                               information from the filetxt jobacct
 *                               storage.
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
#include "src/common/jobacct_common.h"
/* Map field names to positions */

/* slurmd uses "(uint32_t) -2" to track data for batch allocations
 * which have no logical jobsteps. */
#define BATCH_JOB_TIMESTAMP 0
#define EXPIRE_READ_LENGTH 10
#define MAX_RECORD_FIELDS 100

typedef struct expired_rec {  /* table of expired jobs */
	uint32_t job;
	time_t job_submit;
	char *line;
} expired_rec_t;

typedef struct header {
	uint32_t jobnum;
	char	*partition;
	char	*blockid;
	time_t 	job_submit;
	time_t 	timestamp;
	uint32_t uid;
	uint32_t gid;
	uint16_t rec_type;
} filetxt_header_t;

typedef struct {
	uint32_t job_start_seen,		/* useful flags */
		job_step_seen,
		job_terminated_seen,
		jobnum_superseded;	/* older jobnum was reused */
	filetxt_header_t header;
	uint16_t show_full;
	char	*nodes;
	char	*jobname;
	uint16_t track_steps;
	int32_t priority;
	uint32_t ncpus;
	uint32_t ntasks;
	enum job_states	status;
	int32_t	exitcode;
	uint32_t elapsed;
	time_t end;
	uint32_t tot_cpu_sec;
	uint32_t tot_cpu_usec;
	struct rusage rusage;
	sacct_t sacct;
	List    steps;
	char    *account;
	uint32_t requid;
} filetxt_job_rec_t;

typedef struct {
	filetxt_header_t   header;
	uint32_t	stepnum;	/* job's step number */
	char	        *nodes;
	char	        *stepname;
	enum job_states	status;
	int32_t	        exitcode;
	uint32_t	ntasks; 
	uint32_t        ncpus;
	uint32_t	elapsed;
	time_t          end;
	uint32_t	tot_cpu_sec;
	uint32_t        tot_cpu_usec;
	struct rusage   rusage;
	sacct_t         sacct;
	char            *account;
	uint32_t requid;
} filetxt_step_rec_t;

/* Fields common to all records */
enum {	F_JOB =	0,
	F_PARTITION,	
	F_JOB_SUBMIT,	
	F_TIMESTAMP,	
	F_UID,	
	F_GID,	
	F_BLOCKID,	
	F_RESERVED2,	
	F_RECTYPE,	
	HEADER_LENGTH
};

/* JOB_START fields */
enum {	F_JOBNAME = HEADER_LENGTH,
	F_TRACK_STEPS,		
	F_PRIORITY,	
	F_NCPUS,		
	F_NODES,
	F_JOB_ACCOUNT,
	JOB_START_LENGTH
};

/* JOB_STEP fields */
enum {	F_JOBSTEP = HEADER_LENGTH,
	F_STATUS,
	F_EXITCODE,
	F_NTASKS,
	F_STEPNCPUS,
	F_ELAPSED,
	F_CPU_SEC,
	F_CPU_USEC,
	F_USER_SEC,
	F_USER_USEC,
	F_SYS_SEC,
	F_SYS_USEC,
	F_RSS,
	F_IXRSS,
	F_IDRSS,
	F_ISRSS,
	F_MINFLT,
	F_MAJFLT,
	F_NSWAP,
	F_INBLOCKS,
	F_OUBLOCKS,
	F_MSGSND,
	F_MSGRCV,
	F_NSIGNALS,
	F_NVCSW,
	F_NIVCSW,
	F_MAX_VSIZE,
	F_MAX_VSIZE_TASK,
	F_AVE_VSIZE,
	F_MAX_RSS,
	F_MAX_RSS_TASK,
	F_AVE_RSS,
	F_MAX_PAGES,
	F_MAX_PAGES_TASK,
	F_AVE_PAGES,
	F_MIN_CPU,
	F_MIN_CPU_TASK,
	F_AVE_CPU,
	F_STEPNAME,
	F_STEPNODES,
	F_MAX_VSIZE_NODE,
	F_MAX_RSS_NODE,
	F_MAX_PAGES_NODE,
	F_MIN_CPU_NODE,
	F_STEP_ACCOUNT,
	F_STEP_REQUID,
	JOB_STEP_LENGTH
};

/* JOB_TERM / JOB_SUSPEND fields */
enum {	F_TOT_ELAPSED = HEADER_LENGTH,
	F_TERM_STATUS,
	F_JOB_REQUID,
	F_JOB_EXITCODE,
	JOB_TERM_LENGTH
};

static void _destroy_exp(void *object)
{
	expired_rec_t *exp_rec = (expired_rec_t *)object;
	if(exp_rec) {
		xfree(exp_rec->line);
		xfree(exp_rec);
	}
}

static void _free_filetxt_header(void *object)
{
	filetxt_header_t *header = (filetxt_header_t *)object;
	if(header) {
		xfree(header->partition);
#ifdef HAVE_BG
		xfree(header->blockid);
#endif
	}
}

static void _destroy_filetxt_job_rec(void *object)
{
	filetxt_job_rec_t *job = (filetxt_job_rec_t *)object;
	if (job) {
		if(job->steps)
			list_destroy(job->steps);
		_free_filetxt_header(&job->header);
		xfree(job->jobname);
		xfree(job->account);
		xfree(job->nodes);
		xfree(job);
	}
}

static void _destroy_filetxt_step_rec(void *object)
{
	filetxt_step_rec_t *step = (filetxt_step_rec_t *)object;
	if (step) {
		_free_filetxt_header(&step->header);
		xfree(step->stepname);
		xfree(step->nodes);
		xfree(step->account);
		xfree(step);
	}
}

static jobacct_step_rec_t *_create_jobacct_step_rec(
	filetxt_step_rec_t *filetxt_step)
{
	jobacct_step_rec_t *jobacct_step = create_jobacct_step_rec();
	
	jobacct_step->elapsed = filetxt_step->elapsed;
	jobacct_step->end = filetxt_step->header.timestamp;
	jobacct_step->exitcode = filetxt_step->exitcode;
	jobacct_step->ncpus = filetxt_step->ncpus;
	jobacct_step->nodes = xstrdup(filetxt_step->nodes);
	jobacct_step->requid = filetxt_step->requid;
	memcpy(&jobacct_step->sacct, &filetxt_step->sacct, sizeof(sacct_t));
	jobacct_step->start = filetxt_step->header.timestamp -
		jobacct_step->elapsed;
	jobacct_step->state = filetxt_step->status;
	jobacct_step->stepid = filetxt_step->stepnum;
	jobacct_step->stepname = xstrdup(filetxt_step->stepname);
	jobacct_step->sys_cpu_sec = filetxt_step->rusage.ru_stime.tv_sec;
	jobacct_step->sys_cpu_usec = filetxt_step->rusage.ru_stime.tv_usec;
	jobacct_step->tot_cpu_sec = filetxt_step->tot_cpu_sec;
	jobacct_step->tot_cpu_usec = filetxt_step->tot_cpu_usec;
	jobacct_step->user_cpu_sec = filetxt_step->rusage.ru_utime.tv_sec;
	jobacct_step->user_cpu_usec = filetxt_step->rusage.ru_utime.tv_usec;

	return jobacct_step;
}

static jobacct_job_rec_t *_create_jobacct_job_rec(
	filetxt_job_rec_t *filetxt_job)
{
	jobacct_job_rec_t *jobacct_job = create_jobacct_job_rec();
	ListIterator itr = NULL;
	filetxt_step_rec_t *filetxt_step = NULL;

	jobacct_job->associd = 0;
	jobacct_job->account = xstrdup(filetxt_job->account);
	jobacct_job->blockid = xstrdup(filetxt_job->header.blockid);
	jobacct_job->cluster = NULL;
	jobacct_job->elapsed = filetxt_job->elapsed;
	jobacct_job->eligible = filetxt_job->header.job_submit;
	jobacct_job->end = filetxt_job->header.timestamp;
	jobacct_job->exitcode = filetxt_job->exitcode;
	jobacct_job->gid = filetxt_job->header.gid;
	jobacct_job->jobid = filetxt_job->header.jobnum;
	jobacct_job->jobname = xstrdup(filetxt_job->jobname);
	jobacct_job->partition = xstrdup(filetxt_job->header.partition);
	jobacct_job->ncpus = filetxt_job->ncpus;
	jobacct_job->nodes = xstrdup(filetxt_job->nodes);
	jobacct_job->priority = filetxt_job->priority;
	jobacct_job->requid = filetxt_job->requid;
	memcpy(&jobacct_job->sacct, &filetxt_job->sacct, sizeof(sacct_t));
	jobacct_job->start = filetxt_job->header.timestamp -
		jobacct_job->elapsed;
	jobacct_job->state = filetxt_job->status;

	jobacct_job->steps = list_create(destroy_jobacct_step_rec);
	if(filetxt_job->steps) {
		itr = list_iterator_create(filetxt_job->steps);
		while((filetxt_step = list_next(itr))) {
			list_append(jobacct_job->steps,
				    _create_jobacct_step_rec(filetxt_step));
		}
		list_iterator_destroy(itr);
	}
	jobacct_job->submit = filetxt_job->header.job_submit;
	
	jobacct_job->sys_cpu_sec = filetxt_job->rusage.ru_stime.tv_sec;
	jobacct_job->sys_cpu_usec = filetxt_job->rusage.ru_stime.tv_usec;
	jobacct_job->tot_cpu_sec = filetxt_job->tot_cpu_sec;
	jobacct_job->tot_cpu_usec = filetxt_job->tot_cpu_usec;
	jobacct_job->track_steps = filetxt_job->track_steps;
	jobacct_job->uid = filetxt_job->header.uid;
	jobacct_job->user = NULL;
	jobacct_job->user_cpu_sec = filetxt_job->rusage.ru_utime.tv_sec;
	jobacct_job->user_cpu_usec = filetxt_job->rusage.ru_utime.tv_usec;
	
	return jobacct_job;
}

static filetxt_job_rec_t *_create_filetxt_job_rec(filetxt_header_t header)
{
	filetxt_job_rec_t *job = xmalloc(sizeof(filetxt_job_rec_t));
	memcpy(&job->header, &header, sizeof(filetxt_header_t));
	memset(&job->rusage, 0, sizeof(struct rusage));
	memset(&job->sacct, 0, sizeof(sacct_t));
	job->sacct.min_cpu = (float)NO_VAL;
	job->job_start_seen = 0;
	job->job_step_seen = 0;
	job->job_terminated_seen = 0;
	job->jobnum_superseded = 0;
	job->jobname = NULL;
	job->status = JOB_PENDING;
	job->nodes = NULL;
	job->jobname = NULL;
	job->exitcode = 0;
	job->priority = 0;
	job->ntasks = 0;
	job->ncpus = 0;
	job->elapsed = 0;
	job->tot_cpu_sec = 0;
	job->tot_cpu_usec = 0;
	job->steps = list_create(_destroy_filetxt_step_rec);
	job->nodes = NULL;
	job->track_steps = 0;
	job->account = NULL;
	job->requid = -1;

      	return job;
}

static filetxt_step_rec_t *_create_filetxt_step_rec(filetxt_header_t header)
{
	filetxt_step_rec_t *step = xmalloc(sizeof(filetxt_job_rec_t));
	memcpy(&step->header, &header, sizeof(filetxt_header_t));
	memset(&step->rusage, 0, sizeof(struct rusage));
	memset(&step->sacct, 0, sizeof(sacct_t));
	step->stepnum = (uint32_t)NO_VAL;
	step->nodes = NULL;
	step->stepname = NULL;
	step->status = NO_VAL;
	step->exitcode = NO_VAL;
	step->ntasks = (uint32_t)NO_VAL;
	step->ncpus = (uint32_t)NO_VAL;
	step->elapsed = (uint32_t)NO_VAL;
	step->tot_cpu_sec = (uint32_t)NO_VAL;
	step->tot_cpu_usec = (uint32_t)NO_VAL;
	step->account = NULL;
	step->requid = -1;

	return step;
}

/* prefix_filename() -- insert a filename prefix into a path
 *
 * IN:	path = fully-qualified path+file name
 *      prefix = the prefix to insert into the file name
 * RETURNS: pointer to the updated path+file name
 */

static char *_prefix_filename(char *path, char *prefix) {
	char	*out;
	int     i,
		plen;

	plen = strlen(path);
	out = xmalloc(plen+strlen(prefix)+1);
	for (i=plen-1; i>=0; i--)
		if (path[i]=='/') {
			break;
		}
	i++;
	*out = 0;
	strncpy(out, path, i);
	out[i] = 0;
	strcat(out, prefix);
	strcat(out, path+i);
	return(out);
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

static char *_convert_type(int rec_type)
{
	switch(rec_type) {
	case JOB_START:
		return "JOB_START";
	case JOB_STEP:
		return "JOB_STEP";
	case JOB_TERMINATED:
		return "JOB_TERMINATED";
	default:
		return "UNKNOWN";
	}
}

static int _cmp_jrec(const void *a1, const void *a2) {
	expired_rec_t *j1 = (expired_rec_t *) a1;
	expired_rec_t *j2 = (expired_rec_t *) a2;

	if (j1->job <  j2->job)
		return -1;
	else if (j1->job == j2->job) {
		if(j1->job_submit == j2->job_submit)
			return 0;
		else 
			return 1;
	}
	return 1;
}

static void _show_rec(char *f[])
{
	int 	i;
	fprintf(stderr, "rec>");
	for (i=0; f[i]; i++)
		fprintf(stderr, " %s", f[i]);
	fprintf(stderr, "\n");
	return;
}

static void _do_fdump(char* f[], int lc)
{
	int	i=0, j=0;
	char **type;
	char    *header[] = {"job",       /* F_JOB */
			     "partition", /* F_PARTITION */
			     "job_submit", /* F_JOB_SUBMIT */
			     "timestamp", /* F_TIMESTAMP */
			     "uid",	 /* F_UIDGID */
			     "gid",	 /* F_UIDGID */
			     "BlockID",  /* F_BLOCKID */
			     "reserved-2",/* F_RESERVED1 */
			     "recordType",/* F_RECTYPE */
			     NULL};

	char	*start[] = {"jobName",	 /* F_JOBNAME */ 
			    "TrackSteps", /* F_TRACK_STEPS */
			    "priority",	 /* F_PRIORITY */
			    "ncpus",	 /* F_NCPUS */
			    "nodeList", /* F_NODES */
			    "account",   /* F_JOB_ACCOUNT */
			    NULL};
		
	char	*step[] = {"jobStep",	 /* F_JOBSTEP */
			   "status",	 /* F_STATUS */ 
			   "exitcode",	 /* F_EXITCODE */
			   "ntasks",	 /* F_NTASKS */
			   "ncpus",	 /* F_STEPNCPUS */
			   "elapsed",	 /* F_ELAPSED */
			   "cpu_sec",	 /* F_CPU_SEC */
			   "cpu_usec",	 /* F_CPU_USEC */
			   "user_sec",	 /* F_USER_SEC */
			   "user_usec",	 /* F_USER_USEC */
			   "sys_sec",	 /* F_SYS_SEC */
			   "sys_usec",	 /* F_SYS_USEC */
			   "rss",	 /* F_RSS */
			   "ixrss",	 /* F_IXRSS */
			   "idrss",	 /* F_IDRSS */
			   "isrss",	 /* F_ISRSS */
			   "minflt",	 /* F_MINFLT */
			   "majflt",	 /* F_MAJFLT */
			   "nswap",	 /* F_NSWAP */
			   "inblocks",	 /* F_INBLOCKS */
			   "oublocks",	 /* F_OUTBLOCKS */
			   "msgsnd",	 /* F_MSGSND */
			   "msgrcv",	 /* F_MSGRCV */
			   "nsignals",	 /* F_NSIGNALS */
			   "nvcsw",	 /* F_VCSW */
			   "nivcsw",	 /* F_NIVCSW */
			   "max_vsize",	 /* F_MAX_VSIZE */
			   "max_vsize_task",	 /* F_MAX_VSIZE_TASK */
			   "ave_vsize",	 /* F_AVE_VSIZE */
			   "max_rss",	 /* F_MAX_RSS */
			   "max_rss_task",	 /* F_MAX_RSS_TASK */
			   "ave_rss",	 /* F_AVE_RSS */
			   "max_pages",	 /* F_MAX_PAGES */
			   "max_pages_task",	 /* F_MAX_PAGES_TASK */
			   "ave_pages",	 /* F_AVE_PAGES */
			   "min_cputime",	 /* F_MIN_CPU */
			   "min_cputime_task",	 /* F_MIN_CPU_TASK */
			   "ave_cputime",	 /* F_AVE_RSS */
			   "StepName",	 /* F_STEPNAME */
			   "StepNodes",	 /* F_STEPNODES */
			   "max_vsize_node",	 /* F_MAX_VSIZE_NODE */
			   "max_rss_node",	 /* F_MAX_RSS_NODE */
			   "max_pages_node",	 /* F_MAX_PAGES_NODE */
			   "min_cputime_node",	 /* F_MIN_CPU_NODE */
			   "account",    /* F_STEP_ACCOUNT */
			   "requid",     /* F_STEP_REQUID */
			   NULL};
       
	char	*suspend[] = {"Suspend/Run time", /* F_TOT_ELAPSED */
			      "status",	 /* F_STATUS */ 
			      NULL};	 

	char	*term[] = {"totElapsed", /* F_TOT_ELAPSED */
			   "status",	 /* F_STATUS */ 
			   "requid",     /* F_JOB_REQUID */
			   "exitcode",	 /* F_EXITCODE */
			   NULL};	 
		
	i = atoi(f[F_RECTYPE]);
	printf("\n------- Line %d %s -------\n", lc, _convert_type(i));

	for(j=0; j < HEADER_LENGTH; j++) 
		printf("%12s: %s\n", header[j], f[j]);
	switch(i) {
	case JOB_START:
		type = start;
		j = JOB_START_LENGTH;
		break;
	case JOB_STEP:
		type = step;
		j = JOB_STEP_LENGTH;
		break;
	case JOB_SUSPEND:
		type = suspend;
		j = JOB_TERM_LENGTH;
	case JOB_TERMINATED:
		type = term;
		j = JOB_TERM_LENGTH;
		break;
	default:
		while(f[j]) {
			printf("      Field[%02d]: %s\n", j, f[j]); 
			j++;
		}
		return;
	}
	
	for(i=HEADER_LENGTH; i < j; i++)
       		printf("%12s: %s\n", type[i-HEADER_LENGTH], f[i]);	
}

static filetxt_job_rec_t *_find_job_record(List job_list, 
					   filetxt_header_t header,
					   int type)
{
	filetxt_job_rec_t *job = NULL;
	ListIterator itr = list_iterator_create(job_list);

	while((job = (filetxt_job_rec_t *)list_next(itr)) != NULL) {
		if (job->header.jobnum == header.jobnum) {
			if(job->header.job_submit == 0 && type == JOB_START) {
				list_remove(itr);
				_destroy_filetxt_job_rec(job);
				job = NULL;
				break;
			}
		
			if(job->header.job_submit == BATCH_JOB_TIMESTAMP) {
				job->header.job_submit = header.job_submit;
				break;
			}
			
			if(job->header.job_submit == header.job_submit)
				break;
			else {
				/* If we're looking for a later
				 * record with this job number, we
				 * know that this one is an older,
				 * duplicate record.
				 *   We assume that the newer record
				 * will be created if it doesn't
				 * already exist. */
				job->jobnum_superseded = 1;
			}
		}
	}
	list_iterator_destroy(itr);
	return job;
}

static int _remove_job_record(List job_list, uint32_t jobnum)
{
	filetxt_job_rec_t *job = NULL;
	int rc = SLURM_ERROR;
	ListIterator itr = list_iterator_create(job_list);

	while((job = (filetxt_job_rec_t *)list_next(itr)) != NULL) {
		if (job->header.jobnum == jobnum) {
			list_remove(itr);
			_destroy_filetxt_job_rec(job);
			rc = SLURM_SUCCESS;
		}
	}
	list_iterator_destroy(itr);
	return rc;
}

static filetxt_step_rec_t *_find_step_record(filetxt_job_rec_t *job,
					     long stepnum)
{
	filetxt_step_rec_t *step = NULL;
	ListIterator itr = NULL;

	if(!list_count(job->steps))
		return step;
	
	itr = list_iterator_create(job->steps);
	while((step = (filetxt_step_rec_t *)list_next(itr)) != NULL) {
		if (step->stepnum == stepnum)
			break;
	}
	list_iterator_destroy(itr);
	return step;
}

static int _parse_header(char *f[], filetxt_header_t *header)
{
	header->jobnum = atoi(f[F_JOB]);
	header->partition = xstrdup(f[F_PARTITION]);
	header->job_submit = atoi(f[F_JOB_SUBMIT]);
	header->timestamp = atoi(f[F_TIMESTAMP]);
	header->uid = atoi(f[F_UID]);
	header->gid = atoi(f[F_GID]);
	header->blockid = xstrdup(f[F_BLOCKID]);

	return SLURM_SUCCESS;
}

static int _parse_line(char *f[], void **data, int len)
{
	int i = atoi(f[F_RECTYPE]);
	filetxt_job_rec_t **job = (filetxt_job_rec_t **)data;
	filetxt_step_rec_t **step = (filetxt_step_rec_t **)data;
	filetxt_header_t header;
	_parse_header(f, &header);
		
	switch(i) {
	case JOB_START:
		*job = _create_filetxt_job_rec(header);
		(*job)->jobname = xstrdup(f[F_JOBNAME]);
		(*job)->track_steps = atoi(f[F_TRACK_STEPS]);
		(*job)->priority = atoi(f[F_PRIORITY]);
		(*job)->ncpus = atoi(f[F_NCPUS]);
		(*job)->nodes = xstrdup(f[F_NODES]);
		for (i=0; (*job)->nodes[i]; i++) { /* discard trailing <CR> */
			if (isspace((*job)->nodes[i]))
				(*job)->nodes[i] = '\0';
		}
		if (!strcmp((*job)->nodes, "(null)")) {
			xfree((*job)->nodes);
			(*job)->nodes = xstrdup("(unknown)");
		}
		if (len > F_JOB_ACCOUNT) {
			(*job)->account = xstrdup(f[F_JOB_ACCOUNT]);
			for (i=0; (*job)->account[i]; i++) {
				/* discard trailing <CR> */
				if (isspace((*job)->account[i]))
					(*job)->account[i] = '\0';
			}
		}
		break;
	case JOB_STEP:
		*step = _create_filetxt_step_rec(header);
		(*step)->stepnum = atoi(f[F_JOBSTEP]);
		(*step)->status = atoi(f[F_STATUS]);
		(*step)->exitcode = atoi(f[F_EXITCODE]);
		(*step)->ntasks = atoi(f[F_NTASKS]);
		(*step)->ncpus = atoi(f[F_STEPNCPUS]);
		(*step)->elapsed = atoi(f[F_ELAPSED]);
		(*step)->tot_cpu_sec = atoi(f[F_CPU_SEC]);
		(*step)->tot_cpu_usec = atoi(f[F_CPU_USEC]);
		(*step)->rusage.ru_utime.tv_sec = atoi(f[F_USER_SEC]);
		(*step)->rusage.ru_utime.tv_usec = atoi(f[F_USER_USEC]);
		(*step)->rusage.ru_stime.tv_sec = atoi(f[F_SYS_SEC]);
		(*step)->rusage.ru_stime.tv_usec = atoi(f[F_SYS_USEC]);
		(*step)->rusage.ru_maxrss = atoi(f[F_RSS]);
		(*step)->rusage.ru_ixrss = atoi(f[F_IXRSS]);
		(*step)->rusage.ru_idrss = atoi(f[F_IDRSS]);
		(*step)->rusage.ru_isrss = atoi(f[F_ISRSS]);
		(*step)->rusage.ru_minflt = atoi(f[F_MINFLT]);
		(*step)->rusage.ru_majflt = atoi(f[F_MAJFLT]);
		(*step)->rusage.ru_nswap = atoi(f[F_NSWAP]);
		(*step)->rusage.ru_inblock = atoi(f[F_INBLOCKS]);
		(*step)->rusage.ru_oublock = atoi(f[F_OUBLOCKS]);
		(*step)->rusage.ru_msgsnd = atoi(f[F_MSGSND]);
		(*step)->rusage.ru_msgrcv = atoi(f[F_MSGRCV]);
		(*step)->rusage.ru_nsignals = atoi(f[F_NSIGNALS]);
		(*step)->rusage.ru_nvcsw = atoi(f[F_NVCSW]);
		(*step)->rusage.ru_nivcsw = atoi(f[F_NIVCSW]);
		(*step)->sacct.max_vsize = atoi(f[F_MAX_VSIZE]) * 1024;
		if(len > F_STEPNODES) {
			(*step)->sacct.max_vsize_id.taskid = 
				atoi(f[F_MAX_VSIZE_TASK]);
			(*step)->sacct.ave_vsize = atof(f[F_AVE_VSIZE]) * 1024;
			(*step)->sacct.max_rss = atoi(f[F_MAX_RSS]) * 1024;
			(*step)->sacct.max_rss_id.taskid = 
				atoi(f[F_MAX_RSS_TASK]);
			(*step)->sacct.ave_rss = atof(f[F_AVE_RSS]) * 1024;
			(*step)->sacct.max_pages = atoi(f[F_MAX_PAGES]);
			(*step)->sacct.max_pages_id.taskid = 
				atoi(f[F_MAX_PAGES_TASK]);
			(*step)->sacct.ave_pages = atof(f[F_AVE_PAGES]);
			(*step)->sacct.min_cpu = atof(f[F_MIN_CPU]);
			(*step)->sacct.min_cpu_id.taskid = 
				atoi(f[F_MIN_CPU_TASK]);
			(*step)->sacct.ave_cpu = atof(f[F_AVE_CPU]);
			(*step)->stepname = xstrdup(f[F_STEPNAME]);
			(*step)->nodes = xstrdup(f[F_STEPNODES]);
		} else {
			(*step)->sacct.max_vsize_id.taskid = (uint16_t)NO_VAL;
			(*step)->sacct.ave_vsize = (float)NO_VAL;
			(*step)->sacct.max_rss = (uint32_t)NO_VAL;
			(*step)->sacct.max_rss_id.taskid = (uint16_t)NO_VAL;
			(*step)->sacct.ave_rss = (float)NO_VAL;
			(*step)->sacct.max_pages = (uint32_t)NO_VAL;
			(*step)->sacct.max_pages_id.taskid = (uint16_t)NO_VAL;
			(*step)->sacct.ave_pages = (float)NO_VAL;
			(*step)->sacct.min_cpu = (uint32_t)NO_VAL;
			(*step)->sacct.min_cpu_id.taskid = (uint16_t)NO_VAL;
			(*step)->sacct.ave_cpu =  (float)NO_VAL;
			(*step)->stepname = NULL;
			(*step)->nodes = NULL;
		}
		if(len > F_MIN_CPU_NODE) {
			(*step)->sacct.max_vsize_id.nodeid = 
				atoi(f[F_MAX_VSIZE_NODE]);
			(*step)->sacct.max_rss_id.nodeid = 
				atoi(f[F_MAX_RSS_NODE]);
			(*step)->sacct.max_pages_id.nodeid = 
				atoi(f[F_MAX_PAGES_NODE]);
			(*step)->sacct.min_cpu_id.nodeid = 
				atoi(f[F_MIN_CPU_NODE]);
		} else {
			(*step)->sacct.max_vsize_id.nodeid = 
				(uint32_t)NO_VAL;
			(*step)->sacct.max_rss_id.nodeid = 
				(uint32_t)NO_VAL;
			(*step)->sacct.max_pages_id.nodeid = 
				(uint32_t)NO_VAL;
			(*step)->sacct.min_cpu_id.nodeid = 
				(uint32_t)NO_VAL;
		}
		if(len > F_STEP_ACCOUNT)
			(*step)->account = xstrdup(f[F_STEP_ACCOUNT]);
		if(len > F_STEP_REQUID)
			(*step)->requid = atoi(f[F_STEP_REQUID]);
		break;
	case JOB_SUSPEND:
	case JOB_TERMINATED:
		*job = _create_filetxt_job_rec(header);
		(*job)->elapsed = atoi(f[F_TOT_ELAPSED]);
		(*job)->status = atoi(f[F_STATUS]);		
		if(len > F_JOB_REQUID) 
			(*job)->requid = atoi(f[F_JOB_REQUID]);
		if(len > F_JOB_EXITCODE) 
			(*job)->exitcode = atoi(f[F_JOB_EXITCODE]);
		break;
	default:
		printf("UNKOWN TYPE %d",i);
		break;
	}
	return SLURM_SUCCESS;
}

static void _process_start(List job_list, char *f[], int lc,
			   int show_full, int len)
{
	filetxt_job_rec_t *job = NULL;
	filetxt_job_rec_t *temp = NULL;

	_parse_line(f, (void **)&temp, len);
	job = _find_job_record(job_list, temp->header, JOB_START);
	if (job) {	/* Hmmm... that's odd */
		printf("job->header.job_submit = %d",
		       (int)job->header.job_submit);
		if(job->header.job_submit == 0)
			_remove_job_record(job_list, job->header.jobnum);
		else {
			fprintf(stderr,
				"Conflicting JOB_START for job %u at"
				" line %d -- ignoring it\n",
				job->header.jobnum, lc);
			_destroy_filetxt_job_rec(temp);
			return;
		}
	}
	
	job = temp;
	job->show_full = show_full;
	list_append(job_list, job);
	job->job_start_seen = 1;
	
}

static void _process_step(List job_list, char *f[], int lc,
			  int show_full, int len,
			  sacct_parameters_t *params)
{
	filetxt_job_rec_t *job = NULL;
	
	filetxt_step_rec_t *step = NULL;
	filetxt_step_rec_t *temp = NULL;

	_parse_line(f, (void **)&temp, len);

	job = _find_job_record(job_list, temp->header, JOB_STEP);
	
	if (temp->stepnum == -2) {
		_destroy_filetxt_step_rec(temp);
		return;
	}
	if (!job) {	/* fake it for now */
		job = _create_filetxt_job_rec(temp->header);
		job->jobname = xstrdup("(unknown)");
		if (params->opt_verbose > 1) 
			fprintf(stderr, 
				"Note: JOB_STEP record %u.%u preceded "
				"JOB_START record at line %d\n",
				temp->header.jobnum, temp->stepnum, lc);
	}
	job->show_full = show_full;
	
	if ((step = _find_step_record(job, temp->stepnum))) {
		
		if (temp->status == JOB_RUNNING) {
			_destroy_filetxt_step_rec(temp);
			return;/* if "R" record preceded by F or CD; unusual */
		}
		if (step->status != JOB_RUNNING) { /* if not JOB_RUNNING */
			fprintf(stderr,
				"Conflicting JOB_STEP record for "
				"jobstep %u.%u at line %d "
				"-- ignoring it\n",
				step->header.jobnum, 
				step->stepnum, lc);
			_destroy_filetxt_step_rec(temp);
			return;
		}
		step->status = temp->status;
		step->exitcode = temp->exitcode;
		step->ntasks = temp->ntasks;
		step->ncpus = temp->ncpus;
		step->elapsed = temp->elapsed;
		step->tot_cpu_sec = temp->tot_cpu_sec;
		step->tot_cpu_usec = temp->tot_cpu_usec;
		job->requid = temp->requid;
		step->requid = temp->requid;
		memcpy(&step->rusage, &temp->rusage, sizeof(struct rusage));
		memcpy(&step->sacct, &temp->sacct, sizeof(sacct_t));
		xfree(step->stepname);
		step->stepname = xstrdup(temp->stepname);
		step->end = temp->header.timestamp;
		_destroy_filetxt_step_rec(temp);
		goto got_step;
	}
	step = temp;
	temp = NULL;
	list_append(job->steps, step);
	if(list_count(job->steps) > 1)
		job->track_steps = 1;
	if(job->header.timestamp == 0)
		job->header.timestamp = step->header.timestamp;
	job->job_step_seen = 1;
	job->ntasks += step->ntasks;
	if(!job->nodes || !strcmp(job->nodes, "(unknown)")) {
		xfree(job->nodes);
		job->nodes = xstrdup(step->nodes);
	}
	
got_step:
	
		
	if (job->job_terminated_seen == 0) {	/* If the job is still running,
						   this is the most recent
						   status */
		if ( job->exitcode == 0 )
			job->exitcode = step->exitcode;
		job->status = JOB_RUNNING;
		job->elapsed = step->header.timestamp - job->header.timestamp;
	}
}

static void _process_suspend(List job_list, char *f[], int lc,
			     int show_full, int len)
{
	filetxt_job_rec_t *job = NULL;
	filetxt_job_rec_t *temp = NULL;

	_parse_line(f, (void **)&temp, len);
	job = _find_job_record(job_list, temp->header, JOB_SUSPEND);
	if (!job)  {	/* fake it for now */
		job = _create_filetxt_job_rec(temp->header);
		job->jobname = xstrdup("(unknown)");
	} 
			
	job->show_full = show_full;
	if (job->status == JOB_SUSPENDED) 
		job->elapsed -= temp->elapsed;

	//job->header.timestamp = temp->header.timestamp;
	job->status = temp->status;
	_destroy_filetxt_job_rec(temp);
}
	
static void _process_terminated(List job_list, char *f[], int lc,
				int show_full, int len,
				sacct_parameters_t *params)
{
	filetxt_job_rec_t *job = NULL;
	filetxt_job_rec_t *temp = NULL;

	_parse_line(f, (void **)&temp, len);
	job = _find_job_record(job_list, temp->header, JOB_TERMINATED);
	if (!job) {	/* fake it for now */
		job = _create_filetxt_job_rec(temp->header);
		job->jobname = xstrdup("(unknown)");
		if (params->opt_verbose > 1) 
			fprintf(stderr, "Note: JOB_TERMINATED record for job "
				"%u preceded "
				"other job records at line %d\n",
				temp->header.jobnum, lc);
	} else if (job->job_terminated_seen) {
		if (temp->status == JOB_NODE_FAIL) {
			/* multiple node failures - extra TERMINATED records */
			if (params->opt_verbose > 1)
				fprintf(stderr, 
					"Note: Duplicate JOB_TERMINATED "
					"record (nf) for job %u at "
					"line %d\n", 
					temp->header.jobnum, lc);
			/* JOB_TERMINATED/NF records may be preceded
			 * by a JOB_TERMINATED/CA record; NF is much
			 * more interesting.
			 */
			job->status = temp->status;
			goto finished;
		}
		
		fprintf(stderr,
			"Conflicting JOB_TERMINATED record (%s) for "
			"job %u at line %d -- ignoring it\n",
			job_state_string(temp->status), 
			job->header.jobnum, lc);
		goto finished;
	}
	job->job_terminated_seen = 1;
	job->elapsed = temp->elapsed;
	job->end = temp->header.timestamp;
	job->status = temp->status;
	job->requid = temp->requid;
	job->exitcode = temp->exitcode;
	if(list_count(job->steps) > 1)
		job->track_steps = 1;
	job->show_full = show_full;
	
finished:
	_destroy_filetxt_job_rec(temp);
}

extern List filetxt_jobacct_process_get_jobs(List selected_steps,
					     List selected_parts,
					     sacct_parameters_t *params)
{
	char line[BUFFER_SIZE];
	char *f[MAX_RECORD_FIELDS+1];    /* End list with null entry and,
					    possibly, more data than we
					    expected */
	char *fptr;
	int i;
	FILE *fd = NULL;
	int lc = 0;
	int rec_type = -1;
	filetxt_job_rec_t *filetxt_job = NULL;
	jobacct_selected_step_t *selected_step = NULL;
	char *selected_part = NULL;
	ListIterator itr = NULL;
	int show_full = 0;
	List ret_job_list = list_create(destroy_jobacct_job_rec);
	List job_list = list_create(_destroy_filetxt_job_rec);
	fd = _open_log_file(params->opt_filein);
	
	while (fgets(line, BUFFER_SIZE, fd)) {
		lc++;
		fptr = line;	/* break the record into NULL-
				   terminated strings */
				
		for (i = 0; i < MAX_RECORD_FIELDS; i++) {
			f[i] = fptr;
			fptr = strstr(fptr, " ");
			if (fptr == NULL) {
				fptr = strstr(f[i], "\n");
				if (fptr)
					*fptr = 0;
				break; 
			} else
				*fptr++ = 0;
		}
		f[++i] = 0;
		
		if(i < HEADER_LENGTH) {
			continue;
		}
		
		rec_type = atoi(f[F_RECTYPE]);
		
		if (list_count(selected_steps)) {
			itr = list_iterator_create(selected_steps);
			while((selected_step = list_next(itr))) {
				if (strcmp(selected_step->job, f[F_JOB]))
					continue;
				/* job matches; does the step? */
				if(selected_step->step == NULL) {
					show_full = 1;
					list_iterator_destroy(itr);
					goto foundjob;
				} else if (rec_type != JOB_STEP 
					   || !strcmp(f[F_JOBSTEP], 
						      selected_step->step)) {
					list_iterator_destroy(itr);
					goto foundjob;
				} 
			}
			list_iterator_destroy(itr);
			continue;	/* no match */
		} else {
			show_full = 1;
		}
	foundjob:
		
		if (list_count(selected_parts)) {
			itr = list_iterator_create(selected_parts);
			while((selected_part = list_next(itr))) 
				if (!strcasecmp(f[F_PARTITION], 
						selected_part)) {
					list_iterator_destroy(itr);
					goto foundp;
				}
			list_iterator_destroy(itr);
			continue;	/* no match */
		}
	foundp:
		
		if (params->opt_fdump) {
			_do_fdump(f, lc);
			continue;
		}
		
		/* Build suitable tables with all the data */
		switch(rec_type) {
		case JOB_START:
			if(i < F_JOB_ACCOUNT) {
				printf("Bad data on a Job Start\n");
				_show_rec(f);
			} else 
				_process_start(job_list, f, lc, show_full, i);
			break;
		case JOB_STEP:
			if(i < F_MAX_VSIZE) {
				printf("Bad data on a Step entry\n");
				_show_rec(f);
			} else
				_process_step(job_list, f, lc, show_full, i, 
					      params);
			break;
		case JOB_SUSPEND:
			if(i < F_JOB_REQUID) {
				printf("Bad data on a Suspend entry\n");
				_show_rec(f);
			} else
				_process_suspend(job_list, f, lc,
						 show_full, i);
			break;
		case JOB_TERMINATED:
			if(i < F_JOB_REQUID) {
				printf("Bad data on a Job Term\n");
				_show_rec(f);
			} else
				_process_terminated(job_list, f, lc,
						    show_full, i, params);
			break;
		default:
			if (params->opt_verbose > 1)
				fprintf(stderr,
					"Invalid record at line %d of "
					"input file\n",
					lc);
			if (params->opt_verbose > 2)
				_show_rec(f);
			break;
		}
	}
	
	if (ferror(fd)) {
		perror(params->opt_filein);
		exit(1);
	} 
	fclose(fd);

	itr = list_iterator_create(job_list);
	while((filetxt_job = list_next(itr))) {
		list_append(ret_job_list, _create_jobacct_job_rec(filetxt_job));
	}
	list_iterator_destroy(itr);
	list_destroy(job_list);

	return ret_job_list;
}

extern void filetxt_jobacct_process_archive(List selected_parts,
					    sacct_parameters_t *params)
{
	char	line[BUFFER_SIZE],
		*f[EXPIRE_READ_LENGTH],
		*fptr = NULL,
		*logfile_name = NULL,
		*old_logfile_name = NULL;
	int	file_err=0,
		new_file,
		i = 0;
	expired_rec_t *exp_rec = NULL;
	expired_rec_t *exp_rec2 = NULL;
	List keep_list = list_create(_destroy_exp);
	List exp_list = list_create(_destroy_exp);
	List other_list = list_create(_destroy_exp);
	struct	stat statbuf;
	mode_t	prot = 0600;
	uid_t	uid;
	gid_t	gid;
	FILE	*expired_logfile = NULL,
		*new_logfile = NULL;
	FILE *fd = NULL;
	int lc=0;
	int rec_type = -1;
	ListIterator itr = NULL;
	ListIterator itr2 = NULL;
	char *temp = NULL;

	/* Figure out our expiration date */
	time_t		expiry;
	expiry = time(NULL)-params->opt_expire;
	if (params->opt_verbose)
		fprintf(stderr, "Purging jobs completed prior to %d\n",
			(int)expiry);

	/* Open the current or specified logfile, or quit */
	fd = _open_log_file(params->opt_filein);
	if (stat(params->opt_filein, &statbuf)) {
		perror("stat'ing logfile");
		goto finished;
	}
	if ((statbuf.st_mode & S_IFLNK) == S_IFLNK) {
		fprintf(stderr, "%s is a symbolic link; --expire requires "
			"a hard-linked file name\n", params->opt_filein);
		goto finished;
	}
	if (!(statbuf.st_mode & S_IFREG)) {
		fprintf(stderr, "%s is not a regular file; --expire "
			"only works on accounting log files\n",
			params->opt_filein);
		goto finished;
	}
	prot = statbuf.st_mode & 0777;
	gid  = statbuf.st_gid;
	uid  = statbuf.st_uid;
	old_logfile_name = _prefix_filename(params->opt_filein, ".old.");
	if (stat(old_logfile_name, &statbuf)) {
		if (errno != ENOENT) {
			fprintf(stderr,"Error checking for %s: ",
				old_logfile_name);
			perror("");
			goto finished;
		}
	} else {
		fprintf(stderr, "Warning! %s exists -- please remove "
			"or rename it before proceeding\n",
			old_logfile_name);
		goto finished;
	}

	/* create our initial buffer */
	while (fgets(line, BUFFER_SIZE, fd)) {
		lc++;
		fptr = line;	/* break the record into NULL-
				   terminated strings */
		exp_rec = xmalloc(sizeof(expired_rec_t));
		exp_rec->line = xstrdup(line);
	
		for (i = 0; i < EXPIRE_READ_LENGTH; i++) {
			f[i] = fptr;
			fptr = strstr(fptr, " ");
			if (fptr == NULL)
				break; 
			else
				*fptr++ = 0;
		}
		
		exp_rec->job = atoi(f[F_JOB]);
		exp_rec->job_submit = atoi(f[F_JOB_SUBMIT]);
		
		rec_type = atoi(f[F_RECTYPE]);
		/* Odd, but complain some other time */
		if (rec_type == JOB_TERMINATED) {
			if (expiry < atoi(f[F_TIMESTAMP])) {
				list_append(keep_list, exp_rec);
				continue;				
			}
			if (list_count(selected_parts)) {
				itr = list_iterator_create(selected_parts);
				while((temp = list_next(itr))) 
					if(!strcasecmp(f[F_PARTITION], temp)) 
						break;
				list_iterator_destroy(itr);
				if(!temp) {
					list_append(keep_list, exp_rec);
					continue;
				} /* no match */
			}
			list_append(exp_list, exp_rec);
			if (params->opt_verbose > 2)
				fprintf(stderr, "Selected: %8d %d\n",
					exp_rec->job,
					(int)exp_rec->job_submit);
		} else {
			list_append(other_list, exp_rec);
		}
	}
	if (!list_count(exp_list)) {
		printf("No job records were purged.\n");
		goto finished;
	}
	logfile_name = xmalloc(strlen(params->opt_filein)+sizeof(".expired"));
	sprintf(logfile_name, "%s.expired", params->opt_filein);
	new_file = stat(logfile_name, &statbuf);
	if ((expired_logfile = fopen(logfile_name, "a"))==NULL) {
		fprintf(stderr, "Error while opening %s", 
			logfile_name);
		perror("");
		xfree(logfile_name);
		goto finished;
	}
	
	if (new_file) {  /* By default, the expired file looks like the log */
		chmod(logfile_name, prot);
		chown(logfile_name, uid, gid);
	}
	xfree(logfile_name);

	logfile_name = _prefix_filename(params->opt_filein, ".new.");
	if ((new_logfile = fopen(logfile_name, "w"))==NULL) {
		fprintf(stderr, "Error while opening %s",
			logfile_name);
		perror("");
		fclose(expired_logfile);
		goto finished;
	}
	chmod(logfile_name, prot);     /* preserve file protection */
	chown(logfile_name, uid, gid); /* and ownership */
	/* Use line buffering to allow us to safely write
	 * to the log file at the same time as slurmctld. */ 
	if (setvbuf(new_logfile, NULL, _IOLBF, 0)) {
		perror("setvbuf()");
		fclose(expired_logfile);
		goto finished2;
	}

	list_sort(exp_list, (ListCmpF) _cmp_jrec);
	list_sort(keep_list, (ListCmpF) _cmp_jrec);
	
	if (params->opt_verbose > 2) {
		fprintf(stderr, "--- contents of exp_list ---");
		itr = list_iterator_create(exp_list);
		while((exp_rec = list_next(itr))) {
			if (!(i%5))
				fprintf(stderr, "\n");
			else
				fprintf(stderr, "\t");
			fprintf(stderr, "%d", exp_rec->job);
		}
		fprintf(stderr, "\n---- end of exp_list ---\n");
		list_iterator_destroy(itr);
	}
	/* write the expired file */
	itr = list_iterator_create(exp_list);
	while((exp_rec = list_next(itr))) {
		itr2 = list_iterator_create(other_list);
		while((exp_rec2 = list_next(itr2))) {
			if((exp_rec2->job != exp_rec->job) 
			   || (exp_rec2->job_submit != exp_rec->job_submit))
				continue;
			if (fputs(exp_rec2->line, expired_logfile)<0) {
				perror("writing expired_logfile");
				list_iterator_destroy(itr2);
				list_iterator_destroy(itr);
				fclose(expired_logfile);
				goto finished2;
			}
			list_remove(itr2);
			_destroy_exp(exp_rec2);
		}
		list_iterator_destroy(itr2);
		if (fputs(exp_rec->line, expired_logfile)<0) {
			perror("writing expired_logfile");
			list_iterator_destroy(itr);
			fclose(expired_logfile);
			goto finished2;
		}		
	}
	list_iterator_destroy(itr);
	fclose(expired_logfile);
	
	/* write the new log */
	itr = list_iterator_create(keep_list);
	while((exp_rec = list_next(itr))) {
		itr2 = list_iterator_create(other_list);
		while((exp_rec2 = list_next(itr2))) {
			if(exp_rec2->job != exp_rec->job)
				continue;
			if (fputs(exp_rec2->line, new_logfile)<0) {
				perror("writing keep_logfile");
				list_iterator_destroy(itr2);
				list_iterator_destroy(itr);
				goto finished2;
			}
			list_remove(itr2);
			_destroy_exp(exp_rec2);
		}
		list_iterator_destroy(itr2);
		if (fputs(exp_rec->line, new_logfile)<0) {
			perror("writing keep_logfile");
			list_iterator_destroy(itr);
			goto finished2;
		}		
	}
	list_iterator_destroy(itr);
	
	if (rename(params->opt_filein, old_logfile_name)) {
		perror("renaming logfile to .old.");
		goto finished2;
	}
	if (rename(logfile_name, params->opt_filein)) {
		perror("renaming new logfile");
		/* undo it? */
		if (!rename(old_logfile_name, params->opt_filein)) 
			fprintf(stderr, "Please correct the problem "
				"and try again");
		else
			fprintf(stderr, "SEVERE ERROR: Current accounting "
				"log may have been renamed %s;\n"
				"please rename it to \"%s\" if necessary, "
			        "and try again\n",
				old_logfile_name, params->opt_filein);
		goto finished2;
	}
	fflush(new_logfile);	/* Flush the buffers before forking */
	fflush(fd);
	
	file_err = slurm_reconfigure ();
	
	if (file_err) {
		file_err = 1;
		fprintf(stderr, "Error: Attempt to reconfigure "
			"SLURM failed.\n");
		if (rename(old_logfile_name, params->opt_filein)) {
			perror("renaming logfile from .old.");
			goto finished2;
		}

	}
	if (fseek(fd, 0, SEEK_CUR)) {	/* clear EOF */
		perror("looking for late-arriving records");
		goto finished2;
	}
	while (fgets(line, BUFFER_SIZE, fd)) {
		if (fputs(line, new_logfile)<0) {
			perror("writing final records");
			goto finished2;
		}
	}
	
	printf("%d jobs expired.\n", list_count(exp_list));
finished2:
	fclose(new_logfile);
	if (!file_err) {
		if (unlink(old_logfile_name) == -1)
			error("Unable to unlink old logfile %s: %m",
			      old_logfile_name);
	}
finished:
	fclose(fd);
	list_destroy(exp_list);
	list_destroy(keep_list);
	list_destroy(other_list);
	xfree(old_logfile_name);
	xfree(logfile_name);
}
