/*****************************************************************************\
 *  print.c - print functions for sacct
 *
 *  $Id: print.c 7541 2006-03-18 01:44:58Z da $
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
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#include "sacct.h"
#define FORMAT_STRING_SIZE 32

char *_decode_status(int status);
void _elapsed_time(long secs, long usecs, char *str);

void _elapsed_time(long secs, long usecs, char *str)
{
	long	days, hours, minutes, seconds;
	
	while (usecs >= 1E6) {
		secs++;
		usecs -= 1E6;
	}

	seconds =  secs % 60;
	minutes = (secs / 60)   % 60;
	hours   = (secs / 3600) % 24;
	days    =  secs / 86400;

	if (days) 
		snprintf(str, FORMAT_STRING_SIZE,
			 "%ld-%2.2ld:%2.2ld:%2.2ld",
		         days, hours, minutes, seconds);
	else if (hours)
		snprintf(str, FORMAT_STRING_SIZE,
			 "%ld:%2.2ld:%2.2ld",
		         hours, minutes, seconds);
	else
		snprintf(str, FORMAT_STRING_SIZE,
			 "%ld:%2.2ld",
		         minutes, seconds);
}

void print_fields(type_t type, void *object)
{
	int f, pf;
	for (f=0; f<nprintfields; f++) {
		pf = printfields[f];
		if (f)
			printf(" ");
		(fields[pf].print_routine)(type, object);
	}
	printf("\n");
}

/* Field-specific print routines */

void print_cpu(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char str[FORMAT_STRING_SIZE];
	
	switch(type) {
	case HEADLINE:
		printf("%15s", "Cpu");
		break;
	case UNDERSCORE:
		printf("%15s", "---------------");
		break;
	case JOB:
		_elapsed_time(job->tot_cpu_sec, job->tot_cpu_usec, str);
		printf("%15s", str);
		break;
	case JOBSTEP:
		_elapsed_time(step->tot_cpu_sec, step->tot_cpu_usec, str);
		printf("%15s", str);
		break;
	} 
}

void print_elapsed(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char str[FORMAT_STRING_SIZE];

	switch(type) {
	case HEADLINE:
		printf("%15s", "Elapsed");
		break;
	case UNDERSCORE:
		printf("%15s", "---------------");
		break;
	case JOB:
		_elapsed_time(job->elapsed, 0, str);
		printf("%15s", str);
		break;
	case JOBSTEP:
		_elapsed_time(step->elapsed, 0, str);
		printf("%15s", str);
		break;
	} 
}

void print_exitcode(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "ExitCode");
		break;
	case UNDERSCORE:
		printf("%8s", "--------");
		break;
	case JOB:
		printf("%8d", job->exitcode);
		break;
	case JOBSTEP:
		printf("%8d", step->exitcode);
		break;
	} 
}

void print_gid(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%5s", "Gid");
		break;
	case UNDERSCORE:
		printf("%5s", "-----");
		break;
	case JOB:
		printf("%5d", job->header.gid);
		break;
	case JOBSTEP:
		printf("s%5d", step->header.gid);
		break;
	} 
}

void print_group(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	int gid = -1;
	char	*tmp="(unknown)";
	struct	group *gr = NULL;
			
	switch(type) {
	case HEADLINE:
		printf("%-9s", "Group");
		break;
	case UNDERSCORE:
		printf("%-9s", "---------");
		break;
	case JOB:
		gid = job->header.gid;
		break;
	case JOBSTEP:
		gid = step->header.gid;
		break;
	}
	if(gid != -1) {
		if ((gr=getgrgid(gid)))
			tmp=gr->gr_name;
		printf("%-9s", tmp);
	} 
}

void print_idrss(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	struct rusage rusage;
	char outbuf[FORMAT_STRING_SIZE];
	
	switch(type) {
	case HEADLINE:
		printf("%8s", "Idrss");
		break;
	case UNDERSCORE:
		printf("%8s", "------");
		break;
	case JOB:
		rusage = job->rusage;
		break;
	case JOBSTEP:
		rusage = step->rusage;
		break;
	} 
	convert_num((float)rusage.ru_idrss, outbuf);
	printf("%8s", outbuf);
}

void print_inblocks(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Inblocks");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_inblock);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_inblock);
		break;
	} 
}

void print_isrss(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "Isrss");
		break;
	case UNDERSCORE:
		printf("%8s", "------");
		break;
	case JOB:
		printf("%8ld", job->rusage.ru_isrss);
		break;
	case JOBSTEP:
		printf("%8ld", step->rusage.ru_isrss);
		break;
	} 

}

void print_ixrss(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "Ixrss");
		break;
	case UNDERSCORE:
		printf("%8s", "------");
		break;
	case JOB:
		printf("%8ld", job->rusage.ru_ixrss);
		break;
	case JOBSTEP:
		printf("%8ld", step->rusage.ru_ixrss);
		break;
	} 

}

void print_job(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "Job");
		break;
	case UNDERSCORE:
		printf("%8s", "--------");
		break;
	case JOB:
		printf("%8d", job->header.jobnum);
		break;
	case JOBSTEP:
		printf("%8d", step->header.jobnum);
		break;
	} 
}

void print_name(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%-18s", "Jobname");
		break;
	case UNDERSCORE:
		printf("%-18s", "------------------");
		break;
	case JOB:
		if(!job->jobname)
			printf("%-18s", "unknown");			     
		else if(strlen(job->jobname)<19)
			printf("%-18s", job->jobname);
		else
			printf("%-15.15s...", job->jobname);
			
		break;
	case JOBSTEP:
		if(!step->stepname)
			printf("%-18s", "unknown");			     
		else if(strlen(step->stepname)<19)
			printf("%-18s", step->stepname);
		else
			printf("%-15.15s...", step->stepname);
		break;
	} 
}

void print_jobid(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char outbuf[10];

	switch(type) {
	case HEADLINE:
		printf("%-10s", "JobID");
		break;
	case UNDERSCORE:
		printf("%-10s", "----------");
		break;
	case JOB:
		printf("%-10d", job->header.jobnum);
		break;
	case JOBSTEP:
		snprintf(outbuf, sizeof(outbuf), "%u.%u",
			 step->header.jobnum,
			 step->stepnum);
		printf("%-10s", outbuf);
		break;
	} 

}

void print_majflt(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "Majflt");
		break;
	case UNDERSCORE:
		printf("%8s", "------");
		break;
	case JOB:
		printf("%8ld", job->rusage.ru_majflt);
		break;
	case JOBSTEP:
		printf("%8ld", step->rusage.ru_majflt);
		break;
	} 
}

void print_minflt(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "Minflt");
		break;
	case UNDERSCORE:
		printf("%8s", "------");
		break;
	case JOB:
		printf("%8ld", job->rusage.ru_minflt);
		break;
	case JOBSTEP:
		printf("%8ld", step->rusage.ru_minflt);
		break;
	} 
}

void print_msgrcv(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Msgrcv");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_msgrcv);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_msgrcv);
		break;
	} 
}

void print_msgsnd(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Msgsnd");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_msgsnd);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_msgsnd);
		break;
	} 
}

void print_ncpus(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%7s", "Ncpus");
		break;
	case UNDERSCORE:
		printf("%7s", "-------");
		break;
	case JOB:
		printf("%7d", job->ncpus);
		break;
	case JOBSTEP:
		printf("%7d", step->ncpus);
		break;
	} 
}

void print_nivcsw(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Nivcsw");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_nivcsw);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_nivcsw);
		break;
	} 
}

void print_nodes(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%-30s", "Nodes");
		break;
	case UNDERSCORE:
		printf("%-30s", "------------------------------");
		break;
	case JOB:
		printf("%-30s", job->nodes);
		break;
	case JOBSTEP:
		printf("%-30s", "                              ");
		break;
	} 
}

void print_nsignals(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Nsignals");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_nsignals);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_nsignals);
		break;
	} 
}

void print_nswap(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%8s", "Nswap");
		break;
	case UNDERSCORE:
		printf("%8s", "------");
		break;
	case JOB:
		printf("%8ld", job->rusage.ru_nswap);
		break;
	case JOBSTEP:
		printf("%8ld", step->rusage.ru_nswap);
		break;
	} 
}

void print_ntasks(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%7s", "Ntasks");
		break;
	case UNDERSCORE:
		printf("%7s", "-------");
		break;
	case JOB:
		printf("%7d", job->ntasks);
		break;
	case JOBSTEP:
		printf("%7d", step->ntasks);
		break;
	} 
}

void print_nvcsw(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Nvcsw");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_nvcsw);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_nvcsw);
		break;
	} 
}

void print_outblocks(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%9s", "Outblocks");
		break;
	case UNDERSCORE:
		printf("%9s", "---------");
		break;
	case JOB:
		printf("%9ld", job->rusage.ru_oublock);
		break;
	case JOBSTEP:
		printf("%9ld", step->rusage.ru_oublock);
		break;
	} 
}

void print_partition(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%-10s", "Partition");
		break;
	case UNDERSCORE:
		printf("%-10s", "----------");
		break;
	case JOB:
		if(!job->header.partition)
			printf("%-10s", "unknown");			     
		else if(strlen(job->header.partition)<11)
			printf("%-10s", job->header.partition);
		else
			printf("%-7.7s...", job->header.partition);
		
		break;
	case JOBSTEP:
		if(!step->header.partition)
			printf("%-10s", "unknown");			     
		else if(strlen(step->header.partition)<11)
			printf("%-10s", step->header.partition);
		else
			printf("%-7.7s...", step->header.partition);
	
		break;
	} 
}

void print_blockid(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%-16s", "BlockID");
		break;
	case UNDERSCORE:
		printf("%-16s", "----------------");
		break;
	case JOB:
		if(!job->header.blockid)
			printf("%-16s", "unknown");			     
		else if(strlen(job->header.blockid)<17)
			printf("%-16s", job->header.blockid);
		else
			printf("%-13.13s...", job->header.blockid);
		
		break;
	case JOBSTEP:
		if(!step->header.blockid)
			printf("%-16s", "unknown");			     
		else if(strlen(step->header.blockid)<17)
			printf("%-16s", step->header.blockid);
		else
			printf("%-13.13s...", step->header.blockid);
	
		break;
	} 
}

void print_pages(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char outbuf[FORMAT_STRING_SIZE];
	char buf1[FORMAT_STRING_SIZE];
	char buf2[FORMAT_STRING_SIZE];
	sacct_t sacct;
	
	switch(type) {
	case HEADLINE:
		printf("%-22s", "MaxPages/Task - Ave");
		break;
	case UNDERSCORE:
		printf("%-22s", "----------------------");
		break;
	case JOB:
		sacct = job->sacct;
		convert_num((float)sacct.max_pages, buf1);
		if(job->track_steps)
			snprintf(outbuf, FORMAT_STRING_SIZE, "%s/- - -", buf1);
		else {
			convert_num((float)sacct.ave_pages, buf2);
			snprintf(outbuf, FORMAT_STRING_SIZE, "%s/%u - %s", 
				 buf1,
				 sacct.max_pages_task, 
				 buf2);
		}
		printf("%-22s", outbuf);
		break;
	case JOBSTEP:
		sacct = step->sacct;
		convert_num((float)sacct.max_pages, buf1);
		convert_num((float)sacct.ave_pages, buf2);
		snprintf(outbuf, FORMAT_STRING_SIZE, "%s/%u - %s", 
			 buf1,
			 sacct.max_pages_task, 
			 buf2);
		printf("%-22s", outbuf);
		break;
	} 
}

void print_rss(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char outbuf[FORMAT_STRING_SIZE];
	char buf1[FORMAT_STRING_SIZE];
	char buf2[FORMAT_STRING_SIZE];
	sacct_t sacct;
	
	switch(type) {
	case HEADLINE:
		printf("%-22s", "MaxRSS/Task - Ave");
		break;
	case UNDERSCORE:
		printf("%-22s", "----------------------");
		break;
	case JOB:
		sacct = job->sacct;
		convert_num((float)sacct.max_rss, buf1);
		if(job->track_steps)
			snprintf(outbuf, FORMAT_STRING_SIZE, "%s/- - -", buf1);
		else {
			convert_num((float)sacct.ave_rss, buf2);
			snprintf(outbuf, FORMAT_STRING_SIZE, "%s/%u - %s", 
				 buf1,
				 sacct.max_rss_task, 
				 buf2);
		}
		printf("%-22s", outbuf);
		break;
	case JOBSTEP:
		sacct = step->sacct;
		convert_num((float)sacct.max_rss, buf1);
		convert_num((float)sacct.ave_rss, buf2);
		snprintf(outbuf, FORMAT_STRING_SIZE, "%s/%u - %s", 
			 buf1,
			 sacct.max_rss_task, 
			 buf2);
		printf("%-22s", outbuf);
		break;
	} 
}

void print_status(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%-10s", "Status");
		break;
	case UNDERSCORE:
		printf("%-10s", "----------");
		break;
	case JOB:
		printf("%-10s", decode_status_int(job->status));
		break;
	case JOBSTEP:
		printf("%-10s", decode_status_int(step->status));
		break;
	} 
}

void print_submitted(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;

	switch(type) {
	case HEADLINE:
		printf("%-14s", "Submitted");
		break;
	case UNDERSCORE:
		printf("%-14s", "--------------");
		break;
	case JOB:
		printf("%-14d", (int)job->header.job_start);
		break;
	case JOBSTEP:
		printf("%-14d", (int)step->header.job_start);
		break;
	} 
}

void print_systemcpu(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char str[FORMAT_STRING_SIZE];

	switch(type) {
	case HEADLINE:
		printf("%15s", "SystemCpu");
		break;
	case UNDERSCORE:
		printf("%15s", "---------------");
		break;
	case JOB:
		_elapsed_time(job->rusage.ru_stime.tv_sec,
			      job->rusage.ru_stime.tv_usec, str);
		printf("%15s", str);
		break;
	case JOBSTEP:
		_elapsed_time(step->rusage.ru_stime.tv_sec,
			      step->rusage.ru_stime.tv_usec, str);
		printf("%15s", str);
		break;
	} 

}

void print_uid(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	
	switch(type) {
	case HEADLINE:
		printf("%5s", "Uid");
		break;
	case UNDERSCORE:
		printf("%5s", "-----");
		break;
	case JOB:
		printf("%5d", job->header.uid);
		break;
	case JOBSTEP:
		printf("%5d", step->header.uid);
		break;
	} 
}

void print_user(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	int uid = -1;
	char	*tmp="(unknown)";
	struct	passwd *pw = NULL;
	
	switch(type) {
	case HEADLINE:
		printf("%-9s", "User");
		break;
	case UNDERSCORE:
		printf("%-9s", "---------");
		break;
	case JOB:
		uid = job->header.uid;
		break;
	case JOBSTEP:
		uid = step->header.uid;
		break;
	} 
	if(uid != -1) {
		if ((pw=getpwuid(uid)))
			tmp=pw->pw_name;
		printf("%-9s", tmp);
	}
}

void print_usercpu(type_t type, void *object)
{
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char str[FORMAT_STRING_SIZE];
	
	switch(type) {
	case HEADLINE:
		printf("%15s", "UserCpu");
		break;
	case UNDERSCORE:
		printf("%15s", "---------------");
		break;
	case JOB:
		_elapsed_time(job->rusage.ru_utime.tv_sec,
			      job->rusage.ru_utime.tv_usec, str);
		printf("%15s", str);
		break;
	case JOBSTEP:
		_elapsed_time(step->rusage.ru_utime.tv_sec,
			      step->rusage.ru_utime.tv_usec, str);
		printf("%15s", str);
		break;
	} 

}

void print_vsize(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char outbuf[FORMAT_STRING_SIZE];
	char buf1[FORMAT_STRING_SIZE];
	char buf2[FORMAT_STRING_SIZE];
	sacct_t sacct;
	
	switch(type) {
	case HEADLINE:
		printf("%-22s", "MaxVSIZE/Task - Ave");
		break;
	case UNDERSCORE:
		printf("%-22s", "----------------------");
		break;
	case JOB:
		sacct = job->sacct;
		convert_num((float)sacct.max_vsize, buf1);
		if(job->track_steps)
			snprintf(outbuf, FORMAT_STRING_SIZE, "%s/- - -", buf1);
		else {
			convert_num((float)sacct.ave_vsize, buf2);
			snprintf(outbuf, FORMAT_STRING_SIZE, "%s/%u - %s", 
				 buf1,
				 sacct.max_vsize_task, 
				 buf2);
		}
		printf("%-22s", outbuf);
		break;
	case JOBSTEP:
		sacct = step->sacct;
		convert_num((float)sacct.max_vsize, buf1);
		convert_num((float)sacct.ave_vsize, buf2);
		snprintf(outbuf, FORMAT_STRING_SIZE, "%s/%u - %s", 
			 buf1,
			 sacct.max_vsize_task, 
			 buf2);
		printf("%-22s", outbuf);
		break;
	} 
}

void print_cputime(type_t type, void *object)
{ 
	job_rec_t *job = (job_rec_t *)object;
	step_rec_t *step = (step_rec_t *)object;
	char outbuf[FORMAT_STRING_SIZE];
	char buf1[FORMAT_STRING_SIZE];
	char buf2[FORMAT_STRING_SIZE];
	sacct_t sacct;
	
	switch(type) {
	case HEADLINE:
		printf("%-22s", "MinCPUtime/Task - Ave");
		break;
	case UNDERSCORE:
		printf("%-22s", "----------------------");
		break;
	case JOB:
		sacct = job->sacct;
		_elapsed_time((int)sacct.min_cpu, 0, buf1);
		if(job->track_steps)
			snprintf(outbuf, FORMAT_STRING_SIZE, 
				 "%s/- - -", buf1);
		else {
			_elapsed_time((int)sacct.ave_cpu, 0, buf2);
			snprintf(outbuf, FORMAT_STRING_SIZE, 
				 "%s/%u - %s", 
				 buf1,
				 sacct.min_cpu_task, 
				 buf2);
		}
		printf("%-22s", outbuf);
		break;
	case JOBSTEP:
		sacct = step->sacct;
		_elapsed_time((int)sacct.min_cpu, 0, buf1);
		_elapsed_time((int)sacct.ave_cpu, 0, buf2);
		snprintf(outbuf, FORMAT_STRING_SIZE, 
			 "%s/%u - %s", 
			 buf1,
			 sacct.min_cpu_task, 
			 buf2);
		printf("%-22s", outbuf);
		break;
	} 
}


