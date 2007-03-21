/*****************************************************************************\
 *  mysql_jobacct.c - functions the mysql jobacct database.
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

#include "mysql_common.h"
#include "mysql_jobacct.h"

static MYSQL *jobacct_mysql_db = NULL;
static int jobacct_db_init = 0;

static char *job_index = "index_table";
/* jobid partition submittime uid gid blockid */

static char *job_table = "job";
/* index start end name state priority cpus nodes account kill_requid */

static char *step_table = "step";
/* index stepid start end name node_list state kill_requid comp_code num_cpus
   max_vsize max_vsize_task max_vsize_node ave_vsize
   max_rss max_rss_task max_rss_node ave_rss
   max_pages max_pages_task max_pages_node ave_pages
   min_cpu min_cpu_task min_cpu_node ave_cpu */

static char *rusage_table = "step_rusage";
/* index stepid
   cputime_sec cputime_usecs user_sec user_usecs sys_sec sys_usecs
   max_rss max_ixrss max_idrss max_isrss max_minflt max_majflt
   max_nswap inblock outblock msgsnd msgrcv nsignals nvcsw nivcsw */

/* Format of the JOB_STEP record */
static FILE *		LOGFILE;

const char *_jobstep_format = 
"%d "
"%u "	/* stepid */
"%d "	/* completion status */
"%d "	/* completion code */
"%u "	/* nprocs */
"%u "	/* number of cpus */
"%u "	/* elapsed seconds */
"%u "	/* total cputime seconds */
"%u "	/* total cputime microseconds */
"%u "	/* user seconds */
"%u "	/* user microseconds */
"%u "	/* system seconds */
"%u "	/* system microseconds */
"%u "	/* max rss */
"%u "	/* max ixrss */
"%u "	/* max idrss */
"%u "	/* max isrss */
"%u "	/* max minflt */
"%u "	/* max majflt */
"%u "	/* max nswap */
"%u "	/* total inblock */
"%u "	/* total outblock */
"%u "	/* total msgsnd */
"%u "	/* total msgrcv */
"%u "	/* total nsignals */
"%u "	/* total nvcsw */
"%u "	/* total nivcsw */
"%u "	/* max vsize */
"%u "	/* max vsize task */
"%.2f "	/* ave vsize */
"%u "	/* max rss */
"%u "	/* max rss task */
"%.2f "	/* ave rss */
"%u "	/* max pages */
"%u "	/* max pages task */
"%.2f "	/* ave pages */
"%.2f "	/* min cpu */
"%u "	/* min cpu task */
"%.2f "	/* ave cpu */
"%s "	/* step process name */
"%s "	/* step node names */
"%u "	/* max vsize node */
"%u "	/* max rss node */
"%u "	/* max pages node */
"%u "	/* min cpu node */
"%s "   /* account */
"%d";   /* requester user id */
/*
 * Print the record to the log file.
 */

static int _print_record(struct job_record *job_ptr, 
			 time_t time, char *data)
{ 
	static int   rc=SLURM_SUCCESS;
	char *block_id = NULL;
	if(!job_ptr->details) {
		error("job_acct: job=%u doesn't exist", job_ptr->job_id);
		return SLURM_ERROR;
	}
	debug2("_print_record, job=%u, \"%s\"",
	       job_ptr->job_id, data);
#ifdef HAVE_BG
	select_g_get_jobinfo(job_ptr->select_jobinfo, 
			     SELECT_DATA_BLOCK_ID, 
			     &block_id);
		
#endif
	if(!block_id)
		block_id = xstrdup("-");

	slurm_mutex_lock( &mysql_lock );

	if (fprintf(LOGFILE,
		    "%u %s %u %u %d %d %s - %s\n",
		    job_ptr->job_id, job_ptr->partition,
		    (int)job_ptr->details->submit_time, (int)time, 
		    job_ptr->user_id, job_ptr->group_id, block_id, data)
	    < 0)
		rc=SLURM_ERROR;
#ifdef HAVE_FDATSYNC
	fdatasync(LOGFILE_FD);
#endif
	slurm_mutex_unlock( &mysql_lock );
	xfree(block_id);

	return rc;
}

static int _mysql_jobacct_check_tables()
{
	char query[1024];
	snprintf(query, sizeof(query),
		 "create table if not exists %s"
		 "(id int not null auto_increment, "
		 "jobid mediumint unsigned not null, "
		 "partition tinytext not null, "
		 "submit int unsigned not null, "
		 "uid smallint unsigned not null, "
		 "gid smallint unsigned not null, blockid tinytext, "
		 "primary key (id))",
		 job_index);
	if(mysql_db_query(jobacct_mysql_db, jobacct_db_init, query)
	   == SLURM_ERROR) 
		return SLURM_ERROR;
	
	snprintf(query, sizeof(query),
		 "create table if not exists %s(id int not null, "
		 "start int unsigned default 0, end int unsigned default 0, "
		 "suspended int unsigned default 0, "
		 "name tinytext not null, track_steps tinyint, "
		 "state smallint not null, priority int unsigned not null, "
		 "cpus mediumint unsigned not null, nodelist text, "
		 "account tinytext, kill_requid smallint)",
		 job_table);
	if(mysql_db_query(jobacct_mysql_db, jobacct_db_init, query)
	   == SLURM_ERROR) 
		return SLURM_ERROR;
	
	snprintf(query, sizeof(query),
		 "create table if not exists %s(id int not null, "
		 "stepid smallint not null, "
		 "start int unsigned default 0, end int unsigned default 0, "
		 "suspended int unsigned default 0, name text not null, "
		 "nodelist text not null, state smallint not null, "
		 "kill_requid smallint, comp_code smallint, "
		 "cpus mediumint unsigned not null, "
		 "max_vsize mediumint unsigned, "
		 "max_vsize_task smallint unsigned, "
		 "max_vsize_node mediumint unsigned, ave_vsize float, "
		 "max_rss mediumint unsigned, max_rss_task smallint unsigned, "
		 "max_rss_node mediumint unsigned, ave_rss float, "
		 "max_pages mediumint unsigned, "
		 "max_pages_task smallint unsigned, "
		 "max_pages_node mediumint unsigned, ave_pages float, "
		 "min_cpu mediumint unsigned, min_cpu_task smallint unsigned, "
		 "min_cpu_node mediumint unsigned, ave_cpu float)",
		 step_table);
	if(mysql_db_query(jobacct_mysql_db, jobacct_db_init, query)
	   == SLURM_ERROR) 
		return SLURM_ERROR;
	
	snprintf(query, sizeof(query),
		 "create table if not exists %s(id int not null, "
		 "stepid smallint not null, "
		 "cpu_sec int unsigned, cpu_usecs int unsigned, "
		 "user_sec int unsigned, user_usecs int unsigned, "
		 "sys_sec int unsigned, sys_usecs int unsigned, "
		 "max_rss int unsigned, max_ixrss int unsigned, "
		 "max_idrss int unsigned, max_isrss int unsigned, "
		 "max_minflt int unsigned, max_majflt int unsigned, "
		 "max_nswap int unsigned, inblock int unsigned, "
		 "outblock int unsigned, msgsnd int unsigned, "
		 "msgrcv int unsigned, nsignals int unsigned, "
		 "nvcsw int unsigned, invcsw int unsigned)",
		 rusage_table);
	if(mysql_db_query(jobacct_mysql_db, jobacct_db_init, query)
	   == SLURM_ERROR) 
		return SLURM_ERROR;
	
	return SLURM_SUCCESS;
}


extern int mysql_jobacct_init()
{
#ifdef HAVE_MYSQL
	mysql_db_info_t *db_info = create_mysql_db_info();
	int 		rc = SLURM_SUCCESS;
	char *db_name = "slurm_jobacct_db";

	debug2("mysql_connect() called");
	
	mysql_get_db_connection(&jobacct_mysql_db, db_name, db_info,
				&jobacct_db_init);

	_mysql_jobacct_check_tables();

	destroy_mysql_db_info(db_info);

	info("Database init finished");

	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int mysql_jobacct_fini()
{
#ifdef HAVE_MYSQL
	if (jobacct_mysql_db)
		mysql_close(jobacct_mysql_db);
	return SLURM_SUCCESS;
#else
	return SLURM_ERROR;
#endif
}

extern int mysql_jobacct_job_start(struct job_record *job_ptr)
{
	int	i,
		ncpus=0,
		rc=SLURM_SUCCESS;
	char	*jname, *account, *nodes;
	long	priority;
	int track_steps = 0;
	char *block_id = NULL;
	char query[1024];
	
	if(!jobacct_db_init) {
		debug("mysql_jobacct_init was not called or it failed");
		return SLURM_ERROR;
	}

	debug2("mysql_jobacct_job_start() called");
	for (i=0; i < job_ptr->num_cpu_groups; i++)
		ncpus += (job_ptr->cpus_per_node[i])
			* (job_ptr->cpu_count_reps[i]);
	priority = (job_ptr->priority == NO_VAL) ?
		-1L : (long) job_ptr->priority;

	if (job_ptr->name && job_ptr->name[0]) {
		jname = job_ptr->name;
	} else {
		jname = "allocation";
		track_steps = 1;
	}

	if (job_ptr->account && job_ptr->account[0])
		account = job_ptr->account;
	else
		account = "(null)";
	if (job_ptr->nodes && job_ptr->nodes[0])
		nodes = job_ptr->nodes;
	else
		nodes = "(null)";

	if(job_ptr->batch_flag)
		track_steps = 1;
#ifdef HAVE_BG
	select_g_get_jobinfo(job_ptr->select_jobinfo, 
			     SELECT_DATA_BLOCK_ID, 
			     &block_id);
		
#endif
	if(!block_id)
		block_id = xstrdup("-");

	job_ptr->requid = -1; /* force to -1 for sacct to know this
			       * hasn't been set yet */
	snprintf(query, sizeof(query),
		 "insert into %s (jobid, partition, submit, uid, gid, "
		 "blockid) values (%u, '%s', %u, %d, %d, '%s')",
		 job_index, job_ptr->job_id, job_ptr->partition,
		 (int)job_ptr->details->submit_time, job_ptr->user_id,
		 job_ptr->group_id, block_id);
	xfree(block_id);

	if((job_ptr->db_index =
	    mysql_insert_ret_id(jobacct_mysql_db, jobacct_db_init, query))) {
		snprintf(query, sizeof(query),
			 "insert into %s (id, start, name, track_steps, "
			 "priority, cpus, nodelist, account) "
			 "values (%d, %u, '%s', %d, %ld, %u, '%s', '%s')",
			 job_table, job_ptr->db_index, 
			 (int)job_ptr->start_time,
			 jname, track_steps, priority, job_ptr->num_procs,
			 nodes, account);
		rc = mysql_db_query(jobacct_mysql_db, jobacct_db_init, query);
	} else 
		rc = SLURM_ERROR;
	
	return rc;
}

extern int mysql_jobacct_job_complete(struct job_record *job_ptr)
{
	char query[1024];
	char	*account, *nodes;
	int rc=SLURM_SUCCESS;
	
	if(!jobacct_db_init) {
		debug("mysql_jobacct_init was not called or it failed");
		return SLURM_ERROR;
	}
	
	debug2("mysql_jobacct_job_complete() called");
	if (job_ptr->end_time == 0) {
		debug("mysql_jobacct: job %u never started", job_ptr->job_id);
		return SLURM_ERROR;
	}	
	
	if (job_ptr->account && job_ptr->account[0])
		account = job_ptr->account;
	else
		account = "(null)";
	if (job_ptr->nodes && job_ptr->nodes[0])
		nodes = job_ptr->nodes;
	else
		nodes = "(null)";

	if(job_ptr->db_index) {
		snprintf(query, sizeof(query),
			 "update %s set end=%u, state=%d, "
			 "kill_requid=%d where id=%u",
			 job_table, (int)job_ptr->end_time, 
			 job_ptr->job_state & (~JOB_COMPLETING),
			 job_ptr->requid, job_ptr->db_index);
		rc = mysql_db_query(jobacct_mysql_db, jobacct_db_init, query);
	} else 
		rc = SLURM_ERROR;

	return  rc;

}

extern int mysql_jobacct_step_start(struct step_record *step_ptr)
{
	int cpus = 0;
	int rc=SLURM_SUCCESS;
	char node_list[BUFFER_SIZE];
#ifdef HAVE_BG
	char *ionodes = NULL;
#endif
	char query[1024];
	
	if(!jobacct_db_init) {
		debug("jobacct init was not called or it failed");
		return SLURM_ERROR;
	}

#ifdef HAVE_BG
	cpus = step_ptr->job_ptr->num_procs;
	select_g_get_jobinfo(step_ptr->job_ptr->select_jobinfo, 
			     SELECT_DATA_IONODES, 
			     &ionodes);
	if(ionodes) {
		snprintf(node_list, BUFFER_SIZE, 
			 "%s[%s]", step_ptr->job_ptr->nodes, ionodes);
		xfree(ionodes);
	} else
		snprintf(node_list, BUFFER_SIZE, "%s",
			 step_ptr->job_ptr->nodes);
	
#else
	if(!step_ptr->step_layout || !step_ptr->step_layout->task_cnt) {
		cpus = step_ptr->job_ptr->num_procs;
		snprintf(node_list, BUFFER_SIZE, "%s", step_ptr->job_ptr->nodes);
	} else {
		cpus = step_ptr->step_layout->task_cnt;
		snprintf(node_list, BUFFER_SIZE, "%s", 
			 step_ptr->step_layout->node_list);
	}
#endif
	step_ptr->job_ptr->requid = -1; /* force to -1 for sacct to know this
					 * hasn't been set yet  */

	if(step_ptr->job_ptr->db_index) {
		snprintf(query, sizeof(query),
			 "insert into %s (id, stepid, start, name, state, "
			 "cpus, nodelist, kill_requid) "
			 "values (%d, %u, %u, '%s', %d, %u, '%s', %d)",
			 step_table, step_ptr->job_ptr->db_index,
			 step_ptr->step_id, 
			 (int)step_ptr->start_time, step_ptr->name,
			 JOB_RUNNING, cpus, node_list, 
			 step_ptr->job_ptr->requid);
		rc = mysql_db_query(jobacct_mysql_db, jobacct_db_init, query);
	} else 
		rc = SLURM_ERROR;
		 
	return rc;
}

extern int mysql_jobacct_step_complete(struct step_record *step_ptr)
{
	char buf[BUFFER_SIZE];
	time_t now;
	int elapsed;
	int comp_status;
	int cpus = 0;
	char node_list[BUFFER_SIZE];
	struct jobacctinfo *jobacct = (struct jobacctinfo *)step_ptr->jobacct;
#ifdef HAVE_BG
	char *ionodes = NULL;
#endif
	float ave_vsize = 0, ave_rss = 0, ave_pages = 0;
	float ave_cpu = 0, ave_cpu2 = 0;
	char *account;
	char query[1024];
	return SLURM_ERROR;
	if(!jobacct_db_init) {
		debug("jobacct init was not called or it failed");
		return SLURM_ERROR;
	}
	
	now = time(NULL);
	
	if ((elapsed=now-step_ptr->start_time)<0)
		elapsed=0;	/* For *very* short jobs, if clock is wrong */
	if (step_ptr->exit_code)
		comp_status = JOB_FAILED;
	else
		comp_status = JOB_COMPLETE;

#ifdef HAVE_BG
	cpus = step_ptr->job_ptr->num_procs;
	select_g_get_jobinfo(step_ptr->job_ptr->select_jobinfo, 
			     SELECT_DATA_IONODES, 
			     &ionodes);
	if(ionodes) {
		snprintf(node_list, BUFFER_SIZE, 
			 "%s[%s]", step_ptr->job_ptr->nodes, ionodes);
		xfree(ionodes);
	} else
		snprintf(node_list, BUFFER_SIZE, "%s", 
			 step_ptr->job_ptr->nodes);
	
#else
	if(!step_ptr->step_layout || !step_ptr->step_layout->task_cnt) {
		cpus = step_ptr->job_ptr->num_procs;
		snprintf(node_list, BUFFER_SIZE, "%s", step_ptr->job_ptr->nodes);
	
	} else {
		cpus = step_ptr->step_layout->task_cnt;
		snprintf(node_list, BUFFER_SIZE, "%s", 
			 step_ptr->step_layout->node_list);
	}
#endif
	/* figure out the ave of the totals sent */
	if(cpus > 0) {
		ave_vsize = jobacct->tot_vsize;
		ave_vsize /= cpus;
		ave_rss = jobacct->tot_rss;
		ave_rss /= cpus;
		ave_pages = jobacct->tot_pages;
		ave_pages /= cpus;
		ave_cpu = jobacct->tot_cpu;
		ave_cpu /= cpus;	
		ave_cpu /= 100;
	}
 
	if(jobacct->min_cpu != (uint32_t)NO_VAL) {
		ave_cpu2 = jobacct->min_cpu;
		ave_cpu2 /= 100;
	}

	if (step_ptr->job_ptr->account && step_ptr->job_ptr->account[0])
		account = step_ptr->job_ptr->account;
	else
		account = "(null)";

/* 	if(step_ptr->job_ptr->db_index) { */
/* 		snprintf(query, sizeof(query), */
/* 			 "update %s set start=%u, end=%u, state=%d, " */
/* 			 "kill_requid=%d, nodelist='%s', account='%s' " */
/* 			 "where id=%u", */
/* 			 step_table, (int)job_ptr->start_time, */
/* 			 (int)job_ptr->end_time,  */
/* 			 job_ptr->job_state & (~JOB_COMPLETING), */
/* 			 job_ptr->requid, nodes, account,step_ptr-> job_ptr->db_index); */
/* 		rc = mysql_db_query(jobacct_mysql_db, jobacct_db_init, query); */
/* 		if(rc != SLURM_ERROR) { */
/* 			snprintf(query, sizeof(query), */
/* 			 "update %s set start=%u, end=%u, state=%d, " */
/* 			 "kill_requid=%d, nodelist='%s', account='%s' " */
/* 			 "where id=%u", */
/* 			 rusage_table, (int)job_ptr->start_time, */
/* 			 (int)job_ptr->end_time,  */
/* 			 job_ptr->job_state & (~JOB_COMPLETING), */
/* 			 job_ptr->requid, nodes, account, job_ptr->db_index); */
/* 			rc = mysql_db_query(jobacct_mysql_db, jobacct_db_init, */
/* 					    query); */
/* 		} */
/* 	} else  */
/* 		rc = SLURM_ERROR; */


	snprintf(buf, BUFFER_SIZE, _jobstep_format,
		 JOB_STEP,
		 step_ptr->step_id,	/* stepid */
		 comp_status,		/* completion status */
		 step_ptr->exit_code,	/* completion code */
		 cpus,          	/* number of tasks */
		 cpus,                  /* number of cpus */
		 elapsed,	        /* elapsed seconds */
		 /* total cputime seconds */
		 jobacct->rusage.ru_utime.tv_sec	
		 + jobacct->rusage.ru_stime.tv_sec,
		 /* total cputime seconds */
		 jobacct->rusage.ru_utime.tv_usec	
		 + jobacct->rusage.ru_stime.tv_usec,
		 jobacct->rusage.ru_utime.tv_sec,	/* user seconds */
		 jobacct->rusage.ru_utime.tv_usec,/* user microseconds */
		 jobacct->rusage.ru_stime.tv_sec,	/* system seconds */
		 jobacct->rusage.ru_stime.tv_usec,/* system microsecs */
		 jobacct->rusage.ru_maxrss,	/* max rss */
		 jobacct->rusage.ru_ixrss,	/* max ixrss */
		 jobacct->rusage.ru_idrss,	/* max idrss */
		 jobacct->rusage.ru_isrss,	/* max isrss */
		 jobacct->rusage.ru_minflt,	/* max minflt */
		 jobacct->rusage.ru_majflt,	/* max majflt */
		 jobacct->rusage.ru_nswap,	/* max nswap */
		 jobacct->rusage.ru_inblock,	/* total inblock */
		 jobacct->rusage.ru_oublock,	/* total outblock */
		 jobacct->rusage.ru_msgsnd,	/* total msgsnd */
		 jobacct->rusage.ru_msgrcv,	/* total msgrcv */
		 jobacct->rusage.ru_nsignals,	/* total nsignals */
		 jobacct->rusage.ru_nvcsw,	/* total nvcsw */
		 jobacct->rusage.ru_nivcsw,	/* total nivcsw */
		 jobacct->max_vsize,	/* max vsize */
		 jobacct->max_vsize_id.taskid,	/* max vsize node */
		 ave_vsize,	/* ave vsize */
		 jobacct->max_rss,	/* max vsize */
		 jobacct->max_rss_id.taskid,	/* max rss node */
		 ave_rss,	/* ave rss */
		 jobacct->max_pages,	/* max pages */
		 jobacct->max_pages_id.taskid,	/* max pages node */
		 ave_pages,	/* ave pages */
		 ave_cpu2,	/* min cpu */
		 jobacct->min_cpu_id.taskid,	/* min cpu node */
		 ave_cpu,	/* ave cpu */
		 step_ptr->name,      	/* step exe name */
		 node_list, /* name of nodes step running on */
		 jobacct->max_vsize_id.nodeid,	/* max vsize task */
		 jobacct->max_rss_id.nodeid,	/* max rss task */
		 jobacct->max_pages_id.nodeid,	/* max pages task */
		 jobacct->min_cpu_id.nodeid,	/* min cpu task */
		 account,
		 step_ptr->job_ptr->requid); /* requester user id */
		 
	return _print_record(step_ptr->job_ptr, now, buf);
}

extern int mysql_jobacct_suspend(struct job_record *job_ptr)
{
		char buf[BUFFER_SIZE];
	static time_t	now = 0;
	static time_t	temp = 0;
	int elapsed;
	if(!jobacct_db_init) {
		debug("jobacct init was not called or it failed");
		return SLURM_ERROR;
	}
	
	/* tell what time has passed */
	if(!now)
		now = job_ptr->start_time;
	temp = now;
	now = time(NULL);
	
	if ((elapsed=now-temp) < 0)
		elapsed=0;	/* For *very* short jobs, if clock is wrong */
	
	/* here we are really just going for a marker in time to tell when
	 * the process was suspended or resumed (check job state), we don't 
	 * really need to keep track of anything else */
	snprintf(buf, BUFFER_SIZE, "%d %u %d",
		 JOB_SUSPEND,
		 elapsed,
		 job_ptr->job_state & (~JOB_COMPLETING));/* job status */
		
	return _print_record(job_ptr, now, buf);
}

extern List mysql_jobacct_getdata()
{
	return NULL;
}


