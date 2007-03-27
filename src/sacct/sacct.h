/*****************************************************************************\
 *  sacct.h - header file for sacct
 *
 *  $Id: sacct.h 7541 2006-03-18 01:44:58Z da $
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>.
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
#ifndef _SACCT_H
#define _SACCT_H

#include <ctype.h>
#include <errno.h>
#include <grp.h>
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "src/common/getopt.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/list.h"
#include "src/common/hostlist.h"
#include "src/common/slurm_jobacct.h"

#define ERROR 2

#define BRIEF_FIELDS "jobid,status,exitcode"
#define DEFAULT_FIELDS "jobid,jobname,partition,ncpus,status,exitcode"
#define STAT_FIELDS "jobid,vsize,rss,pages,cputime,ntasks,status"
#define LONG_FIELDS "jobid,jobname,partition,vsize,rss,pages,cputime,ntasks,ncpus,elapsed,status,exitcode"

#define BUFFER_SIZE 4096
#define STATUS_COUNT 10

#define MAX_PRINTFIELDS 100

#define SECONDS_IN_MINUTE 60
#define SECONDS_IN_HOUR (60*SECONDS_IN_MINUTE)
#define SECONDS_IN_DAY (24*SECONDS_IN_HOUR)

/* On output, use fields 12-37 from JOB_STEP */

typedef enum {	HEADLINE,
		UNDERSCORE,
		JOB,
		JOBSTEP
} type_t;


typedef struct fields {
	char *name;		/* Specified in --fields= */
	void (*print_routine) ();	/* Who gets to print it? */
} fields_t;

extern fields_t fields[];
extern sacct_parameters_t params;

extern List jobs;

extern int printfields[MAX_PRINTFIELDS],	/* Indexed into fields[] */
	nprintfields;

/* process.c */
void find_hostname(uint32_t pos, char *hosts, char *host);
void aggregate_sacct(sacct_t *dest, sacct_t *from);

/* print.c */
void print_fields(type_t type, void *object);
void print_cpu(type_t type, void *object);
void print_elapsed(type_t type, void *object);
void print_exitcode(type_t type, void *object);
void print_gid(type_t type, void *object);
void print_group(type_t type, void *object);
void print_idrss(type_t type, void *object);
void print_inblocks(type_t type, void *object);
void print_isrss(type_t type, void *object);
void print_ixrss(type_t type, void *object);
void print_job(type_t type, void *object);
void print_name(type_t type, void *object);
void print_jobid(type_t type, void *object);
void print_majflt(type_t type, void *object);
void print_minflt(type_t type, void *object);
void print_msgrcv(type_t type, void *object);
void print_msgsnd(type_t type, void *object);
void print_ncpus(type_t type, void *object);
void print_nivcsw(type_t type, void *object);
void print_nodes(type_t type, void *object);
void print_nsignals(type_t type, void *object);
void print_nswap(type_t type, void *object);
void print_ntasks(type_t type, void *object);
void print_nvcsw(type_t type, void *object);
void print_outblocks(type_t type, void *object);
void print_partition(type_t type, void *object);
void print_blockid(type_t type, void *object);
void print_pages(type_t type, void *object);
void print_rss(type_t type, void *object);
void print_status(type_t type, void *object);
void print_submit(type_t type, void *object);
void print_start(type_t type, void *object);
void print_end(type_t type, void *object);
void print_systemcpu(type_t type, void *object);
void print_uid(type_t type, void *object);
void print_user(type_t type, void *object);
void print_usercpu(type_t type, void *object);
void print_vsize(type_t type, void *object);
void print_cputime(type_t type, void *object);
void print_account(type_t type, void *object);

/* options.c */
int decode_status_char(char *status);
char *decode_status_int(int status);
int get_data(void);
void parse_command_line(int argc, char **argv);
void do_dump(void);
void do_expire(void);
void do_help(void);
void do_list(void);
void do_stat(void);
void sacct_init();
void sacct_fini();

/* sacct_stat.c */
extern int sacct_stat(uint32_t jobid, uint32_t stepid);

#endif /* !_SACCT_H */
