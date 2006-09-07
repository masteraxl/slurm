/*****************************************************************************\
 *  opt.c - options processing for sattach
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona1@llnl.gov>, et. al.
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <string.h>		/* strcpy, strncasecmp */

#ifdef HAVE_STRINGS_H
#  include <strings.h>
#endif

#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#if HAVE_GETOPT_H
#  include <getopt.h>
#else
#  include "src/common/getopt.h"
#endif

#include <fcntl.h>
#include <stdarg.h>		/* va_start   */
#include <stdio.h>
#include <stdlib.h>		/* getenv     */
#include <pwd.h>		/* getpwuid   */
#include <ctype.h>		/* isdigit    */
#include <sys/param.h>		/* MAXPATHLEN */
#include <sys/stat.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>

#include "src/common/list.h"
#include "src/common/log.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/uid.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/slurm_rlimits_info.h"
#include "src/common/read_config.h" /* contains getnodename() */

#include "src/sattach/opt.h"

#include "src/common/mpi.h"

/* generic getopt_long flags, integers and *not* valid characters */
#define LONG_OPT_INFO_ONLY 0x100

/*---- global variables, defined in opt.h ----*/
opt_t opt;

/*---- forward declarations of static functions  ----*/

typedef struct env_vars env_vars_t;

static void  _help(void);

/* fill in default options  */
static void _opt_default(void);

/* set options based upon env vars  */
static void _opt_env(void);

static void _opt_args(int argc, char **argv);

/* list known options and their settings  */
static void  _opt_list(void);

/* verify options sanity  */
static bool _opt_verify(void);

static void  _print_version(void);

static void _process_env_var(env_vars_t *e, const char *val);

static void  _usage(void);

/*---[ end forward declarations of static functions ]---------------------*/

int initialize_and_process_args(int argc, char *argv[])
{
	/* initialize option defaults */
	_opt_default();

	/* initialize options with env vars */
	_opt_env();

	/* initialize options with argv */
	_opt_args(argc, argv);

	if (opt.verbose > 1)
		_opt_list();

	return 1;

}

static void _print_version(void)
{
	printf("%s %s\n", PACKAGE, SLURM_VERSION);
}

/*
 * print error message to stderr with opt.progname prepended
 */
#undef USE_ARGERROR
#if USE_ARGERROR
static void argerror(const char *msg, ...)
{
	va_list ap;
	char buf[256];

	va_start(ap, msg);
	vsnprintf(buf, sizeof(buf), msg, ap);

	fprintf(stderr, "%s: %s\n",
		opt.progname ? opt.progname : "sbatch", buf);
	va_end(ap);
}
#else
#  define argerror error
#endif				/* USE_ARGERROR */

/*
 * _opt_default(): used by initialize_and_process_args to set defaults
 */
static void _opt_default()
{
	struct passwd *pw;
	static slurm_step_io_fds_t fds = SLURM_STEP_IO_FDS_INITIALIZER;

	if ((pw = getpwuid(getuid())) != NULL) {
		strncpy(opt.user, pw->pw_name, MAX_USERNAME);
		opt.uid = pw->pw_uid;
	} else
		error("who are you?");

	opt.gid = getgid();

	opt.progname = NULL;

	opt.jobid = NO_VAL;
	opt.jobid_set = false;

	opt.quiet = 0;
	opt.verbose = 0;

	opt.euid = (uid_t) -1;
	opt.egid = (gid_t) -1;
	
	opt.labelio = false;
	opt.ctrl_comm_ifhn  = xshort_hostname();
	memcpy(&opt.fds, &fds, sizeof(fds));
	opt.info_only = false;
}

/*---[ env var processing ]-----------------------------------------------*/

/*
 * try to use a similar scheme as popt. 
 * 
 * in order to add a new env var (to be processed like an option):
 *
 * define a new entry into env_vars[], if the option is a simple int
 * or string you may be able to get away with adding a pointer to the
 * option to set. Otherwise, process var based on "type" in _opt_env.
 */
struct env_vars {
	const char *var;
	int type;
	void *arg;
	void *set_flag;
};

env_vars_t env_vars[] = {
  {NULL, 0, NULL, NULL}
};


/*
 * _opt_env(): used by initialize_and_process_args to set options via
 *            environment variables. See comments above for how to
 *            extend srun to process different vars
 */
static void _opt_env()
{
	char       *val = NULL;
	env_vars_t *e   = env_vars;

	while (e->var) {
		if ((val = getenv(e->var)) != NULL) 
			_process_env_var(e, val);
		e++;
	}
}


static void
_process_env_var(env_vars_t *e, const char *val)
{
	debug2("now processing env var %s=%s", e->var, val);

	if (e->set_flag) {
		*((bool *) e->set_flag) = true;
	}

	switch (e->type) {
	default:
		/* do nothing */
		break;
	}
}

void set_options(const int argc, char **argv)
{
	int opt_char, option_index = 0;
	static struct option long_options[] = {
		{"help", 	no_argument,   0, 'h'},
		{"label",       no_argument,   0, 'l'},
		{"quiet",       no_argument,   0, 'q'},
		{"usage",       no_argument,   0, 'u'},
		{"verbose",     no_argument,   0, 'v'},
		{"version",     no_argument,   0, 'V'},
		{"info",        no_argument,   0, LONG_OPT_INFO_ONLY},
		{NULL}
	};
	char *opt_string = "+hlquvV";

	opt.progname = xbasename(argv[0]);
	optind = 0;		
	while((opt_char = getopt_long(argc, argv, opt_string,
				      long_options, &option_index)) != -1) {
		switch (opt_char) {
			
		case '?':
			fprintf(stderr, "Try \"sbatch --help\" for more "
				"information\n");
			exit(1);
			break;
		case 'h':
			_help();
			exit(0);
		case 'l':
			opt.labelio = true;
			break;
		case 'q':
			opt.quiet++;
			break;
		case 'u':
			_usage();
			exit(0);
		case 'v':
			opt.verbose++;
			break;
		case 'V':
			_print_version();
			exit(0);
			break;
		case LONG_OPT_INFO_ONLY:
			opt.info_only = true;
			break;
		default:
			fatal("Unrecognized command line parameter %c",
			      opt_char);
		}
	}
}

static void _parse_jobid_stepid(char *jobid_str)
{
	char *ptr, *job, *step;
	long jobid, stepid;

	verbose("jobid/stepid string = %s\n", jobid_str);
	job = xstrdup(jobid_str);
	ptr = index(job, '.');
	if (ptr == NULL) {
		error("Did not find a period in the step ID string");
		_usage();
		xfree(job);
		exit(1);
	} else {
		*ptr = '\0';
		step = ptr + 1;
	}

	jobid = strtol(job, &ptr, 10);
	if (!xstring_is_whitespace(ptr)) {
		error("\"%s\" does not look like a jobid", job);
		_usage();
		xfree(job);
		exit(1);
	}

	stepid = strtol(step, &ptr, 10);
	if (!xstring_is_whitespace(ptr)) {
		error("\"%s\" does not look like a stepid", step);
		_usage();
		xfree(job);
		exit(1);
	}

	opt.jobid = (uint32_t) jobid;
	opt.stepid = (uint32_t) stepid;

	xfree(job);
}

/*
 * _opt_args() : set options via commandline args and popt
 */
static void _opt_args(int argc, char **argv)
{
	char **rest = NULL;
	int leftover;

	set_options(argc, argv);

	leftover = 0;
	if (optind < argc) {
		rest = argv + optind;
		while (rest[leftover] != NULL)
			leftover++;
	}
	if (leftover != 1) {
		error("too many parameters");
		_usage();
		exit(1);
	}

	_parse_jobid_stepid(*(argv + optind));

	if (!_opt_verify())
		exit(1);
}

/* 
 * _opt_verify : perform some post option processing verification
 *
 */
static bool _opt_verify(void)
{
	bool verified = true;

	if (opt.quiet && opt.verbose) {
		error ("don't specify both --verbose (-v) and --quiet (-Q)");
		verified = false;
	}

	return verified;
}

#define tf_(b) (b == true) ? "true" : "false"

static void _opt_list()
{
	info("defined options for program `%s'", opt.progname);
	info("--------------- ---------------------");
	info("job ID         : %u", opt.jobid);
	info("step ID        : %u", opt.stepid);
	info("user           : `%s'", opt.user);
	info("uid            : %ld", (long) opt.uid);
	info("gid            : %ld", (long) opt.gid);
	info("verbose        : %d", opt.verbose);
}

static void _usage(void)
{
 	printf("Usage: sattach [options] <jobid.stepid>\n");
}

static void _help(void)
{
        printf("Usage: sattach <jobid.stepid>\n");
	printf(
"  -v, --verbose               verbose mode (multiple -v's increase verbosity)\n"
"  -Q, --quiet                 quiet mode (suppress informational messages)\n"
		);
	printf("\n");

	printf(
"Help options:\n"
"  -h, --help                  show this help message\n"
"  -u, --usage                 display brief usage message\n"
"\n"
"Other options:\n"
"  -V, --version               output version information and exit\n"
"\n"
);

}
