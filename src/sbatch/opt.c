/*****************************************************************************\
 *  opt.c - options processing for sbatch
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona1@llnl.gov>, et. al.
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

#include "src/sbatch/opt.h"

/* generic OPT_ definitions -- mainly for use with env vars  */
#define OPT_NONE        0x00
#define OPT_INT         0x01
#define OPT_STRING      0x02
#define OPT_DEBUG       0x03
#define OPT_NODES       0x05
#define OPT_BOOL        0x06
#define OPT_CORE        0x07
#define OPT_CONN_TYPE	0x08
#define OPT_NO_ROTATE	0x0a
#define OPT_GEOMETRY	0x0b
#define OPT_MULTI	0x0f
#define OPT_EXCLUSIVE	0x10

/* generic getopt_long flags, integers and *not* valid characters */
#define LONG_OPT_JOBID       0x105
#define LONG_OPT_TMP         0x106
#define LONG_OPT_MEM         0x107
#define LONG_OPT_MINCPU      0x108
#define LONG_OPT_CONT        0x109
#define LONG_OPT_UID         0x10a
#define LONG_OPT_GID         0x10b
#define LONG_OPT_MINSOCKETS  0x10c
#define LONG_OPT_MINCORES    0x10d
#define LONG_OPT_MINTHREADS  0x10e
#define LONG_OPT_CORE	     0x10f
#define LONG_OPT_CONNTYPE    0x110
#define LONG_OPT_EXCLUSIVE   0x111
#define LONG_OPT_BEGIN       0x112
#define LONG_OPT_MAIL_TYPE   0x113
#define LONG_OPT_MAIL_USER   0x114
#define LONG_OPT_NICE        0x115
#define LONG_OPT_NO_REQUEUE  0x116
#define LONG_OPT_COMMENT     0x117
#define LONG_OPT_WRAP        0x118
#define LONG_OPT_BLRTS_IMAGE     0x140
#define LONG_OPT_LINUX_IMAGE     0x141
#define LONG_OPT_MLOADER_IMAGE   0x142
#define LONG_OPT_RAMDISK_IMAGE   0x143
#define LONG_OPT_REBOOT          0x144
#define LONG_OPT_TASKSPERNODE    0x145
#define LONG_OPT_GET_USER_ENV    0x146

/*---- global variables, defined in opt.h ----*/
opt_t opt;

/*---- forward declarations of static functions  ----*/

typedef struct env_vars env_vars_t;

/* return command name from its full path name */
static char * _base_name(char* command);

static List  _create_path_list(void);

/* Get a decimal integer from arg */
static int  _get_int(const char *arg, const char *what);

static void  _help(void);

/* fill in default options  */
static void _opt_default(void);

/* set options from batch script */
static void _opt_batch_script(const void *body, int size);

/* set options based upon env vars  */
static void _opt_env(void);

/* list known options and their settings  */
static void  _opt_list(void);

/* verify options sanity  */
static bool _opt_verify(void);

static void  _print_version(void);

static void _process_env_var(env_vars_t *e, const char *val);

static uint16_t _parse_mail_type(const char *arg);
static char *_print_mail_type(const uint16_t type);

/* search PATH for command returns full path */
static char *_search_path(char *, bool, int);

static long  _to_bytes(const char *arg);

static void  _usage(void);
static bool  _verify_node_count(const char *arg, int *min, int *max);
static int   _verify_geometry(const char *arg, uint16_t *geometry);
static int   _verify_conn_type(const char *arg);
static char *_fullpath(const char *filename);
static void _set_options(int argc, char **argv);

/*---[ end forward declarations of static functions ]---------------------*/

static void _print_version(void)
{
	printf("%s %s\n", PACKAGE, SLURM_VERSION);
}

/*
 * verify that a connection type in arg is of known form
 * returns the connection_type or -1 if not recognized
 */
static int _verify_conn_type(const char *arg)
{
	int len = strlen(arg);

	if (!strncasecmp(arg, "MESH", len))
		return SELECT_MESH;
	else if (!strncasecmp(arg, "TORUS", len))
		return SELECT_TORUS;
	else if (!strncasecmp(arg, "NAV", len))
		return SELECT_NAV;

	error("invalid --conn-type argument %s ignored.", arg);
	return -1;
}

/*
 * verify geometry arguments, must have proper count
 * returns -1 on error, 0 otherwise
 */
static int _verify_geometry(const char *arg, uint16_t *geometry)
{
	char* token, *delimiter = ",x", *next_ptr;
	int i, rc = 0;
	char* geometry_tmp = xstrdup(arg);
	char* original_ptr = geometry_tmp;

	token = strtok_r(geometry_tmp, delimiter, &next_ptr);
	for (i=0; i<SYSTEM_DIMENSIONS; i++) {
		if (token == NULL) {
			error("insufficient dimensions in --geometry");
			rc = -1;
			break;
		}
		geometry[i] = (uint16_t)atoi(token);
		if (geometry[i] == 0 || geometry[i] == (uint16_t)NO_VAL) {
			error("invalid --geometry argument");
			rc = -1;
			break;
		}
		geometry_tmp = next_ptr;
		token = strtok_r(geometry_tmp, delimiter, &next_ptr);
	}
	if (token != NULL) {
		error("too many dimensions in --geometry");
		rc = -1;
	}

	if (original_ptr)
		xfree(original_ptr);

	return rc;
}

/* 
 * verify that a node count in arg is of a known form (count or min-max)
 * OUT min, max specified minimum and maximum node counts
 * RET true if valid
 */
static bool 
_verify_node_count(const char *arg, int *min_nodes, int *max_nodes)
{
	char *end_ptr;
	double val1, val2;
	
	val1 = strtod(arg, &end_ptr);
	if (end_ptr[0] == 'k' || end_ptr[0] == 'K') {
		val1 *= 1024;
		end_ptr++;
	}

 	if (end_ptr[0] == '\0') {
		*min_nodes = val1;
		return true;
	}
	
	if (end_ptr[0] != '-')
		return false;

	val2 = strtod(&end_ptr[1], &end_ptr);
	if (end_ptr[0] == 'k' || end_ptr[0] == 'K') {
		val2 *= 1024;
		end_ptr++;
	}

	if (end_ptr[0] == '\0') {
		*min_nodes = val1;
		*max_nodes = val2;
		return true;
	} else
		return false;

}

/* return command name from its full path name */
static char * _base_name(char* command)
{
	char *char_ptr, *name;
	int i;

	if (command == NULL)
		return NULL;

	char_ptr = strrchr(command, (int)'/');
	if (char_ptr == NULL)
		char_ptr = command;
	else
		char_ptr++;

	i = strlen(char_ptr);
	name = xmalloc(i+1);
	strcpy(name, char_ptr);
	return name;
}

/*
 * _to_bytes(): verify that arg is numeric with optional "G" or "M" at end
 * if "G" or "M" is there, multiply by proper power of 2 and return
 * number in bytes
 */
static long _to_bytes(const char *arg)
{
	char *buf;
	char *endptr;
	int end;
	int multiplier = 1;
	long result;

	buf = xstrdup(arg);

	end = strlen(buf) - 1;

	if (isdigit(buf[end])) {
		result = strtol(buf, &endptr, 10);

		if (*endptr != '\0')
			result = -result;

	} else {

		switch (toupper(buf[end])) {

		case 'G':
			multiplier = 1024;
			break;

		case 'M':
			/* do nothing */
			break;

		default:
			multiplier = -1;
		}

		buf[end] = '\0';

		result = multiplier * strtol(buf, &endptr, 10);

		if (*endptr != '\0')
			result = -result;
	}

	return result;
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
	char buf[MAXPATHLEN + 1];
	struct passwd *pw;
	int i;

	if ((pw = getpwuid(getuid())) != NULL) {
		strncpy(opt.user, pw->pw_name, MAX_USERNAME);
		opt.uid = pw->pw_uid;
	} else
		error("who are you?");

	opt.script_argc = 0;
	opt.script_argv = NULL;

	opt.gid = getgid();

	if ((getcwd(buf, MAXPATHLEN)) == NULL) 
		fatal("getcwd failed: %m");
	opt.cwd = xstrdup(buf);

	opt.progname = NULL;

	opt.nprocs = 1;
	opt.nprocs_set = false;
	opt.cpus_per_task = 1; 
	opt.cpus_set = false;
	opt.min_nodes = 1;
	opt.max_nodes = 0;
	opt.tasks_per_node   = -1;
	opt.nodes_set = false;
	opt.time_limit = 0;
	opt.partition = NULL;

	opt.job_name = NULL;
	opt.jobid    = NO_VAL;
	opt.jobid_set = false;
	opt.dependency = NO_VAL;
	opt.account  = NULL;
	opt.comment  = NULL;

	opt.shared = (uint16_t)NO_VAL;
	opt.no_kill = false;

	opt.immediate	= false;
	opt.no_requeue	= false;

	opt.quiet = 0;
	opt.verbose = 0;

	/* constraint default (-1 is no constraint) */
	opt.mincpus	    = -1;
	opt.minsockets      = -1;
	opt.mincores        = -1;
	opt.minthreads      = -1;
	opt.realmem	    = -1;
	opt.tmpdisk	    = -1;

	opt.hold	    = false;
	opt.constraints	    = NULL;
	opt.contiguous	    = false;
	opt.nodelist	    = NULL;
	opt.exc_nodes	    = NULL;

	for (i=0; i<SYSTEM_DIMENSIONS; i++)
		opt.geometry[i]	    = (uint16_t) NO_VAL;
	opt.reboot          = false;
	opt.no_rotate	    = false;
	opt.conn_type	    = (uint16_t) NO_VAL;

	opt.euid	    = (uid_t) -1;
	opt.egid	    = (gid_t) -1;
	
	opt.ifname = xstrdup("/dev/null");
	opt.ofname = NULL;
	opt.efname = NULL;
	opt.get_user_env = false;
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
  {"SBATCH_ACCOUNT",       OPT_STRING,     &opt.account,       NULL           },
  {"SBATCH_BLRTS_IMAGE",   OPT_STRING,     &opt.blrtsimage,    NULL           },
  {"SBATCH_CONN_TYPE",     OPT_CONN_TYPE,  NULL,               NULL           },
  {"SBATCH_DEBUG",         OPT_DEBUG,      NULL,               NULL           },
  {"SBATCH_GEOMETRY",      OPT_GEOMETRY,   NULL,               NULL           },
  {"SBATCH_IMMEDIATE",     OPT_BOOL,       &opt.immediate,     NULL           },
  {"SBATCH_JOBID",         OPT_INT,        &opt.jobid,         NULL           },
  {"SBATCH_JOB_NAME",      OPT_STRING,     &opt.job_name,      NULL           },
  {"SBATCH_LINUX_IMAGE",   OPT_STRING,     &opt.linuximage,    NULL           },
  {"SBATCH_MLOADER_IMAGE", OPT_STRING,     &opt.mloaderimage,  NULL           },
  {"SBATCH_NO_REQUEUE",    OPT_BOOL,       &opt.no_requeue,    NULL           },
  {"SBATCH_NO_ROTATE",     OPT_BOOL,       &opt.no_rotate,     NULL           },
  {"SBATCH_PARTITION",     OPT_STRING,     &opt.partition,     NULL           },
  {"SBATCH_RAMDISK_IMAGE", OPT_STRING,     &opt.ramdiskimage,  NULL           },
  {"SBATCH_TIMELIMIT",     OPT_STRING,     &opt.time_limit_str,NULL           },
  {"SBATCH_EXCLUSIVE",     OPT_EXCLUSIVE,  NULL,               NULL           },
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
	char *end = NULL;

	debug2("now processing env var %s=%s", e->var, val);

	if (e->set_flag) {
		*((bool *) e->set_flag) = true;
	}

	switch (e->type) {
	case OPT_STRING:
		*((char **) e->arg) = xstrdup(val);
		break;
	case OPT_INT:
		if (val != NULL) {
			*((int *) e->arg) = (int) strtol(val, &end, 10);
			if (!(end && *end == '\0')) 
				error("%s=%s invalid. ignoring...", e->var, val);
		}
		break;

	case OPT_BOOL:
		/* A boolean env variable is true if:
		 *  - set, but no argument
		 *  - argument is "yes"
		 *  - argument is a non-zero number
		 */
		if (val == NULL || strcmp(val, "") == 0) {
			*((bool *)e->arg) = true;
		} else if (strcasecmp(val, "yes") == 0) {
			*((bool *)e->arg) = true;
		} else if ((strtol(val, &end, 10) != 0)
			   && end != val) {
			*((bool *)e->arg) = true;
		} else {
			*((bool *)e->arg) = false;
		}
		break;

	case OPT_DEBUG:
		if (val != NULL) {
			opt.verbose = (int) strtol(val, &end, 10);
			if (!(end && *end == '\0')) 
				error("%s=%s invalid", e->var, val);
		}
		break;

	case OPT_NODES:
		opt.nodes_set = _verify_node_count( val, 
						    &opt.min_nodes, 
						    &opt.max_nodes );
		if (opt.nodes_set == false) {
			error("\"%s=%s\" -- invalid node count. ignoring...",
			      e->var, val);
		}
		break;

	case OPT_CONN_TYPE:
		opt.conn_type = _verify_conn_type(val);
		break;
	
	case OPT_NO_ROTATE:
		opt.no_rotate = true;
		break;

	case OPT_GEOMETRY:
		if (_verify_geometry(val, opt.geometry)) {
			error("\"%s=%s\" -- invalid geometry, ignoring...",
			      e->var, val);
		}
		break;

	case OPT_EXCLUSIVE:
		opt.shared = 0;
		break;

	default:
		/* do nothing */
		break;
	}
}


/*---[ command line option processing ]-----------------------------------*/

static struct option long_options[] = {
	{"batch",         no_argument,       0, 'b'}, /* batch option
							 is only here for
							 moab tansition
							 doesn't do anything */
	{"cpus-per-task", required_argument, 0, 'c'},
	{"constraint",    required_argument, 0, 'C'},
	{"dependency",    required_argument, 0, 'd'},
	{"workdir",       required_argument, 0, 'D'},
	{"error",         required_argument, 0, 'e'},
	{"nodefile",      required_argument, 0, 'F'},
	{"geometry",      required_argument, 0, 'g'},
	{"help",          no_argument,       0, 'h'},
	{"hold",          no_argument,       0, 'H'}, /* undocumented */
	{"input",         required_argument, 0, 'i'},
	{"immediate",     no_argument,       0, 'I'},
	{"job-name",      required_argument, 0, 'J'},
	{"no-kill",       no_argument,       0, 'k'},
	{"tasks",         required_argument, 0, 'n'},
	{"nodes",         required_argument, 0, 'N'},
	{"output",        required_argument, 0, 'o'},
	{"partition",     required_argument, 0, 'p'},
	{"quiet",         no_argument,       0, 'q'},
	{"no-rotate",     no_argument,       0, 'R'},
	{"share",         no_argument,       0, 's'},
	{"time",          required_argument, 0, 't'},
	{"usage",         no_argument,       0, 'u'},
	{"account",       required_argument, 0, 'U'},
	{"verbose",       no_argument,       0, 'v'},
	{"version",       no_argument,       0, 'V'},
	{"nodelist",      required_argument, 0, 'w'},
	{"exclude",       required_argument, 0, 'x'},
	{"contiguous",    no_argument,       0, LONG_OPT_CONT},
	{"exclusive",     no_argument,       0, LONG_OPT_EXCLUSIVE},
	{"mincpus",       required_argument, 0, LONG_OPT_MINCPU},
	{"minsockets",    required_argument, 0, LONG_OPT_MINSOCKETS},
	{"mincores",      required_argument, 0, LONG_OPT_MINCORES},
	{"minthreads",    required_argument, 0, LONG_OPT_MINTHREADS},
	{"mem",           required_argument, 0, LONG_OPT_MEM},
	{"tmp",           required_argument, 0, LONG_OPT_TMP},
	{"jobid",         required_argument, 0, LONG_OPT_JOBID},
	{"uid",           required_argument, 0, LONG_OPT_UID},
	{"gid",           required_argument, 0, LONG_OPT_GID},
	{"conn-type",     required_argument, 0, LONG_OPT_CONNTYPE},
	{"begin",         required_argument, 0, LONG_OPT_BEGIN},
	{"mail-type",     required_argument, 0, LONG_OPT_MAIL_TYPE},
	{"mail-user",     required_argument, 0, LONG_OPT_MAIL_USER},
	{"nice",          optional_argument, 0, LONG_OPT_NICE},
	{"no-requeue",    no_argument,       0, LONG_OPT_NO_REQUEUE},
	{"comment",       required_argument, 0, LONG_OPT_COMMENT},
	{"blrts-image",   required_argument, 0, LONG_OPT_BLRTS_IMAGE},
	{"linux-image",   required_argument, 0, LONG_OPT_LINUX_IMAGE},
	{"mloader-image", required_argument, 0, LONG_OPT_MLOADER_IMAGE},
	{"ramdisk-image", required_argument, 0, LONG_OPT_RAMDISK_IMAGE},
	{"reboot",        no_argument,       0, LONG_OPT_REBOOT},
	{"tasks-per-node",  required_argument,0,LONG_OPT_TASKSPERNODE},
	{"ntasks-per-node", required_argument,0,LONG_OPT_TASKSPERNODE}, 
	{"wrap",          required_argument, 0, LONG_OPT_WRAP},
	{"get-user-env",  no_argument,       0, LONG_OPT_GET_USER_ENV},
	{NULL,            0,                 0, 0}
};

static char *opt_string =
	"+a:bc:C:d:D:e:F:g:hHi:IJ:kn:N:o:Op:qR:st:uU:vVw:x:";


/*
 * process_options_first_pass()
 *
 * In this first pass we only look at the command line options, and we
 * will only handle a few options (help, usage, quiet, verbose, version),
 * and look for the script name and arguments (if provided).
 *
 * We will parse the environment variable options, batch script options,
 * and all of the rest of the command line options in
 * process_options_second_pass().
 *
 * Return a pointer to the batch script file name is provided on the command
 * line, otherwise return NULL, and the script will need to be read from
 * standard input.
 */
char *process_options_first_pass(int argc, char **argv)
{
	int opt_char, option_index = 0;
	char *str = NULL;

	/* initialize option defaults */
	_opt_default();

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
		case LONG_OPT_WRAP:
			opt.wrap = xstrdup(optarg);
			break;
		default:
			/* will be parsed in second pass function */
			break;
		}
	}
	xfree(str);

	if (argc > optind && opt.wrap != NULL) {
		fatal("Script arguments are not permitted with the"
		      " --wrap option.");
	}
	if (argc > optind) {
		int i;
		char **leftover;

		opt.script_argc = argc - optind;
		leftover = argv + optind;
		opt.script_argv = (char **) xmalloc((opt.script_argc + 1)
						    * sizeof(char *));
		for (i = 0; i < opt.script_argc; i++)
			opt.script_argv[i] = xstrdup(leftover[i]);
		opt.script_argv[i] = NULL;
	}
	if (opt.script_argc > 0) {
		char *fullpath;
		char *cmd       = opt.script_argv[0];
		int  mode       = R_OK;

		if ((fullpath = _search_path(cmd, true, mode))) {
			xfree(opt.script_argv[0]);
			opt.script_argv[0] = fullpath;
		} 

		return opt.script_argv[0];
	} else {
		return NULL;
	}
}

/* process options:
 * 1. update options with option set in the script
 * 2. update options with env vars
 * 3. update options with commandline args
 * 4. perform some verification that options are reasonable
 */
int process_options_second_pass(int argc, char *argv[],
				const void *script_body, int script_size)
{
	/* set options from batch script */
	_opt_batch_script(script_body, script_size);

	/* set options from env vars */
	_opt_env();

	/* set options from command line */
	_set_options(argc, argv);

	if (!_opt_verify())
		exit(1);

	if (opt.verbose > 3)
		_opt_list();

	return 1;

}

/*
 * _next_line - Interpret the contents of a byte buffer as characters in
 *	a file.  _next_line will find and return the next line in the buffer.
 *
 *	If "state" is NULL, it will start at the beginning of the buffer.
 *	_next_line will update the "state" pointer to point at the
 *	spot in the buffer where it left off.
 *
 * IN buf - buffer containing file contents
 * IN size - size of buffer "buf"
 * IN/OUT state - used by _next_line to determine where the last line ended
 *
 * RET - xmalloc'ed character string, or NULL if no lines remaining in buf.
 */
static char *_next_line(const void *buf, int size, void **state)
{
	char *line;
	char *current, *ptr;

	if (*state == NULL) /* initial state */
		*state = (void *)buf;

	if ((*state - buf) >= size) /* final state */
		return NULL;

	ptr = current = (char *)*state;
	while ((*ptr != '\n') && (ptr < ((char *)buf+size)))
		ptr++;
	if (*ptr == '\n')
		ptr++;
	
	line = xstrndup(current, (ptr-current));

	*state = (void *)ptr;
	return line;
}

/*
 * _get_argument - scans a line for something that looks like a command line
 *	argument, and return an xmalloc'ed string containing the argument.
 *	Quotes can be used to group characters, including whitespace.
 *	Quotes can be included in an argument be escaping the quotes,
 *	preceding the quote with a backslash (\").
 *
 * IN - line
 * OUT - skipped - number of characters parsed from line
 * RET - xmalloc'ed argument string (may be shorter than "skipped")
 *       or NULL if no arguments remaining
 */
static char *_get_argument(const char *line, int *skipped)
{
	char *ptr;
	char argument[BUFSIZ];
	bool escape_flag = false;
	bool no_isspace_check = false;
	int i;

	ptr = (char *)line;
	*skipped = 0;

	/* skip whitespace */
	while (isspace(*ptr) && *ptr != '\0') {
		ptr++;
	}

	if (*ptr == '\0')
		return NULL;

	/* copy argument into "argument" buffer, */
	i = 0;
	while ((no_isspace_check || !isspace(*ptr))
	       && *ptr != '\n'
	       && *ptr != '\0') {

		if (escape_flag) {
			escape_flag = false;
			argument[i] = *ptr;
			ptr++;
			i++;
		} else if (*ptr == '\\') {
			escape_flag = true;
			ptr++;
		} else if (*ptr == '"') {
			/* toggle the no_isspace_check flag */
			no_isspace_check = no_isspace_check? false : true;
			ptr++;
		} else if (*ptr == '#') {
			/* found an un-escaped #, rest of line is a comment */
			break;
		} else {
			argument[i] = *ptr;
			ptr++;
			i++;
		}
	}

	*skipped = ptr - line;
	if (i > 0) {
		return xstrndup(argument, i);
	} else {
		return NULL;
	}
}

/*
 * set options from batch script
 *
 * Build an argv-style array of options from the script "body",
 * then pass the array to _set_options for() further parsing.
 */
static void _opt_batch_script(const void *body, int size)
{
	char *magic_word = "#SBATCH";
	int magic_word_len;
	int argc;
	char **argv;
	void *state = NULL;
	char *line;
	char *option;
	char *ptr;
	int skipped = 0;
	int i;

	magic_word_len = strlen(magic_word);
	/* getopt_long skips over the first argument, so fill it in */
	argc = 1;
	argv = xmalloc(sizeof(char *));
	argv[0] = "sbatch";

	while((line = _next_line(body, size, &state)) != NULL) {
		if (strncmp(line, magic_word, magic_word_len) != 0) {
			xfree(line);
			continue;
		}

		/* this line starts with the magic word */
		ptr = line + magic_word_len;
		while ((option = _get_argument(ptr, &skipped)) != NULL) {
			debug2("Found in script, argument \"%s\"", option);
			argc += 1;
			xrealloc(argv, sizeof(char*) * argc);
			argv[argc-1] = option;
			ptr += skipped;
		}
		xfree(line);
	}

	if (argc > 0)
		_set_options(argc, argv);

	for (i = 1; i < argc; i++)
		xfree(argv[i]);
	xfree(argv);
}

static void _set_options(int argc, char **argv)
{
	int opt_char, option_index = 0;
	char *tmp;

	optind = 0;
	while((opt_char = getopt_long(argc, argv, opt_string,
				      long_options, &option_index)) != -1) {
		switch (opt_char) {
		case '?':
			fatal("Try \"sbatch --help\" for more information");
			break;
		case 'b':
			/* Only here for Moab transition not suppose
			   to do anything */
			break;
		case 'c':
			opt.cpus_set = true;
			opt.cpus_per_task = 
				_get_int(optarg, "cpus-per-task");
			break;
		case 'C':
			xfree(opt.constraints);
			opt.constraints = xstrdup(optarg);
			break;
		case 'd':
			opt.dependency = _get_int(optarg, "dependency");
			break;
		case 'D':
			xfree(opt.cwd);
			opt.cwd = xstrdup(optarg);
			break;
		case 'e':
			xfree(opt.efname);
			if (strncasecmp(optarg, "none", (size_t)4) == 0)
				opt.efname = xstrdup("/dev/null");
			else
				opt.efname = _fullpath(optarg);
			break;
		case 'F':
			xfree(opt.nodelist);
			tmp = slurm_read_hostfile(optarg, 0);
			if (tmp != NULL) {
				opt.nodelist = xstrdup(tmp);
				free(tmp);
			} else {
				error("\"%s\" is not a valid node file");
				exit(1);
			}
			break;
		case 'g':
			if (_verify_geometry(optarg, opt.geometry))
				exit(1);
			break;
		case 'h':
			_help();
			exit(0);
		case 'H':
			opt.hold = true;
			break;
		case 'i':
			xfree(opt.ifname);
			if (strncasecmp(optarg, "none", (size_t)4) == 0)
				opt.ifname = xstrdup("/dev/null");
			else
				opt.ifname = _fullpath(optarg);
			break;
		case 'I':
			opt.immediate = true;
			break;
		case 'J':
			xfree(opt.job_name);
			opt.job_name = xstrdup(optarg);
			break;
		case 'k':
			opt.no_kill = true;
			break;
		case 'n':
			opt.nprocs_set = true;
			opt.nprocs = 
				_get_int(optarg, "number of tasks");
			break;
		case 'N':
			opt.nodes_set = 
				_verify_node_count(optarg, 
						   &opt.min_nodes,
						   &opt.max_nodes);
			if (opt.nodes_set == false) {
				error("invalid node count `%s'", 
				      optarg);
				exit(1);
			}
			break;
		case 'o':
			xfree(opt.ofname);
			if (strncasecmp(optarg, "none", (size_t)4) == 0)
				opt.ofname = xstrdup("/dev/null");
			else
				opt.ofname = _fullpath(optarg);
			break;
		case 'p':
			xfree(opt.partition);
			opt.partition = xstrdup(optarg);
			break;
		case 'q':
			opt.quiet++;
			break;
		case 'R':
			opt.no_rotate = true;
			break;
		case 's':
			opt.shared = 1;
			break;
		case 't':
			xfree(opt.time_limit_str);
			opt.time_limit_str = xstrdup(optarg);
			break;
		case 'u':
			_usage();
			exit(0);
		case 'U':
			xfree(opt.account);
			opt.account = xstrdup(optarg);
			break;
		case 'v':
			opt.verbose++;
			break;
		case 'V':
			_print_version();
			exit(0);
			break;
		case 'w':
			xfree(opt.nodelist);
			opt.nodelist = xstrdup(optarg);
#ifdef HAVE_BG
			info("\tThe nodelist option should only be used if\n"
			     "\tthe block you are asking for can be created.\n"
			     "\tPlease consult smap before using this option\n"
			     "\tor your job may be stuck with no way to run.");
#endif
			break;
		case 'x':
			xfree(opt.exc_nodes);
			opt.exc_nodes = xstrdup(optarg);
			break;
		case LONG_OPT_CONT:
			opt.contiguous = true;
			break;
                case LONG_OPT_EXCLUSIVE:
                        opt.shared = 0;
                        break;
		case LONG_OPT_MINCPU:
			opt.mincpus = _get_int(optarg, "mincpus");
			if (opt.mincpus < 0) {
				error("invalid mincpus constraint %s", 
				      optarg);
				exit(1);
			}
			break;
		case LONG_OPT_MINSOCKETS:
			opt.minsockets = _get_int(optarg, "minsockets");
			if (opt.minsockets < 0) {
				error("invalid minsockets constraint %s", 
				      optarg);
				exit(1);
			}
			break;
		case LONG_OPT_MINCORES:
			opt.mincores = _get_int(optarg, "mincores");
			if (opt.mincores < 0) {
				error("invalid mincores constraint %s", 
				      optarg);
				exit(1);
			}
			break;
		case LONG_OPT_MINTHREADS:
			opt.minthreads = _get_int(optarg, "minthreads");
			if (opt.minthreads < 0) {
				error("invalid minthreads constraint %s", 
				      optarg);
				exit(1);
			}
			break;
		case LONG_OPT_MEM:
			opt.realmem = (int) _to_bytes(optarg);
			if (opt.realmem < 0) {
				error("invalid memory constraint %s", 
				      optarg);
				exit(1);
			}
			break;
		case LONG_OPT_TMP:
			opt.tmpdisk = _to_bytes(optarg);
			if (opt.tmpdisk < 0) {
				error("invalid tmp value %s", optarg);
				exit(1);
			}
			break;
		case LONG_OPT_JOBID:
			opt.jobid = _get_int(optarg, "jobid");
			opt.jobid_set = true;
			break;
		case LONG_OPT_UID:
			opt.euid = uid_from_string (optarg);
			if (opt.euid == (uid_t) -1)
				fatal ("--uid=\"%s\" invalid", optarg);
			break;
		case LONG_OPT_GID:
			opt.egid = gid_from_string (optarg);
			if (opt.egid == (gid_t) -1)
				fatal ("--gid=\"%s\" invalid", optarg);
			break;
		case LONG_OPT_CONNTYPE:
			opt.conn_type = _verify_conn_type(optarg);
			break;
		case LONG_OPT_BEGIN:
			opt.begin = parse_time(optarg);
			break;
		case LONG_OPT_MAIL_TYPE:
			opt.mail_type |= _parse_mail_type(optarg);
			if (opt.mail_type == 0)
				fatal("--mail-type=%s invalid", optarg);
			break;
		case LONG_OPT_MAIL_USER:
			xfree(opt.mail_user);
			opt.mail_user = xstrdup(optarg);
			break;
		case LONG_OPT_NICE:
			if (optarg)
				opt.nice = strtol(optarg, NULL, 10);
			else
				opt.nice = 100;
			if (abs(opt.nice) > NICE_OFFSET) {
				error("Invalid nice value, must be between "
				      "-%d and %d", NICE_OFFSET, NICE_OFFSET);
				exit(1);
			}
			break;
		case LONG_OPT_NO_REQUEUE:
			opt.no_requeue = true;
			break;
		case LONG_OPT_COMMENT:
			xfree(opt.comment);
			opt.comment = xstrdup(optarg);
			break;
		case LONG_OPT_BLRTS_IMAGE:
			xfree(opt.blrtsimage);
			opt.blrtsimage = xstrdup(optarg);
			break;
		case LONG_OPT_LINUX_IMAGE:
			xfree(opt.linuximage);
			opt.linuximage = xstrdup(optarg);
			break;
		case LONG_OPT_MLOADER_IMAGE:
			xfree(opt.mloaderimage);
			opt.mloaderimage = xstrdup(optarg);
			break;
		case LONG_OPT_RAMDISK_IMAGE:
			xfree(opt.ramdiskimage);
			opt.ramdiskimage = xstrdup(optarg);
			break;
		case LONG_OPT_REBOOT:
			opt.reboot = true;
			break;
		case LONG_OPT_TASKSPERNODE:
			opt.tasks_per_node = _get_int(optarg, "ntasks-per-node");
			break;
		case LONG_OPT_WRAP:
			/* handled in process_options_first_pass() */
			break;
		case LONG_OPT_GET_USER_ENV:
			opt.get_user_env = true;
			break;
		default:
			fatal("Unrecognized command line parameter %c",
			      opt_char);
		}
	}

	if (optind < argc) {
		fatal("Invalid argument: %s", argv[optind]);
	}
}

/* 
 * _opt_verify : perform some post option processing verification
 *
 */
static bool _opt_verify(void)
{
	bool verified = true;

	if (opt.quiet && opt.verbose) {
		error ("don't specify both --verbose (-v) and --quiet (-q)");
		verified = false;
	}

	if (opt.mincpus < opt.cpus_per_task)
		opt.mincpus = opt.cpus_per_task;

	if ((opt.job_name == NULL) && (opt.script_argc > 0))
		opt.job_name = _base_name(opt.script_argv[0]);

	/* check for realistic arguments */
	if (opt.nprocs <= 0) {
		error("%s: invalid number of processes (-n %d)",
		      opt.progname, opt.nprocs);
		verified = false;
	}

	if (opt.cpus_per_task <= 0) {
		error("%s: invalid number of cpus per task (-c %d)\n",
		      opt.progname, opt.cpus_per_task);
		verified = false;
	}

	if ((opt.min_nodes <= 0) || (opt.max_nodes < 0) || 
	    (opt.max_nodes && (opt.min_nodes > opt.max_nodes))) {
		error("%s: invalid number of nodes (-N %d-%d)\n",
		      opt.progname, opt.min_nodes, opt.max_nodes);
		verified = false;
	}

	/* massage the numbers */
	if (opt.nodes_set && !opt.nprocs_set) {
		/* 1 proc / node default */
		opt.nprocs = opt.min_nodes;

	} else if (opt.nodes_set && opt.nprocs_set) {

		/* 
		 *  make sure # of procs >= min_nodes 
		 */
		if (opt.nprocs < opt.min_nodes) {

			info ("Warning: can't run %d processes on %d " 
			      "nodes, setting nnodes to %d", 
			      opt.nprocs, opt.min_nodes, opt.nprocs);

			opt.min_nodes = opt.nprocs;
			if (   opt.max_nodes 
			       && (opt.min_nodes > opt.max_nodes) )
				opt.max_nodes = opt.min_nodes;
		}

	} /* else if (opt.nprocs_set && !opt.nodes_set) */

	if (opt.time_limit_str) {
		opt.time_limit = time_str2mins(opt.time_limit_str);
		if (opt.time_limit < 0) {
			error("Invalid time limit specification");
			exit(1);
		}
	} else
		opt.time_limit = INFINITE;

	if ((opt.euid != (uid_t) -1) && (opt.euid != opt.uid)) 
		opt.uid = opt.euid;

	if ((opt.egid != (gid_t) -1) && (opt.egid != opt.gid)) 
		opt.gid = opt.egid;

        if ((opt.egid != (gid_t) -1) && (opt.egid != opt.gid))
	        opt.gid = opt.egid;

	if (opt.immediate) {
		char *sched_name = slurm_get_sched_type();
		if (strcmp(sched_name, "sched/wiki") == 0) {
			info("WARNING: Ignoring the -I/--immediate option "
				"(not supported by Maui)");
			opt.immediate = false;
		}
		xfree(sched_name);
	}

	return verified;
}

static uint16_t _parse_mail_type(const char *arg)
{
	uint16_t rc;

	if (strcasecmp(arg, "BEGIN") == 0)
		rc = MAIL_JOB_BEGIN;
	else if  (strcasecmp(arg, "END") == 0)
		rc = MAIL_JOB_END;
	else if (strcasecmp(arg, "FAIL") == 0)
		rc = MAIL_JOB_FAIL;
	else if (strcasecmp(arg, "ALL") == 0)
		rc = MAIL_JOB_BEGIN |  MAIL_JOB_END |  MAIL_JOB_FAIL;
	else
		rc = 0;		/* failure */

	return rc;
}
static char *_print_mail_type(const uint16_t type)
{
	if (type == 0)
		return "NONE";

	if (type == MAIL_JOB_BEGIN)
		return "BEGIN";
	if (type == MAIL_JOB_END)
		return "END";
	if (type == MAIL_JOB_FAIL)
		return "FAIL";
	if (type == (MAIL_JOB_BEGIN |  MAIL_JOB_END |  MAIL_JOB_FAIL))
		return "ALL";

	return "MULTIPLE";
}

static void
_freeF(void *data)
{
	xfree(data);
}

static List
_create_path_list(void)
{
	List l = list_create(_freeF);
	char *path = xstrdup(getenv("PATH"));
	char *c, *lc;

	if (!path) {
		error("Error in PATH environment variable");
		list_destroy(l);
		return NULL;
	}

	c = lc = path;

	while (*c != '\0') {
		if (*c == ':') {
			/* nullify and push token onto list */
			*c = '\0';
			if (lc != NULL && strlen(lc) > 0)
				list_append(l, xstrdup(lc));
			lc = ++c;
		} else
			c++;
	}

	if (strlen(lc) > 0)
		list_append(l, xstrdup(lc));

	xfree(path);

	return l;
}

static char *
_search_path(char *cmd, bool check_current_dir, int access_mode)
{
	List         l        = _create_path_list();
	ListIterator i        = NULL;
	char *path, *fullpath = NULL;

	if (  (cmd[0] == '.' || cmd[0] == '/') 
           && (access(cmd, access_mode) == 0 ) ) {
		if (cmd[0] == '.')
			xstrfmtcat(fullpath, "%s/", opt.cwd);
		xstrcat(fullpath, cmd);
		goto done;
	}

	if (check_current_dir) 
		list_prepend(l, xstrdup(opt.cwd));

	i = list_iterator_create(l);
	while ((path = list_next(i))) {
		xstrfmtcat(fullpath, "%s/%s", path, cmd);

		if (access(fullpath, access_mode) == 0)
			goto done;

		xfree(fullpath);
		fullpath = NULL;
	}
  done:
	list_destroy(l);
	return fullpath;
}


/* helper function for printing options
 * 
 * warning: returns pointer to memory allocated on the stack.
 */
static char *print_constraints()
{
	char *buf = xstrdup("");

	if (opt.mincpus > 0)
		xstrfmtcat(buf, "mincpus=%d ", opt.mincpus);

	if (opt.realmem > 0)
		xstrfmtcat(buf, "mem=%dM ", opt.realmem);

	if (opt.tmpdisk > 0)
		xstrfmtcat(buf, "tmp=%ld ", opt.tmpdisk);

	if (opt.contiguous == true)
		xstrcat(buf, "contiguous ");
 
	if (opt.nodelist != NULL)
		xstrfmtcat(buf, "nodelist=%s ", opt.nodelist);

	if (opt.exc_nodes != NULL)
		xstrfmtcat(buf, "exclude=%s ", opt.exc_nodes);

	if (opt.constraints != NULL)
		xstrfmtcat(buf, "constraints=`%s' ", opt.constraints);

	return buf;
}

static char * 
print_commandline()
{
	int i;
	char buf[256];

	buf[0] = '\0';
	for (i = 0; i < opt.script_argc; i++)
		snprintf(buf, 256,  "%s", opt.script_argv[i]);
	return xstrdup(buf);
}

static char *
print_geometry()
{
	int i;
	char buf[32], *rc = NULL;

	if ((SYSTEM_DIMENSIONS == 0)
	||  (opt.geometry[0] == (uint16_t)NO_VAL))
		return NULL;

	for (i=0; i<SYSTEM_DIMENSIONS; i++) {
		if (i > 0)
			snprintf(buf, sizeof(buf), "x%u", opt.geometry[i]);
		else
			snprintf(buf, sizeof(buf), "%u", opt.geometry[i]);
		xstrcat(rc, buf);
	}

	return rc;
}


/*
 *  Get a decimal integer from arg.
 *
 *  Returns the integer on success, exits program on failure.
 * 
 */
static int
_get_int(const char *arg, const char *what)
{
	char *p;
	long int result = strtol(arg, &p, 10);

	if ((*p != '\0') || (result < 0L)) {
		error ("Invalid numeric value \"%s\" for %s.", arg, what);
		exit(1);
	}

	if (result > INT_MAX) {
		error ("Numeric argument (%ld) to big for %s.", result, what);
	}

	return (int) result;
}


/*
 * Return an absolute path for the "filename".  If "filename" is already
 * an absolute path, it returns a copy.  Free the returned with xfree().
 */
static char *_fullpath(const char *filename)
{
	char cwd[BUFSIZ];
	char *ptr = NULL;

	if (filename[0] == '/') {
		return xstrdup(filename);
	} else {
		if (getcwd(cwd, BUFSIZ) == NULL) {
			error("could not get current working directory");
			return NULL;
		}
		ptr = xstrdup(cwd);
		xstrcat(ptr, "/");
		xstrcat(ptr, filename);
		return ptr;
	}
}

#define tf_(b) (b == true) ? "true" : "false"

static void _opt_list()
{
	char *str;

	info("defined options for program `%s'", opt.progname);
	info("--------------- ---------------------");

	info("user           : `%s'", opt.user);
	info("uid            : %ld", (long) opt.uid);
	info("gid            : %ld", (long) opt.gid);
	info("cwd            : %s", opt.cwd);
	info("nprocs         : %d %s", opt.nprocs,
		opt.nprocs_set ? "(set)" : "(default)");
	info("cpus_per_task  : %d %s", opt.cpus_per_task,
		opt.cpus_set ? "(set)" : "(default)");
	if (opt.max_nodes)
		info("nodes          : %d-%d", opt.min_nodes, opt.max_nodes);
	else {
		info("nodes          : %d %s", opt.min_nodes,
			opt.nodes_set ? "(set)" : "(default)");
	}
	info("jobid          : %u %s", opt.jobid, 
		opt.jobid_set ? "(set)" : "(default)");
	info("partition      : %s",
		opt.partition == NULL ? "default" : opt.partition);
	info("job name       : `%s'", opt.job_name);
	info("verbose        : %d", opt.verbose);
	info("immediate      : %s", tf_(opt.immediate));
	info("no-requeue     : %s", tf_(opt.no_requeue));
	if (opt.time_limit == INFINITE)
		info("time_limit     : INFINITE");
	else
		info("time_limit     : %d", opt.time_limit);
	if (opt.nice)
		info("nice           : %d", opt.nice);
	info("account        : %s", opt.account);
	info("comment        : %s", opt.comment);
	if (opt.dependency == NO_VAL)
		info("dependency     : none");
	else
		info("dependency     : %u", opt.dependency);
	str = print_constraints();
	info("constraints    : %s", str);
	xfree(str);
	if (opt.conn_type != (uint16_t) NO_VAL)
		info("conn_type      : %u", opt.conn_type);
	str = print_geometry();
	info("geometry       : %s", str);
	xfree(str);
	info("reboot         : %s", opt.reboot ? "no" : "yes");
	info("rotate         : %s", opt.no_rotate ? "yes" : "no");

	if (opt.blrtsimage)
		info("BlrtsImage     : %s", opt.blrtsimage);
	if (opt.linuximage)
		info("LinuxImage     : %s", opt.linuximage);
	if (opt.mloaderimage)
		info("MloaderImage   : %s", opt.mloaderimage);
	if (opt.ramdiskimage)
		info("RamDiskImage   : %s", opt.ramdiskimage);

	if (opt.begin) {
		char time_str[32];
		slurm_make_time_str(&opt.begin, time_str, sizeof(time_str));
		info("begin          : %s", time_str);
	}
	info("mail_type      : %s", _print_mail_type(opt.mail_type));
	info("mail_user      : %s", opt.mail_user);
	info("tasks-per-node : %d", opt.tasks_per_node);
	str = print_commandline();
	info("remote command : `%s'", str);
	xfree(str);

}

static void _usage(void)
{
 	printf(
"Usage: sbatch [-N nnodes] [-n ntasks]\n"
"              [-c ncpus] [-r n] [-p partition] [--hold] [-t minutes]\n"
"              [-D path] [--immediate] [--no-kill]\n"
"              [--input file] [--output file] [--error file]\n"
"              [--workdir=directory] [--share] [-m dist] [-J jobname]\n"
"              [--jobid=id] [--verbose] [--gid=group] [--uid=user]\n"
"              [-W sec] [--minsockets=n] [--mincores=n] [--minthreads=n]\n"
"              [--contiguous] [--mincpus=n] [--mem=MB] [--tmp=MB] [-C list]\n"
"              [--account=name] [--dependency=jobid] [--comment=name]\n"
#ifdef HAVE_BG		/* Blue gene specific options */
"              [--geometry=XxYxZ] [--conn-type=type] [--no-rotate] [ --reboot]\n"
"              [--blrts-image=path] [--linux-image=path]\n"
"              [--mloader-image=path] [--ramdisk-image=path]\n"
#endif
"              [--mail-type=type] [--mail-user=user][--nice[=value]]\n"
"              [--no-requeue] [--ntasks-per-node=n]\n"
"              [--nodefile=file] [--nodelist=hosts] [--exclude=hosts]\n"
"              executable [args...]\n");
}

static void _help(void)
{
        printf (
"Usage: sbatch [OPTIONS...] executable [args...]\n"
"\n"
"Parallel run options:\n"
"  -n, --ntasks=ntasks         number of tasks to run\n"
"  -N, --nodes=N               number of nodes on which to run (N = min[-max])\n"
"  -c, --cpus-per-task=ncpus   number of cpus required per task\n"
"      --ntasks-per-node=n     number of tasks to invoke on each node\n"
"  -i, --input=in              file for batch script's standard input\n"
"  -o, --output=out            file for batch script's standard output\n"
"  -e, --error=err             file for batch script's standard error\n"
"  -p, --partition=partition   partition requested\n"
"  -H, --hold                  submit job in held state\n"
"  -t, --time=minutes          time limit\n"
"  -D, --chdir=path            change remote current working directory\n"
"  -I, --immediate             exit if resources are not immediately available\n"
"  -k, --no-kill               do not kill job on node failure\n"
"  -s, --share                 share nodes with other jobs\n"
"  -J, --job-name=jobname      name of job\n"
"      --jobid=id              run under already allocated job\n"
"  -v, --verbose               verbose mode (multiple -v's increase verbosity)\n"
"  -q, --quiet                 quiet mode (suppress informational messages)\n"
"  -d, --dependency=jobid      defer job until specified jobid completes\n"
"  -D, --workdir=directory     set working directory for batch script\n"
"      --nice[=value]          decrease secheduling priority by value\n"
"  -U, --account=name          charge job to specified account\n"
"      --begin=time            defer job until HH:MM DD/MM/YY\n"
"      --comment=name          arbitrary comment\n"
"      --mail-type=type        notify on state change: BEGIN, END, FAIL or ALL\n"
"      --mail-user=user        who to send email notification for job state changes\n"
"      --gid=group_id          group ID to run job as (user root only)\n"
"      --uid=user_id           user ID to run job as (user root only)\n"
"      --get-user-env          used by Moab.  See srun man page.\n"
"      --no-requeue            if set, do not permit the job to be requeued\n"
"\n"
"Constraint options:\n"
"      --mincpus=n             minimum number of cpus per node\n"
"      --minsockets=n          minimum number of sockets per node\n"
"      --mincores=n            minimum number of cores per cpu\n"
"      --minthreads=n          minimum number of threads per core\n"
"      --mem=MB                minimum amount of real memory\n"
"      --tmp=MB                minimum amount of temporary disk\n"
"      --contiguous            demand a contiguous range of nodes\n"
"  -C, --constraint=list       specify a list of constraints\n"
"  -F, --nodefile=filename     request a specific list of hosts\n"
"  -w, --nodelist=hosts...     request a specific list of hosts\n"
"  -x, --exclude=hosts...      exclude a specific list of hosts\n"
"\n"
"Consumable resources related options:\n" 
"      --exclusive             allocate nodes in exclusive mode when\n" 
"                              cpu consumable resource is enabled\n"
"\n"
#ifdef HAVE_BG				/* Blue gene specific options */
"Blue Gene related options:\n"
"  -g, --geometry=XxYxZ        geometry constraints of the job\n"
"  -R, --no-rotate             disable geometry rotation\n"
"      --reboot                reboot block before starting job\n"
"      --conn-type=type        constraint on type of connection, MESH or TORUS\n"
"                              if not set, then tries to fit TORUS else MESH\n"
"      --blrts-image=path      path to blrts image for bluegene block.  Default if not set\n"
"      --linux-image=path      path to linux image for bluegene block.  Default if not set\n"
"      --mloader-image=path    path to mloader image for bluegene block.  Default if not set\n"
"      --ramdisk-image=path    path to ramdisk image for bluegene block.  Default if not set\n"
"\n"
#endif
"Help options:\n"
"  -h, --help                  show this help message\n"
"  -u, --usage                 display brief usage message\n"
"\n"
"Other options:\n"
"  -V, --version               output version information and exit\n"
"\n"
);

}
