/*****************************************************************************\
 *  scontrol.c - administration tool for slurm.
 *	provides interface to read, write, update, and configurations.
 *****************************************************************************
 *  Copyright (C) 2002-2007 The Regents of the University of California.
 *  Copyright (C) 2008-2010 Lawrence Livermore National Security.
 *  Portions Copyright (C) 2008 Vijay Ramasubramanian.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
 *  CODE-OCEC-09-009. All rights reserved.
 *
 *  This file is part of SLURM, a resource management program.
 *  For details, see <https://computing.llnl.gov/linux/slurm/>.
 *  Please also read the included file: DISCLAIMER.
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

#include "scontrol.h"
#include "src/plugins/select/bluegene/wrap_rm_api.h"
#include "src/common/proc_args.h"

#define OPT_LONG_HIDE   0x102

char *command_name;
int all_flag;		/* display even hidden partitions */
int detail_flag;	/* display additional details */
int exit_code;		/* scontrol's exit code, =1 on any error at any time */
int exit_flag;		/* program to terminate if =1 */
int input_words;	/* number of words of input permitted */
int one_liner;		/* one record per line if =1 */
int quiet_flag;		/* quiet=1, verbose=-1, normal=0 */
int verbosity;		/* count of "-v" options */

static void	_create_it (int argc, char *argv[]);
static void	_delete_it (int argc, char *argv[]);
static void     _show_it (int argc, char *argv[]);
static int	_get_command (int *argc, char *argv[]);
static void     _ping_slurmctld(char *control_machine,
				char *backup_controller);
static void	_print_config (char *config_param);
static void     _print_daemons (void);
static void	_print_ping (void);
static void	_print_slurmd(char *hostlist);
static void     _print_version( void );
static int	_process_command (int argc, char *argv[]);
static void	_update_it (int argc, char *argv[]);
static int	_update_bluegene_block (int argc, char *argv[]);
static int      _update_bluegene_subbp (int argc, char *argv[]);
static int	_update_slurmctld_debug(char *val);
static void	_usage ();

int
main (int argc, char *argv[])
{
	int error_code = SLURM_SUCCESS, i, opt_char, input_field_count;
	char **input_fields;
	log_options_t opts = LOG_OPTS_STDERR_ONLY ;

	int option_index;
	static struct option long_options[] = {
		{"all",      0, 0, 'a'},
		{"details",  0, 0, 'd'},
		{"help",     0, 0, 'h'},
		{"hide",     0, 0, OPT_LONG_HIDE},
		{"oneliner", 0, 0, 'o'},
		{"quiet",    0, 0, 'Q'},
		{"usage",    0, 0, 'h'},
		{"verbose",  0, 0, 'v'},
		{"version",  0, 0, 'V'},
		{NULL,       0, 0, 0}
	};

	command_name      = argv[0];
	all_flag          = 0;
	detail_flag       = 0;
	exit_code         = 0;
	exit_flag         = 0;
	input_field_count = 0;
	quiet_flag        = 0;
	verbosity         = 0;
	log_init("scontrol", opts, SYSLOG_FACILITY_DAEMON, NULL);

	if (getenv ("SCONTROL_ALL"))
		all_flag= 1;

	while((opt_char = getopt_long(argc, argv, "adhoQvV",
			long_options, &option_index)) != -1) {
		switch (opt_char) {
		case (int)'?':
			fprintf(stderr, "Try \"scontrol --help\" for "
				"more information\n");
			exit(1);
			break;
		case (int)'a':
			all_flag = 1;
			break;
		case (int)'d':
			detail_flag = 1;
			break;
		case (int)'h':
			_usage ();
			exit(exit_code);
			break;
		case OPT_LONG_HIDE:
			all_flag = 0;
			detail_flag = 0;
			break;
		case (int)'o':
			one_liner = 1;
			break;
		case (int)'Q':
			quiet_flag = 1;
			break;
		case (int)'v':
			quiet_flag = -1;
			verbosity++;
			break;
		case (int)'V':
			_print_version();
			exit(exit_code);
			break;
		default:
			exit_code = 1;
			fprintf(stderr, "getopt error, returned %c\n",
				opt_char);
			exit(exit_code);
		}
	}

	if (verbosity) {
		opts.stderr_level += verbosity;
		log_alter(opts, SYSLOG_FACILITY_USER, NULL);
	}

	if (argc > MAX_INPUT_FIELDS)	/* bogus input, but continue anyway */
		input_words = argc;
	else
		input_words = 128;
	input_fields = (char **) xmalloc (sizeof (char *) * input_words);
	if (optind < argc) {
		for (i = optind; i < argc; i++) {
			input_fields[input_field_count++] = argv[i];
		}
	}

	if (input_field_count)
		exit_flag = 1;
	else
		error_code = _get_command (&input_field_count, input_fields);

	while (error_code == SLURM_SUCCESS) {
		error_code = _process_command (input_field_count,
					       input_fields);
		if (error_code || exit_flag)
			break;
		error_code = _get_command (&input_field_count, input_fields);
	}

	exit(exit_code);
}

static void _print_version(void)
{
	print_slurm_version();
	if (quiet_flag == -1) {
		long version = slurm_api_version();
		printf("slurm_api_version: %ld, %ld.%ld.%ld\n", version,
			SLURM_VERSION_MAJOR(version),
			SLURM_VERSION_MINOR(version),
			SLURM_VERSION_MICRO(version));
	}
}


#if !HAVE_READLINE
/*
 * Alternative to readline if readline is not available
 */
static char *_getline(const char *prompt)
{
	char buf[4096];
	char *line;
	int len;

	printf("%s", prompt);

	/* we only set this here to avoid a warning.  We throw it away
	   later. */
	line = fgets(buf, 4096, stdin);

	len = strlen(buf);
	if ((len > 0) && (buf[len-1] == '\n'))
		buf[len-1] = '\0';
	else
		len++;
	line = malloc (len * sizeof(char));
	return strncpy(line, buf, len);
}
#endif

/*
 * _get_command - get a command from the user
 * OUT argc - location to store count of arguments
 * OUT argv - location to store the argument list
 */
static int
_get_command (int *argc, char **argv)
{
	char *in_line;
	static char *last_in_line = NULL;
	int i, in_line_size;
	static int last_in_line_size = 0;

	*argc = 0;

#if HAVE_READLINE
	in_line = readline ("scontrol: ");
#else
	in_line = _getline("scontrol: ");
#endif
	if (in_line == NULL)
		return 0;
	else if (strcmp (in_line, "!!") == 0) {
		free (in_line);
		in_line = last_in_line;
		in_line_size = last_in_line_size;
	} else {
		if (last_in_line)
			free (last_in_line);
		last_in_line = in_line;
		last_in_line_size = in_line_size = strlen (in_line);
	}

#if HAVE_READLINE
	add_history(in_line);
#endif

	/* break in_line into tokens */
	for (i = 0; i < in_line_size; i++) {
		bool double_quote = false, single_quote = false;
		if (in_line[i] == '\0')
			break;
		if (isspace ((int) in_line[i]))
			continue;
		if (((*argc) + 1) > MAX_INPUT_FIELDS) {	/* bogus input line */
			exit_code = 1;
			fprintf (stderr,
				 "%s: can not process over %d words\n",
				 command_name, input_words);
			return E2BIG;
		}
		argv[(*argc)++] = &in_line[i];
		for (i++; i < in_line_size; i++) {
			if (in_line[i] == '\042') {
				double_quote = !double_quote;
				continue;
			}
			if (in_line[i] == '\047') {
				single_quote = !single_quote;
				continue;
			}
			if (in_line[i] == '\0')
				break;
			if (double_quote || single_quote)
				continue;
			if (isspace ((int) in_line[i])) {
				in_line[i] = '\0';
				break;
			}
		}
	}
	return 0;
}

/*
 * _print_config - print the specified configuration parameter and value
 * IN config_param - NULL to print all parameters and values
 */
static void
_print_config (char *config_param)
{
	int error_code;
	static slurm_ctl_conf_info_msg_t *old_slurm_ctl_conf_ptr = NULL;
	slurm_ctl_conf_info_msg_t  *slurm_ctl_conf_ptr = NULL;

	if (old_slurm_ctl_conf_ptr) {
		error_code = slurm_load_ctl_conf (
				old_slurm_ctl_conf_ptr->last_update,
			 &slurm_ctl_conf_ptr);
		if (error_code == SLURM_SUCCESS)
			slurm_free_ctl_conf(old_slurm_ctl_conf_ptr);
		else if (slurm_get_errno () == SLURM_NO_CHANGE_IN_DATA) {
			slurm_ctl_conf_ptr = old_slurm_ctl_conf_ptr;
			error_code = SLURM_SUCCESS;
			if (quiet_flag == -1) {
				printf ("slurm_load_ctl_conf no change "
					"in data\n");
			}
		}
	}
	else
		error_code = slurm_load_ctl_conf ((time_t) NULL,
						  &slurm_ctl_conf_ptr);

	if (error_code) {
		exit_code = 1;
		if (quiet_flag != 1)
			slurm_perror ("slurm_load_ctl_conf error");
	}
	else
		old_slurm_ctl_conf_ptr = slurm_ctl_conf_ptr;


	if (error_code == SLURM_SUCCESS) {
		slurm_print_ctl_conf (stdout, slurm_ctl_conf_ptr) ;
		fprintf(stdout, "\n");
	}
	if (slurm_ctl_conf_ptr)
		_ping_slurmctld (slurm_ctl_conf_ptr->control_machine,
				 slurm_ctl_conf_ptr->backup_controller);
}

/* Print slurmd status on localhost.
 * Parse hostlist in the future */
static void _print_slurmd(char *hostlist)
{
	slurmd_status_t *slurmd_status;

	if (slurm_load_slurmd_status(&slurmd_status)) {
		exit_code = 1;
		if (quiet_flag != 1)
			slurm_perror("slurm_load_slurmd_status");
	} else {
		slurm_print_slurmd_status(stdout, slurmd_status);
		slurm_free_slurmd_status(slurmd_status);
	}
}

/* Print state of controllers only */
static void
_print_ping (void)
{
	slurm_ctl_conf_info_msg_t *conf;
	char *primary, *secondary;

	slurm_conf_init(NULL);

	conf = slurm_conf_lock();
	primary = xstrdup(conf->control_machine);
	secondary = xstrdup(conf->backup_controller);
	slurm_conf_unlock();

	_ping_slurmctld (primary, secondary);

	xfree(primary);
	xfree(secondary);
}

/* Report if slurmctld daemons are responding */
static void
_ping_slurmctld(char *control_machine, char *backup_controller)
{
	static char *state[2] = { "UP", "DOWN" };
	int primary = 1, secondary = 1;
	int down_msg = 0;

	if (slurm_ping(1) == SLURM_SUCCESS)
		primary = 0;
	if (slurm_ping(2) == SLURM_SUCCESS)
		secondary = 0;
	fprintf(stdout, "Slurmctld(primary/backup) ");
	if (control_machine || backup_controller) {
		fprintf(stdout, "at ");
		if (control_machine) {
			fprintf(stdout, "%s/", control_machine);
			if (primary)
				down_msg = 1;
		} else
			fprintf(stdout, "(NULL)/");
		if (backup_controller) {
			fprintf(stdout, "%s ", backup_controller);
			if (secondary)
				down_msg = 1;
		} else
			fprintf(stdout, "(NULL) ");
	}
	fprintf(stdout, "are %s/%s\n",
		state[primary], state[secondary]);

	if (down_msg && (getuid() == 0)) {
		fprintf(stdout, "*****************************************\n");
		fprintf(stdout, "** RESTORE SLURMCTLD DAEMON TO SERVICE **\n");
		fprintf(stdout, "*****************************************\n");
	}
}

/*
 * _print_daemons - report what daemons should be running on this node
 */
static void
_print_daemons (void)
{
	slurm_ctl_conf_info_msg_t *conf;
	char me[MAX_SLURM_NAME], *b, *c, *n;
	int actld = 0, ctld = 0, d = 0;
	char daemon_list[] = "slurmctld slurmd";

	slurm_conf_init(NULL);
	conf = slurm_conf_lock();

	gethostname_short(me, MAX_SLURM_NAME);
	if ((b = conf->backup_controller)) {
		if ((strcmp(b, me) == 0) ||
		    (strcasecmp(b, "localhost") == 0))
			ctld = 1;
	}
	if ((c = conf->control_machine)) {
		actld = 1;
		if ((strcmp(c, me) == 0) ||
		    (strcasecmp(c, "localhost") == 0))
			ctld = 1;
	}
	slurm_conf_unlock();

	if ((n = slurm_conf_get_nodename(me))) {
		d = 1;
		xfree(n);
	} else if ((n = slurm_conf_get_aliased_nodename())) {
		d = 1;
		xfree(n);
	} else if ((n = slurm_conf_get_nodename("localhost"))) {
		d = 1;
		xfree(n);
	}

	strcpy(daemon_list, "");
	if (actld && ctld)
		strcat(daemon_list, "slurmctld ");
	if (actld && d)
		strcat(daemon_list, "slurmd");
	fprintf (stdout, "%s\n", daemon_list) ;
}

/*
 * _process_command - process the user's command
 * IN argc - count of arguments
 * IN argv - the arguments
 * RET 0 or errno (only for errors fatal to scontrol)
 */
static int
_process_command (int argc, char *argv[])
{
	int error_code = 0;
	char *tag = argv[0];
	int taglen = 0;

	if (argc < 1) {
		exit_code = 1;
		if (quiet_flag == -1)
			fprintf(stderr, "no input");
		return 0;
	} else if(tag)
		taglen = strlen(tag);
	else {
		if (quiet_flag == -1)
			fprintf(stderr, "input problem");
		return 0;
	}

	if (strncasecmp (tag, "abort", MAX(taglen, 5)) == 0) {
		/* require full command name */
		if (argc > 2) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		error_code = slurm_shutdown (1);
		if (error_code) {
			exit_code = 1;
			if (quiet_flag != 1)
				slurm_perror ("slurm_shutdown error");
		}
	}
	else if (strncasecmp (tag, "all", MAX(taglen, 2)) == 0)
		all_flag = 1;
	else if (strncasecmp (tag, "completing", MAX(taglen, 2)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		scontrol_print_completing();
	}
	else if (strncasecmp (tag, "create", MAX(taglen, 2)) == 0) {
		if (argc < 2) {
			exit_code = 1;
			fprintf (stderr, "too few arguments for %s keyword\n",
				 tag);
			return 0;
		}
		_create_it ((argc - 1), &argv[1]);
	}
	else if (strncasecmp (tag, "details", MAX(taglen, 1)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
			return 0;
		}
		detail_flag = 1;
	}
	else if (strncasecmp (tag, "exit", MAX(taglen, 1)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		exit_flag = 1;
	}
	else if (strncasecmp (tag, "help", MAX(taglen, 2)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		_usage ();
	}
	else if (strncasecmp (tag, "hide", MAX(taglen, 2)) == 0) {
		all_flag = 0;
		detail_flag = 0;
	}
	else if (strncasecmp (tag, "oneliner", MAX(taglen, 1)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		one_liner = 1;
	}
	else if (strncasecmp (tag, "pidinfo", MAX(taglen, 3)) == 0) {
		if (argc > 2) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		} else if (argc < 2) {
			exit_code = 1;
			fprintf (stderr,
				 "missing argument for keyword:%s\n",
				 tag);
		} else
			scontrol_pid_info ((pid_t) atol (argv[1]) );
	}
	else if (strncasecmp (tag, "ping", MAX(taglen, 3)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		_print_ping ();
	}
	else if ((strncasecmp (tag, "\\q", 2) == 0) ||
		 (strncasecmp (tag, "quiet", MAX(taglen, 4)) == 0)) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, "too many arguments for keyword:%s\n",
				 tag);
		}
		quiet_flag = 1;
	}
	else if (strncasecmp (tag, "quit", MAX(taglen, 4)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		exit_flag = 1;
	}
	else if (strncasecmp (tag, "reconfigure", MAX(taglen, 3)) == 0) {
		if (argc > 2) {
			exit_code = 1;
			fprintf (stderr, "too many arguments for keyword:%s\n",
			         tag);
		}
		error_code = slurm_reconfigure ();
		if (error_code) {
			exit_code = 1;
			if (quiet_flag != 1)
				slurm_perror ("slurm_reconfigure error");
		}
	}
	else if (strncasecmp (tag, "checkpoint", MAX(taglen, 2)) == 0) {
		if (argc > 5) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
				        "too many arguments for keyword:%s\n",
				        tag);
		}
		else if (argc < 3) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
				        "too few arguments for keyword:%s\n",
				        tag);
		}
		else {
			error_code = scontrol_checkpoint(argv[1], argv[2],
							 argc - 3, &argv[3]);
			if (error_code) {
				exit_code = 1;
				if (quiet_flag != 1)
					slurm_perror(
						"scontrol_checkpoint error");
			}
		}
	}
	else if (strncasecmp (tag, "requeue", MAX(taglen, 3)) == 0) {
		if (argc > 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too many arguments for keyword:%s\n",
					tag);
		} else if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too few arguments for keyword:%s\n",
					tag);
		} else {
			error_code = scontrol_requeue(argv[1]);
			if (error_code) {
				exit_code = 1;
				if (quiet_flag != 1)
					slurm_perror ("slurm_requeue error");
			}
		}

	}
	else if ((strncasecmp (tag, "suspend", MAX(taglen, 2)) == 0)
	||       (strncasecmp (tag, "resume", MAX(taglen, 3)) == 0)) {
		if (argc > 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too many arguments for keyword:%s\n",
					tag);
		}
		else if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too few arguments for keyword:%s\n",
					tag);
		} else {
			error_code =scontrol_suspend(argv[0], argv[1]);
			if (error_code) {
				exit_code = 1;
				if (quiet_flag != 1)
					slurm_perror ("slurm_suspend error");
			}
		}
	}
	else if (strncasecmp (tag, "setdebug", MAX(taglen, 2)) == 0) {
		if (argc > 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too many arguments for keyword:%s\n",
					tag);
		} else if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too few arguments for keyword:%s\n",
					tag);
		} else {
			int level = -1;
			char *endptr;
			char *levels[] = {
				"quiet", "fatal", "error", "info", "verbose",
				"debug", "debug2", "debug3", "debug4",
				"debug5", NULL};
			int index = 0;
			while (levels[index]) {
				if (strcasecmp(argv[1], levels[index]) == 0) {
					level = index;
					break;
				}
				index ++;
			}
			if (level == -1) {
				/* effective levels: 0 - 9 */
				level = (int)strtoul (argv[1], &endptr, 10);
				if (*endptr != '\0' || level > 9) {
					level = -1;
					exit_code = 1;
					if (quiet_flag != 1)
						fprintf(stderr, "invalid "
							"debug level: %s\n", 
							argv[1]);
				}
			}
			if (level != -1) {
				error_code = slurm_set_debug_level(level);
				if (error_code) {
					exit_code = 1;
					if (quiet_flag != 1)
						slurm_perror(
							"slurm_set_debug_level "
							"error");
				}
			}
		}
	}
	else if (strncasecmp (tag, "schedloglevel", MAX(taglen, 2)) == 0) {
		if (argc > 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too many arguments for keyword:%s\n",
					tag);
		} else if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
					"too few arguments for keyword:%s\n",
					tag);
		} else {
			int level = -1;
			char *endptr;
			char *levels[] = {
				"disable", "enable", NULL};
			int index = 0;
			while (levels[index]) {
				if (strcasecmp(argv[1], levels[index]) == 0) {
					level = index;
					break;
				}
				index ++;
			}
			if (level == -1) {
				/* effective levels: 0 - 1 */
				level = (int)strtoul (argv[1], &endptr, 10);
				if (*endptr != '\0' || level > 1) {
					level = -1;
					exit_code = 1;
					if (quiet_flag != 1)
						fprintf(stderr, "invalid schedlog "
							"level: %s\n", argv[1]);
				}
			}
			if (level != -1) {
				error_code = slurm_set_schedlog_level(level);
				if (error_code) {
					exit_code = 1;
					if (quiet_flag != 1)
						slurm_perror(
							"slurm_set_schedlog_level"
							" error");
				}
			}
		}
	}
	else if (strncasecmp (tag, "show", MAX(taglen, 3)) == 0) {
		_show_it (argc, argv);
	}
	else if (strncasecmp (tag, "takeover", MAX(taglen, 8)) == 0) {
		char *secondary = NULL;
		slurm_ctl_conf_info_msg_t  *slurm_ctl_conf_ptr = NULL;

		slurm_ctl_conf_ptr = slurm_conf_lock();
		secondary = xstrdup(slurm_ctl_conf_ptr->backup_controller);
		slurm_conf_unlock();

		if ( secondary && secondary[0] != '\0' ) {
			error_code = slurm_takeover();
			if (error_code) {
				exit_code = 1;
				if (quiet_flag != 1)
					slurm_perror("slurm_takeover error");
			}
		} else {
			fprintf(stderr, "slurm_takeover error: no backup "
				"controller defined\n");
		}
		xfree(secondary);
	}
	else if (strncasecmp (tag, "shutdown", MAX(taglen, 8)) == 0) {
		/* require full command name */
		uint16_t options = 0;
		if (argc == 2) {
			if (strcmp(argv[1], "slurmctld") &&
			    strcmp(argv[1], "controller")) {
				error_code = 1;
				exit_code = 1;
				fprintf (stderr,
					 "invalid shutdown argument:%s\n",
					 argv[1]);
			} else
				options= 2;
		} else if (argc > 2) {
			error_code = 1;
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		}
		if (error_code == 0) {
			error_code = slurm_shutdown(options);
			if (error_code) {
				exit_code = 1;
				if (quiet_flag != 1)
					slurm_perror ("slurm_shutdown error");
			}
		}
	}
	else if (strncasecmp (tag, "update", MAX(taglen, 1)) == 0) {
		if (argc < 2) {
			exit_code = 1;
			fprintf (stderr, "too few arguments for %s keyword\n",
				 tag);
			return 0;
		}
		_update_it ((argc - 1), &argv[1]);
	}
	else if (strncasecmp (tag, "delete", MAX(taglen, 1)) == 0) {
		if (argc < 2) {
			exit_code = 1;
			fprintf (stderr, "too few arguments for %s keyword\n",
				 tag);
			return 0;
		}
		_delete_it ((argc - 1), &argv[1]);
	}
	else if (strncasecmp (tag, "verbose", MAX(taglen, 4)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for %s keyword\n",
				 tag);
		}
		quiet_flag = -1;
	}
	else if (strncasecmp (tag, "version", MAX(taglen, 4)) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for %s keyword\n",
				 tag);
		}
		_print_version();
	}
	else if (strncasecmp (tag, "listpids", MAX(taglen, 1)) == 0) {
		if (argc > 3) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for keyword:%s\n",
				 tag);
		} else {
			scontrol_list_pids (argc == 1 ? NULL : argv[1],
					    argc <= 2 ? NULL : argv[2]);
		}
	}
	else if (strncasecmp (tag, "notify", MAX(taglen, 1)) == 0) {
		if (argc < 3) {
			exit_code = 1;
			fprintf (stderr,
				 "too few arguments for keyword:%s\n",
				 tag);
		} else if (scontrol_job_notify(argc-1, &argv[1])) {
			exit_code = 1;
			slurm_perror("job notify failure");
		}
	}
	else {
		exit_code = 1;
		fprintf (stderr, "invalid keyword: %s\n", tag);
	}

	return 0;
}


/*
 * _create_it - create a slurm configuration per the supplied arguments
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void
_create_it (int argc, char *argv[])
{
	/* Scan for "res" first, anywhere in the args.  When creating
	   a reservation there is a partition= option, which we don't
	   want to mistake for a requestion to create a partition. */
	int i, error_code = SLURM_SUCCESS;
	for (i=0; i<argc; i++) {
		char *tag = argv[i];
		char *val = strchr(argv[i], '=');
		int taglen;

		if (val) {
			taglen = val - argv[i];
			val++;
		} else {
			taglen = strlen(tag);
		}
		if (!strncasecmp(tag, "ReservationName", MAX(taglen, 3))) {
			error_code = scontrol_create_res(argc, argv);
			break;
		} else if (!strncasecmp(tag, "PartitionName", MAX(taglen, 3))) {
			error_code = scontrol_create_part(argc, argv);
			break;
		}
	}

	if (i >= argc) {
		exit_code = 1;
		error("Invalid creation entity: %s\n", argv[0]);
	} else if (error_code)
		exit_code = 1;
}




/*
 * _delete_it - delete the specified slurm entity
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void
_delete_it (int argc, char *argv[])
{
	char *tag = NULL, *val = NULL;
	int taglen = 0;

	if (argc != 1) {
		error("Only one option follows delete.  %d given.\n", argc);
		exit_code = 1;
		return;
	}

	tag = argv[0];
	val = strchr(argv[0], '=');
	if (val) {
		taglen = val - argv[0];
		val++;
	} else {
		error("Proper format is 'delete Partition=p'"
		      " or 'delete Reservation=r'\n");
		exit_code = 1;
		return;
	}

	/* First identify the entity type to delete */
	if (strncasecmp (tag, "PartitionName", MAX(taglen, 3)) == 0) {
		delete_part_msg_t part_msg;
		part_msg.name = val;
		if (slurm_delete_partition(&part_msg)) {
			char errmsg[64];
			snprintf(errmsg, 64, "delete_partition %s", argv[0]);
			slurm_perror(errmsg);
		}
	} else if (strncasecmp (tag, "ReservationName", MAX(taglen, 3)) == 0) {
		reservation_name_msg_t   res_msg;
		res_msg.name = val;
		if (slurm_delete_reservation(&res_msg)) {
			char errmsg[64];
			snprintf(errmsg, 64, "delete_reservation %s", argv[0]);
			slurm_perror(errmsg);
		}
	} else if (strncasecmp (tag, "BlockName", MAX(taglen, 3)) == 0) {
#ifdef HAVE_BG
		update_block_msg_t   block_msg;
		slurm_init_update_block_msg ( &block_msg );

		block_msg.bg_block_id = val;
		block_msg.state = RM_PARTITION_NAV;
		if (slurm_update_block(&block_msg)) {
			char errmsg[64];
			snprintf(errmsg, 64, "delete_block %s", argv[0]);
			slurm_perror(errmsg);
		}
#else
		exit_code = 1;
		fprintf(stderr, "This only works on a bluegene system.\n");
#endif
	} else {
		exit_code = 1;
		fprintf(stderr, "Invalid deletion entity: %s\n", argv[0]);
	}
}


/*
 * _show_it - print a description of the specified slurm entity
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void
_show_it (int argc, char *argv[])
{
	char *tag = NULL, *val = NULL;
	int taglen = 0;

	if (argc > 3) {
		exit_code = 1;
		if (quiet_flag != 1)
			fprintf(stderr,
			        "too many arguments for keyword:%s\n",
			        argv[0]);
		return;
	}
	else if (argc < 2) {
		exit_code = 1;
		if (quiet_flag != 1)
			fprintf(stderr,
			        "too few arguments for keyword:%s\n", argv[0]);
		return;
	}

	tag = argv[1];
	taglen = strlen(tag);
	val = strchr(argv[1], '=');
	if (val) {
		taglen = val - argv[1];
		val++;
	} else if (argc == 3) {
		val = argv[2];
	} else {
		val = NULL;
	}

	if (strncasecmp (tag, "blocks", MAX(taglen, 1)) == 0) {
		scontrol_print_block (val);
	} else if (strncasecmp (tag, "config", MAX(taglen, 1)) == 0) {
		_print_config (val);
	} else if (strncasecmp (tag, "daemons", MAX(taglen, 1)) == 0) {
		if (val) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr,
				        "too many arguments for keyword:%s\n",
				        argv[0]);
		}
		_print_daemons ();
	} else if (strncasecmp (tag, "jobs", MAX(taglen, 1)) == 0 ||
		 strncasecmp (tag, "jobid", MAX(taglen, 1)) == 0 ) {
		scontrol_print_job (val);
	} else if (strncasecmp (tag, "hostnames", MAX(taglen, 5)) == 0) {
		if (val)
			scontrol_print_hosts(val);
		else
			scontrol_print_hosts(getenv("SLURM_NODELIST"));
	} else if (strncasecmp (tag, "hostlist", MAX(taglen, 5)) == 0) {
		if (!val) {
			exit_code = 1;
			fprintf(stderr, "invalid encode argument\n");
			_usage();
		} else if (scontrol_encode_hostlist(val))
			exit_code = 1;
	} else if (strncasecmp (tag, "nodes", MAX(taglen, 1)) == 0) {
		scontrol_print_node_list (val);
	} else if (strncasecmp (tag, "partitions", MAX(taglen, 1)) == 0 ||
		 strncasecmp (tag, "partitionname", MAX(taglen, 1)) == 0) {
		scontrol_print_part (val);
	} else if (strncasecmp (tag, "reservations", MAX(taglen, 1)) == 0 ||
		 strncasecmp (tag, "reservationname", MAX(taglen, 1)) == 0) {
		scontrol_print_res (val);
	} else if (strncasecmp (tag, "slurmd", MAX(taglen, 2)) == 0) {
		_print_slurmd (val);
	} else if (strncasecmp (tag, "steps", MAX(taglen, 2)) == 0) {
		scontrol_print_step (val);
	} else if (strncasecmp (tag, "topology", MAX(taglen, 1)) == 0) {
		scontrol_print_topo (val);
	} else {
		exit_code = 1;
		if (quiet_flag != 1)
			fprintf (stderr,
				 "invalid entity:%s for keyword:%s \n",
				 tag, argv[0]);
	}

}



/*
 * _update_it - update the slurm configuration per the supplied arguments
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void
_update_it (int argc, char *argv[])
{
	char *val = NULL;
	int i, error_code = SLURM_SUCCESS;
	int nodetag=0, partag=0, jobtag=0;
	int blocktag=0, subtag=0, restag=0;
	int debugtag=0, steptag=0;

	/* First identify the entity to update */
	for (i=0; i<argc; i++) {
		char *tag = argv[i];
		int taglen = 0;
		val = strchr(argv[i], '=');
		if (!val)
			continue;
		taglen = val - argv[i];
		val++;

		if (!strncasecmp(tag, "NodeName", MAX(taglen, 3))) {
			nodetag=1;
		} else if (!strncasecmp(tag, "PartitionName", MAX(taglen, 3))) {
			partag=1;
		} else if (!strncasecmp(tag, "JobId", MAX(taglen, 3))) {
			jobtag=1;
		} else if (!strncasecmp(tag, "StepId", MAX(taglen, 4))) {
			steptag=1;
		} else if (!strncasecmp(tag, "BlockName", MAX(taglen, 3))) {
			blocktag=1;
		} else if (!strncasecmp(tag, "SubBPName", MAX(taglen, 3))) {
			subtag=1;
		} else if (!strncasecmp(tag, "ReservationName",
					MAX(taglen, 3))) {
			restag=1;
		} else if (!strncasecmp(tag, "SlurmctldDebug",
					MAX(taglen, 2))) {
			debugtag=1;
		}
	}

	/* The order of tests matters here.  An update job request can include
	 * partition and reservation tags, possibly before the jobid tag, but
	 * none of the other updates have a jobid tag, so check jobtag first.
	 * Likewise, check restag next, because reservations can have a
	 * partition tag.  The order of the rest doesn't matter because there
	 * aren't any other duplicate tags.  */

	if (jobtag)
		error_code = scontrol_update_job (argc, argv);
	else if (steptag)
		error_code = scontrol_update_step (argc, argv);
	else if (restag)
		error_code = scontrol_update_res (argc, argv);
	else if (nodetag)
		error_code = scontrol_update_node (argc, argv);
	else if (partag)
		error_code = scontrol_update_part (argc, argv);
	else if (blocktag)
		error_code = _update_bluegene_block (argc, argv);
	else if (subtag)
		error_code = _update_bluegene_subbp (argc, argv);
	else if (debugtag)
		error_code = _update_slurmctld_debug(val);
	else {
		exit_code = 1;
		fprintf(stderr, "No valid entity in update command\n");
		fprintf(stderr, "Input line must include \"NodeName\", ");
#ifdef HAVE_BG
		fprintf(stderr, "\"BlockName\", \"SubBPName\" "
			"(i.e. bgl000[0-3]),");
#endif
		fprintf(stderr, "\"PartitionName\", \"Reservation\", "
				"\"JobId\", or \"SlurmctldDebug\" \n");
	}

	if (error_code) {
		exit_code = 1;
		slurm_perror ("slurm_update error");
	}
}

/*
 * _update_bluegene_block - update the bluegene block per the
 *	supplied arguments
 * IN argc - count of arguments
 * IN argv - list of arguments
 * RET 0 if no slurm error, errno otherwise. parsing error prints
 *			error message and returns 0
 */
static int
_update_bluegene_block (int argc, char *argv[])
{
#ifdef HAVE_BG
	int i, update_cnt = 0;
	update_block_msg_t block_msg;

	slurm_init_update_block_msg ( &block_msg );

	for (i=0; i<argc; i++) {
		char *tag = argv[i];
		char *val = strchr(argv[i], '=');
		int taglen = 0, vallen = 0;

		if (val) {
			taglen = val - argv[i];
			val++;
			vallen = strlen(val);
		} else {
			exit_code = 1;
			error("Invalid input for BlueGene block update %s",
			      argv[i]);
			return 0;
		}

		if (strncasecmp(tag, "BlockName", MAX(taglen, 2)) == 0) {
			block_msg.bg_block_id = val;
		} else if (strncasecmp(tag, "State", MAX(taglen, 2)) == 0) {
			if (strncasecmp(val, "ERROR", MAX(vallen, 1)) == 0)
				block_msg.state = RM_PARTITION_ERROR;
			else if (strncasecmp(val, "FREE", MAX(vallen, 1)) == 0)
				block_msg.state = RM_PARTITION_FREE;
			else if (strncasecmp(val, "REMOVE",
					     MAX(vallen, 1)) == 0)
				block_msg.state = RM_PARTITION_NAV;
			else {
				exit_code = 1;
				fprintf (stderr, "Invalid input: %s\n",
					 argv[i]);
				fprintf (stderr, "Acceptable State values "
					"are FREE, ERROR, REMOVE\n");
				return 0;
			}
			update_cnt++;
		} else {
			exit_code = 1;
			error("Invalid input for BlueGene block update %s",
			      argv[i]);
			return 0;
		}
	}

	if(!block_msg.bg_block_id) {
		error("You didn't supply a block name.");
		return 0;
	} else if (block_msg.state == (uint16_t)NO_VAL) {
		error("You didn't give me a state to set %s to "
		      "(i.e. FREE, ERROR).", block_msg.nodes);
		return 0;
	}

	if (slurm_update_block(&block_msg)) {
		exit_code = 1;
		return slurm_get_errno ();
	} else
		return 0;
#else
	exit_code = 1;
	fprintf(stderr, "This only works on a bluegene system.\n");
	return 0;
#endif
}

/*
 * _update_bluegene_subbp - update the bluegene nodecards per the
 *	supplied arguments
 * IN argc - count of arguments
 * IN argv - list of arguments
 * RET 0 if no slurm error, errno otherwise. parsing error prints
 *			error message and returns 0
 */
static int
_update_bluegene_subbp (int argc, char *argv[])
{
#ifdef HAVE_BG
	int i, update_cnt = 0;
	update_block_msg_t block_msg;

	slurm_init_update_block_msg ( &block_msg );

	for (i=0; i<argc; i++) {
		char *tag = argv[i];
		char *val = strchr(argv[i], '=');
		int taglen = 0, vallen = 0;

		if (val) {
			taglen = val - argv[i];
			val++;
			vallen = strlen(val);
		} else {
			exit_code = 1;
			error("Invalid input for BlueGene SubBPName update %s",
			      argv[i]);
			return 0;
		}

		if (strncasecmp(tag, "SubBPName", MAX(taglen, 2)) == 0)
			block_msg.nodes = val;
		else if (strncasecmp(tag, "State", MAX(taglen, 2)) == 0) {
			if (strncasecmp(val, "ERROR", MAX(vallen, 1)) == 0)
				block_msg.state = RM_PARTITION_ERROR;
			else if (strncasecmp(val, "FREE", MAX(vallen, 1)) == 0)
				block_msg.state = RM_PARTITION_FREE;
			else {
				exit_code = 1;
				fprintf (stderr, "Invalid input: %s\n",
					 argv[i]);
				fprintf (stderr, "Acceptable State values "
					"are FREE and ERROR\n");
				return 0;
			}
			update_cnt++;
		} else {
			exit_code = 1;
			error("Invalid input for BlueGene SubBPName update %s",
			      argv[i]);
			return 0;
		}
	}

	if(!block_msg.nodes) {
		error("You didn't supply an ionode list.");
		return 0;
	} else if (block_msg.state == (uint16_t)NO_VAL) {
		error("You didn't give me a state to set %s to "
		      "(i.e. FREE, ERROR).", block_msg.nodes);
		return 0;
	}

	if (slurm_update_block(&block_msg)) {
		exit_code = 1;
		return slurm_get_errno ();
	} else
		return 0;
#else
	exit_code = 1;
	fprintf(stderr, "This only works on a bluegene system.\n");
	return 0;
#endif
}

/*
 * _update_slurmctld_debug - update the slurmctld debug level
 * IN  val - new value
 * RET 0 if no slurm error, errno otherwise. parsing error prints
 *			error message and returns 0
 */
static int _update_slurmctld_debug(char *val)
{
	char *endptr;
	int error_code = SLURM_SUCCESS;
	uint32_t level = (uint32_t)strtoul(val, &endptr, 10);

	if (*endptr != '\0' || level > 9) {
		error_code = 1;
		if (quiet_flag != 1)
			fprintf(stderr, "invalid debug level: %s\n",
				val);
	} else {
		error_code = slurm_set_debug_level(level);
	}

	return error_code;
}

/* _usage - show the valid scontrol commands */
void
_usage () {
	printf ("\
scontrol [<OPTION>] [<COMMAND>]                                            \n\
    Valid <OPTION> values are:                                             \n\
     -a or --all: equivalent to \"all\" command                            \n\
     -d or --detail: equivalent to \"detail\" command                      \n\
     -h or --help: equivalent to \"help\" command                          \n\
     --hide: equivalent to \"hide\" command                                \n\
     -o or --oneliner: equivalent to \"oneliner\" command                  \n\
     -Q or --quiet: equivalent to \"quiet\" command                        \n\
     -v or --verbose: equivalent to \"verbose\" command                    \n\
     -V or --version: equivalent to \"version\" command                    \n\
                                                                           \n\
  <keyword> may be omitted from the execute line and scontrol will execute \n\
  in interactive mode. It will process commands as entered until explicitly\n\
  terminated.                                                              \n\
                                                                           \n\
    Valid <COMMAND> values are:                                            \n\
     abort                    shutdown slurm controller immediately        \n\
                              generating a core file.                      \n\
     all                      display information about all partitions,    \n\
                              including hidden partitions.                 \n\
     checkpoint <CH_OP><ID>   perform a checkpoint operation on identified \n\
                              job or job step \n\
     completing               display jobs in completing state along with  \n\
                              their completing or down nodes               \n\
     create <SPECIFICATIONS>  create a new partition or reservation        \n\
     detail                   evokes additional details from the \"show\"  \n\
                              command                                      \n\
     delete <SPECIFICATIONS>  delete the specified partition or reservation\n\
                              On Dynamic layout Bluegene systems you can also\n\
                              delete blocks.                               \n\
     exit                     terminate scontrol                           \n\
     help                     print this description of use.               \n\
     hide                     do not display information about hidden      \n\
                              partitions                                   \n\
     listpids <job_id<.step>> List pids associated with the given jobid, or\n\
                              all jobs if no id is given (This will only   \n\
                              display the processes on the node which the  \n\
                              scontrol is ran on, and only for those       \n\
                              processes spawned by SLURM and their         \n\
                              descendants)                                 \n\
     notify <job_id> msg      send message to specified job                \n\
     oneliner                 report output one record per line.           \n\
     pidinfo <pid>            return slurm job information for given pid.  \n\
     ping                     print status of slurmctld daemons.           \n\
     quiet                    print no messages other than error messages. \n\
     quit                     terminate this command.                      \n\
     reconfigure              re-read configuration files.                 \n\
     requeue <job_id>         re-queue a batch job                         \n\
     setdebug <level>         set slurmctld debug level                    \n\
     schedloglevel <slevel>   set scheduler log level                      \n\
     show <ENTITY> [<ID>]     display state of identified entity, default  \n\
                              is all records.                              \n\
     shutdown <OPTS>          shutdown slurm daemons                       \n\
     takeover                 ask slurm backup controller to take over     \n\
                              (the primary controller will be stopped)     \n\
     suspend <job_id>         susend specified job                         \n\
     resume <job_id>          resume previously suspended job              \n\
     update <SPECIFICATIONS>  update job, node, partition, reservation,    \n\
                              step or bluegene block/subbp configuration   \n\
     verbose                  enable detailed logging.                     \n\
     version                  display tool version number.                 \n\
     !!                       Repeat the last command entered.             \n\
                                                                           \n\
  <ENTITY> may be \"config\", \"daemons\", \"job\", \"node\", \"partition\"\n\
       \"reservation\", \"hostlist\", \"hostnames\", \"slurmd\",           \n\
       \"topology\", or \"step\"                                           \n\
       (also for BlueGene only: \"block\" or \"subbp\").                  \n\
                                                                           \n\
  <ID> may be a configuration parameter name, job id, node name, partition \n\
       name, reservation name, job step id, or hostlist or pathname to a   \n\
       list of host names.                                                 \n\
                                                                           \n\
  <HOSTLIST> may either be a comma separated list of host names or the     \n\
       absolute pathname of a file (with leading '/' containing host names \n\
       either separated by commas or new-lines                             \n\
                                                                           \n\
  <LEVEL> may be an integer value like SlurmctldDebug in the slurm.conf    \n\
       file or the name of the most detailed errors to report (e.g. \"info\",\n\
       \"verbose\", \"debug\", \"debug2\", etc.).                          \n\
                                                                           \n\
  <SLEVEL> may be an integer value like SlurmSchedLogLevel in the          \n\
       slurm.conf file or \"enable\" or \"disable\".                       \n\
                                                                           \n\
  <OPTS> may be \"slurmctld\" to shutdown just the slurmctld daemon,       \n\
       otherwise all slurm daemons are shutdown                            \n\
                                                                           \n\
  Node names may be specified using simple range expressions,              \n\
  (e.g. \"lx[10-20]\" corresponds to lx10, lx11, lx12, ...)                \n\
  The job step id is the job id followed by a period and the step id.      \n\
                                                                           \n\
  <SPECIFICATIONS> are specified in the same format as the configuration   \n\
  file. You may wish to use the \"show\" keyword then use its output as    \n\
  input for the update keyword, editing as needed.  Bluegene blocks/subbps \n\
  are only able to be set to an error or free state.  You can also remove  \n\
  blocks by specifying 'remove' as the state.  The remove option is only   \n\
  valid on Dynamic layout systems.                                         \n\
  (Bluegene systems only)                                                  \n\
                                                                           \n\
  <CH_OP> identify checkpoint operations and may be \"able\", \"disable\", \n\
  \"enable\", \"create\", \"vacate\", \"restart\", or \"error\".           \n\
  Additional options include \"ImageDir=<dir>\", \"MaxWait=<seconds>\" and \n\
  \"StickToNodes\"   \n\
                                                                           \n\
  All commands and options are case-insensitive, although node names and   \n\
  partition names tests are case-sensitive (node names \"LX\" and \"lx\"   \n\
  are distinct).                                                       \n\n");

}
