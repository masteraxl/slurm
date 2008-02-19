/*****************************************************************************\
 *  sacctmgr.c - administration tool for slurm's accounting. 
 *	         provides interface to read, write, update, and configure
 *               accounting.
 *****************************************************************************
 *  Copyright (C) 2002-2008 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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

#include "sacctmgr.h"

#define OPT_LONG_HIDE   0x102
#define BUFFER_SIZE 4096

FILE *history_fp = NULL;
int history_pos = 0;
int history_count = 0;
char **history = NULL;

char *command_name;
int all_flag;		/* display even hidden partitions */
int exit_code;		/* scontrol's exit code, =1 on any error at any time */
int exit_flag;		/* program to terminate if =1 */
int input_words;	/* number of words of input permitted */
int one_liner;		/* one record per line if =1 */
int quiet_flag;		/* quiet=1, verbose=-1, normal=0 */

static void	_show_it (int argc, char *argv[]);
static void	_create_it (int argc, char *argv[]);
static void	_update_it (int argc, char *argv[]);
static void	_delete_it (int argc, char *argv[]);
static int	_get_command (int *argc, char *argv[]);
static void     _print_version( void );
static int	_process_command (int argc, char *argv[]);
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
		{"help",     0, 0, 'h'},
		{"hide",     0, 0, OPT_LONG_HIDE},
		{"oneliner", 0, 0, 'o'},
		{"quiet",    0, 0, 'q'},
		{"usage",    0, 0, 'h'},
		{"verbose",  0, 0, 'v'},
		{"version",  0, 0, 'V'},
		{NULL,       0, 0, 0}
	};

	command_name      = argv[0];
	all_flag          = 0;
	exit_code         = 0;
	exit_flag         = 0;
	input_field_count = 0;
	quiet_flag        = 0;
	log_init("sacctmgr", opts, SYSLOG_FACILITY_DAEMON, NULL);

	if (getenv ("SACCTMGR_ALL"))
		all_flag= 1;

	while((opt_char = getopt_long(argc, argv, "ahoqvV",
			long_options, &option_index)) != -1) {
		switch (opt_char) {
		case (int)'?':
			fprintf(stderr, "Try \"sacctmgr --help\" for more information\n");
			exit(1);
			break;
		case (int)'a':
			all_flag = 1;
			break;
		case (int)'h':
			_usage ();
			exit(exit_code);
			break;
		case OPT_LONG_HIDE:
			all_flag = 0;
			break;
		case (int)'o':
			one_liner = 1;
			break;
		case (int)'q':
			quiet_flag = 1;
			break;
		case (int)'v':
			quiet_flag = -1;
			break;
		case (int)'V':
			_print_version();
			exit(exit_code);
			break;
		default:
			exit_code = 1;
			fprintf(stderr, "getopt error, returned %c\n", opt_char);
			exit(exit_code);
		}
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

#if !HAVE_READLINE
/*
 * Alternative to readline if readline is not available
 */
static char *
getline(const char *prompt)
{
	char buf[4096];
	char *line;
	int len;
	printf("%s", prompt);

	fgets(buf, 4096, stdin);
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
	in_line = readline ("sacctmgr: ");
#else
	in_line = getline("sacctmgr: ");
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


static void _print_version(void)
{
	printf("%s %s\n", PACKAGE, SLURM_VERSION);
	if (quiet_flag == -1) {
		long version = slurm_api_version();
		printf("slurm_api_version: %ld, %ld.%ld.%ld\n", version,
			SLURM_VERSION_MAJOR(version), 
			SLURM_VERSION_MINOR(version),
			SLURM_VERSION_MICRO(version));
	}
}

/*
 * _process_command - process the user's command
 * IN argc - count of arguments
 * IN argv - the arguments
 * RET 0 or errno (only for errors fatal to sacctmgr)
 */
static int
_process_command (int argc, char *argv[]) 
{
	if (argc < 1) {
		exit_code = 1;
		if (quiet_flag == -1)
			fprintf(stderr, "no input");
	} else if (strncasecmp (argv[0], "all", 3) == 0) {
		all_flag = 1;
	} else if (strncasecmp (argv[0], "exit", 1) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, 
				 "too many arguments for keyword:%s\n", 
				 argv[0]);
		}
		exit_flag = 1;
	} else if (strncasecmp (argv[0], "help", 2) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, 
				 "too many arguments for keyword:%s\n",
				 argv[0]);
		}
		_usage ();
	} else if (strncasecmp (argv[0], "hide", 2) == 0) {
		all_flag = 0;
	} else if (strncasecmp (argv[0], "oneliner", 1) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, 
				 "too many arguments for keyword:%s\n",
				 argv[0]);
		}
		one_liner = 1;
	} else if (strncasecmp (argv[0], "quiet", 4) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, "too many arguments for keyword:%s\n",
				 argv[0]);
		}
		quiet_flag = 1;
	}
	else if (strncasecmp (argv[0], "quit", 4) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, 
				 "too many arguments for keyword:%s\n", 
				 argv[0]);
		}
		exit_flag = 1;
	} else if (strncasecmp (argv[0], "create", 6) == 0) {
		if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr, 
				        "too few arguments for keyword:%s\n", 
				        argv[0]);
		}
		_create_it((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "show", 3) == 0) {
		if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr, 
				        "too few arguments for keyword:%s\n", 
				        argv[0]);
		}
		_show_it((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "update", 1) == 0) {
		if (argc < 2) {
			exit_code = 1;
			fprintf (stderr, "too few arguments for %s keyword\n",
				 argv[0]);
			return 0;
		}		
		_update_it((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "delete", 3) == 0) {
		if (argc < 2) {
			exit_code = 1;
			fprintf (stderr, "too few arguments for %s keyword\n",
				 argv[0]);
			return 0;
		}
		_delete_it((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "verbose", 4) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for %s keyword\n",
				 argv[0]);
		}		
		quiet_flag = -1;
	} else if (strncasecmp (argv[0], "version", 4) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr,
				 "too many arguments for %s keyword\n",
				 argv[0]);
		}		
		_print_version();
	} else {
		exit_code = 1;
		fprintf (stderr, "invalid keyword: %s\n", argv[0]);
	}

	return 0;
}


/* 
 * _create_it - create the slurm configuration per the supplied arguments 
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void _create_it (int argc, char *argv[]) 
{
	int i, error_code = SLURM_SUCCESS;

	/* First identify the entity to create */
	for (i=0; i<argc; i++) {
		if (strncasecmp (argv[i], "Association", 11) == 0) {
			error_code = sacctmgr_create_association(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "User", 4) == 0) {
			error_code = sacctmgr_create_user(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Account", 7) == 0) {
			error_code = sacctmgr_create_account(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Cluster", 7) == 0) {
			error_code = sacctmgr_create_cluster(
				(argc - 1), &argv[1]);
			break;
		}		
	}
	
	if (i >= argc) {
		exit_code = 1;
		fprintf(stderr, "No valid entity in create command\n");
		fprintf(stderr, "Input line must include \"Association\", ");
		fprintf(stderr, "\"UserName\", \"AccountName\", ");
		fprintf(stderr, "or \"ClusterName\"\n");
	} else if (error_code) {
		exit_code = 1;
		slurm_perror ("sacctmgr_create error");
	}
}

/* 
 * _show_it - list the slurm configuration per the supplied arguments 
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void _show_it (int argc, char *argv[]) 
{
	int error_code = SLURM_SUCCESS;
		
	/* First identify the entity to list */
	if (strncasecmp (argv[0], "Association", 11) == 0) {
		error_code = sacctmgr_list_association((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "User", 4) == 0) {
		error_code = sacctmgr_list_user((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "Account", 7) == 0) {
		error_code = sacctmgr_list_account((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "Cluster", 7) == 0) {
		error_code = sacctmgr_list_cluster((argc - 1), &argv[1]);
	} else {
		exit_code = 1;
		fprintf(stderr, "No valid entity in list command\n");
		fprintf(stderr, "Input line must include \"Association\", ");
		fprintf(stderr, "\"UserName\", \"AccountName\", ");
		fprintf(stderr, "or \"ClusterName\"\n");
	} 
	
	if (error_code) {
		exit_code = 1;
		slurm_perror ("sacctmgr_list error");
	}
}


/* 
 * _update_it - update the slurm configuration per the supplied arguments 
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void _update_it (int argc, char *argv[]) 
{
	int i, error_code = SLURM_SUCCESS;

	/* First identify the entity to update */
	for (i=0; i<argc; i++) {
		if (strncasecmp (argv[i], "Association", 11) == 0) {
			error_code = sacctmgr_update_association(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "User", 4) == 0) {
			error_code = sacctmgr_update_user((argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Account", 7) == 0) {
			error_code = sacctmgr_update_account(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Cluster", 7) == 0) {
			error_code = sacctmgr_update_cluster(
				(argc - 1), &argv[1]);
			break;
		}		
	}
	
	if (i >= argc) {
		exit_code = 1;
		fprintf(stderr, "No valid entity in update command\n");
		fprintf(stderr, "Input line must include \"Association\", ");
		fprintf(stderr, "\"UserName\", \"AccountName\", ");
		fprintf(stderr, "or \"ClusterName\"\n");
	} else if (error_code) {
		exit_code = 1;
		slurm_perror ("sacctmgr_update error");
	}
}

/* 
 * _delete_it - delete the slurm configuration per the supplied arguments 
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void _delete_it (int argc, char *argv[]) 
{
	int i, error_code = SLURM_SUCCESS;

	/* First identify the entity to delete */
	for (i=0; i<argc; i++) {
		if (strncasecmp (argv[i], "Association", 11) == 0) {
			error_code = sacctmgr_delete_association(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "User", 4) == 0) {
			error_code = sacctmgr_delete_user((argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Account", 7) == 0) {
			error_code = sacctmgr_delete_account(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Cluster", 7) == 0) {
			error_code = sacctmgr_delete_cluster(
				(argc - 1), &argv[1]);
			break;
		}		
	}
	
	if (i >= argc) {
		exit_code = 1;
		fprintf(stderr, "No valid entity in delete command\n");
		fprintf(stderr, "Input line must include \"Association\", ");
		fprintf(stderr, "\"UserName\", \"AccountName\", ");
		fprintf(stderr, "or \"ClusterName\"\n");
	} else if (error_code) {
		exit_code = 1;
		slurm_perror ("sacctmgr_delete error");
	}
}

/* _usage - show the valid sacctmgr commands */
void _usage () {
	printf ("\
sacctmgr [<OPTION>] [<COMMAND>]                                            \n\
    Valid <OPTION> values are:                                             \n\
     -a or --all: equivalent to \"all\" command                            \n\
     -h or --help: equivalent to \"help\" command                          \n\
     --hide: equivalent to \"hide\" command                                \n\
     -o or --oneliner: equivalent to \"oneliner\" command                  \n\
     -q or --quiet: equivalent to \"quiet\" command                        \n\
     -v or --verbose: equivalent to \"verbose\" command                    \n\
     -V or --version: equivalent to \"version\" command                    \n\
                                                                           \n\
  <keyword> may be omitted from the execute line and sacctmgr will execute \n\
  in interactive mode. It will process commands as entered until explicitly\n\
  terminated.                                                              \n\
                                                                           \n\
    Valid <COMMAND> values are:                                            \n\
     abort                    shutdown slurm controller immediately        \n\
                              generating a core file.                      \n\
     all                      display information about all partitions,    \n\
                              including hidden partitions.                 \n\
     checkpoint <CH_OP><step> perform a checkpoint operation on identified \n\
                              job step \n\
     completing               display jobs in completing state along with  \n\
                              their completing or down nodes               \n\
     delete <SPECIFICATIONS>  delete the specified partition, kill its jobs\n\
     exit                     terminate sacctmgr                           \n\
     help                     print this description of use.               \n\
     hide                     do not display information about hidden partitions.\n\
     listpids <job_id<.step>> List pids associated with the given jobid, or\n\
                              all jobs if no id is given (This will only   \n\
                              display the processes on the node which the  \n\
                              sacctmgr is ran on, and only for those       \n\
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
     setdebug <LEVEL>         reset slurmctld debug level                  \n\
     show <ENTITY> [<ID>]     display state of identified entity, default  \n\
                              is all records.                              \n\
     shutdown                 shutdown slurm controller.                   \n\
     suspend <job_id>         susend specified job                         \n\
     resume <job_id>          resume previously suspended job              \n\
     setdebug <level>         set slurmctld debug level                    \n\
     update <SPECIFICATIONS>  update job, node, partition, or bluegene     \n\
                              block/subbp configuration                    \n\
     verbose                  enable detailed logging.                     \n\
     version                  display tool version number.                 \n\
     !!                       Repeat the last command entered.             \n\
                                                                           \n\
  <ENTITY> may be \"config\", \"daemons\", \"job\", \"node\", \"partition\"\n\
           \"hostlist\", \"hostnames\", \"slurmd\",                        \n\
           (for BlueGene only: \"block\", \"subbp\" or \"step\").          \n\
                                                                           \n\
  <ID> may be a configuration parameter name, job id, node name, partition \n\
       name, job step id, or hostlist or pathname to a list of host names. \n\
                                                                           \n\
  <HOSTLIST> may either be a comma separated list of host names or the     \n\
       absolute pathname of a file (with leading '/' containing host names \n\
       either separated by commas or new-lines                             \n\
                                                                           \n\
  <LEVEL> may be an integer value like SlurmctldDebug in the slurm.conf    \n\
       file or the name of the most detailed errors to report (e.g. \"info\",\n\
       \"verbose\", \"debug\", \"debug2\", etc.).                          \n\
                                                                           \n\
  Node names may be specified using simple range expressions,              \n\
  (e.g. \"lx[10-20]\" corresponsds to lx10, lx11, lx12, ...)               \n\
  The job step id is the job id followed by a period and the step id.      \n\
                                                                           \n\
  <SPECIFICATIONS> are specified in the same format as the configuration   \n\
  file. You may wish to use the \"show\" keyword then use its output as    \n\
  input for the update keyword, editing as needed.  Bluegene blocks/subbps \n\
  are only able to be set to an error or free state.                       \n\
  (Bluegene systems only)                                                  \n\
                                                                           \n\
  <CH_OP> identify checkpoint operations and may be \"able\", \"disable\", \n\
  \"enable\", \"create\", \"vacate\", \"restart\", or \"error\".           \n\
                                                                           \n\
  All commands and options are case-insensitive, although node names and   \n\
  partition names tests are case-sensitive (node names \"LX\" and \"lx\"   \n\
  are distinct).                                                       \n\n");

}

