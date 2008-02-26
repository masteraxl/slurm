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
int execute_flag;       /* immediate execute=1, else = 0 */
List sacctmgr_action_list = NULL;
List sacctmgr_user_list = NULL;
List sacctmgr_association_list = NULL;
List sacctmgr_account_list = NULL;
List sacctmgr_cluster_list = NULL;

static void	_show_it (int argc, char *argv[]);
static void	_add_it (int argc, char *argv[]);
static void	_modify_it (int argc, char *argv[]);
static void	_delete_it (int argc, char *argv[]);
static int	_get_command (int *argc, char *argv[]);
static void     _print_version( void );
static int	_process_command (int argc, char *argv[]);
static void     _commit ();
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
		{"immediate",0, 0, 'i'},
		{"oneliner", 0, 0, 'o'},
		{"quiet",    0, 0, 'q'},
		{"usage",    0, 0, 'h'},
		{"verbose",  0, 0, 'v'},
		{"version",  0, 0, 'V'},
		{NULL,       0, 0, 0}
	};

	command_name      = argv[0];
	all_flag          = 0;
	execute_flag      = 0;
	exit_code         = 0;
	exit_flag         = 0;
	input_field_count = 0;
	quiet_flag        = 0;
	log_init("sacctmgr", opts, SYSLOG_FACILITY_DAEMON, NULL);

	if (getenv ("SACCTMGR_ALL"))
		all_flag= 1;

	while((opt_char = getopt_long(argc, argv, "ahioqvV",
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
		case (int)'i':
			execute_flag = 1;
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

	if(sacctmgr_action_list) {
		if(list_count(sacctmgr_action_list) > 0) {
			if(commit_check("Would you like to commit "
					"these changes?")) 
				_commit();			
			else 
				printf("Changes discarded.\n");
		}
		list_destroy(sacctmgr_action_list);
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
	} else if (strncasecmp (argv[0], "commit", 6) == 0) {
		_commit();
	} else if (strncasecmp (argv[0], "exit", 1) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, 
				 "too many arguments for keyword:%s\n", 
				 argv[0]);
		}
		if(list_count(sacctmgr_action_list) > 0) {
			char tmp_char[255];

			snprintf(tmp_char, sizeof(tmp_char),
				 "There are %d action(s) that haven't been "
				 "committed yet, would you like to commit "
				 "before exit?", 
				 list_count(sacctmgr_action_list));
			if(commit_check(tmp_char)) 
				_commit();			
			else 
				printf("Changes discarded.\n");
			list_destroy(sacctmgr_action_list);
			sacctmgr_action_list = NULL;
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
	} else if (strncasecmp (argv[0], "immediate", 9) == 0) {
		execute_flag = 1;
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
	} else if (strncasecmp (argv[0], "quit", 4) == 0) {
		if (argc > 1) {
			exit_code = 1;
			fprintf (stderr, 
				 "too many arguments for keyword:%s\n", 
				 argv[0]);
		}
		exit_flag = 1;
	} else if (strncasecmp (argv[0], "add", 3) == 0) {
		if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr, 
				        "too few arguments for keyword:%s\n", 
				        argv[0]);
		}
		_add_it((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "show", 3) == 0) {
		if (argc < 2) {
			exit_code = 1;
			if (quiet_flag != 1)
				fprintf(stderr, 
				        "too few arguments for keyword:%s\n", 
				        argv[0]);
		}
		_show_it((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "modify", 1) == 0) {
		if (argc < 2) {
			exit_code = 1;
			fprintf (stderr, "too few arguments for %s keyword\n",
				 argv[0]);
			return 0;
		}		
		_modify_it((argc - 1), &argv[1]);
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
 * _add_it - add the entity per the supplied arguments 
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void _add_it (int argc, char *argv[]) 
{
	int i, error_code = SLURM_SUCCESS;
	sacctmgr_init();

	/* First identify the entity to add */
	for (i=0; i<argc; i++) {
		if (strncasecmp (argv[i], "User", 4) == 0) {
			error_code = sacctmgr_add_user(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Account", 7) == 0) {
			error_code = sacctmgr_add_account(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Cluster", 7) == 0) {
			error_code = sacctmgr_add_cluster(
				(argc - 1), &argv[1]);
			break;
		}		
	}
	
	if (i >= argc) {
		exit_code = 1;
		fprintf(stderr, "No valid entity in add command\n");
		fprintf(stderr, "Input line must include \"Association\", ");
		fprintf(stderr, "\"UserName\", \"AccountName\", ");
		fprintf(stderr, "or \"ClusterName\"\n");
	} else if (error_code) {
		exit_code = 1;
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
	if (strncasecmp (argv[0], "User", 4) == 0) {
		error_code = sacctmgr_list_user((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "Account", 7) == 0) {
		error_code = sacctmgr_list_account((argc - 1), &argv[1]);
	} else if (strncasecmp (argv[0], "Cluster", 7) == 0) {
		error_code = sacctmgr_list_cluster((argc - 1), &argv[1]);
	} else {
		exit_code = 1;
		fprintf(stderr, "No valid entity in list command\n");
		fprintf(stderr, "Input line must include ");
		fprintf(stderr, "\"User\", \"Account\", ");
		fprintf(stderr, "or \"Cluster\"\n");
	} 
	
	if (error_code) {
		exit_code = 1;
	}
}


/* 
 * _modify_it - modify the slurm configuration per the supplied arguments 
 * IN argc - count of arguments
 * IN argv - list of arguments
 */
static void _modify_it (int argc, char *argv[]) 
{
	int i, error_code = SLURM_SUCCESS;

	sacctmgr_init();

	/* First identify the entity to modify */
	for (i=0; i<argc; i++) {
		if (strncasecmp (argv[i], "User", 4) == 0) {
			error_code = sacctmgr_modify_user((argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Account", 7) == 0) {
			error_code = sacctmgr_modify_account(
				(argc - 1), &argv[1]);
			break;
		} else if (strncasecmp (argv[i], "Cluster", 7) == 0) {
			error_code = sacctmgr_modify_cluster(
				(argc - 1), &argv[1]);
			break;
		}		
	}
	
	if (i >= argc) {
		exit_code = 1;
		fprintf(stderr, "No valid entity in modify command\n");
		fprintf(stderr, "Input line must include ");
		fprintf(stderr, "\"User\", \"Account\", ");
		fprintf(stderr, "or \"Cluster\"\n");
	} else if (error_code) {
		exit_code = 1;
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
		if (strncasecmp (argv[i], "User", 4) == 0) {
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
		fprintf(stderr, "Input line must include ");
		fprintf(stderr, "\"User\", \"Account\", ");
		fprintf(stderr, "or \"Cluster\"\n");
	} else if (error_code) {
		exit_code = 1;
	}
}

static void _commit ()
{
	int rc = SLURM_SUCCESS;
	ListIterator itr = NULL;
	sacctmgr_action_t *action = NULL;

	if(!sacctmgr_action_list) {
		error("No actions to commit");
		return;
	}
	
	itr = list_iterator_create(sacctmgr_action_list);
	while((action = list_next(itr))) {
		/* if(rc != SLURM_SUCCESS) { */
/* 			error("_commit: last command returned error."); */
/* 			break; */
/* 		} */
		switch(action->type) {
		case SACCTMGR_ACTION_NOTSET:
			error("This action does not have a type.");
			break;
		case SACCTMGR_USER_CREATE:
			rc = account_storage_g_add_users(action->list);		
			break;
		case SACCTMGR_ACCOUNT_CREATE:
			rc = account_storage_g_add_accounts(action->list);
			break;
		case SACCTMGR_CLUSTER_CREATE:
			rc = account_storage_g_add_clusters(action->list);
			break;
		case SACCTMGR_ASSOCIATION_CREATE:
			rc = account_storage_g_add_associations(action->list);
			break;
		case SACCTMGR_USER_MODIFY:
			rc = account_storage_g_modify_users(action->cond,
							    action->rec);
			break;
		case SACCTMGR_USER_DELETE:
			rc = account_storage_g_remove_users(action->cond);
			break;
		case SACCTMGR_ACCOUNT_MODIFY:
			rc = account_storage_g_modify_accounts(action->cond,
							       action->rec);
			break;
		case SACCTMGR_ACCOUNT_DELETE:
			rc = account_storage_g_remove_accounts(action->cond);
			break;
		case SACCTMGR_CLUSTER_MODIFY:
			rc = account_storage_g_modify_clusters(action->cond,
							       action->rec);
			break;
		case SACCTMGR_CLUSTER_DELETE:
			rc = account_storage_g_remove_clusters(action->cond);
			break;
		case SACCTMGR_ASSOCIATION_MODIFY:
			rc = account_storage_g_modify_associations(action->cond,
								   action->rec);
			break;
		case SACCTMGR_ASSOCIATION_DELETE:
			rc = account_storage_g_remove_associations(
				action->cond);
			break;
		case SACCTMGR_ADMIN_MODIFY:
			rc = account_storage_g_modify_user_admin_level(
				action->cond);
			break;
		case SACCTMGR_COORD_CREATE:
			rc = account_storage_g_add_coord(action->rec,
							 action->cond);
			break;
		case SACCTMGR_COORD_DELETE:
			rc = account_storage_g_remove_coord(action->rec,
							    action->cond);
			break;	
		default:
			error("unknown action %d", action->type);
			break;
		}
	}
	list_iterator_destroy(itr);
}

/* _usage - show the valid sacctmgr commands */
void _usage () {
	printf ("\
sacctmgr [<OPTION>] [<COMMAND>]                                            \n\
    Valid <OPTION> values are:                                             \n\
     -a or --all: equivalent to \"all\" command                            \n\
     -h or --help: equivalent to \"help\" command                          \n\
     --hide: equivalent to \"hide\" command                                \n\
     -i or --immediate: equivalent to \"immediate\" command                \n\
     -o or --oneliner: equivalent to \"oneliner\" command                  \n\
     -q or --quiet: equivalent to \"quiet\" command                        \n\
     -s or --associations: equivalent to \"associations\" command          \n\
     -v or --verbose: equivalent to \"verbose\" command                    \n\
     -V or --version: equivalent to \"version\" command                    \n\
                                                                           \n\
  <keyword> may be omitted from the execute line and sacctmgr will execute \n\
  in interactive mode. It will process commands as entered until explicitly\n\
  terminated.                                                              \n\
                                                                           \n\
    Valid <COMMAND> values are:                                            \n\
     all                      display information about all entities,      \n\
                              including hidden/deleted ones.               \n\
     add <ENTITY> <SPECS>     add entity                                   \n\
     associations             when using show/list will list the           \n\
                              associations asspciated with the entity.     \n\
     commit                   commit changes done with create, modify,     \n\
                              or delete                                    \n\
     delete <ENTITY> <SPECS>  delete the specified entity(s)               \n\
     exit                     terminate sacctmgr                           \n\
     help                     print this description of use.               \n\
     hide                     do not display information about             \n\
                              hidden/deleted entities.                     \n\
     immediate                commit changes immediately                   \n\
     list <ENTITY> [<SPECS>]  display info of identified entity, default   \n\
                              is display all.                              \n\
     modify <ENTITY> <SPECS>  modify entity                                \n\
     oneliner                 report output one record per line.           \n\
     quiet                    print no messages other than error messages. \n\
     quit                     terminate this command.                      \n\
     show                     same as list                                 \n\
     verbose                  enable detailed logging.                     \n\
     version                  display tool version number.                 \n\
     !!                       Repeat the last command entered.             \n\
                                                                           \n\
  <ENTITY> may be \"user\", \"cluster\", \"account\", or \"association\".  \n\
                                                                           \n\
  <SPECS> are different for each command entity pair.                      \n\
       list user          - Names=, DefaultAccounts=, ExpediteLevel=,      \n\
                            and AdminLevel=                                \n\
       add user           - Names=, DefaultAccount=, ExpediteLevel=,       \n\
                            and AdminLevel=                                \n\
       modify user        - Names=, DefaultAccounts=, ExpediteLevel=,      \n\
                            and AdminLevel=                                \n\
       delete user        - Names=, DefaultAccounts=, ExpediteLevel=,      \n\
                            and AdminLevel=                                \n\
                                                                           \n\
       list account       - Names=, Descriptions=, ExpediteLevel=,         \n\
                            and Organizations=                             \n\
       add account        - Names=, Descriptions=, ExpediteLevel=,         \n\
                            and Organizations=                             \n\
       modify account     - Names=, Descriptions=, ExpediteLevel=,         \n\
                            and Organizations=                             \n\
       delete account     - Names=, Descriptions=, ExpediteLevel=,         \n\
                            and Organizations=                             \n\
                                                                           \n\
       list cluster       - Names=                                         \n\
       add cluster        - Name=, and InterfaceNode=                      \n\
       modify cluster     - Name=, and InterfaceNode=                      \n\
       delete cluster     - Names=                                         \n\
                                                                           \n\
                                                                           \n\
  All commands entitys, and options are case-insensitive.               \n\n");
	
}

