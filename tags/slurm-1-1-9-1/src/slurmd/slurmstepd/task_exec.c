/*****************************************************************************\
 *  task_exec.c - execute program according to task rank
 *
 *  NOTE: This code could be moved to srun if desired. That would mean the 
 *  logic would be executed once per job instead of once per task. This would
 *  require substantial modifications to the srun, slurmd, slurmstepd, and
 *  communications logic; so we'll stick with the simple solution for now. 
 *****************************************************************************
 *  Produced at National University of Defense Technology (China)
 *  Written by Hongjia Cao <hjcao@nudt.edu.cn>
 *  and
 *  Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>.
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

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "src/common/log.h"
#include "src/common/xassert.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"

#define BUF_SIZE 256

/*
 * Test if the specified rank is included in the supplied task range 
 * IN rank    - this task's rank
 * IN spec    - a line from the configuration file
 * OUT offset - the task's offset within rank range of the configuration file
 * RET 1 if within range, 0 otherwise
 */
static int
_in_range(int rank, char* spec, int *offset)
{
	char* range;
	char* p;
	char* upper;
	int high_num, low_num, passed = 0;

	xassert(offset);
	
	if (spec[0] == '*' && spec[1] == '\0') {
		*offset = rank;
		return 1;
	}

	for (range = strtok (spec, ","); range != NULL; 
			range = strtok (NULL, ",")) {
		p = range;
		while (*p != '\0' && isdigit (*p))
			p ++;
		if (*p == '\0') { /* single rank */
			if (rank == atoi (range)) {
				*offset = passed;
				return 1;
			}
			passed ++;

		} else if (*p == '-') { /* lower-upper */
			upper = ++ p;
			while (isdigit (*p))
				p ++;
			if (*p != '\0') {
				error ("Invalid task range specification (%s) "
					"ignored.", range);
				continue;
			};
			low_num  = atoi (range);
			high_num = atoi (upper);
			if ((rank >= low_num) && (rank <= high_num)) {
				*offset = passed + (rank - low_num);
				return 1;
			} else
				passed += (1 + high_num - low_num);

		} else {
			error ("Invalid task range specification (%s) ignored.",
				range);
		}
	}
	return 0;
}

/* substitute "%t" or "%o" in argument with task number or range offset */
static void
_sub_expression(char *args_spec, int task_rank, int task_offset)
{
	char tmp[BUF_SIZE];

	if (args_spec[1] == 't') {
		/* task rank */
		strcpy(tmp, &args_spec[2]);
		sprintf(args_spec, "%d%s", task_rank, tmp);

	} else if (args_spec[1] == 'o') {
		/* task offset */
		strcpy(tmp, &args_spec[2]);
		sprintf(args_spec, "%d%s", task_offset, tmp);
	}
}

/* Given a program name, translate it to a fully qualified pathname
 * as needed based upon the PATH environment variable */
static char *
_build_path(char* fname, char **prog_env)
{
	int i;
	char *path_env = NULL, *dir;
	static char file_name[256], file_path[256];	/* return values */
	struct stat buf;

	/* make copy of file name (end at white space) */
	snprintf(file_name, sizeof(file_name), "%s", fname);
	for (i=0; i<sizeof(file_name); i++) {
		if (file_name[i] == '\0')
			break;
		if (!isspace(file_name[i]))
			continue;
		file_name[i] = '\0';
		break;
	}

	/* check if already absolute path */
	if (file_name[0] == '/')
		return file_name;

	/* search for the file using PATH environment variable */
	for (i=0; ; i++) {
		if (prog_env[i] == NULL)
			return file_name;
		if (strncmp(prog_env[i], "PATH=", 5))
			continue;
		path_env = xstrdup(&prog_env[i][5]);
		break;
	}

	dir = strtok(path_env, ":");
	while (dir) {
		snprintf(file_path, sizeof(file_path), "%s/%s", dir, file_name);
		if (stat(file_path, &buf) == 0)
			break;
		dir = strtok(NULL, ":");
	}
	if (dir == NULL)	/* not found */
		snprintf(file_path, sizeof(file_path), "%s", file_name);

	xfree(path_env);
	return file_path;
}

extern int
task_exec(char *config_data, char **prog_env, int task_rank)
{
	char *line;
	int line_num = 0;
	int task_offset;
	char* p, *s, *ptrptr;
	char* rank_spec, *prog_spec = NULL, *args_spec;
	int prog_argc = 0;
	char* prog_argv[(BUF_SIZE - 4)/ 2];

	if (task_rank < 0)
		return -1;

	line = strtok_r(config_data, "\n", &ptrptr);
	while (line) {
		if (line_num > 0)
			line = strtok_r(NULL, "\n", &ptrptr);
		if (line == NULL) {
			error("Could not identify executable program for this task");
			return -1;
		}
		line_num ++;
		if (strlen (line) >= (BUF_SIZE - 1)) {
			error ("Line %d of configuration file too long", 
				line_num);
			return -1;
		}
		
		p = line;
		while (*p != '\0' && isspace (*p)) /* remove leading spaces */
			p ++;
		
		if (*p == '#') /* only whole-line comments handled */
			continue;

		if (*p == '\0') /* blank line ignored */
			continue;
		
		rank_spec = p;

		while (*p != '\0' && !isspace (*p))
			p ++;
		if (*p == '\0') {
			error("Invalid configuration line: %s", line);
			return -1;
		}
		*p ++ = '\0';

		if (!_in_range (task_rank, rank_spec, &task_offset))
			continue;

		while(*p != '\0' && isspace (*p))
			p ++;
		prog_spec = _build_path(p, prog_env);

		if (prog_spec[0] == '\0') {
			error("Program for task rank %d not specified.", 
				task_rank);
			return -1;
		}
		
		prog_argv[0] = prog_spec; 
		prog_argv[1] = NULL;
		prog_argc = 1;

		while (*p != '\0' && !isspace (*p))
			p ++;

		/* If *p is already \0, then we are at the end of line;
		   therre are no command line parameters. */
		if (*p != '\0')
			*p++ = '\0';

		while (*p != '\0' && isspace (*p))
			p ++;

		args_spec = p;
		while (*args_spec != '\0') { 
			/* Only simple quote and escape supported */
			prog_argv[prog_argc ++] = args_spec;
		CONT:	while (*args_spec != '\0' && *args_spec != '\\'
			&&     *args_spec != '%'
			&&     *args_spec != '\'' && !isspace (*args_spec)) {
			        args_spec ++;
		        }
			if (*args_spec == '\0') {
				/* the last argument */
				break;

			} else if (*args_spec == '%') {
				_sub_expression(args_spec, task_rank, 
					task_offset);
				args_spec ++;
				goto CONT;

			} else if (*args_spec == '\\') {
				/* escape, just remove the backslash */
				s = args_spec ++;
				p = args_spec;
				do {
					*s ++ = *p;
				} while (*p ++ != '\0');
				goto CONT;
				
			} else if (*args_spec == '\'') {
				/* single quote, 
				 * preserve all characters quoted. */
				p = args_spec + 1;
				while (*p != '\0' && *p != '\'') {
					/* remove quote */
					*args_spec ++ = *p ++;
				}
				if (*p == '\0') {
					/* closing quote not found */
					error("Program arguments specification"
						" format invalid: %s.", 
						prog_argv[prog_argc -1]);
					return -1;
				}
				p ++; /* skip closing quote */
				s = args_spec;
				do {
					*s ++ = *p;
				} while (*p ++ != '\0');
				goto CONT;
				
			} else {
				/* space */
				*args_spec ++ = '\0';
				while (*args_spec != '\0' 
				&& isspace (*args_spec))
					args_spec ++;
			}
		}

		prog_argv[prog_argc] = NULL;
		if (execve (prog_spec, prog_argv, prog_env) == -1) {
			error("Error executing program \"%s\": %m", prog_spec);
			return -1;
		}
	}

	error("Program for task rank %d not specified.", task_rank);
	return -1;
}