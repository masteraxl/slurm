/****************************************************************************\
 *  opts.c - sfree command line option processing functions
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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

#include "sfree.h"

/* FUNCTIONS */
static void _help(void);
static void _print_version(void);
static void _usage(void);

/*
 * parse_command_line, fill in params data structure with data
 */
void parse_command_line(int argc, char *argv[])
{
	int opt_char;
	int option_index;
	
	static struct option long_options[] = {
		{"all",       no_argument,       0, 'a'},
		{"bgblock",  required_argument, 0, 'b'},
		{"partition", required_argument, 0, 'p'},
		{"version",   no_argument,       0, 'V'},
		{"help",      no_argument,       0, 'h'},
		{"usage",     no_argument,       0, 'u'},
		{NULL, 0, 0, 0}
	};

	while ((opt_char =
		getopt_long(argc, argv, "ab:hup:V",
			    long_options, &option_index)) != -1) {
		switch (opt_char) {
		case (int) '?':
			fprintf(stderr,
				"Try \"sfree --help\" for more information\n");
			exit(1);
			break;
		case (int) 'a':
			all_blocks = 1;
			break;
		case (int) 'V':
			_print_version();
			exit(0);
		case (int) 'b':
		case (int) 'p':
			bg_block_id = optarg;
			break;
		case (int) 'h':
		case (int) OPT_LONG_HELP:
			_help();
			exit(0);
		case (int) 'u':
		case (int) OPT_LONG_USAGE:
			_usage();
			exit(0);
		}
	}

}

void snprint_time(char *buf, size_t buf_size, time_t time)
{
	if (time == INFINITE) {
		snprintf(buf, buf_size, "UNLIMITED");
	} else {
		long days, hours, minutes, seconds;
		seconds = time % 60;
		minutes = (time / 60) % 60;
		hours = (time / 3600) % 24;
		days = time / 86400;

		if (days)
			snprintf(buf, buf_size,
				"%ld:%2.2ld:%2.2ld:%2.2ld",
				days, hours, minutes, seconds);
		else if (hours)
			snprintf(buf, buf_size,
				"%ld:%2.2ld:%2.2ld", 
				hours, minutes, seconds);
		else
			snprintf(buf, buf_size,
				"%ld:%2.2ld", minutes,seconds);
	}
}

static void _print_version(void)
{
	printf("%s %s\n", PACKAGE, SLURM_VERSION);
}

static void _usage(void)
{
	printf("Usage: sfree [-huVa] [-b]\n");
}

static void _help(void)
{
	/* We still honor -p and --partition, 
	 * but don't tell users about them here */

	printf("\
Usage: sfree [OPTIONS]\n\
  -b, --bgblock             free specific bgblock named\n\
  -a, --all                  free all bgblocks\n\
  -V, --version              output version information and exit\n\
\nHelp options:\n\
  --help                     show this help message\n\
  --usage                    display brief usage message\n");
}
