 /*****************************************************************************\
 *  test12.2.prog.c - Simple test program for SLURM regression test12.2.
 *  Usage: test12.2.prog <exit_code> <sleep_secs> <mem_kb>
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

main (int argc, char **argv) 
{
	int exit_code, sleep_time, mem_kb;
	char *mem;

	if (argc != 4) {
		fprintf(stderr, 
			"Usage: %s <exit_code> <sleep_time> <mem_kb>\n",
			argv[0]);
		exit(1);
	}

	exit_code  = atoi(argv[1]);
	sleep_time = atoi(argv[2]);
	mem_kb     = atoi(argv[3]);

	mem = malloc(mem_kb * 1024);
	sleep(sleep_time);
	exit(exit_code);
}
