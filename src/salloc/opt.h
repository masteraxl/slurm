/*****************************************************************************\
 *  opt.h - definitions for salloc option processing
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona1@llnl.gov>,
 *    Christopher J. Morrone <morrone2@llnl.gov>, et. al.
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

#ifndef _HAVE_OPT_H
#define _HAVE_OPT_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "src/common/macros.h" /* true and false */
#include "src/common/env.h"

#define MAX_USERNAME	9
#define DEFAULT_BELL_DELAY 10

typedef enum {BELL_NEVER, BELL_AFTER_DELAY, BELL_ALWAYS} bell_flag_t;

typedef struct salloc_options {

	char *progname;		/* argv[0] of this program or 
				 * configuration file if multi_prog */
	char user[MAX_USERNAME];/* local username		*/
	uid_t uid;		/* local uid			*/
	gid_t gid;		/* local gid			*/
	uid_t euid;		/* effective user --uid=user	*/
	gid_t egid;		/* effective group --gid=group	*/

	int  nprocs;		/* --nprocs=n,      -n n	*/
	bool nprocs_set;	/* true if nprocs explicitly set */
	int  cpus_per_task;	/* --cpus-per-task=n, -c n	*/
	bool cpus_set;		/* true if cpus_per_task explicitly set */
	int  min_nodes;		/* --nodes=n,       -N n	*/ 
	int  max_nodes;		/* --nodes=x-n,       -N x-n	*/ 
	bool nodes_set;		/* true if nodes explicitly set */
	int  time_limit;	/* --time,   -t			*/
	char *partition;	/* --partition=n,   -p n   	*/
	enum task_dist_states
		distribution;	/* --distribution=, -m dist	*/
	char *job_name;		/* --job-name=,     -J name	*/
	unsigned int jobid;	/* --jobid=jobid		*/
	unsigned int dependency;/* --dependency, -P jobid	*/
	int nice;		/* --nice			*/
	char *account;		/* --account, -U acct_name	*/
	char *comment;		/* --comment			*/

	int immediate;		/* -i, --immediate      	*/

	bool hold;		/* --hold, -H			*/
	bool no_kill;		/* --no-kill, -k		*/
	int kill_command_signal;/* --kill-command, -K           */
	bool kill_command_signal_set;
	uint16_t shared;	/* --share,   -s		*/
	int  max_wait;		/* --wait,    -W		*/
	int  quiet;
	int  verbose;

	/* constraint options */
	int mincpus;		/* --mincpus=n			*/
	int minsockets;		/* --minsockets=n		*/
	int mincores;		/* --mincores=n			*/
	int minthreads;		/* --minthreads=n		*/
	int realmem;		/* --mem=n			*/
	long tmpdisk;		/* --tmp=n			*/
	char *constraints;	/* --constraints=, -C constraint*/
	bool contiguous;	/* --contiguous			*/
	char *nodelist;		/* --nodelist=node1,node2,...	*/
	char *exc_nodes;	/* --exclude=node1,node2,... -x	*/

	uint16_t geometry[SYSTEM_DIMENSIONS]; /* --geometry, -g	*/
	bool reboot;		/* --reboot			*/
	bool no_rotate;		/* --no_rotate, -R		*/
	int16_t conn_type;	/* --conn-type 			*/
	time_t begin;		/* --begin			*/
	uint16_t mail_type;	/* --mail-type			*/
	char *mail_user;	/* --mail-user			*/
	bell_flag_t bell;       /* --bell, --no-bell            */
} opt_t;

extern opt_t opt;

/* process options:
 * 1. set defaults
 * 2. update options with env vars
 * 3. update options with commandline args
 * 4. perform some verification that options are reasonable
 */
int initialize_and_process_args(int argc, char *argv[]);

/* set options based upon commandline args */
void set_options(const int argc, char **argv);


#endif	/* _HAVE_OPT_H */
