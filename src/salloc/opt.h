/*****************************************************************************\
 *  opt.h - definitions for srun option processing
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona1@llnl.gov>,
 *    Christopher J. Morrone <morrone2@llnl.gov>, et. al.
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
//#include "src/common/mpi.h"

#define MAX_USERNAME	9


/* global variables relating to user options */
extern char **remote_argv;
extern int remote_argc;
extern int _verbose;

#define format_task_dist_states(t) (t == SLURM_DIST_BLOCK) ? "block" :   \
		                 (t == SLURM_DIST_CYCLIC) ? "cyclic" : \
			         (t == SLURM_DIST_ARBITRARY) ? "arbitrary" : \
			         "unknown"

typedef struct salloc_options {

	char *progname;		/* argv[0] of this program or 
				 * configuration file if multi_prog */
	bool multi_prog;	/* multiple programs to execute */
	char user[MAX_USERNAME];/* local username		*/
	uid_t uid;		/* local uid			*/
	gid_t gid;		/* local gid			*/
	uid_t euid;		/* effective user --uid=user	*/
	gid_t egid;		/* effective group --gid=group	*/
	char *cwd;		/* current working directory	*/

	int  nprocs;		/* --nprocs=n,      -n n	*/
	bool nprocs_set;	/* true if nprocs explicitly set */
	int  cpus_per_task;	/* --cpus-per-task=n, -c n	*/
	bool cpus_set;		/* true if cpus_per_task explicitly set */
	int  min_nodes;		/* --nodes=n,       -N n	*/ 
	int  max_nodes;		/* --nodes=x-n,       -N x-n	*/ 
	cpu_bind_type_t cpu_bind_type; /* --cpu_bind=           */
	char *cpu_bind;		/* binding map for map/mask_cpu */
	mem_bind_type_t mem_bind_type; /* --mem_bind=		*/
	char *mem_bind;		/* binding map for map/mask_mem	*/
	bool nodes_set;		/* true if nodes explicitly set */
	int  time_limit;	/* --time,   -t			*/
	char *partition;	/* --partition=n,   -p n   	*/
	enum task_dist_states
		distribution;	/* --distribution=, -m dist	*/
	char *job_name;		/* --job-name=,     -J name	*/
	unsigned int jobid;     /* --jobid=jobid                */
	bool jobid_set;		/* true of jobid explicitly set */
	char *mpi_type;		/* --mpi=type			*/
	unsigned int dependency;/* --dependency, -P jobid	*/
	int nice;		/* --nice			*/
	char *account;		/* --account, -U acct_name	*/

	int immediate;		/* -i, --immediate      	*/

	bool hold;		/* --hold, -H			*/
	bool noshell;		/* --noshell                    */
	bool overcommit;	/* --overcommit,   -O		*/
	bool no_kill;		/* --no-kill, -k		*/
	bool kill_bad_exit;	/* --kill-on-bad-exit, -K	*/
	bool no_requeue;	/* --no-requeue			*/
	bool share;		/* --share,   -s		*/
	int  max_wait;		/* --wait,    -W		*/
	bool quit_on_intr;      /* --quit-on-interrupt, -q      */
	bool disable_status;    /* --disable-status, -X         */
	int  quiet;
	bool test_only;		/* --test-only			*/
	char *propagate;	/* --propagate[=RLIMIT_CORE,...]*/
	char *task_epilog;	/* --task-epilog=		*/
	char *task_prolog;	/* --task-prolog=		*/

	/* constraint options */
	int mincpus;		/* --mincpus=n			*/
	int realmem;		/* --mem=n			*/
	long tmpdisk;		/* --tmp=n			*/
	char *constraints;	/* --constraints=, -C constraint*/
	bool contiguous;	/* --contiguous			*/
	char *nodelist;		/* --nodelist=node1,node2,...	*/
	char *exc_nodes;	/* --exclude=node1,node2,... -x	*/
	char *relative;		/* --relative -r N              */
	bool no_alloc;		/* --no-allocate, -Z		*/
	int  max_launch_time;   /* Undocumented                 */
	int  max_exit_timeout;  /* Undocumented                 */
	int  msg_timeout;       /* Undocumented                 */
	char *network;		/* --network=			*/
        bool exclusive;         /* --exclusive                  */

	uint16_t geometry[SYSTEM_DIMENSIONS]; /* --geometry, -g	*/
	bool no_rotate;		/* --no_rotate, -R		*/
	int16_t conn_type;	/* --conn-type 			*/
	time_t begin;		/* --begin			*/
	uint16_t mail_type;	/* --mail-type			*/
	char *mail_user;	/* --mail-user			*/
	char *ctrl_comm_ifhn;	/* --ctrl-comm-ifhn		*/
} opt_t;

extern opt_t opt;

/* return whether any constraints were specified by the user 
 * (if new constraints are added above, might want to add them to this
 *  macro or move this to a function if it gets a little complicated)
 */
#define constraints_given() opt.mincpus != -1 || opt.realmem != -1 ||\
                            opt.tmpdisk != -1 || opt.contiguous   

/* process options:
 * 1. set defaults
 * 2. update options with env vars
 * 3. update options with commandline args
 * 4. perform some verification that options are reasonable
 */
int initialize_and_process_args(int argc, char *argv[]);

/* set options based upon commandline args */
void set_options(const int argc, char **argv, int first);


#endif	/* _HAVE_OPT_H */
