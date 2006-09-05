/*****************************************************************************\
 *  opt.h - definitions for srun option processing
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

#ifndef _HAVE_OPT_H
#define _HAVE_OPT_H

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "src/common/macros.h" /* true and false */
#include "src/srun/core-format.h"
#include "src/common/env.h"
#include "src/srun/fname.h"

#define MAX_THREADS	32
#define MAX_USERNAME	9

#define INT_UNASSIGNED ((int)-1)

/* global variables relating to user options */
extern char **remote_argv;
extern int remote_argc;
extern int _verbose;

/* mutually exclusive modes for srun */
enum modes {
	MODE_UNKNOWN	= 0,
	MODE_NORMAL	= 1,
	MODE_IMMEDIATE	= 2,
	MODE_ATTACH	= 3,
	MODE_ALLOCATE	= 4,
	MODE_BATCH	= 5
};

extern enum modes mode;

#define format_task_dist_states(t) (t == SLURM_DIST_BLOCK) ? "block" :   \
		                 (t == SLURM_DIST_CYCLIC) ? "cyclic" : \
		                 (t == SLURM_DIST_PLANE) ? "plane" : \
		                 (t == SLURM_DIST_CYCLIC_CYCLIC) ? "cyclic:cyclic" : \
		                 (t == SLURM_DIST_CYCLIC_BLOCK) ? "cyclic:block" : \
		                 (t == SLURM_DIST_BLOCK_CYCLIC) ? "block:cyclic" : \
		                 (t == SLURM_DIST_BLOCK_BLOCK) ? "block:block" : \
			         (t == SLURM_DIST_ARBITRARY) ? "arbitrary" : \
			         "unknown"

typedef struct srun_options {

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
	int  max_threads;	/* --threads, -T (threads in srun) */
	int  min_nodes;		/* --nodes=n,       -N n	*/ 
	int  max_nodes;		/* --nodes=x-n,       -N x-n	*/ 
        int  min_sockets_per_node; /* --sockets-per-node=n      */
        int  max_sockets_per_node; /* --sockets-per-node=x-n    */
        int  min_cores_per_socket; /* --cores-per-socket=n      */
        int  max_cores_per_socket; /* --cores-per-socket=x-n    */
        int  min_threads_per_core; /* --threads-per-core=n      */
        int  max_threads_per_core; /* --threads-per-core=x-n    */
        int  ntasks_per_node;   /* --ntasks-per-node=n		*/
        int  ntasks_per_socket; /* --ntasks-per-socket=n	*/
        int  ntasks_per_core;   /* --ntasks-per-core=n		*/
	cpu_bind_type_t cpu_bind_type; /* --cpu_bind=           */
	char *cpu_bind;		/* binding map for map/mask_cpu */
	mem_bind_type_t mem_bind_type; /* --mem_bind=		*/
	char *mem_bind;		/* binding map for map/mask_mem	*/
	bool nodes_set;		/* true if nodes explicitly set */
	int  time_limit;	/* --time,   -t			*/
	char *partition;	/* --partition=n,   -p n   	*/
	enum task_dist_states
	        distribution;	/* --distribution=, -m dist	*/
        uint32_t plane_size;    /* lllp distribution -> plane_size for
				 * when -m plane=<# of lllp per
				 * plane> */      
	char *job_name;		/* --job-name=,     -J name	*/
	unsigned int jobid;     /* --jobid=jobid                */
	bool jobid_set;		/* true if jobid explicitly set */
	char *mpi_type;		/* --mpi=type			*/
	unsigned int dependency;/* --dependency, -P jobid	*/
	int nice;		/* --nice			*/
	char *account;		/* --account, -U acct_name	*/

	char *ofname;		/* --output -o filename         */
	char *ifname;		/* --input  -i filename         */
	char *efname;		/* --error, -e filename         */

	int  slurmd_debug;	/* --slurmd-debug, -D           */
	core_format_t core_type;/* --core= 	        	*/
	char *attach;		/* --attach=id	    -a id	*/ 
	bool join;		/* --join, 	    -j		*/

	/* no longer need these, they are set globally : 	*/
	/*int verbose;*/	/* -v, --verbose		*/	
	/*int debug;*/		/* -d, --debug			*/

	int immediate;		/* -i, --immediate      	*/

	bool hold;		/* --hold, -H			*/
	bool labelio;		/* --label-output, -l		*/
	bool unbuffered;        /* --unbuffered,   -u           */
	bool allocate;		/* --allocate, 	   -A		*/
	bool noshell;		/* --noshell                    */
	bool overcommit;	/* --overcommit,   -O		*/
	bool batch;		/* --batch,   -b		*/
	bool no_kill;		/* --no-kill, -k		*/
	bool kill_bad_exit;	/* --kill-on-bad-exit, -K	*/
	bool no_requeue;	/* --no-requeue			*/
	uint16_t shared;	/* --share,   -s		*/
	int  max_wait;		/* --wait,    -W		*/
	bool quit_on_intr;      /* --quit-on-interrupt, -q      */
	bool disable_status;    /* --disable-status, -X         */
	int  quiet;
	bool parallel_debug;	/* srun controlled by debugger	*/
	bool debugger_test;	/* --debugger-test		*/
	bool test_only;		/* --test-only			*/
	char *propagate;	/* --propagate[=RLIMIT_CORE,...]*/
	char *task_epilog;	/* --task-epilog=		*/
	char *task_prolog;	/* --task-prolog=		*/
        bool printreq;          /* --print-request              */

	/* constraint options */
	int job_min_cpus;	/* --mincpus=n			*/
	int job_min_sockets;	/* --minsockets=n		*/
	int job_min_cores;	/* --mincores=n			*/
	int job_min_threads;	/* --minthreads=n		*/
	int job_min_memory;	/* --mem=n			*/
	long job_min_tmp_disk;	/* --tmp=n			*/
	char *constraints;	/* --constraints=, -C constraint*/
	bool contiguous;	/* --contiguous			*/
	char *nodelist;		/* --nodelist=node1,node2,...	*/
	char *alloc_nodelist;   /* grabbed from the environment */
	char *exc_nodes;	/* --exclude=node1,node2,... -x	*/
	int  relative;		/* --relative -r N              */
	bool relative_set;
	bool no_alloc;		/* --no-allocate, -Z		*/
	int  max_launch_time;   /* Undocumented                 */
	int  max_exit_timeout;  /* Undocumented                 */
	int  msg_timeout;       /* Undocumented                 */
	char *network;		/* --network=			*/

	uint16_t geometry[SYSTEM_DIMENSIONS]; /* --geometry, -g	*/
	bool no_rotate;		/* --no_rotate, -R		*/
	int16_t conn_type;	/* --conn-type 			*/
	char *prolog;           /* --prolog                     */
	char *epilog;           /* --epilog                     */
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
#define constraints_given() opt.job_min_cpus != INT_UNASSIGNED ||\
			    opt.job_min_memory != INT_UNASSIGNED ||\
			    opt.job_min_tmp_disk != INT_UNASSIGNED ||\
			    opt.job_min_sockets != INT_UNASSIGNED ||\
			    opt.job_min_cores != INT_UNASSIGNED ||\
			    opt.job_min_threads != INT_UNASSIGNED ||\
			    opt.contiguous   

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
