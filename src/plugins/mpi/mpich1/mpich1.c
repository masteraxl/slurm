/*****************************************************************************\
 *  mpich1.c - srun support for MPICH1
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).  
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif

#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <strings.h>
#include <sys/poll.h>
#include <sys/time.h>

#include "src/common/slurm_xlator.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/net.h"
#include "src/common/fd.h"

#include "src/plugins/mpi/mpich1/mpich1.h"

/* NOTE: AIX lacks timersub */
#ifndef timersub
#  define timersub(a, b, result)					\
	do {								\
		(result)->tv_sec = (a)->tv_sec - (b)->tv_sec;		\
		(result)->tv_usec = (a)->tv_usec - (b)->tv_usec;	\
		if ((result)->tv_usec < 0) {				\
			--(result)->tv_sec;				\
			(result)->tv_usec += 1000000;			\
		}							\
	} while (0)
#endif

/*
 *  Information read from each MPICH1 process
 */
struct mpich1_info
{
	int do_poll;          
	int fd;             /* fd for socket connection to MPI task  */
	int rank;           /* This process' MPI rank                */
	int pidlen;         /* length of pid buffer                  */
	char *pid;          /* This rank's local pid (V3 only)       */
	int hostidlen;      /* Host id length                        */
	int hostid;         /* Separate hostid (for protocol v5)     */
	int addrlen;        /* Length of addr array in bytes         */

	int *addr;          /* This process' address array, which for
	                     *  process rank N in an M process job 
	                     *  looks like:
	                     *
	                     *   qp0,qp1,..,lid,qpN+1,..,qpM-1, hostid
	                     *
	                     *  Where position N is this rank's lid,
	                     *  and the hostid is tacked onto the end
	                     *  of the array (for protocol version 3)
	                     */
};

/*  Globals for the mpich1 thread.
 */
int mpich1_verbose = 0;
static time_t first_abort_time = 0;

/*  Per-job step state information.  The MPI plugin may be called
 *  multiple times from the SLURM API's slurm_step_launch() in the
 *  same process.
 */
struct mpich1_state {
	pthread_t tid;
	struct mpich1_info **mvarray;
	int fd;
	int nprocs;
	int protocol_version;
	int protocol_phase;
	int connect_once;
	int do_timing;

	mpi_plugin_client_info_t job[1];
};

#define mpich1_debug(args...) \
	do { \
		if (mpich1_verbose) \
			info ("mpich1: " args); \
	} while (0);

#define mpich1_debug2(args...) \
	do { \
		if (mpich1_verbose > 1) \
			info ("mpich1: " args); \
	} while (0);

static struct mpich1_info * mpich1_info_create (void)
{
	struct mpich1_info *mvi = xmalloc (sizeof (*mvi));
	memset (mvi, 0, sizeof (*mvi));
	mvi->fd = -1;
	mvi->rank = -1;
	return (mvi);
}

static void mpich1_info_destroy (struct mpich1_info *mvi)
{
	xfree (mvi->addr);
	xfree (mvi->pid);
	xfree (mvi);
	return;
}

static int mpich1_requires_pids (mpich1_state_t *st)
{
	if ( st->protocol_version == 5
	  || st->protocol_version == 6 )
		return (1);
	return (0);
}

/*
 *  Return non-zero if protocol version has two phases.
 */
static int mpich1_dual_phase (mpich1_state_t *st)
{
	return (st->protocol_version == 5 || st->protocol_version == 6);
}

static int mpich1_abort_sends_rank (mpich1_state_t *st)
{
	if (st->protocol_version >= 3)
		return (1);
	return (0);
}

/*
 *  Create an mpich1_info object by reading information from
 *   file descriptor `fd'
 */
static int mpich1_get_task_info (mpich1_state_t *st,
				  struct mpich1_info *mvi)
{
	int fd = mvi->fd;

	if (fd_read_n (fd, &mvi->addrlen, sizeof (int)) <= 0)
		return error ("mpich1: Unable to read addrlen for rank %d: %m", 
				mvi->rank);

	mvi->addr = xmalloc (mvi->addrlen);

	if (fd_read_n (fd, mvi->addr, mvi->addrlen) <= 0)
		return error ("mpich1: Unable to read addr info for rank %d: %m", 
				mvi->rank);

	if (!mpich1_requires_pids (st))
		return (0);

	if (fd_read_n (fd, &mvi->pidlen, sizeof (int)) <= 0) {
		return error ("mpich1: Unable to read pidlen for rank %d: %m", 
				mvi->rank);
	}

	mvi->pid = xmalloc (mvi->pidlen);

	if (fd_read_n (fd, mvi->pid, mvi->pidlen) <= 0) {
		return error ("mpich1: Unable to read pid for rank %d: %m", 
				mvi->rank);
	}

	mvi->do_poll = 0;

	return (0);
}

static int mpich1_get_hostid (struct mpich1_info *mvi)
{
	if (fd_read_n (mvi->fd, &mvi->hostidlen, sizeof (int)) < 0) {
		return error ("mpich1: Unable to read hostidlen for rank %d: %m",
				mvi->rank);
	}
	if (mvi->hostidlen != sizeof (int)) {
		return error ("mpich1: Unexpected size for hostidlen (%d)", 
				mvi->hostidlen);
	}
	if (fd_read_n (mvi->fd, &mvi->hostid, sizeof (int)) < 0) {
		return error ("mpich1: unable to read hostid from rank %d", 
				mvi->rank);
	}

	return (0);
}

static int mpich1_get_task_header (mpich1_state_t *st,
				    int fd, int *version, int *rank)
{
	/*
	 *  dual phase only sends version on first pass
	 */
	if (!mpich1_dual_phase (st) || st->protocol_phase == 0) {
		if (fd_read_n (fd, version, sizeof (int)) < 0) 
			return error ("mpich1: Unable to read version from task: %m");
	} 

	if (fd_read_n (fd, rank, sizeof (int)) < 0) 
		return error ("mpich1: Unable to read task rank: %m");

	if (mpich1_dual_phase (st) && st->protocol_phase > 0)
		return (0);

	if (st->protocol_version == -1)
		st->protocol_version = *version;
	else if (st->protocol_version != *version) {
		return error ("mpich1: rank %d version %d != %d",
			      *rank, *version, st->protocol_version);
	}

	return (0);

}

static int mpich1_handle_task (mpich1_state_t *st,
				int fd, struct mpich1_info *mvi)
{
	mvi->fd = fd;

	switch (st->protocol_version) {
		case 1:
		case 2:
		case 3:
			return mpich1_get_task_info (st, mvi);
		case 5:
		case 6:
			if (st->protocol_phase == 0)
				return mpich1_get_hostid (mvi);
			else
				return mpich1_get_task_info (st, mvi);
		case 8:
			return (0);
		default:
			return (error ("mpich1: Unsupported protocol version %d", 
				       st->protocol_version));
	}

	return (0);
}

/*
 *  Broadcast addr information to all connected mpich1 processes.
 *   The format of the information sent back to each process is:
 *
 *   for rank N in M process job:
 *   
 *    lid info :  lid0,lid1,...lidM-1
 *    qp info  :  qp0, qp1, ..., -1, qpN+1, ...,qpM-1
 *    hostids  :  hostid0,hostid1,...,hostidM-1
 *
 *   total of 3*nprocs ints.
 *
 */   
static void mpich1_bcast_addrs (mpich1_state_t *st)
{
	struct mpich1_info *m;
	int out_addrs_len = 3 * st->nprocs * sizeof (int);
	int *out_addrs = xmalloc (out_addrs_len);
	int i = 0;
	int j = 0;

	for (i = 0; i < st->nprocs; i++) {
		m = st->mvarray[i];
		/*
		 * lids are found in addrs[rank] for each process
		 */
		out_addrs[i] = m->addr[m->rank];

		/*
		 * hostids are the last entry in addrs
		 */
		out_addrs[2 * st->nprocs + i] =
			m->addr[(m->addrlen/sizeof (int)) - 1];
	}

	for (i = 0; i < st->nprocs; i++) {
		m = st->mvarray[i];

		/*
		 * qp array is tailored to each process.
		 */
		for (j = 0; j < st->nprocs; j++)  
			out_addrs[st->nprocs + j] = 
				(i == j) ? -1 : st->mvarray[j]->addr[i];

		fd_write_n (m->fd, out_addrs, out_addrs_len);

		/*
		 * Protocol version 3 requires pid list to be sent next
		 */
		if (mpich1_requires_pids (st)) {
			for (j = 0; j < st->nprocs; j++)
				fd_write_n (m->fd, &st->mvarray[j]->pid,
					    st->mvarray[j]->pidlen);
		}

	}

	xfree (out_addrs);
	return;
}

static void mpich1_bcast_hostids (mpich1_state_t *st)
{
	int *  hostids;
	int    i   = 0;
	size_t len = st->nprocs * sizeof (int);

	hostids = xmalloc (len);

	for (i = 0; i < st->nprocs; i++)
		hostids [i] = st->mvarray[i]->hostid;

	for (i = 0; i < st->nprocs; i++) {
		struct mpich1_info *mvi = st->mvarray[i];
		int co, rc;
		if (fd_write_n (mvi->fd, hostids, len) < 0)
			error ("mpich1: write hostid rank %d: %m", mvi->rank);

		if ((rc = fd_read_n (mvi->fd, &co, sizeof (int))) <= 0) {
			close (mvi->fd);
			st->connect_once = 0;
		} else
			mvi->do_poll = 1;
	}

	xfree (hostids);
}

/* Write size bytes from buf into socket for rank */
static void mpich1_send (mpich1_state_t *st, void* buf, int size, int rank)
{
	struct mpich1_info *mvi = st->mvarray [rank];
	if (fd_write_n (mvi->fd, buf, size) < 0)
		error ("mpich1: write hostid rank %d: %m", mvi->rank);
}

/* Read size bytes from socket for rank into buf */
static void mpich1_recv (mpich1_state_t *st, void* buf, int size, int rank)
{
	struct mpich1_info *mvi = st->mvarray [rank];

	int rc;
	if ((rc = fd_read_n (mvi->fd, buf, size)) <= 0) {
		error("mpich1 reading from %d", mvi->rank);
	}
}

/* Read an integer from socket for rank */
static int mpich1_recv_int (mpich1_state_t *st, int rank)
{
	int buf;
	mpich1_recv(st, &buf, sizeof(buf), rank);
	return buf;
}

/* Scatter data in buf to ranks using chunks of size bytes */
static void mpich1_scatterbcast (mpich1_state_t *st, void* buf, int size)
{
	int i;
	for (i = 0; i < st->nprocs; i++)
		mpich1_send(st, buf + i*size, size, i);
}

/* Broadcast buf to each rank, which is size bytes big */
static void mpich1_allgatherbcast (mpich1_state_t *st, void* buf, int size)
{
	int i;
	for (i = 0; i < st->nprocs; i++)
		mpich1_send(st, buf, size, i);
}

/* Perform alltoall using data in buf with elements of size bytes */
static void mpich1_alltoallbcast (mpich1_state_t *st, void* buf, int size)
{
	int pbufsize = size * st->nprocs;
	void* pbuf = xmalloc(pbufsize);	

	int i, src;
	for (i = 0; i < st->nprocs; i++) {
		for (src = 0; src < st->nprocs; src++) {
			memcpy( pbuf + size*src,
				buf + size*(src*st->nprocs + i),
				size
				);
		}
		mpich1_send(st, pbuf, pbufsize, i);
	}
	
	xfree(pbuf);
}

/* Check that new == curr value if curr has been initialized */
static int set_current (int curr, int new)
{
	if (curr == -1)
		curr = new;
	if (new != curr) {
		error("PMGR unexpected value: received %d, expecting %d", 
			new, curr);
	}
	return curr;
}

/* 
 * This function carries out pmgr_collective operations to
 * bootstrap MPI.  These collective operations are modeled after
 * MPI collectives -- all tasks must call them in the same order
 * and with consistent parameters.
 *
 * Until a 'CLOSE' or 'ABORT' message is seen, we continuously loop
 * processing ops
 *   For each op, we read one packet from each rank (socket)
 *     A packet consists of an integer OP CODE, followed by variable
 *     length data depending on the operation
 *   After reading a packet from each rank, srun completes the
 *   operation by broadcasting data back to any destinations,
 *   depending on the operation being performed
 *
 * Note: Although there are op codes available for PMGR_OPEN and
 * PMGR_ABORT, neither is fully implemented and should not be used.
 */
static void mpich1_processops (mpich1_state_t *st)
{
	/* Until a 'CLOSE' or 'ABORT' message is seen, we continuously 
	 *  loop processing ops
	 */
	int exit = 0;
	while (!exit) {
	int opcode = -1;
	int root   = -1;
	int size   = -1;
	void* buf = NULL;

	mpich1_debug ("Processing PMGR opcodes");

	// for each process, read in one opcode and its associated data
	int i;
	for (i = 0; i < st->nprocs; i++) {
		struct mpich1_info *mvi = st->mvarray [i];

		// read in opcode
		opcode = set_current(opcode, mpich1_recv_int(st, i));

		// read in additional data depending on current opcode
		int rank, code;
		switch(opcode) {
		case 0: // PMGR_OPEN (followed by rank)
			rank = mpich1_recv_int(st, i);
			break;
		case 1: // PMGR_CLOSE (no data, close the socket)
			close(mvi->fd);
			break;
		case 2: // PMGR_ABORT (followed by exit code)
			code = mpich1_recv_int(st, i);
			error("mpich1 abort with code %d from rank %d", 
				code, i);
			break;
		case 3: // PMGR_BARRIER (no data)
			break;
		case 4: // PMGR_BCAST (root, size of message, 
			// then message data (from root only))
			root = set_current(root, mpich1_recv_int(st, i));
			size = set_current(size, mpich1_recv_int(st, i));
			if (!buf) buf = (void*) xmalloc(size);
			if (i == root) mpich1_recv(st, buf, size, i);
			break;
		case 5: // PMGR_GATHER (root, size of message, 
			// then message data)
			root = set_current(root, mpich1_recv_int(st, i));
			size = set_current(size, mpich1_recv_int(st, i));
			if (!buf) buf = (void*) xmalloc(size * st->nprocs);
			mpich1_recv(st, buf + size*i, size, i);
			break;
		case 6: // PMGR_SCATTER (root, size of message, 
			// then message data)
			root = set_current(root, mpich1_recv_int(st, i));
			size = set_current(size, mpich1_recv_int(st, i));
			if (!buf) buf = (void*) xmalloc(size * st->nprocs);
			if (i == root) mpich1_recv(st, buf, size * st->nprocs, i);
			break;
		case 7: // PMGR_ALLGATHER (size of message, then message data)
			size = set_current(size, mpich1_recv_int(st, i));
			if (!buf) buf = (void*) xmalloc(size * st->nprocs);
			mpich1_recv(st, buf + size*i, size, i);
			break;
		case 8: // PMGR_ALLTOALL (size of message, then message data)
			size = set_current(size, mpich1_recv_int(st, i));
			if (!buf) buf = (void*) xmalloc(size * st->nprocs * st->nprocs);
			mpich1_recv(st, buf + (size*st->nprocs)*i, size * st->nprocs, i);
			break;
		default:
			error("Unrecognized PMGR opcode: %d", opcode);
		}
	}

	// Complete any operations
	switch(opcode) {
		case 0: // PMGR_OPEN
			mpich1_debug ("Completed PMGR_OPEN");
			break;
		case 1: // PMGR_CLOSE
			mpich1_debug ("Completed PMGR_CLOSE");
			exit = 1;
			break;
		case 2: // PMGR_ABORT
			mpich1_debug ("Completed PMGR_ABORT");
			exit = 1;
			break;
		case 3: // PMGR_BARRIER (just echo the opcode back)
			mpich1_debug ("Completing PMGR_BARRIER");
			mpich1_allgatherbcast (st, &opcode, sizeof(opcode));
			mpich1_debug ("Completed PMGR_BARRIER");
			break;
		case 4: // PMGR_BCAST
			mpich1_debug ("Completing PMGR_BCAST");
			mpich1_allgatherbcast (st, buf, size);
			mpich1_debug ("Completed PMGR_BCAST");
			break;
		case 5: // PMGR_GATHER
			mpich1_debug ("Completing PMGR_GATHER");
			mpich1_send (st, buf, size * st->nprocs, root);
			mpich1_debug ("Completed PMGR_GATHER");
			break;
		case 6: // PMGR_SCATTER
			mpich1_debug ("Completing PMGR_SCATTER");
			mpich1_scatterbcast (st, buf, size);
			mpich1_debug ("Completed PMGR_SCATTER");
			break;
		case 7: // PMGR_ALLGATHER
			mpich1_debug ("Completing PMGR_ALLGATHER");
			mpich1_allgatherbcast (st, buf, size * st->nprocs);
			mpich1_debug ("Completed PMGR_ALLGATHER");
			break;
		case 8: // PMGR_ALLTOALL
			mpich1_debug ("Completing PMGR_ALLTOALL");
			mpich1_alltoallbcast (st, buf, size);
			mpich1_debug ("Completed PMGR_ALLTOALL");
			break;
		default:
			error("Unrecognized PMGR opcode: %d", opcode);
	}

	xfree(buf);
  } // while(!exit)
  mpich1_debug ("Completed processing PMGR opcodes");
}

static void mpich1_bcast (mpich1_state_t *st)
{
	if (!mpich1_dual_phase (st) || st->protocol_phase > 0)
		return mpich1_bcast_addrs (st);
	else
		return mpich1_bcast_hostids (st);
}

static void mpich1_barrier (mpich1_state_t *st)
{
	int i;
	struct mpich1_info *m;
	/*
	 *  Simple barrier to wait for qp's to come up. 
	 *   Once all processes have written their rank over the socket,
	 *   simply write their rank right back to them.
	 */

	debug ("mpich1: starting barrier");

	for (i = 0; i < st->nprocs; i++) {
		int j;
		m = st->mvarray[i];
		if (fd_read_n (m->fd, &j, sizeof (j)) == -1)
			error("mpich1 read on barrier");
	}

	debug ("mpich1: completed barrier for all tasks");

	for (i = 0; i < st->nprocs; i++) {
		m = st->mvarray[i];
		if (fd_write_n (m->fd, &i, sizeof (i)) == -1)
			error("mpich1: write on barrier: %m");
		close (m->fd);
		m->fd = -1;
	}

	return;
}

static void 
mpich1_print_abort_message (mpich1_state_t *st, int rank,
			     int dest, char *msg, int msglen)
{
	slurm_step_layout_t *sl = st->job->step_layout;
	char *host;
	char *msgstr;

	if (!mpich1_abort_sends_rank (st)) {
		info ("mpich1: Received ABORT message from an MPI process.");
		return;
	}

	if (msg && (msglen > 0)) {
		/* 
		 *  Remove trailing newline if it exists (syslog will add newline)
		 */
		if (msg [msglen - 1] == '\n')
			msg [msglen - 1] = '\0';

		msgstr = msg;
	} 
	else {
		msgstr = "";
		msglen = 0;
	}

	host = slurm_step_layout_host_name(
		sl, slurm_step_layout_host_id(sl, rank));

	if (dest >= 0) {
		const char *dsthost = slurm_step_layout_host_name (sl, dest);

		info ("mpich1: %M: ABORT from MPI rank %d [on %s] dest rank %d [on %s]",
		      rank, host, dest, dsthost);

		/*
		 *  Log the abort event to syslog
		 *   so that system administrators know about possible HW events.
		 */
		openlog ("srun", 0, LOG_USER);
		syslog (LOG_WARNING, 
				"MVAPICH ABORT [jobid=%u.%u src=%d(%s) dst=%d(%s)]: %s",
				st->job->jobid, st->job->stepid, 
				rank, host, dest, dsthost, msgstr);
		closelog();
	}
	else {
		info ("mpich1: %M: ABORT from MPI rank %d [on %s]", 
				rank, host);
		/*
		 *  Log the abort event to syslog
		 *   so that system administrators know about possible HW events.
		 */
		openlog ("srun", 0, LOG_USER);
		syslog (LOG_WARNING, 
				"MVAPICH ABORT [jobid=%u.%u src=%d(%s) dst=-1()]: %s",
				st->job->jobid, st->job->stepid, 
				rank, host, msgstr);
		closelog();

	}
	return;
}


static int mpich1_abort_timeout (void)
{
	int timeout;

	if (first_abort_time == 0)
		return (-1);

	timeout = 60 - (time (NULL) - first_abort_time);

	if (timeout < 0)
		return (0);

	return (timeout * 1000);
}

static int mpich1_accept (uint32_t jobid, uint32_t stepid, int fd)
{
	slurm_addr addr;
	int rc;
	struct pollfd pfds[1];

	pfds->fd = fd;
	pfds->events = POLLIN;

	while ((rc = poll (pfds, 1, mpich1_abort_timeout ())) < 0) {
		if (errno != EINTR)
			return (-1);
	}

	/* 
	 *  If poll() timed out, forcibly kill job and exit instead of
	 *   waiting longer for remote IO, process exit, etc.
	 */
	if (rc == 0) {
		error("Timeout waiting for all tasks after MVAPICH ABORT. Exiting.");
		slurm_signal_job_step(jobid, stepid, SIGKILL);
		exit(1);
		/* NORETURN */
	}

	return (slurm_accept_msg_conn (fd, &addr));
}


static void mpich1_wait_for_abort(mpich1_state_t *st)
{
	int src, dst;
	int ranks[2];
	int n;
	char msg [1024] = "";
	int msglen = 0;

	/*
	 *  Wait for abort notification from any process.
	 *  For mpich 0.9.4, it appears that an MPI_Abort is registered
	 *   simply by connecting to this socket and immediately closing
	 *   the connection. In other versions, the process may write
	 *   its rank.
	 */
	while (1) {
		int newfd = mpich1_accept (st->job->jobid, st->job->stepid,
					    st->fd);

		if (newfd == -1) {
			fatal("MPI master failed to accept (abort-wait)");
		}

		fd_set_blocking (newfd);

		ranks[1] = -1;
		if ((n = fd_read_n (newfd, &ranks, sizeof (ranks))) < 0) {
			error("mpich1: MPI recv (abort-wait) failed");
			close (newfd);
			continue;
		}

		/*
		 *  If we read both src/dest rank, then also try to 
		 *   read an error message. If this fails, msglen will
		 *   stay zero and no message will be printed.
		 */
		if (ranks[1] >= 0) {
			dst = ranks[0];
			src = ranks[1];
			fd_read_n (newfd, &msglen, sizeof (int));
			if (msglen)
				fd_read_n (newfd, msg, msglen);
		} else {
			src = ranks[0];
			dst = -1;
		}

		close(newfd);

		mpich1_print_abort_message (st, src, dst, msg, msglen);
		slurm_signal_job_step(st->job->jobid, st->job->stepid, SIGKILL);
		if (!first_abort_time)
			first_abort_time = time (NULL);
	}

	return; /* but not reached */
}

static void mpich1_mvarray_create (mpich1_state_t *st)
{
	int i;
	st->mvarray = xmalloc (st->nprocs * sizeof (*(st->mvarray)));
	for (i = 0; i < st->nprocs; i++) {
		st->mvarray [i] = mpich1_info_create ();
		st->mvarray [i]->rank = i;
	}
}

static void mpich1_mvarray_destroy (mpich1_state_t *st)
{
	int i;
	for (i = 0; i < st->nprocs; i++)
		mpich1_info_destroy (st->mvarray[i]);
	xfree (st->mvarray);
}

static int mpich1_rank_from_fd (mpich1_state_t *st, int fd)
{
	int rank = 0;
	while (st->mvarray[rank]->fd != fd)
		rank++;
	return (rank);
}

static int mpich1_handle_connection (mpich1_state_t *st, int fd)
{
	int version, rank;

	if (st->protocol_phase == 0 || !st->connect_once) {
		if (mpich1_get_task_header (st, fd, &version, &rank) < 0)
			return (-1);

		st->mvarray[rank]->rank = rank;

		if (rank > st->nprocs - 1) { 
			return (error ("mpich1: task reported invalid rank (%d)", 
					rank));
		}
	} else {
		rank = mpich1_rank_from_fd (st, fd);
	}

	if (mpich1_handle_task (st, fd, st->mvarray[rank]) < 0) 
		return (-1);

	return (0);
}

static int poll_mpich1_fds (mpich1_state_t *st)
{
	int i = 0;
	int j = 0;
	int rc;
	int fd;
	int nfds = 0;
	struct pollfd *fds = xmalloc (st->nprocs * sizeof (struct pollfd));

	for (i = 0; i < st->nprocs; i++) {
		if (st->mvarray[i]->do_poll) {
			fds[j].fd = st->mvarray[i]->fd;
			fds[j].events = POLLIN;
			j++;
			nfds++;
		}
	}

	mpich1_debug2 ("Going to poll %d fds", nfds);
	if ((rc = poll (fds, nfds, -1)) < 0) {
		error ("mpich1: poll: %m");
		xfree (fds);
		return SLURM_ERROR;
	}

	i = 0;
	while (fds[i].revents != POLLIN)
		i++;

	fd = fds[i].fd;
	xfree (fds);

	return (fd);
}

static int mpich1_get_next_connection (mpich1_state_t *st)
{
	slurm_addr addr;
	int fd;

	if (st->connect_once && st->protocol_phase > 0) {
		return (poll_mpich1_fds (st));
	} 
		
	if ((fd = slurm_accept_msg_conn (st->fd, &addr)) < 0) {
		error ("mpich1: accept: %m");
		return (-1);
	}
	mpich1_debug2 ("accept() = %d", fd);

	return (fd);
}

static void do_timings (mpich1_state_t *st)
{
	static int initialized = 0;
	static struct timeval initv = { 0, 0 };
	struct timeval tv;
	struct timeval result;

	if (!st->do_timing)
		return;

	if (!initialized) {
		if (gettimeofday (&initv, NULL) < 0)
			error ("mpich1: do_timings(): gettimeofday(): %m\n");
		initialized = 1;
		return;
	}

	if (gettimeofday (&tv, NULL) < 0) {
		error ("mpich1: do_timings(): gettimeofday(): %m\n");
		return;
	}

	timersub (&tv, &initv, &result);

	info ("mpich1: Intialization took %d.%03d seconds", result.tv_sec,
			result.tv_usec/1000);

	return;
}

static void *mpich1_thr(void *arg)
{
	mpich1_state_t *st = arg;
	mpi_plugin_client_info_t *job = st->job;
	int i = 0;
	int first = 1;

	debug ("mpich1: thread started: %ld", pthread_self ());

	mpich1_mvarray_create (st);

again:
	i = 0;
	while (i < st->nprocs) {
		int fd;
		
		mpich1_debug ("Waiting to accept remote connection %d of %d\n", 
				i, st->nprocs);

		if ((fd = mpich1_get_next_connection (st)) < 0) {
			error ("mpich1: accept: %m");
			goto fail;
		}

		if (first) {
			mpich1_debug ("first task checked in");
			do_timings (st);
			first = 0;
		}

		if (mpich1_handle_connection (st, fd) < 0) 
			goto fail;

		i++;
	}

	if (st->protocol_version == 8) {
		mpich1_processops(st);
	} else {
		mpich1_debug ("bcasting mpich1 info to %d tasks", st->nprocs);
		mpich1_bcast (st);

		if (mpich1_dual_phase (st) && st->protocol_phase == 0) {
			st->protocol_phase = 1;
			goto again;
		}

		mpich1_debug ("calling mpich1_barrier");
		mpich1_barrier (st);
		mpich1_debug ("all tasks have checked in");
	}

	do_timings (st);

	mpich1_wait_for_abort (st);

	mpich1_mvarray_destroy (st);

	return (NULL);

fail:
	error ("mpich1: fatal error, killing job");
	slurm_signal_job_step(job->jobid, job->stepid, SIGKILL);
	return (void *)0;
}

static int process_environment (mpich1_state_t *st)
{
	char *val;

	if (getenv ("MVAPICH_CONNECT_TWICE"))
		st->connect_once = 0;

	if ((val = getenv ("SLURM_MVAPICH_DEBUG"))) {
		int level = atoi (val);
		if (level > 0)
			mpich1_verbose = level;
	}

	if (getenv ("SLURM_MVAPICH_TIMING"))
		st->do_timing = 1;

	return (0);
}

static mpich1_state_t * mpich1_state_create(const mpi_plugin_client_info_t *job)
{
	mpich1_state_t *state;

	state = (mpich1_state_t *)xmalloc(sizeof(mpich1_state_t));

	state->tid		= (pthread_t)-1;
	state->mvarray          = NULL;
	state->fd               = -1;
	state->nprocs           = job->step_layout->task_cnt;
	state->protocol_version = -1;
	state->protocol_phase   = 0;
	state->connect_once     = 1;
	state->do_timing        = 0;

	*(state->job) = *job;

	return state;
}

static void mpich1_state_destroy(mpich1_state_t *st)
{
	xfree(st);
}

extern mpich1_state_t *mpich1_thr_create(const mpi_plugin_client_info_t *job,
					   char ***env)
{
	short port;
	pthread_attr_t attr;
	mpich1_state_t *st = NULL;

	st = mpich1_state_create(job);
	if (process_environment (st) < 0) {
		error ("mpich1: Failed to read environment settings\n");
		mpich1_state_destroy(st);
		return NULL;
	}
	if (net_stream_listen(&st->fd, &port) < 0) {
		error ("Unable to create ib listen port: %m");
		mpich1_state_destroy(st);
		return NULL;
	}

	/*
	 * Accept in a separate thread.
	 */
	slurm_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
	if (pthread_create(&st->tid, &attr, &mpich1_thr, (void *)st)) {
		slurm_attr_destroy(&attr);
		mpich1_state_destroy(st);
		return NULL;
	}
	slurm_attr_destroy(&attr);

	/*
	 *  Set some environment variables in current env so they'll get
	 *   passed to all remote tasks
	 */
	env_array_overwrite_fmt(env, "MPIRUN_PORT",   "%hu", port);
	env_array_overwrite_fmt(env, "MPIRUN_NPROCS", "%d", st->nprocs);
	env_array_overwrite_fmt(env, "MPIRUN_ID",     "%d", st->job->jobid);
	if (st->connect_once) {
		env_array_overwrite_fmt(env, "MPIRUN_CONNECT_ONCE", "1");
	}

	verbose ("mpich1 master listening on port %d", port);

	return st;
}

extern int mpich1_thr_destroy(mpich1_state_t *st)
{
	if (st != NULL) {
		if (st->tid != (pthread_t)-1) {
			pthread_cancel(st->tid);
			pthread_join(st->tid, NULL);
		}
		mpich1_state_destroy(st);
	}
	return SLURM_SUCCESS;
}
