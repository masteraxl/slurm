/****************************************************************************\
 *  io.c - process stdin, stdout, and stderr for parallel jobs.
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona@llnl.gov>, et. al.
 *  UCRL-CODE-2002-040.
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

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/select.h>
#include <sys/poll.h>
#include <arpa/inet.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#include "src/common/fd.h"
#include "src/common/hostlist.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/pack.h"
#include "src/common/slurm_protocol_defs.h"
#include "src/common/slurm_protocol_pack.h"
#include "src/common/slurm_cred.h"
#include "src/common/xassert.h"
#include "src/common/xmalloc.h"
#include "src/common/xsignal.h"
#include "src/common/io_hdr.h"
#include "src/common/net.h"

#include "src/srun/io.h"
#include "src/srun/srun_job.h"
#include "src/srun/opt.h"

static int    fmt_width       = 0;

/* fd_info struct used in poll() loop to map fds back to task number,
 * appropriate output type (stdout/err), and original fd
 */
typedef struct fd_info {
	int taskid;	/* corresponding task id		*/
	int *fd; 	/* pointer to fd in job->out/err array 	*/
	FILE *fp;	/* fp on which to write output		*/
	cbuf_t buf;
} fd_info_t;

static void	_handle_io_init_msg(int fd, srun_job_t *job);
static int	_close_stream(int *fd, FILE *out, int tasknum);
static int	_handle_pollerr(fd_info_t *info);
static ssize_t	_readx(int fd, char *buf, size_t maxbytes);
static int      _read_io_init_msg(int fd, srun_job_t *job, char *host);
static int      _wid(int n);

static bool _listening_socket_readable(eio_obj_t *obj);
static int _listening_socket_read(eio_obj_t *obj, List objs);

/* True if an EOF needs to be broadcast to all tasks
 */
static bool stdin_got_eof = false;
static bool stdin_open    = true;
static uint32_t nbytes    = 0;
static uint32_t nwritten  = 0;

#if 0
struct io_operations server_ops = {
        readable:	&_server_readable,
	writable:	&_server_writable,
	handle_read:	&_server_read,
	handle_write:	&_server_write,
	handle_error:	&_server_error,
	handle_close:   &_obj_close
};

struct io_operations file_read_ops = {
        readable:	&_file_readable,
	handle_read:	&_file_read,
        handle_error:	&_file_error,
	handle_close:   &_obj_close
};

struct io_operations file_write_ops = {
	writable:	&_file_writable,
	handle_write:	&_file_write,
	handle_error:	&_file_error,
	handle_close:   &_obj_close
};
#endif

struct io_operations listening_socket_ops = {
	readable:	&_listening_socket_readable,
	handle_read:	&_listening_socket_read
};

static bool 
_listening_socket_readable(eio_obj_t *obj)
{
	debug3("Called _listening_socket_readable");
	return true;
}

static int
_listening_socket_read(eio_obj_t *obj, List objs)
{
	srun_job_t *job = (srun_job_t *)obj->arg;

	debug3("Called _listening_socket_read");
	_handle_io_init_msg(obj->fd, job);
}

static void
_set_listensocks_nonblocking(srun_job_t *job)
{
	int i;
	for (i = 0; i < job->num_listen; i++) 
		fd_set_nonblocking(job->listensock[i]);
}

static void *
_io_thr_internal(void *job_arg)
{
	srun_job_t *job  = (srun_job_t *) job_arg;
	sigset_t set;

	xassert(job != NULL);

	debug3("IO thread pid = %lu", (unsigned long) getpid());

	/* Block SIGHUP because it is interrupting file stream functions
	 * (fprintf, fflush, etc.) and causing data loss on stdout.
	 */
	sigemptyset(&set);
	sigaddset(&set, SIGHUP);
 	pthread_sigmask(SIG_BLOCK, &set, NULL);

	_set_listensocks_nonblocking(job);

	/* start the eio engine */
	io_handle_events(job->eio, job->eio_objs);

	debug("IO thread exiting");

	return NULL;
}

static FILE *
_fopen(char *filename)
{
	FILE *fp;

	xassert(filename != NULL);

	if (!(fp = fopen(filename, "w"))) 
		error ("Unable to open `%s' for writing: %m", filename);

	return fp;
}


static struct io_operations *
_ops_copy(struct io_operations *ops)
{
	struct io_operations *ret = xmalloc(sizeof(*ops));
	/* 
	 * Copy initial client_ops 
	 */
	*ret = *ops;
	return ret;
}


static eio_obj_t *
_create_listensock_eio(int fd, srun_job_t *job)
{
	eio_obj_t *eio = NULL;

	eio = (eio_obj_t *)xmalloc(sizeof(eio_obj_t));
	eio->fd = fd;
	eio->arg = (void *)job;
	eio->ops = _ops_copy(&listening_socket_ops);
	return eio;
}

int
io_thr_create(srun_job_t *job)
{
	int i;
	pthread_attr_t attr;

	if (opt.labelio)
		fmt_width = _wid(opt.nprocs);

	for (i = 0; i < job->num_listen; i++) {
		eio_obj_t *obj;

		if (net_stream_listen(&job->listensock[i],
				      &job->listenport[i]) < 0)
			fatal("unable to initialize stdio listen socket: %m");
		debug("initialized stdio listening socket, port %d\n",
		      ntohs(job->listenport[i]));
		net_set_low_water(job->listensock[i], 140);
		obj = _create_listensock_eio(job->listensock[i], job);
		list_enqueue(job->eio_objs, obj);
	}

	/* FIXME - Need to open files here (or perhaps earlier) */

	xsignal(SIGTTIN, SIG_IGN);

	slurm_attr_init(&attr);
	if (errno = pthread_create(&job->ioid, &attr,
				   &_io_thr_internal, (void *) job))
		return SLURM_ERROR;

	debug("Started IO server thread (%lu)", (unsigned long) job->ioid);

	return SLURM_SUCCESS;
}

static int
_read_io_init_msg(int fd, srun_job_t *job, char *host)
{
	struct slurm_io_init_msg msg;
	char *sig;
	int siglen;

	if (io_init_msg_read_from_fd(fd, &msg) != SLURM_SUCCESS) {
		error("failed reading io init message");
		goto fail;
	}

	if (slurm_cred_get_signature(job->cred, &sig, &siglen) < 0) {
		error ("Couldn't get existing cred signature");
		goto fail;
	}

	if (io_init_msg_validate(&msg, sig) < 0)
		goto fail;

	if (msg.nodeid >= job->nhosts) {
		error ("Invalid nodeid %d from %s", msg.nodeid, host);
		goto fail;
	}

	debug2("Validated IO connection from %s, node rank %u, sd=%d",
	       host, msg.nodeid, fd);
	return SLURM_SUCCESS;

    fail:
	close(fd);
	return SLURM_ERROR;
}


static bool 
_is_fd_ready(int fd)
{
	struct pollfd pfd[1];
	int    rc;

	pfd[0].fd     = fd;
	pfd[0].events = POLLIN;

	rc = poll(pfd, 1, 10);

	return ((rc == 1) && (pfd[0].revents & POLLIN));
}


static void
_handle_io_init_msg(int fd, srun_job_t *job)
{
	int j;
	debug2("Activity on IO server socket %d", fd);

	for (j = 0; j < 15; j++) {
		int sd;
		struct sockaddr addr;
		struct sockaddr_in *sin;
		socklen_t size = sizeof(addr);
		char buf[INET_ADDRSTRLEN];
		
		/* 
		 * Return early if fd is not now ready
		 */
		if (!_is_fd_ready(fd))
			return;

		while ((sd = accept(fd, &addr, &size)) < 0) {
			if (errno == EINTR)
				continue;
			if (errno == EAGAIN)	/* No more connections */
				return;
			if ((errno == ECONNABORTED) || 
			    (errno == EWOULDBLOCK)) {
				return;
			}
			error("Unable to accept new connection: %m\n");
			return;
		}

		sin = (struct sockaddr_in *) &addr;
		inet_ntop(AF_INET, &sin->sin_addr, buf, INET_ADDRSTRLEN);

		debug3("Accepted IO connection: ip=%s sd=%d", buf, sd); 

		/*
		 * On AIX the new socket [sd] seems to inherit the O_NONBLOCK
		 * flag from the listening socket [fd], so we need to explicitly
		 * set it back to blocking mode.
		 * (XXX: This should eventually be fixed by making
		 *  reads of IO headers nonblocking)
		 */
		fd_set_blocking(sd);

		/*
		 * Read IO header and update job structure appropriately
		 */
		if (_read_io_init_msg(sd, job, buf) < 0)
			continue;

		fd_set_nonblocking(sd);
	}

}

static ssize_t 
_readx(int fd, char *buf, size_t maxbytes)
{
	size_t n;

	if ((n = read(fd, (void *) buf, maxbytes)) < 0) {
		if (errno == EINTR)
			return -1;
		if ((errno == EAGAIN) || 
		    (errno == EWOULDBLOCK))
			return -1;
		error("readx fd %d: %m", fd, n);
		return -1; /* shutdown socket, cleanup. */
	}
	return n;
}	


/*
 * io_node_fail - Some nodes have failed.  Identify affected I/O streams.
 * Flag them as done and signal the I/O thread.
 */
extern int 
io_node_fail(char *nodelist, srun_job_t *job)
{
	hostlist_t fail_list = hostlist_create(nodelist);
	char *node_name;
	int node_inx;

	if (!fail_list) {
		error("Invalid node list `%s' specified", nodelist);
		return SLURM_ERROR;
 	}

	while ( (node_name = hostlist_shift(fail_list)) ) {
		for (node_inx=0; node_inx<job->nhosts; node_inx++) {
			if (strcmp(node_name, job->host[node_inx]))
				continue;
			break;
		}
	}

	eio_handle_signal(job->eio);
	hostlist_destroy(fail_list);
	return SLURM_SUCCESS;
}

static int
_wid(int n)
{
	int width = 1;
	n--;    /* For zero origin */
	while (n /= 10)
		width++;
	return width;
}
