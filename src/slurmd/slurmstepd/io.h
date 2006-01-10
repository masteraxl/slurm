/*****************************************************************************\
 * src/slurmd/slurmstepd/io.h - slurmstepd standard IO routines
 * $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <mgrondona@llnl.gov>.
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

#ifndef _IO_H
#define _IO_H

#include "src/common/eio.h"

#include "src/slurmd/slurmstepd/slurmstepd_job.h"

/*
 * The message cache uses up free message buffers, so STDIO_MAX_MSG_CACHE
 * must be a number smaller than STDIO_MAX_FREE_BUF.
 */
#define STDIO_MAX_FREE_BUF 1024
#define STDIO_MAX_MSG_CACHE 128

struct io_buf {
	int ref_count;
	uint32_t length;
	void *data;
};

struct io_buf *alloc_io_buf(void);
void free_io_buf(struct io_buf *buf);

/* 
 * Create a TCP connection back the initial client (e.g. srun).
 *
 * Since this is the first client connection and the IO engine has not
 * yet started, we initialize the msg_queue as an empty list and
 * directly add the eio_obj_t to the eio handle with eio_new_initial_handle.
 */
int io_initial_client_connect(srun_info_t *srun, slurmd_job_t *job);

/* 
 * Initiate a TCP connection back to a waiting client (e.g. srun).
 *
 * Create a new eio client object and wake up the eio engine so that
 * it can see the new object.
 */
int io_client_connect(srun_info_t *srun, slurmd_job_t *job);

/*
 * Initialize each task's standard I/O file descriptors.  The file descriptors
 * may be files, or may be the end of a pipe which is handled by an eio_obj_t.
 */
int io_init_tasks_stdio(slurmd_job_t *job);

/*
 * Start IO handling thread.
 * Initializes IO pipes, creates IO objects and appends them to job->objs,
 * and opens 2*ntask initial connections for stdout/err, also appending these
 * to job->objs list.
 */
int io_thread_start(slurmd_job_t *job);

int io_dup_stdio(slurmd_task_info_t *t);

/*
 *  Close the tasks' ends of the stdio pipes.
 *  Presumably the tasks have already been started, and
 *  have their copies of these file descriptors.
 */
void io_close_task_fds(slurmd_job_t *job);

void io_close_all(slurmd_job_t *job);

#endif /* !_IO_H */
