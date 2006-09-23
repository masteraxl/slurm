/*****************************************************************************\
 *  step_launch.h - launch a parallel job step
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher J. Morrone <morrone2@llnl.gov>
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
#ifndef _STEP_LAUNCH_H
#define _STEP_LAUNCH_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <unistd.h>
#include <stdint.h>
#include <pthread.h>

#include <slurm/slurm.h>

#include "src/common/slurm_step_layout.h"
#include "src/common/eio.h"
#include "src/common/bitstring.h"

#include "src/api/step_io.h"

typedef struct {
	int connected;
	int *sockets; /* array of socket file descriptors */
} user_managed_io_t;

struct step_launch_state {
	pthread_mutex_t lock;
	pthread_cond_t cond;
	int tasks_requested;
	bitstr_t *tasks_started; /* or attempted to start, but failed */
	bitstr_t *tasks_exited;  /* or never started correctly */
	bool abort;
	bool abort_action_taken;

	/* message thread variables */
	eio_handle_t *msg_handle;
	pthread_t msg_thread;
	/* set to -1 if slaunch message handler should not attempt to handle */
	int slurmctld_socket_fd;
	uint16_t num_resp_port;
	uint16_t *resp_port; /* array of message response ports */

	/* io variables */
	bool user_managed_io;
	union {
		client_io_t *normal;
		user_managed_io_t *user;
	} io;
	slurm_step_layout_t *layout; /* a pointer into the ctx
					step_resp, do not free */

	/* user registered callbacks */
	slurm_job_step_launch_callbacks_t callback;
};

/*
 * Create a launch state structure for a specified step context, "ctx".
 */
struct step_launch_state * step_launch_state_create(slurm_step_ctx ctx);

/*
 * Free the memory associated with the a launch state structure.
 */
void step_launch_state_destroy(struct step_launch_state *sls);

#endif /* _STEP_LAUNCH_H */
