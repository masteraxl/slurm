/*****************************************************************************\
 *  jobacct_storage_slurmdbd.c - SlurmDBD slurm job accounting plugin.
 *****************************************************************************
 *  Copyright (C) 2002-2008 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Aubke <da@llnl.gov>.
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

#if HAVE_STDINT_H
#  include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <slurm/slurm_errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "src/common/pack.h"
#include "src/common/slurmdbd_defs.h"
#include "src/common/xstring.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmd/slurmd/slurmd.h"

static char *   slurmdbd_addr = NULL;
static char *   slurmdbd_host = NULL;
static uint16_t slurmdbd_port = 0;
static slurm_fd slurmdbd_fd   = -1;
static pthread_mutex_t slurmdbd_lock = PTHREAD_MUTEX_INITIALIZER;

static void   _close_slurmdbd_fd(void);
static char * _get_conf_path(void);
static void   _open_slurmdbd_fd(void);
static int    _read_slurmdbd_conf(void);

/*
 * These variables are required by the generic plugin interface.  If they
 * are not found in the plugin, the plugin loader will ignore it.
 *
 * plugin_name - a string giving a human-readable description of the
 * plugin.  There is no maximum length, but the symbol must refer to
 * a valid string.
 *
 * plugin_type - a string suggesting the type of the plugin or its
 * applicability to a particular form of data or method of data handling.
 * If the low-level plugin API is used, the contents of this string are
 * unimportant and may be anything.  SLURM uses the higher-level plugin
 * interface which requires this string to be of the form
 *
 *	<application>/<method>
 *
 * where <application> is a description of the intended application of
 * the plugin (e.g., "jobacct" for SLURM job completion logging) and <method>
 * is a description of how this plugin satisfies that application.  SLURM will
 * only load job completion logging plugins if the plugin_type string has a 
 * prefix of "jobacct/".
 *
 * plugin_version - an unsigned 32-bit integer giving the version number
 * of the plugin.  If major and minor revisions are desired, the major
 * version number may be multiplied by a suitable magnitude constant such
 * as 100 or 1000.  Various SLURM versions will likely require a certain
 * minimum versions for their plugins as the job accounting API 
 * matures.
 */
const char plugin_name[] = "Job accounting storage SLURMDBD plugin";
const char plugin_type[] = "jobacct_storage/slurmdbd";
const uint32_t plugin_version = 100;

/*
 * init() is called when the plugin is loaded, before any other functions
 * are called.  Put global initialization here.
 */
extern int init ( void )
{
	slurm_mutex_lock(&slurmdbd_lock);
	if (slurmdbd_fd < 0) {
		/* since this can be loaded from many different places
		   only tell us once. */
		verbose("%s loaded", plugin_name);
		_read_slurmdbd_conf();
		_open_slurmdbd_fd();
	} else {
		debug4("%s loaded", plugin_name);
	}
	slurm_mutex_unlock(&slurmdbd_lock);

	return SLURM_SUCCESS;
}

extern int fini ( void )
{
	slurm_mutex_lock(&slurmdbd_lock);
	_close_slurmdbd_fd();
	xfree(slurmdbd_addr);
	xfree(slurmdbd_host);
	slurm_mutex_unlock(&slurmdbd_lock);
	return SLURM_SUCCESS;
}

/* Open a connection to the Slurm DBD and set slurm_fd */
static void _open_slurmdbd_fd(void)
{
	slurm_addr dbd_addr;

	if (slurmdbd_fd < 0) {
		slurm_set_addr(&dbd_addr, slurmdbd_port, slurmdbd_addr);
		if (dbd_addr.sin_port == 0)
			error("Unable to locate SlurmDBD host %s:%s", 
			      slurmdbd_host, slurmdbd_addr);
		else {
			slurmdbd_fd = slurm_open_msg_conn(&dbd_addr);
			if (slurmdbd_fd < 0)
				error("slurmdbd: slurm_open_msg_conn: %m");
		}
	}
	if (slurmdbd_fd >= 0) {
		/* FIXME: Send an authentication message now */
		write(slurmdbd_fd, "OPEN", 5);
	}
}

/* Send termination message to Slurm DBD and close the connection */
static void _close_slurmdbd_fd(void)
{
	if (slurmdbd_fd >= 0) {
		/* FIXME: Send a termination message now */
		write(slurmdbd_fd, "FINI", 5);
		close(slurmdbd_fd);
		slurmdbd_fd = -1;
	}
}

/* Read the slurmdbd.conf file to get the DbdPort value */
static int _read_slurmdbd_conf(void)
{
	s_p_options_t options[] = {
		{"DbdAddr", S_P_STRING},
		{"DbdHost", S_P_STRING},
		{"DbdPort", S_P_UINT16},
		{"DebugLevel", S_P_UINT16},
		{"LogFile", S_P_STRING},
		{"PidFile", S_P_STRING},
		{"StoragePassword", S_P_STRING},
		{"StorageUser", S_P_STRING},
		{NULL} };
	s_p_hashtbl_t *tbl;
	char *conf_path;
	struct stat buf;

	/* Get the slurmdbd.conf path and validate the file */
	conf_path = _get_conf_path();
	if ((conf_path == NULL) || (stat(conf_path, &buf) == -1)) {
		info("No slurmdbd.conf file (%s)", conf_path);
	} else {
		debug("Reading slurmdbd.conf file %s", conf_path);

		tbl = s_p_hashtbl_create(options);
		if (s_p_parse_file(tbl, conf_path) == SLURM_ERROR) {
			fatal("Could not open/read/parse slurmdbd.conf file %s",
		 	     conf_path);
		}

		if (!s_p_get_string(&slurmdbd_host,"DbdHost", tbl)) {
			error("slurmdbd.conf lacks DbdHost parameter");
			slurmdbd_host = xstrdup("localhost");
		}
		if (!s_p_get_string(&slurmdbd_addr, "DbdAddr", tbl)) {
			slurmdbd_addr = xstrdup(slurmdbd_host);
		}
		if (!s_p_get_uint16(&slurmdbd_port, "DbdPort", tbl))
			slurmdbd_port = SLURMDBD_PORT;

		s_p_hashtbl_destroy(tbl);
	}

	xfree(conf_path);

	return SLURM_SUCCESS;
}

/* Return the pathname of the slurmdbd.conf file.
 * xfree() the value returned */
static char * _get_conf_path(void)
{
	char *val = getenv("SLURM_CONF");
	char *path = NULL;
	int i;

	if (!val)
		val = default_slurm_config_file;

	/* Replace file name on end of path */
	i = strlen(val) + 15;
	path = xmalloc(i);
	strcpy(path, val);
	val = strrchr(path, (int)'/');
	if (val)	/* absolute path */
		val++;
	else		/* not absolute path */
		val = path;
	strcpy(val, "slurmdbd.conf");

	return path;
}

static int _send_msg(Buf buffer)
{
	slurm_mutex_lock(&slurmdbd_lock);
	if (slurmdbd_fd < 0)
		_open_slurmdbd_fd();
	write(slurmdbd_fd, "SEND", 5);
	slurm_mutex_unlock(&slurmdbd_lock);
	return SLURM_SUCCESS;
}

/* 
 * Initialize the storage make sure tables are created and in working
 * order
 */
extern int jobacct_storage_p_init(char *location)
{
	return SLURM_SUCCESS;
}

/*
 * finish up storage connection
 */
extern int jobacct_storage_p_fini()
{
	return SLURM_SUCCESS;
}

/* 
 * load into the storage the start of a job
 */
extern int jobacct_storage_p_job_start(struct job_record *job_ptr)
{
	int rc;
	dbd_job_start_msg_t msg;
	Buf buffer = init_buf(1024);

	msg.job_id  = job_ptr->job_id;
	slurm_dbd_pack_job_start_msg(&msg, buffer);
	rc = _send_msg(buffer);
	free_buf(buffer);
	return rc;
}

/* 
 * load into the storage the end of a job
 */
extern int jobacct_storage_p_job_complete(struct job_record *job_ptr)
{
	int rc;
	dbd_job_comp_msg_t msg;
	Buf buffer = init_buf(1024);

	msg.job_id  = job_ptr->job_id;
	slurm_dbd_pack_job_complete_msg(&msg, buffer);
	rc = _send_msg(buffer);
	free_buf(buffer);
	return rc;
}

/* 
 * load into the storage the start of a job step
 */
extern int jobacct_storage_p_step_start(struct step_record *step_ptr)
{
	int rc;
	dbd_step_start_msg_t msg;
	Buf buffer = init_buf(1024);

	msg.job_id  = step_ptr->job_ptr->job_id;
	msg.step_id = step_ptr->step_id;
	slurm_dbd_pack_step_start_msg(&msg, buffer);
	rc = _send_msg(buffer);
	free_buf(buffer);
	return rc;
}

/* 
 * load into the storage the end of a job step
 */
extern int jobacct_storage_p_step_complete(struct step_record *step_ptr)
{
	int rc;
	dbd_step_comp_msg_t msg;
	Buf buffer = init_buf(1024);

	msg.job_id  = step_ptr->job_ptr->job_id;
	msg.step_id = step_ptr->step_id;
	slurm_dbd_pack_step_complete_msg(&msg, buffer);
	rc = _send_msg(buffer);
	free_buf(buffer);
	return rc;
}

/* 
 * load into the storage a suspention of a job
 */
extern int jobacct_storage_p_suspend(struct job_record *job_ptr)
{
	int rc;
	dbd_job_suspend_msg_t msg;
	Buf buffer = init_buf(1024);

	msg.job_id = job_ptr->job_id;
	slurm_dbd_pack_job_suspend_msg(&msg, buffer);
	rc = _send_msg(buffer);
	free_buf(buffer);
	return rc;
}

/* 
 * get info from the storage 
 * returns List of job_rec_t *
 * note List needs to be freed when called
 */
extern void jobacct_storage_p_get_jobs(List job_list,
					List selected_steps,
					List selected_parts,
					void *params)
{
	dbd_get_jobs_msg_t msg;
	Buf buffer = init_buf(1024);

	msg.job_id = NO_VAL;
	slurm_dbd_pack_get_jobs_msg(&msg, buffer);
	_send_msg(buffer);
	free_buf(buffer);
	return;
}

/* 
 * expire old info from the storage 
 */
extern void jobacct_storage_p_archive(List selected_parts,
				       void *params)
{
	return;
}
