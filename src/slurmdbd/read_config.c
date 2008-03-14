/*****************************************************************************\
 *  read_config.c - functions for reading slurmdbd.conf
 *****************************************************************************
 *  Copyright (C) 2003-2007 The Regents of the University of California.
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>
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

#include <pwd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <slurm/slurm_errno.h>

#include "src/common/macros.h"
#include "src/common/log.h"
#include "src/common/parse_config.h"
#include "src/common/read_config.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/slurmdbd/read_config.h"

/* Global variables */
pthread_mutex_t conf_mutex = PTHREAD_MUTEX_INITIALIZER;
//slurm_dbd_conf_t *slurmdbd_conf = NULL;

/* Local functions */
static void _clear_slurmdbd_conf(void);
static char * _get_conf_path(void);

/*
 * free_slurmdbd_conf - free storage associated with the global variable 
 *	slurmdbd_conf
 */
extern void free_slurmdbd_conf(void)
{
	slurm_mutex_lock(&conf_mutex);
	_clear_slurmdbd_conf();
	xfree(slurmdbd_conf);
	slurm_mutex_unlock(&conf_mutex);
}

static void _clear_slurmdbd_conf(void)
{
	if (slurmdbd_conf) {
		xfree(slurmdbd_conf->auth_info);
		xfree(slurmdbd_conf->auth_type);
		xfree(slurmdbd_conf->dbd_addr);
		xfree(slurmdbd_conf->dbd_host);
		slurmdbd_conf->dbd_port = 0;
		xfree(slurmdbd_conf->log_file);
		xfree(slurmdbd_conf->pid_file);
		xfree(slurmdbd_conf->plugindir);
		xfree(slurmdbd_conf->slurm_user_name);
		xfree(slurmdbd_conf->storage_host);
		xfree(slurmdbd_conf->storage_loc);
		xfree(slurmdbd_conf->storage_pass);
		slurmdbd_conf->storage_port = 0;
		xfree(slurmdbd_conf->storage_type);
		xfree(slurmdbd_conf->storage_user);
	}
}

/*
 * read_slurmdbd_conf - load the SlurmDBD configuration from the slurmdbd.conf  
 *	file. Store result into global variable slurmdbd_conf. 
 *	This function can be called more than once.
 * RET SLURM_SUCCESS if no error, otherwise an error code
 */
extern int read_slurmdbd_conf(void)
{
	s_p_options_t options[] = {
		{"AuthInfo", S_P_STRING},
		{"AuthType", S_P_STRING},
		{"DbdAddr", S_P_STRING},
		{"DbdHost", S_P_STRING},
		{"DbdPort", S_P_UINT16},
		{"DebugLevel", S_P_UINT16},
		{"LogFile", S_P_STRING},
		{"MessageTimeout", S_P_UINT16},
		{"PidFile", S_P_STRING},
		{"PluginDir", S_P_STRING},
		{"SlurmUser", S_P_STRING},
		{"StorageHost", S_P_STRING},
		{"StorageLoc", S_P_STRING},
		{"StoragePass", S_P_STRING},
		{"StoragePort", S_P_UINT16},
		{"StorageType", S_P_STRING},
		{"StorageUser", S_P_STRING},
		{NULL} };
	s_p_hashtbl_t *tbl;
	char *conf_path;
	struct stat buf;

	/* Set initial values */
	slurm_mutex_lock(&conf_mutex);
	if (slurmdbd_conf == NULL)
		slurmdbd_conf = xmalloc(sizeof(slurm_dbd_conf_t));
	slurmdbd_conf->debug_level = LOG_LEVEL_INFO;
	_clear_slurmdbd_conf();

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

		s_p_get_string(&slurmdbd_conf->auth_info, "AuthInfo", tbl);
		s_p_get_string(&slurmdbd_conf->auth_type, "AuthType", tbl);
		s_p_get_string(&slurmdbd_conf->dbd_host, "DbdHost", tbl);
		s_p_get_string(&slurmdbd_conf->dbd_addr, "DbdAddr", tbl);
		s_p_get_uint16(&slurmdbd_conf->dbd_port, "DbdPort", tbl);
		s_p_get_uint16(&slurmdbd_conf->debug_level, "DebugLevel", tbl);
		s_p_get_string(&slurmdbd_conf->log_file, "LogFile", tbl);
		if (!s_p_get_uint16(&slurmdbd_conf->msg_timeout,
				    "MessageTimeout", tbl))
			slurmdbd_conf->msg_timeout = DEFAULT_MSG_TIMEOUT;
		else if (slurmdbd_conf->msg_timeout > 100) {
			info("WARNING: MessageTimeout is too high for "
			     "effective fault-tolerance");
		}
		s_p_get_string(&slurmdbd_conf->pid_file, "PidFile", tbl);
		s_p_get_string(&slurmdbd_conf->plugindir, "PluginDir", tbl);
		s_p_get_string(&slurmdbd_conf->slurm_user_name, "SlurmUser",
			       tbl);
		s_p_get_string(&slurmdbd_conf->storage_host,
				"StorageHost", tbl);
		s_p_get_string(&slurmdbd_conf->storage_loc,
				"StorageLoc", tbl);
		s_p_get_string(&slurmdbd_conf->storage_pass,
				"StoragePass", tbl);
		s_p_get_uint16(&slurmdbd_conf->storage_port,
			       "StoragePort", tbl);
		s_p_get_string(&slurmdbd_conf->storage_type,
			       "StorageType", tbl);
		s_p_get_string(&slurmdbd_conf->storage_user,
				"StorageUser", tbl);

		s_p_hashtbl_destroy(tbl);
	}

	xfree(conf_path);
	if (slurmdbd_conf->auth_type == NULL)
		slurmdbd_conf->auth_type = xstrdup(DEFAULT_SLURMDBD_AUTHTYPE);
	if (slurmdbd_conf->dbd_host == NULL) {
		error("slurmdbd.conf lacks DbdHost parameter, using 'localhost'");
		slurmdbd_conf->dbd_host = xstrdup("localhost");
	}
	if (slurmdbd_conf->dbd_addr == NULL)
		slurmdbd_conf->dbd_addr = xstrdup(slurmdbd_conf->dbd_host);
	if (slurmdbd_conf->pid_file == NULL)
		slurmdbd_conf->pid_file = xstrdup(DEFAULT_SLURMDBD_PIDFILE);
	if (slurmdbd_conf->dbd_port == 0)
		slurmdbd_conf->dbd_port = SLURMDBD_PORT;
	if(slurmdbd_conf->plugindir == NULL)
		slurmdbd_conf->plugindir = xstrdup(default_plugin_path);
	if (slurmdbd_conf->slurm_user_name) {
		struct passwd *slurm_passwd;
		slurm_passwd = getpwnam(slurmdbd_conf->slurm_user_name);
		if (slurm_passwd == NULL) {
			fatal("Invalid user for SlurmUser %s, ignored",
			      slurmdbd_conf->slurm_user_name);
		} else
			slurmdbd_conf->slurm_user_id = slurm_passwd->pw_uid;
	} else {
		slurmdbd_conf->slurm_user_name = xstrdup("root");
		slurmdbd_conf->slurm_user_id = 0;
	}
	if (slurmdbd_conf->storage_type == NULL)
		fatal("StorageType must be specified");
				
	slurm_mutex_unlock(&conf_mutex);
	return SLURM_SUCCESS;
}

/* Log the current configuration using verbose() */
extern void log_config(void)
{
	debug2("AuthInfo          = %s", slurmdbd_conf->auth_info);
	debug2("AuthType          = %s", slurmdbd_conf->auth_type);
	debug2("DbdAddr           = %s", slurmdbd_conf->dbd_addr);
	debug2("DbdHost           = %s", slurmdbd_conf->dbd_host);
	debug2("DbdPort           = %u", slurmdbd_conf->dbd_port);
	debug2("DebugLevel        = %u", slurmdbd_conf->debug_level);
	debug2("LogFile           = %s", slurmdbd_conf->log_file);
	debug2("MessageTimeout    = %u", slurmdbd_conf->msg_timeout);
	debug2("PidFile           = %s", slurmdbd_conf->pid_file);
	debug2("PluginDir         = %s", slurmdbd_conf->plugindir);
	debug2("SlurmUser         = %s(%u)", 
		slurmdbd_conf->slurm_user_name, slurmdbd_conf->slurm_user_id); 
	debug2("StorageHost       = %s", slurmdbd_conf->storage_host);
	debug2("StorageLoc        = %s", slurmdbd_conf->storage_loc);
	debug2("StoragePass       = %s", slurmdbd_conf->storage_pass);
	debug2("StoragePort       = %u", slurmdbd_conf->storage_port);
	debug2("StorageType       = %s", slurmdbd_conf->storage_type);
	debug2("StorageUser       = %s", slurmdbd_conf->storage_user);
}

/* Return the DbdPort value */
extern uint16_t get_dbd_port(void)
{
	uint16_t port;

	slurm_mutex_lock(&conf_mutex);
	port = slurmdbd_conf->dbd_port;
	slurm_mutex_unlock(&conf_mutex);
	return port;
}

extern void slurmdbd_conf_lock(void)
{
	slurm_mutex_lock(&conf_mutex);
}

extern void slurmdbd_conf_unlock(void)
{
	slurm_mutex_unlock(&conf_mutex);
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
