/*****************************************************************************\
 *  read_config.c - read the overall slurm configuration file
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov>.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include <slurm/slurm.h>

#include "src/common/hostlist.h"
#include "src/common/slurm_protocol_defs.h"
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/parse_spec.h"
#include "src/common/read_config.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/slurm_rlimits_info.h"
#include "src/common/parse_config.h"
#include "src/common/slurm_config.h"

static pthread_mutex_t conf_lock = PTHREAD_MUTEX_INITIALIZER;
static s_p_hashtbl_t *conf_hashtbl;
static slurm_ctl_conf_t *conf_ptr;
static bool conf_initialized = false;

static int defunct_option(void **dest, slurm_parser_enum_t type,
			  const char *key, const char *value,
			  const char *line);
static void defunct_destroy(void *ptr);
static void validate_and_set_defaults(slurm_ctl_conf_t *conf,
				      s_p_hashtbl_t *hashtbl);

s_p_options_t slurm_conf_options[] = {
	{"AuthType", S_P_STRING},
	{"CheckpointType", S_P_STRING},
	{"CacheGroups", S_P_UINT16},
	{"BackupAddr", S_P_STRING},
	{"BackupController", S_P_STRING},
	{"ControlAddr", S_P_STRING},
	{"ControlMachine", S_P_STRING},
	{"Epilog", S_P_STRING},
	{"FastSchedule", S_P_UINT16},
	{"FirstJobId", S_P_UINT32},
	{"HashBase", S_P_LONG, defunct_option, defunct_destroy},
	{"HeartbeatInterval", S_P_LONG, defunct_option, defunct_destroy},
	{"InactiveLimit", S_P_UINT16},
	{"JobAcctloc", S_P_STRING},
	{"JobAcctParameters", S_P_STRING},
	{"JobAcctType", S_P_STRING},
	{"JobCompLoc", S_P_STRING},
	{"JobCompType", S_P_STRING},
	{"JobCredentialPrivateKey", S_P_STRING},
	{"JobCredentialPublicCertificate", S_P_STRING},
	{"KillTree", S_P_UINT16, defunct_option, defunct_destroy},
	{"KillWait", S_P_UINT16},
	{"MaxJobCount", S_P_UINT16},
	{"MinJobAge", S_P_UINT16},
	{"MpichGmDirectSupport", S_P_LONG},
	{"MpiDefault", S_P_STRING},
	{"NodeName", S_P_ARRAY,
	 parse_nodename, destroy_nodename},
	{"PartitionName", S_P_ARRAY,
	 parse_partitionname, destroy_partitionname},
	{"PluginDir", S_P_STRING},
	{"ProctrackType", S_P_STRING},
	{"Prolog", S_P_STRING},
	{"PropagateResourceLimitsExcept", S_P_STRING},
	{"PropagateResourceLimits", S_P_STRING},
	{"ReturnToService", S_P_UINT16},
	{"SchedulerAuth", S_P_STRING},
	{"SchedulerPort", S_P_UINT16},
	{"SchedulerRootFilter", S_P_UINT16},
	{"SchedulerType", S_P_STRING},
	{"SelectType", S_P_STRING},
	{"SlurmUser", S_P_STRING},
	{"SlurmctldDebug", S_P_UINT16},
	{"SlurmctldLogFile", S_P_STRING},
	{"SlurmctldPidFile", S_P_STRING},
	{"SlurmctldPort", S_P_UINT32},
	{"SlurmctldTimeout", S_P_UINT16},
	{"SlurmdDebug", S_P_UINT16},
	{"SlurmdLogFile", S_P_STRING},
	{"SlurmdPidFile",  S_P_STRING},
	{"SlurmdPort", S_P_UINT32},
	{"SlurmdSpoolDir", S_P_STRING},
	{"SlurmdTimeout", S_P_UINT16},
	{"SrunEpilog", S_P_STRING},
	{"SrunProlog", S_P_STRING},
	{"StateSaveLocation", S_P_STRING},
	{"SwitchType", S_P_STRING},
	{"TaskEpilog", S_P_STRING},
	{"TaskProlog", S_P_STRING},
	{"TaskPlugin", S_P_STRING},
	{"TmpFS", S_P_STRING},
	{"TreeWidth", S_P_UINT16},
	{"WaitTime", S_P_UINT16},
	{NULL}
};

s_p_options_t slurm_nodename_options[] = {
	{"NodeName", S_P_STRING},
	{"NodeHostname", S_P_STRING},
	{"NodeAddr", S_P_STRING},
	{"Feature", S_P_STRING},
	{"Port", S_P_LONG},
	{"Procs", S_P_LONG},
	{"RealMemory", S_P_LONG},
	{"Reason", S_P_STRING},
	{"State", S_P_STRING},
	{"TmpDisk", S_P_LONG},
	{"Weight", S_P_LONG},
	{NULL}
};

s_p_options_t slurm_partition_options[] = {
	{"PartitionName", S_P_STRING},
	{"AllowGroups", S_P_STRING},
	{"Default", S_P_STRING},
	{"Hidden", S_P_STRING},
	{"RootOnly", S_P_STRING},
	{"MaxTime", S_P_STRING},
	{"MaxNodes", S_P_LONG},
	{"MinNodes", S_P_LONG},
	{"Nodes", S_P_STRING},
	{"Shared", S_P_STRING},
	{"State", S_P_STRING},
	{NULL}
};

inline static void _normalize_debug_level(uint16_t *level);

typedef struct names_ll_s {
	char *node_hostname;
	char *node_name;
	struct names_ll_s *next;
} names_ll_t;
bool all_slurmd_hosts = false;
#define NAME_HASH_LEN 512
static names_ll_t *host_to_node_hashtbl[NAME_HASH_LEN] = {NULL};
static names_ll_t *node_to_host_hashtbl[NAME_HASH_LEN] = {NULL};
static char *this_hostname = NULL;

static int defunct_option(void **dest, slurm_parser_enum_t type,
			  const char *key, const char *value,
			  const char *line)
{
	error("The option \"%s\" is defunct, see man slurm.conf.", key);
	return -1;
}

static void defunct_destroy(void *ptr)
{
	/* do nothing */
}


static void _free_name_hashtbl()
{
	int i;
	names_ll_t *p, *q;

	for (i=0; i<NAME_HASH_LEN; i++) {
		p = host_to_node_hashtbl[i];
		while (p) {
			xfree(p->node_hostname);
			xfree(p->node_name);
			q = p->next;
			xfree(p);
			p = q;
		}
		host_to_node_hashtbl[i] = NULL;
		p = node_to_host_hashtbl[i];
		while (p) {
			xfree(p->node_hostname);
			xfree(p->node_name);
			q = p->next;
			xfree(p);
			p = q;
		}
		node_to_host_hashtbl[i] = NULL;
	}
	xfree(this_hostname);
}

static void _init_name_hashtbl()
{
	return;
}

static int _get_hash_idx(char *s)
{
	int i;

	i = 0;
	while (*s) i += (int)*s++;
	return i % NAME_HASH_LEN;
}

static void _push_to_hashtbl(char *node, char *host)
{
	int idx;
	names_ll_t *p, *new;
	char *hh;

	hh = host ? host : node;
	idx = _get_hash_idx(hh);
#ifndef HAVE_FRONT_END		/* Operate only on front-end */
	p = host_to_node_hashtbl[idx];
	while (p) {
		if (strcmp(p->node_hostname, hh)==0) {
			fatal("Duplicated NodeHostname %s in the config file",
				hh);
			return;
		}
		p = p->next;
	}
#endif
	new = (names_ll_t *)xmalloc(sizeof(*new));
	new->node_hostname = xstrdup(hh);
	new->node_name = xstrdup(node);
	new->next = host_to_node_hashtbl[idx];
	host_to_node_hashtbl[idx] = new;

	idx = _get_hash_idx(node);
	p = node_to_host_hashtbl[idx];
	while (p) {
		if (strcmp(p->node_name, node)==0) {
			fatal("Duplicated NodeName %s in the config file",
				node);
			return;
		}
		p = p->next;
	}
	new = (names_ll_t *)xmalloc(sizeof(*new));
	new->node_name = xstrdup(node);
	new->node_hostname = xstrdup(hh);
	new->next = node_to_host_hashtbl[idx];
	node_to_host_hashtbl[idx] = new;
}

/*
 * Register the given NodeName in the alias table.
 * If node_hostname is NULL, only node_name will be used and 
 * no lookup table record is created.
 */
static void _register_conf_node_aliases(char *node_name, char *node_hostname)
{
	hostlist_t node_list = NULL, host_list = NULL;
	char *hn = NULL, *nn;

	if (node_name == NULL || *node_name == '\0')
		return;
	if (strcasecmp(node_name, "DEFAULT") == 0) {
		if (node_hostname) {
			fatal("NodeHostname for NodeName=DEFAULT is illegal");
		}
		return;
	}
	if (!this_hostname) {
		this_hostname = xmalloc(MAX_NAME_LEN);
		getnodename(this_hostname, MAX_NAME_LEN);
	}
	if (strcasecmp(node_name, "localhost") == 0)
		node_name = this_hostname;
	if (node_hostname == NULL)
		node_hostname = node_name;
	if (strcasecmp(node_hostname, "localhost") == 0)
		node_hostname = this_hostname;

	node_list = hostlist_create(node_name);
#ifdef HAVE_FRONT_END	/* Common NodeHostname for all NodeName values */
	/* Expect one common node_hostname for all back-end nodes */
	hn = node_hostname;
#else
	host_list = hostlist_create(node_hostname);
	if (hostlist_count(node_list) != hostlist_count(host_list))
		fatal("NodeName and NodeHostname have different "
			"number of records");
#endif
	while ((nn = hostlist_shift(node_list))) {
		if (host_list)
			hn = hostlist_shift(host_list);
		_push_to_hashtbl(nn, hn);
		if (host_list)
			free(hn);
		free(nn);
	}
	hostlist_destroy(node_list);
	if (host_list)
		hostlist_destroy(host_list);

	return;
}

/*
 * get_conf_node_hostname - Return the NodeHostname for given NodeName
 */
extern char *get_conf_node_hostname(char *node_name)
{
	int idx;
	names_ll_t *p;

	idx = _get_hash_idx(node_name);
	p = node_to_host_hashtbl[idx];
	while (p) {
		if (strcmp(p->node_name, node_name) == 0) {
			return xstrdup(p->node_hostname);
		}
		p = p->next;
	}

	if (all_slurmd_hosts)
		return NULL;
	else {
		 /* Assume identical if we didn't explicitly save all pairs */
		return xstrdup(node_name);
	}
}

/*
 * get_conf_node_name - Return the NodeName for given NodeHostname
 */
extern char *get_conf_node_name(char *node_hostname)
{
	int idx;
	names_ll_t *p;

	idx = _get_hash_idx(node_hostname);
	p = host_to_node_hashtbl[idx];
	while (p) {
		if (strcmp(p->node_hostname, node_hostname) == 0) {
			return xstrdup(p->node_name);
		}
		p = p->next;
	}

	if (all_slurmd_hosts)
		return NULL;
	else {
		/* Assume identical if we didn't explicitly save all pairs */
		return xstrdup(node_hostname);
	}
}




/* getnodename - equivalent to gethostname, but return only the first 
 * component of the fully qualified name 
 * (e.g. "linux123.foo.bar" becomes "linux123") 
 * OUT name
 */
int
getnodename (char *name, size_t len)
{
	int error_code, name_len;
	char *dot_ptr, path_name[1024];

	error_code = gethostname (path_name, sizeof(path_name));
	if (error_code)
		return error_code;

	dot_ptr = strchr (path_name, '.');
	if (dot_ptr == NULL)
		dot_ptr = path_name + strlen(path_name);
	else
		dot_ptr[0] = '\0';

	name_len = (dot_ptr - path_name);
	if (name_len > len)
		return ENAMETOOLONG;

	strcpy (name, path_name);
	return 0;
}

/* 
 * free_slurm_conf - free all storage associated with a slurm_ctl_conf_t.   
 * IN/OUT ctl_conf_ptr - pointer to data structure to be freed
 */
void
free_slurm_conf (slurm_ctl_conf_t *ctl_conf_ptr)
{
	xfree (ctl_conf_ptr->authtype);
	xfree (ctl_conf_ptr->checkpoint_type);
	xfree (ctl_conf_ptr->backup_addr);
	xfree (ctl_conf_ptr->backup_controller);
	xfree (ctl_conf_ptr->control_addr);
	xfree (ctl_conf_ptr->control_machine);
	xfree (ctl_conf_ptr->epilog);
	xfree (ctl_conf_ptr->job_acct_loc);
	xfree (ctl_conf_ptr->job_acct_parameters);
	xfree (ctl_conf_ptr->job_acct_type);
	xfree (ctl_conf_ptr->job_comp_loc);
	xfree (ctl_conf_ptr->job_comp_type);
	xfree (ctl_conf_ptr->job_credential_private_key);
	xfree (ctl_conf_ptr->job_credential_public_certificate);
	xfree (ctl_conf_ptr->mpi_default);
	xfree (ctl_conf_ptr->plugindir);
	xfree (ctl_conf_ptr->proctrack_type);
	xfree (ctl_conf_ptr->prolog);
	xfree (ctl_conf_ptr->propagate_rlimits_except);
	xfree (ctl_conf_ptr->propagate_rlimits);
	xfree (ctl_conf_ptr->schedauth);
	xfree (ctl_conf_ptr->schedtype);
	xfree (ctl_conf_ptr->select_type);
	xfree (ctl_conf_ptr->slurm_conf);
	xfree (ctl_conf_ptr->slurm_user_name);
	xfree (ctl_conf_ptr->slurmctld_logfile);
	xfree (ctl_conf_ptr->slurmctld_pidfile);
	xfree (ctl_conf_ptr->slurmd_logfile);
	xfree (ctl_conf_ptr->slurmd_pidfile);
	xfree (ctl_conf_ptr->slurmd_spooldir);
	xfree (ctl_conf_ptr->state_save_location);
	xfree (ctl_conf_ptr->switch_type);
	xfree (ctl_conf_ptr->tmp_fs);
	xfree (ctl_conf_ptr->task_epilog);
	xfree (ctl_conf_ptr->task_prolog);
	xfree (ctl_conf_ptr->task_plugin);
	xfree (ctl_conf_ptr->tmp_fs);
	xfree (ctl_conf_ptr->srun_prolog);
	xfree (ctl_conf_ptr->srun_epilog);
	xfree (ctl_conf_ptr->node_prefix);
	
	_free_name_hashtbl();
}

/* 
 * init_slurm_conf - initialize or re-initialize the slurm configuration 
 *	values to defaults (NULL or NO_VAL). Note that the configuration
 *	file pathname (slurm_conf) is not changed.    
 * IN/OUT ctl_conf_ptr - pointer to data structure to be initialized
 */
void
init_slurm_conf (slurm_ctl_conf_t *ctl_conf_ptr)
{
	ctl_conf_ptr->last_update		= time(NULL);
	xfree (ctl_conf_ptr->authtype);
	ctl_conf_ptr->cache_groups		= (uint16_t) NO_VAL;
	xfree (ctl_conf_ptr->checkpoint_type);
	xfree (ctl_conf_ptr->backup_addr);
	xfree (ctl_conf_ptr->backup_controller);
	xfree (ctl_conf_ptr->control_addr);
	xfree (ctl_conf_ptr->control_machine);
	xfree (ctl_conf_ptr->epilog);
	ctl_conf_ptr->fast_schedule		= (uint16_t) NO_VAL;
	ctl_conf_ptr->first_job_id		= (uint32_t) NO_VAL;
	ctl_conf_ptr->inactive_limit		= (uint16_t) NO_VAL;
	xfree (ctl_conf_ptr->job_acct_loc);
	xfree (ctl_conf_ptr->job_acct_parameters);
	xfree (ctl_conf_ptr->job_acct_type);
	xfree (ctl_conf_ptr->job_comp_loc);
	xfree (ctl_conf_ptr->job_comp_type);
	xfree (ctl_conf_ptr->job_credential_private_key);
	xfree (ctl_conf_ptr->job_credential_public_certificate);
	ctl_conf_ptr->kill_wait			= (uint16_t) NO_VAL;
	ctl_conf_ptr->max_job_cnt		= (uint16_t) NO_VAL;
	ctl_conf_ptr->min_job_age		= (uint16_t) NO_VAL;
	xfree (ctl_conf_ptr->mpi_default);
	xfree (ctl_conf_ptr->plugindir);
	xfree (ctl_conf_ptr->proctrack_type);
	xfree (ctl_conf_ptr->prolog);
	xfree (ctl_conf_ptr->propagate_rlimits_except);
	xfree (ctl_conf_ptr->propagate_rlimits);
	ctl_conf_ptr->ret2service		= (uint16_t) NO_VAL;
	xfree( ctl_conf_ptr->schedauth );
	ctl_conf_ptr->schedport			= (uint16_t) NO_VAL;
	ctl_conf_ptr->schedrootfltr		= (uint16_t) NO_VAL;
	xfree( ctl_conf_ptr->schedtype );
	xfree( ctl_conf_ptr->select_type );
	ctl_conf_ptr->slurm_user_id		= (uint16_t) NO_VAL; 
	xfree (ctl_conf_ptr->slurm_user_name);
	ctl_conf_ptr->slurmctld_debug		= (uint16_t) NO_VAL; 
	xfree (ctl_conf_ptr->slurmctld_logfile);
	xfree (ctl_conf_ptr->slurmctld_pidfile);
	ctl_conf_ptr->slurmctld_port		= (uint32_t) NO_VAL;
	ctl_conf_ptr->slurmctld_timeout		= (uint16_t) NO_VAL;
	ctl_conf_ptr->slurmd_debug		= (uint16_t) NO_VAL; 
	xfree (ctl_conf_ptr->slurmd_logfile);
	xfree (ctl_conf_ptr->slurmd_pidfile);
 	ctl_conf_ptr->slurmd_port		= (uint32_t) NO_VAL;
	xfree (ctl_conf_ptr->slurmd_spooldir);
	ctl_conf_ptr->slurmd_timeout		= (uint16_t) NO_VAL;
	xfree (ctl_conf_ptr->state_save_location);
	xfree (ctl_conf_ptr->switch_type);
	xfree (ctl_conf_ptr->task_epilog);
	xfree (ctl_conf_ptr->task_prolog);
	xfree (ctl_conf_ptr->task_plugin);
	xfree (ctl_conf_ptr->tmp_fs);
	ctl_conf_ptr->wait_time			= (uint16_t) NO_VAL;
	xfree (ctl_conf_ptr->srun_prolog);
	xfree (ctl_conf_ptr->srun_epilog);
	xfree (ctl_conf_ptr->node_prefix);
	ctl_conf_ptr->tree_width       		= (uint16_t) NO_VAL;
	
	_free_name_hashtbl();
	_init_name_hashtbl();

	return;
}

/*
 * read_slurm_conf_ctl - load the slurm configuration from the configured 
 *	file. 
 * OUT ctl_conf_ptr - pointer to data structure to be filled
 * IN  slurmd_hosts - if true then build a list of hosts on which slurmd runs
 *	(only useful for "scontrol show daemons" command). Otherwise only
 *	record nodes in which NodeName and NodeHostname differ.
 * RET 0 if no error, otherwise an error code
 */
extern int 
read_slurm_conf_ctl (slurm_ctl_conf_t *ctl_conf_ptr, bool slurmd_hosts) 
{
	s_p_hashtbl_t *hashtbl;

	assert(ctl_conf_ptr);
	/* zero the conf structure */
	init_slurm_conf (ctl_conf_ptr);
	/* memset(ctl_conf_ptr, 0, sizeof(slurm_ctl_conf_t)); */

	if (ctl_conf_ptr->slurm_conf == NULL) {
		char *val = getenv("SLURM_CONF");

		if (val == NULL) {
			val = SLURM_CONFIG_FILE;
		}
		ctl_conf_ptr->slurm_conf = xstrdup (val);
	}

	hashtbl = s_p_hashtbl_create(slurm_conf_options);
	s_p_parse_file(hashtbl, ctl_conf_ptr->slurm_conf);

	s_p_dump_values(hashtbl, slurm_conf_options);

	validate_and_set_defaults(ctl_conf_ptr, hashtbl);
	s_p_hashtbl_destroy(hashtbl);
	
	return SLURM_SUCCESS;
}


/* caller must lock conf_lock */
static void
_init_slurm_conf(char *file_name)
{
	conf_ptr = (slurm_ctl_conf_t *)xmalloc(sizeof(slurm_ctl_conf_t));

	if (file_name == NULL) {
		file_name = getenv("SLURM_CONF");
		if (file_name == NULL)
			file_name = SLURM_CONFIG_FILE;
	}

	conf_hashtbl = s_p_hashtbl_create(slurm_conf_options);
	s_p_parse_file(conf_hashtbl, file_name);
	/* s_p_dump_values(conf_hashtbl, slurm_conf_options); */
	validate_and_set_defaults(conf_ptr, conf_hashtbl);
	conf_ptr->slurm_conf = xstrdup(file_name);
}

/* caller must lock conf_lock */
static void
_destroy_slurm_conf()
{
	s_p_hashtbl_destroy(conf_hashtbl);
	free_slurm_conf(conf_ptr);
	xfree(conf_ptr);
}

/*
 * slurm_conf_init - load the slurm configuration from the configured file.
 * IN file_name - name of the slurm configuration file to be read
 *                If file_name is NULL, then the compiled-in default
 *                file name is used.
 * Note: If the conf structures have already been initialized by a call to
 *       slurm_conf_read, any subsequent calls will do nothing until
 *       slurm_conf_destroy is called.
 * RET 0 if no error, otherwise an error code.
 */
extern int
slurm_conf_init(char *file_name)
{
	pthread_mutex_lock(&conf_lock);

	if (conf_initialized) {
		pthread_mutex_unlock(&conf_lock);
		return SLURM_SUCCESS;
	}

	_init_slurm_conf(file_name);
	conf_initialized = true;

	pthread_mutex_unlock(&conf_lock);
	return SLURM_SUCCESS;
}

extern int 
slurm_conf_reinit(char *file_name)
{
	pthread_mutex_lock(&conf_lock);

	if (conf_initialized) {
		_destroy_slurm_conf();
	}

	_init_slurm_conf(file_name);
	conf_initialized = true;

	pthread_mutex_unlock(&conf_lock);
	return SLURM_SUCCESS;

}

extern int
slurm_conf_destroy(void)
{
	pthread_mutex_lock(&conf_lock);

	if (!conf_initialized) {
		pthread_mutex_unlock(&conf_lock);
		return SLURM_SUCCESS;
	}

	_destroy_slurm_conf();

	conf_initialized = false;
	pthread_mutex_unlock(&conf_lock);

	return SLURM_SUCCESS;
}

extern slurm_ctl_conf_t *
slurm_conf_lock(void)
{
	pthread_mutex_lock(&conf_lock);

	if (!conf_initialized) {
		_init_slurm_conf(NULL);
		conf_initialized = true;
	}

	return conf_ptr;
}

extern void
slurm_conf_unlock(void)
{
	pthread_mutex_unlock(&conf_lock);
}

/* 
 * report_leftover - report any un-parsed (non-whitespace) characters on the
 * configuration input line (we over-write parsed characters with whitespace).
 * IN in_line - what is left of the configuration input line.
 * IN line_num - line number of the configuration file.
 */
void
report_leftover (char *in_line, int line_num)
{
	int i;

	for (i = 0; i < strlen (in_line); i++) {
		if (isspace ((int) in_line[i]) || (in_line[i] == '\n'))
			continue;
		error ("Ignored input on line %d of configuration: %s",
			line_num, &in_line[i]);
		break;
	}
}

/* validate configuration
 *
 * IN/OUT ctl_conf_ptr - a configuration as loaded by read_slurm_conf_ctl
 *
 * NOTE: a backup_controller or control_machine of "localhost" are over-written
 *	with this machine's name.
 * NOTE: if backup_addr is NULL, it is over-written by backup_controller
 * NOTE: if control_addr is NULL, it is over-written by control_machine
 */
void
validate_config (slurm_ctl_conf_t *ctl_conf_ptr)
{
	if ((ctl_conf_ptr->backup_controller != NULL) &&
	    (strcasecmp("localhost", ctl_conf_ptr->backup_controller) == 0)) {
		xfree (ctl_conf_ptr->backup_controller);
		ctl_conf_ptr->backup_controller = xmalloc (MAX_NAME_LEN);
		if ( getnodename (ctl_conf_ptr->backup_controller, 
		                  MAX_NAME_LEN) ) 
			fatal ("getnodename: %m");
	}

	if ((ctl_conf_ptr->backup_addr == NULL) && 
	    (ctl_conf_ptr->backup_controller != NULL))
		ctl_conf_ptr->backup_addr = 
				xstrdup (ctl_conf_ptr->backup_controller);

	if ((ctl_conf_ptr->backup_controller == NULL) && 
	    (ctl_conf_ptr->backup_addr != NULL)) {
		error ("BackupAddr specified without BackupController");
		xfree (ctl_conf_ptr->backup_addr);
	}

	if (ctl_conf_ptr->control_machine == NULL)
		fatal ("validate_config: ControlMachine not specified.");
	else if (strcasecmp("localhost", ctl_conf_ptr->control_machine) == 0) {
		xfree (ctl_conf_ptr->control_machine);
		ctl_conf_ptr->control_machine = xmalloc (MAX_NAME_LEN);
		if ( getnodename (ctl_conf_ptr->control_machine, 
		                  MAX_NAME_LEN) ) 
			fatal ("getnodename: %m");
	}

	if ((ctl_conf_ptr->control_addr == NULL) && 
	    (ctl_conf_ptr->control_machine != NULL))
		ctl_conf_ptr->control_addr = 
				xstrdup (ctl_conf_ptr->control_machine);

	if ((ctl_conf_ptr->backup_controller != NULL) && 
	    (strcmp (ctl_conf_ptr->backup_controller, 
	             ctl_conf_ptr->control_machine) == 0)) {
		error ("ControlMachine and BackupController identical");
		xfree (ctl_conf_ptr->backup_addr);
		xfree (ctl_conf_ptr->backup_controller);
	}

	if (ctl_conf_ptr->job_credential_private_key == NULL)
		fatal ("JobCredentialPrivateKey not set");
	if (ctl_conf_ptr->job_credential_public_certificate == NULL)
		fatal ("JobCredentialPublicCertificate not set");

	if (ctl_conf_ptr->max_job_cnt < 1)
		fatal ("MaxJobCount=%u, No jobs permitted",
		       ctl_conf_ptr->max_job_cnt);

	if (ctl_conf_ptr->authtype == NULL)
		ctl_conf_ptr->authtype = xstrdup(DEFAULT_AUTH_TYPE);

	if (ctl_conf_ptr->cache_groups == (uint16_t) NO_VAL)
		ctl_conf_ptr->cache_groups = DEFAULT_CACHE_GROUPS;

	if (ctl_conf_ptr->checkpoint_type == NULL)
		 ctl_conf_ptr->checkpoint_type = 
			xstrdup(DEFAULT_CHECKPOINT_TYPE);

	if (ctl_conf_ptr->fast_schedule == (uint16_t) NO_VAL)
		ctl_conf_ptr->fast_schedule = DEFAULT_FAST_SCHEDULE;

	if (ctl_conf_ptr->first_job_id == (uint32_t) NO_VAL)
		ctl_conf_ptr->first_job_id = DEFAULT_FIRST_JOB_ID;

	if (ctl_conf_ptr->inactive_limit == (uint16_t) NO_VAL)
		ctl_conf_ptr->inactive_limit = DEFAULT_INACTIVE_LIMIT;

	if (ctl_conf_ptr->job_acct_loc == NULL)
		ctl_conf_ptr->job_acct_loc = xstrdup(DEFAULT_JOB_ACCT_LOC);

	if (ctl_conf_ptr->job_acct_parameters == NULL)
		ctl_conf_ptr->job_acct_parameters =
				xstrdup(DEFAULT_JOB_ACCT_PARAMETERS);

	if (ctl_conf_ptr->job_acct_type == NULL)
		ctl_conf_ptr->job_acct_type = xstrdup(DEFAULT_JOB_ACCT_TYPE);

	if (ctl_conf_ptr->job_comp_type == NULL)
		ctl_conf_ptr->job_comp_type = xstrdup(DEFAULT_JOB_COMP_TYPE);

	if (ctl_conf_ptr->kill_wait == (uint16_t) NO_VAL)
		ctl_conf_ptr->kill_wait = DEFAULT_KILL_WAIT;

	if (ctl_conf_ptr->max_job_cnt == (uint16_t) NO_VAL)
		ctl_conf_ptr->max_job_cnt = DEFAULT_MAX_JOB_COUNT;

	if (ctl_conf_ptr->min_job_age == (uint16_t) NO_VAL)
		ctl_conf_ptr->min_job_age = DEFAULT_MIN_JOB_AGE;

	if (ctl_conf_ptr->mpi_default == NULL)
		ctl_conf_ptr->mpi_default = xstrdup(DEFAULT_MPI_DEFAULT);
	if (ctl_conf_ptr->plugindir == NULL)
		ctl_conf_ptr->plugindir = xstrdup(SLURM_PLUGIN_PATH);

	if (ctl_conf_ptr->switch_type == NULL)
		ctl_conf_ptr->switch_type = xstrdup(DEFAULT_SWITCH_TYPE);

	if (ctl_conf_ptr->proctrack_type == NULL) {
		if (!strcmp(ctl_conf_ptr->switch_type,"switch/elan"))
			ctl_conf_ptr->proctrack_type = 
					xstrdup("proctrack/rms");
		else
			ctl_conf_ptr->proctrack_type = 
					xstrdup(DEFAULT_PROCTRACK_TYPE);
	}
	if ((!strcmp(ctl_conf_ptr->switch_type,   "switch/elan"))
	&&  (!strcmp(ctl_conf_ptr->proctrack_type,"proctrack/linuxproc")))
		fatal("proctrack/linuxproc is incompatable with switch/elan");

        if (ctl_conf_ptr->propagate_rlimits_except) {
                if ((parse_rlimits( ctl_conf_ptr->propagate_rlimits_except,
                                   NO_PROPAGATE_RLIMITS )) < 0)
                        fatal( "Bad PropagateResourceLimitsExcept: %s",
                                ctl_conf_ptr->propagate_rlimits_except );
        }
        else {
                if (ctl_conf_ptr->propagate_rlimits == NULL)
                        ctl_conf_ptr->propagate_rlimits = xstrdup( "ALL" );
                if ((parse_rlimits( ctl_conf_ptr->propagate_rlimits,
                                   PROPAGATE_RLIMITS )) < 0)
                        fatal( "Bad PropagateResourceLimits: %s",
                                ctl_conf_ptr->propagate_rlimits );
        }

	if (ctl_conf_ptr->ret2service == (uint16_t) NO_VAL)
		ctl_conf_ptr->ret2service = DEFAULT_RETURN_TO_SERVICE;

	if (ctl_conf_ptr->schedrootfltr == (uint16_t) NO_VAL)
		ctl_conf_ptr->schedrootfltr = DEFAULT_SCHEDROOTFILTER;

	if (ctl_conf_ptr->schedtype == NULL)
		ctl_conf_ptr->schedtype = xstrdup(DEFAULT_SCHEDTYPE);

	if (ctl_conf_ptr->select_type == NULL)
		ctl_conf_ptr->select_type = xstrdup(DEFAULT_SELECT_TYPE);

	if (ctl_conf_ptr->slurm_user_name == NULL) {
		ctl_conf_ptr->slurm_user_name = xstrdup("root");
		ctl_conf_ptr->slurm_user_id   = 0;
	}

	if (ctl_conf_ptr->slurmctld_debug != (uint16_t) NO_VAL)
		_normalize_debug_level(&ctl_conf_ptr->slurmctld_debug);
	else
		ctl_conf_ptr->slurmctld_debug = LOG_LEVEL_INFO;

	if (ctl_conf_ptr->slurmctld_pidfile == NULL)
		ctl_conf_ptr->slurmctld_pidfile =
			xstrdup(DEFAULT_SLURMCTLD_PIDFILE);

	if (ctl_conf_ptr->slurmctld_port == (uint32_t) NO_VAL) 
		ctl_conf_ptr->slurmctld_port = SLURMCTLD_PORT;

	if (ctl_conf_ptr->slurmctld_timeout == (uint16_t) NO_VAL)
		ctl_conf_ptr->slurmctld_timeout = DEFAULT_SLURMCTLD_TIMEOUT;

	if (ctl_conf_ptr->slurmd_debug != (uint16_t) NO_VAL)
		_normalize_debug_level(&ctl_conf_ptr->slurmd_debug);
	else
		ctl_conf_ptr->slurmd_debug = LOG_LEVEL_INFO;

	if (ctl_conf_ptr->slurmd_pidfile == NULL)
		ctl_conf_ptr->slurmd_pidfile = xstrdup(DEFAULT_SLURMD_PIDFILE);

#ifndef MULTIPLE_SLURMD
	if (ctl_conf_ptr->slurmd_port == (uint32_t) NO_VAL) 
		ctl_conf_ptr->slurmd_port = SLURMD_PORT;
#endif

	if (ctl_conf_ptr->slurmd_spooldir == NULL)
		ctl_conf_ptr->slurmd_spooldir = xstrdup(DEFAULT_SPOOLDIR);

	if (ctl_conf_ptr->slurmd_timeout == (uint16_t) NO_VAL)
		ctl_conf_ptr->slurmd_timeout = DEFAULT_SLURMD_TIMEOUT;

	if (ctl_conf_ptr->state_save_location == NULL)
		ctl_conf_ptr->state_save_location = xstrdup(
				DEFAULT_SAVE_STATE_LOC);

	/* see above for switch_type, order dependent */

	if (ctl_conf_ptr->task_plugin == NULL)
		ctl_conf_ptr->task_plugin = xstrdup(DEFAULT_TASK_PLUGIN);

	if (ctl_conf_ptr->tmp_fs == NULL)
		ctl_conf_ptr->tmp_fs = xstrdup(DEFAULT_TMP_FS);

	if (ctl_conf_ptr->wait_time == (uint16_t) NO_VAL)
		ctl_conf_ptr->wait_time = DEFAULT_WAIT_TIME;
	
	if (ctl_conf_ptr->tree_width == (uint16_t) NO_VAL) 
		ctl_conf_ptr->tree_width = DEFAULT_TREE_WIDTH;

}

/* Normalize supplied debug level to be in range per log.h definitions */
static void _normalize_debug_level(uint16_t *level)
{
	if (*level > LOG_LEVEL_DEBUG3) {
		error("Normalizing debug level from %u to %d", 
		      *level, LOG_LEVEL_DEBUG3);
		*level = LOG_LEVEL_DEBUG3;
	}
	/* level is uint16, always > LOG_LEVEL_QUIET(0), can't underflow */
}

/* 
 *
 * IN/OUT ctl_conf_ptr - a configuration as loaded by read_slurm_conf_ctl
 *
 * NOTE: a backup_controller or control_machine of "localhost" are over-written
 *	with this machine's name.
 * NOTE: if backup_addr is NULL, it is over-written by backup_controller
 * NOTE: if control_addr is NULL, it is over-written by control_machine
 */
static void
validate_and_set_defaults(slurm_ctl_conf_t *conf, s_p_hashtbl_t *hashtbl)
{
	if (s_p_get_string(hashtbl, "BackupController",
			   &conf->backup_controller)
	    && strcasecmp("localhost", conf->backup_controller) == 0) {
		xfree(conf->backup_controller);
		conf->backup_controller = xmalloc (MAX_NAME_LEN);
		if (getnodename(conf->backup_controller, MAX_NAME_LEN)) 
			fatal("getnodename: %m");
	}
	if (s_p_get_string(hashtbl, "BackupAddr", &conf->backup_addr)) {
		if (conf->backup_controller == NULL) {
			error("BackupAddr specified without BackupController");
			xfree(conf->backup_addr);
		}
	} else {
		if (conf->backup_controller != NULL)
			conf->backup_addr = xstrdup(conf->backup_controller);
	}

	if (!s_p_get_string(hashtbl, "ControlMachine", &conf->control_machine))
		fatal ("validate_config: ControlMachine not specified.");
	else if (strcasecmp("localhost", conf->control_machine) == 0) {
		xfree (conf->control_machine);
		conf->control_machine = xmalloc(MAX_NAME_LEN);
		if (getnodename(conf->control_machine, MAX_NAME_LEN))
			fatal("getnodename: %m");
	}

	if (!s_p_get_string(hashtbl, "ControlAddr", &conf->control_addr)
	    && conf->control_machine != NULL)
		conf->control_addr = xstrdup (conf->control_machine);

	if ((conf->backup_controller != NULL)
	    && (strcmp(conf->backup_controller, conf->control_machine) == 0)) {
		error("ControlMachine and BackupController identical");
		xfree(conf->backup_addr);
		xfree(conf->backup_controller);
	}

	if (!s_p_get_string(hashtbl, "JobCredentialPrivateKey",
			    &conf->job_credential_private_key))
		fatal("JobCredentialPrivateKey not set");
	if (!s_p_get_string(hashtbl, "JobCredentialPublicCertificate",
			    &conf->job_credential_public_certificate))
		fatal("JobCredentialPublicCertificate not set");

	if (s_p_get_uint16(hashtbl, "MaxJobCount", &conf->max_job_cnt)
	    && conf->max_job_cnt < 1)
		fatal("MaxJobCount=%u, No jobs permitted", conf->max_job_cnt);

	if (!s_p_get_string(hashtbl, "AuthType", &conf->authtype))
		conf->authtype = xstrdup(DEFAULT_AUTH_TYPE);

	if (!s_p_get_uint16(hashtbl, "CacheGroups", &conf->cache_groups))
		conf->cache_groups = DEFAULT_CACHE_GROUPS;

	if (!s_p_get_string(hashtbl, "CheckpointType", &conf->checkpoint_type))
		conf->checkpoint_type = xstrdup(DEFAULT_CHECKPOINT_TYPE);

	s_p_get_string(hashtbl, "Epilog", &conf->epilog);

	if (!s_p_get_uint16(hashtbl, "FastSchedule", &conf->fast_schedule))
		conf->fast_schedule = DEFAULT_FAST_SCHEDULE;

	if (!s_p_get_uint32(hashtbl, "FirstJobId", &conf->first_job_id))
		conf->first_job_id = DEFAULT_FIRST_JOB_ID;

	if (s_p_get_uint16(hashtbl, "InactiveLimit", &conf->inactive_limit)) {
#ifdef HAVE_BG
		/* Inactive limit must be zero on Blue Gene */
		error("InactiveLimit=%ld is invalid on Blue Gene",
		      cont->inactive_limit);
		conf->inactive_limit = 0; /* default value too */
#endif
	} else {
		conf->inactive_limit = DEFAULT_INACTIVE_LIMIT;
	}

	if (!s_p_get_string(hashtbl, "JobAcctLoc", &conf->job_acct_loc))
		conf->job_acct_loc = xstrdup(DEFAULT_JOB_ACCT_LOC);

	if (!s_p_get_string(hashtbl, "JobAcctParameters",
			    &conf->job_acct_parameters))
		conf->job_acct_parameters =
			xstrdup(DEFAULT_JOB_ACCT_PARAMETERS);

	if (!s_p_get_string(hashtbl, "JobAcctType", &conf->job_acct_type))
		conf->job_acct_type = xstrdup(DEFAULT_JOB_ACCT_TYPE);

	s_p_get_string(hashtbl, "JobCompLoc", &conf->job_comp_loc);

	if (!s_p_get_string(hashtbl, "JobCompType", &conf->job_comp_type))
		conf->job_comp_type = xstrdup(DEFAULT_JOB_COMP_TYPE);

	if (!s_p_get_uint16(hashtbl, "KillWait", &conf->kill_wait))
		conf->kill_wait = DEFAULT_KILL_WAIT;

	if (!s_p_get_uint16(hashtbl, "MaxJobCount", &conf->max_job_cnt))
		conf->max_job_cnt = DEFAULT_MAX_JOB_COUNT;

	if (!s_p_get_uint16(hashtbl, "MinJobAge", &conf->min_job_age))
		conf->min_job_age = DEFAULT_MIN_JOB_AGE;

	if (!s_p_get_string(hashtbl, "MpiDefault", &conf->mpi_default))
		conf->mpi_default = xstrdup(DEFAULT_MPI_DEFAULT);

	if (!s_p_get_string(hashtbl, "PluginDir", &conf->plugindir))
		conf->plugindir = xstrdup(SLURM_PLUGIN_PATH);

	if (!s_p_get_string(hashtbl, "SwitchType", &conf->switch_type))
		conf->switch_type = xstrdup(DEFAULT_SWITCH_TYPE);

	if (!s_p_get_string(hashtbl, "ProctrackType", &conf->proctrack_type)) {
		if (!strcmp(conf->switch_type,"switch/elan"))
			conf->proctrack_type = xstrdup("proctrack/rms");
		else
			conf->proctrack_type = 
				xstrdup(DEFAULT_PROCTRACK_TYPE);
	}
	if ((!strcmp(conf->switch_type,      "switch/elan"))
	    && (!strcmp(conf->proctrack_type,"proctrack/linuxproc")))
		fatal("proctrack/linuxproc is incompatable with switch/elan");

	s_p_get_string(hashtbl, "Prolog", &conf->prolog);

	/* FIXME - figure out how to convert to s_p_get_* */
        if (conf->propagate_rlimits_except) {
                if ((parse_rlimits( conf->propagate_rlimits_except,
                                   NO_PROPAGATE_RLIMITS )) < 0)
                        fatal( "Bad PropagateResourceLimitsExcept: %s",
                                conf->propagate_rlimits_except );
        }
        else {
                if (conf->propagate_rlimits == NULL)
                        conf->propagate_rlimits = xstrdup( "ALL" );
                if ((parse_rlimits( conf->propagate_rlimits,
                                   PROPAGATE_RLIMITS )) < 0)
                        fatal( "Bad PropagateResourceLimits: %s",
                                conf->propagate_rlimits );
        }

	if (!s_p_get_uint16(hashtbl, "ReturnToService", &conf->ret2service))
		conf->ret2service = DEFAULT_RETURN_TO_SERVICE;

	s_p_get_string(hashtbl, "SchedulerAuth", &conf->schedauth);

	if (s_p_get_uint16(hashtbl, "SchedulerPort", &conf->schedport)) {
		if (conf->schedport == 0) {
			error("SchedulerPort=0 is invalid");
			conf->schedport = (uint16_t)NO_VAL;
		}
	}

	if (!s_p_get_uint16(hashtbl, "SchedulerRootFilter",
			    &conf->schedrootfltr))
		conf->schedrootfltr = DEFAULT_SCHEDROOTFILTER;

	if (!s_p_get_string(hashtbl, "SchedulerType", &conf->schedtype))
		conf->schedtype = xstrdup(DEFAULT_SCHEDTYPE);

	if (!s_p_get_string(hashtbl, "SelectType", &conf->select_type))
		conf->select_type = xstrdup(DEFAULT_SELECT_TYPE);

	if (!s_p_get_string(hashtbl, "SlurmUser", &conf->slurm_user_name)) {
		conf->slurm_user_name = xstrdup("root");
		conf->slurm_user_id   = 0;
	} else {
		struct passwd *slurm_passwd;
		slurm_passwd = getpwnam(conf->slurm_user_name);
		if (slurm_passwd == NULL) {
			error ("Invalid user for SlurmUser %s, ignored",
			       conf->slurm_user_name);
			xfree(conf->slurm_user_name);
		} else {
			if (slurm_passwd->pw_uid > 0xffff)
				error("SlurmUser numeric overflow, "
				      "will be fixed soon");
			else
				conf->slurm_user_id = slurm_passwd->pw_uid;
		}
	}

	if (s_p_get_uint16(hashtbl, "SlurmctldDebug", &conf->slurmctld_debug))
		_normalize_debug_level(&conf->slurmctld_debug);
	else
		conf->slurmctld_debug = LOG_LEVEL_INFO;

	if (!s_p_get_string(hashtbl, "SlurmctldPidFile",
			    &conf->slurmctld_pidfile))
		conf->slurmctld_pidfile = xstrdup(DEFAULT_SLURMCTLD_PIDFILE);

	s_p_get_string(hashtbl, "SlurmctldLogFile", &conf->slurmctld_logfile);

	if (!s_p_get_uint32(hashtbl, "SlurmctldPort", &conf->slurmctld_port))
		conf->slurmctld_port = SLURMCTLD_PORT;

	if (!s_p_get_uint16(hashtbl, "SlurmctldTimeout",
			    &conf->slurmctld_timeout))
		conf->slurmctld_timeout = DEFAULT_SLURMCTLD_TIMEOUT;

	if (s_p_get_uint16(hashtbl, "SlurmdDebug", &conf->slurmd_debug))
		_normalize_debug_level(&conf->slurmd_debug);
	else
		conf->slurmd_debug = LOG_LEVEL_INFO;

	s_p_get_string(hashtbl, "SlurmdLogFile", &conf->slurmd_logfile);

	if (!s_p_get_string(hashtbl, "SlurmdPidFile", &conf->slurmd_pidfile))
		conf->slurmd_pidfile = xstrdup(DEFAULT_SLURMD_PIDFILE);

	if (!s_p_get_uint32(hashtbl, "SlurmdPort", &conf->slurmd_port))
		conf->slurmd_port = SLURMD_PORT;

	if (!s_p_get_string(hashtbl, "SlurmdSpoolDir", &conf->slurmd_spooldir))
		conf->slurmd_spooldir = xstrdup(DEFAULT_SPOOLDIR);

	if (!s_p_get_uint16(hashtbl, "SlurmdTimeout", &conf->slurmd_timeout))
		conf->slurmd_timeout = DEFAULT_SLURMD_TIMEOUT;

	s_p_get_string(hashtbl, "SrunProlog", &conf->srun_prolog);
	s_p_get_string(hashtbl, "SrunEpilog", &conf->srun_epilog);

	if (!s_p_get_string(hashtbl, "StateSaveLocation",
			    &conf->state_save_location))
		conf->state_save_location = xstrdup(DEFAULT_SAVE_STATE_LOC);

	/* see above for switch_type, order dependent */

	if (!s_p_get_string(hashtbl, "TaskPlugin", &conf->task_plugin))
		conf->task_plugin = xstrdup(DEFAULT_TASK_PLUGIN);

	s_p_get_string(hashtbl, "TaskEpilog", &conf->task_epilog);
	s_p_get_string(hashtbl, "TaskProlog", &conf->task_prolog);

	if (!s_p_get_string(hashtbl, "TmpFS", &conf->tmp_fs))
		conf->tmp_fs = xstrdup(DEFAULT_TMP_FS);

	if (!s_p_get_uint16(hashtbl, "WaitTime", &conf->wait_time))
		conf->wait_time = DEFAULT_WAIT_TIME;
	
	if (s_p_get_uint16(hashtbl, "TreeWidth", &conf->schedport)) {
		if (conf->tree_width == 0) {
			error("TreeWidth=0 is invalid");
			conf->tree_width = 50; /* default? */
		}
	} else {
		conf->tree_width = DEFAULT_TREE_WIDTH;
	}
}
