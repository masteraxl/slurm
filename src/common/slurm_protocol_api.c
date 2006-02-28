/*****************************************************************************\
 *  slurm_protocol_api.c - high-level slurm communication functions
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Kevin Tew <tew1@llnl.gov>, et. al.
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

#if HAVE_CONFIG_H
#  include "config.h"
#endif

/* GLOBAL INCLUDES */

#ifdef WITH_PTHREADS
#  include <pthread.h>
#endif /* WITH_PTHREADS */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/* PROJECT INCLUDES */
#include "src/common/macros.h"
#include "src/common/pack.h"
#include "src/common/parse_spec.h"
#include "src/common/read_config.h"
#include "src/common/slurm_auth.h"
#include "src/common/slurm_protocol_interface.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/slurm_protocol_common.h"
#include "src/common/slurm_protocol_pack.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/log.h"
#include "src/common/forward.h"

/* EXTERNAL VARIABLES */

/* #DEFINES */
#define _DEBUG	0
#define MAX_SHUTDOWN_RETRY 5
#define MAX_RETRIES 3

/* STATIC VARIABLES */
static pthread_mutex_t config_lock = PTHREAD_MUTEX_INITIALIZER;
static slurm_protocol_config_t proto_conf_default;
static slurm_protocol_config_t *proto_conf = &proto_conf_default;
static slurm_ctl_conf_t slurmctld_conf;

/* STATIC FUNCTIONS */
static void _remap_slurmctld_errno(void);

/**********************************************************************\
 * protocol configuration functions
\**********************************************************************/
/* slurm_set_api_config
 * sets the slurm_protocol_config object
 * NOT THREAD SAFE
 * IN protocol_conf		-  slurm_protocol_config object
 *
 * XXX: Why isn't the "config_lock" mutex used here?
 */
int slurm_set_api_config(slurm_protocol_config_t * protocol_conf)
{
	proto_conf = protocol_conf;
	return SLURM_SUCCESS;
}

/* slurm_get_api_config
 * returns a pointer to the current slurm_protocol_config object
 * RET slurm_protocol_config_t  - current slurm_protocol_config object
 */
slurm_protocol_config_t *slurm_get_api_config()
{
	return proto_conf;
}

/* slurm_api_set_conf_file
 *      set slurm configuration file to a non-default value
 * pathname IN - pathname of slurm configuration file to be used
 */
extern void  slurm_api_set_conf_file(char *pathname)
{
	if (pathname == NULL)
		return;
	xfree(slurmctld_conf.slurm_conf);
	slurmctld_conf.slurm_conf = xstrdup(pathname);
	return;
}

/* slurm_api_set_default_config
 *      called by the send_controller_msg function to insure that at least 
 *	the compiled in default slurm_protocol_config object is initialized
 * RET int		 - return code
 */
int slurm_api_set_default_config()
{
	int rc = SLURM_SUCCESS;
	struct stat config_stat;
	static time_t last_config_update = (time_t) 0;

	slurm_mutex_lock(&config_lock);
	config_stat.st_mtime = 0;
	if (slurmctld_conf.slurm_conf
	&&  (stat(slurmctld_conf.slurm_conf, &config_stat) < 0)) {
		error("Can't stat %s: %m", slurmctld_conf.slurm_conf);
		rc = SLURM_ERROR;
		goto cleanup;
	}

	if (last_config_update
	&&  (slurmctld_conf.slurm_conf
	&&   (last_config_update == config_stat.st_mtime))
	&&  slurmctld_conf.control_addr
	&&  slurmctld_conf.slurmctld_port)
  		goto cleanup;

	init_slurm_conf(&slurmctld_conf);
	read_slurm_conf_ctl(&slurmctld_conf, false);
	if (!config_stat.st_mtime)
		stat(slurmctld_conf.slurm_conf, &config_stat);
	last_config_update = config_stat.st_mtime;

	if ((slurmctld_conf.control_addr == NULL) ||
	    (slurmctld_conf.slurmctld_port == 0)) {
		error("Unable to establish control machine or port");
		rc = SLURM_ERROR;
		goto cleanup;
	}

	slurm_set_addr(&proto_conf_default.primary_controller,
		       slurmctld_conf.slurmctld_port,
		       slurmctld_conf.control_addr);
	if (proto_conf_default.primary_controller.sin_port == 0) {
		error("Unable to establish control machine address");
		rc = SLURM_ERROR;
		goto cleanup;
	}

	if (slurmctld_conf.backup_addr) {
		slurm_set_addr(&proto_conf_default.secondary_controller,
			       slurmctld_conf.slurmctld_port,
			       slurmctld_conf.backup_addr);
	}
	proto_conf = &proto_conf_default;

      cleanup:
	slurm_mutex_unlock(&config_lock);
	return rc;
}

/* slurm_api_clear_config
 * execute this only at program termination to free all memory */
void slurm_api_clear_config(void)
{
	slurm_mutex_lock(&config_lock);
	slurmctld_conf.slurmd_port = 0;
	free_slurm_conf(&slurmctld_conf);
	slurm_mutex_unlock(&config_lock);
}

/* update internal configuration data structure as needed.
 *	exit with lock set */
static inline void _lock_update_config()
{
	slurm_api_set_default_config();
	slurm_mutex_lock(&config_lock);
}

/* slurm_get_mpi_default
 * get default mpi value from slurmctld_conf object
 * RET char *   - mpi default value from slurm.conf,  MUST be xfreed by caller
 */
char *slurm_get_mpi_default(void)
{
	char *mpi_default;

	_lock_update_config();
	mpi_default = xstrdup(slurmctld_conf.mpi_default);
	slurm_mutex_unlock(&config_lock);
	return mpi_default;
}

/* slurm_get_plugin_dir
 * get plugin directory from slurmctld_conf object
 * RET char *   - plugin directory, MUST be xfreed by caller
 */
char *slurm_get_plugin_dir(void)
{
	char *plugin_dir;

	_lock_update_config();
	plugin_dir = xstrdup(slurmctld_conf.plugindir);
	slurm_mutex_unlock(&config_lock);
	return plugin_dir;
}

/* slurm_get_auth_type
 * returns the authentication type from slurmctld_conf object
 * RET char *    - auth type, MUST be xfreed by caller
 */
char *slurm_get_auth_type(void)
{
	char *auth_type;

	_lock_update_config();
	auth_type = xstrdup(slurmctld_conf.authtype);
	slurm_mutex_unlock(&config_lock);
	return auth_type;
}

/* slurm_get_fast_schedule
 * returns the value of fast_schedule in slurmctld_conf object
 */
extern uint16_t slurm_get_fast_schedule(void)
{
	uint16_t fast_val;

	_lock_update_config();
	fast_val = slurmctld_conf.fast_schedule;
	slurm_mutex_unlock(&config_lock);
	return fast_val;
}

/* slurm_set_tree_width
 * sets the value of tree_width in slurmctld_conf object
 * RET 0 or error code
 */
extern int slurm_set_tree_width(uint16_t tree_width)
{
	if (tree_width == 0) {
		error("can't have span count of 0");
		return SLURM_ERROR;
	}
	slurmctld_conf.tree_width = tree_width;
	return SLURM_SUCCESS;
}
/* slurm_get_tree_width
 * returns the value of tree_width in slurmctld_conf object
 */
extern uint16_t slurm_get_tree_width(void)
{
	uint16_t tree_width;

	_lock_update_config();
	tree_width = slurmctld_conf.tree_width;
	slurm_mutex_unlock(&config_lock);
	return tree_width;
}

/* slurm_set_auth_type
 * set the authentication type in slurmctld_conf object
 * used for security testing purposes
 * RET 0 or error code
 */
extern int slurm_set_auth_type(char *auth_type)
{
	_lock_update_config();
	xfree(slurmctld_conf.authtype);
	slurmctld_conf.authtype = xstrdup(auth_type);
	slurm_mutex_unlock(&config_lock);
	return 0;
}

/* slurm_get_jobacct_loc
 * returns the job accounting loc from the slurmctld_conf object
 * RET char *    - job accounting loc,  MUST be xfreed by caller
 */
char *slurm_get_jobacct_loc(void)
{
	char *jobacct_loc;

	_lock_update_config();
	jobacct_loc = xstrdup(slurmctld_conf.job_acct_loc);
	slurm_mutex_unlock(&config_lock);
	return jobacct_loc;
}

/* slurm_get_jobacct_parameters
 * returns the job accounting parameters from the slurmctld_conf object
 * RET char *    - job accounting parameters,  MUST be xfreed by caller
 */
char *slurm_get_jobacct_parameters(void)
{
	char *jobacct_parameters;

	_lock_update_config();
	jobacct_parameters = xstrdup(slurmctld_conf.job_acct_parameters);
	slurm_mutex_unlock(&config_lock);
	return jobacct_parameters;
}

/* slurm_get_jobacct_type
 * returns the job accounting type from the slurmctld_conf object
 * RET char *    - job accounting type,  MUST be xfreed by caller
 */
char *slurm_get_jobacct_type(void)
{
	char *jobacct_type;

	_lock_update_config();
	jobacct_type = xstrdup(slurmctld_conf.job_acct_type);
	slurm_mutex_unlock(&config_lock);
	return jobacct_type;
}

/* slurm_get_jobcomp_type
 * returns the job completion logger type from slurmctld_conf object
 * RET char *    - job completion type,  MUST be xfreed by caller
 */
char *slurm_get_jobcomp_type(void)
{
	char *jobcomp_type;

	_lock_update_config();
	jobcomp_type = xstrdup(slurmctld_conf.job_comp_type);
	slurm_mutex_unlock(&config_lock);
	return jobcomp_type;
}

/* slurm_get_proctrack_type
 * get ProctrackType from slurmctld_conf object
 * RET char *   - proctrack type, MUST be xfreed by caller
 */
char *slurm_get_proctrack_type(void)
{
	char *proctrack_type;

	_lock_update_config();
	proctrack_type = xstrdup(slurmctld_conf.proctrack_type);
	slurm_mutex_unlock(&config_lock);
	return proctrack_type;
}

/* slurm_get_slurmd_port
 * returns slurmd port from slurmctld_conf object
 * RET uint16_t	- slurmd port
 */
uint16_t slurm_get_slurmd_port(void)
{
	uint16_t slurmd_port;

	_lock_update_config();
	slurmd_port = slurmctld_conf.slurmd_port;
	slurm_mutex_unlock(&config_lock);
	return slurmd_port;
}

/* slurm_get_slurm_user_id
 * returns slurmd uid from slurmctld_conf object
 * RET uint32_t	- slurm user id
 */
uint32_t slurm_get_slurm_user_id(void)
{
	uint32_t slurm_uid;

	_lock_update_config();
	slurm_uid = slurmctld_conf.slurm_user_id;
	slurm_mutex_unlock(&config_lock);
	return slurm_uid;
}

/* slurm_get_sched_type
 * get sched type from slurmctld_conf object
 * RET char *   - sched type, MUST be xfreed by caller
 */
char *slurm_get_sched_type(void)
{
	char *sched_type;

	_lock_update_config();
	sched_type = xstrdup(slurmctld_conf.schedtype);
	slurm_mutex_unlock(&config_lock);
	return sched_type;
}

/* slurm_get_select_type
 * get select_type from slurmctld_conf object
 * RET char *   - select_type, MUST be xfreed by caller
 */
char *slurm_get_select_type(void)
{
	char *select_type;

	_lock_update_config();
	select_type = xstrdup(slurmctld_conf.select_type);
	slurm_mutex_unlock(&config_lock);
	return select_type;
}

/* slurm_get_switch_type
 * get switch type from slurmctld_conf object
 * RET char *   - switch type, MUST be xfreed by caller
 */
char *slurm_get_switch_type(void)
{
	char *switch_type;

	_lock_update_config();
	switch_type = xstrdup(slurmctld_conf.switch_type);
	slurm_mutex_unlock(&config_lock);
	return switch_type;
}

/* slurm_get_wait_time
 * returns wait_time from slurmctld_conf object
 * RET uint16_t	- wait_time
 */
uint16_t slurm_get_wait_time(void)
{
	uint16_t wait_time;

	_lock_update_config();
	wait_time = slurmctld_conf.wait_time;
	slurm_mutex_unlock(&config_lock);
	return wait_time;
}

/* slurm_get_srun_prolog
 * return the name of the srun prolog program
 * RET char *   - name of prolog program, must be xfreed by caller
 */
char *slurm_get_srun_prolog(void)
{
	char *prolog;

	_lock_update_config();
	prolog = xstrdup(slurmctld_conf.srun_prolog);
	slurm_mutex_unlock(&config_lock);
	return prolog;
}

/* slurm_get_srun_epilog
 * return the name of the srun epilog program
 * RET char *   - name of epilog program, must be xfreed by caller
 */
char *slurm_get_srun_epilog(void)
{
	char *epilog;

	_lock_update_config();
	epilog = xstrdup(slurmctld_conf.srun_epilog);
	slurm_mutex_unlock(&config_lock);
	return epilog;
}

/* slurm_get_task_epilog
 * RET task_epilog name, must be xfreed by caller */
char *slurm_get_task_epilog(void)
{
        char *task_epilog;

        _lock_update_config();
        task_epilog = xstrdup(slurmctld_conf.task_epilog);
        slurm_mutex_unlock(&config_lock);
        return task_epilog;
}

/* slurm_get_task_prolog
 * RET task_prolog name, must be xfreed by caller */
char *slurm_get_task_prolog(void)
{
        char *task_prolog;
        
	_lock_update_config();
        task_prolog = xstrdup(slurmctld_conf.task_prolog);
        slurm_mutex_unlock(&config_lock);
        return task_prolog;
}

/* slurm_get_task_plugin
 * RET task_plugin name, must be xfreed by caller */
char *slurm_get_task_plugin(void)
{
        char *task_plugin;

        _lock_update_config();
        task_plugin = xstrdup(slurmctld_conf.task_plugin);
        slurm_mutex_unlock(&config_lock);
        return task_plugin;
}
/* Change general slurm communication errors to slurmctld specific errors */
static void _remap_slurmctld_errno(void)
{
	int err = slurm_get_errno();

	if (err == SLURM_COMMUNICATIONS_CONNECTION_ERROR)
		slurm_seterrno(SLURMCTLD_COMMUNICATIONS_CONNECTION_ERROR);
	else if (err ==  SLURM_COMMUNICATIONS_SEND_ERROR)
		slurm_seterrno(SLURMCTLD_COMMUNICATIONS_SEND_ERROR);
	else if (err == SLURM_COMMUNICATIONS_RECEIVE_ERROR)
		slurm_seterrno(SLURMCTLD_COMMUNICATIONS_RECEIVE_ERROR);
	else if (err == SLURM_COMMUNICATIONS_SHUTDOWN_ERROR)
		slurm_seterrno(SLURMCTLD_COMMUNICATIONS_SHUTDOWN_ERROR);
}

/**********************************************************************\
 * general message management functions used by slurmctld, slurmd
\**********************************************************************/

/* 
 *  Initialize a slurm server at port "port"
 * 
 * IN  port     - port to bind the msg server to
 * RET slurm_fd - file descriptor of the connection created
 */
slurm_fd slurm_init_msg_engine_port(uint16_t port)
{
	slurm_addr addr;

	slurm_set_addr_any(&addr, port);
	return _slurm_init_msg_engine(&addr);
}

/* 
 *  Same as above, but initialize using a slurm address "addr"
 *
 * IN  addr     - slurm_addr to bind the msg server to 
 * RET slurm_fd - file descriptor of the connection created
 */
slurm_fd slurm_init_msg_engine(slurm_addr *addr)
{
	return _slurm_init_msg_engine(addr);
}

/* 
 *  Close an established message engine.
 *    Returns SLURM_SUCCESS or SLURM_FAILURE.
 *
 * IN  fd  - an open file descriptor to close
 * RET int - the return code
 */
int slurm_shutdown_msg_engine(slurm_fd fd)
{
	int rc = _slurm_close(fd);
	if (rc)
		slurm_seterrno(SLURM_COMMUNICATIONS_SHUTDOWN_ERROR);
	return rc;
}

/* 
 *   Close an established message connection.
 *     Returns SLURM_SUCCESS or SLURM_FAILURE.
 *
 * IN  fd  - an open file descriptor to close
 * RET int - the return code
 */
int slurm_shutdown_msg_conn(slurm_fd fd)
{
	return _slurm_close(fd);
}

/**********************************************************************\
 * msg connection establishment functions used by msg clients
\**********************************************************************/

/* In the bsd socket implementation it creates a SOCK_STREAM socket  
 *	and calls connect on it a SOCK_DGRAM socket called with connect   
 *	is defined to only receive messages from the address/port pair  
 *	argument of the connect call slurm_address - for now it is  
 *	really just a sockaddr_in
 * IN slurm_address     - slurm_addr of the connection destination
 * RET slurm_fd         - file descriptor of the connection created
 */
slurm_fd slurm_open_msg_conn(slurm_addr * slurm_address)
{
	return _slurm_open_msg_conn(slurm_address);
}

/* calls connect to make a connection-less datagram connection to the 
 *	primary or secondary slurmctld message engine
 * RET slurm_fd	- file descriptor of the connection created
 */
slurm_fd slurm_open_controller_conn()
{
	slurm_fd fd;

	if (slurm_api_set_default_config() < 0)
		return SLURM_FAILURE;

	if ((fd = slurm_open_msg_conn(&proto_conf->primary_controller)) >= 0)
		return fd;
	
	debug("Failed to contact primary controller: %m");

	if (!slurmctld_conf.backup_controller) 
		goto fail;

	if ((fd = slurm_open_msg_conn(&proto_conf->secondary_controller)) >= 0)
		return fd;

	debug("Failed to contact secondary controller: %m");

    fail:
	slurm_seterrno_ret(SLURMCTLD_COMMUNICATIONS_CONNECTION_ERROR);
}

/* calls connect to make a connection-less datagram connection to the 
 *	primary or secondary slurmctld message engine
 * RET slurm_fd - file descriptor of the connection created
 * IN dest      - controller to contact, primary or secondary
 */
slurm_fd slurm_open_controller_conn_spec(enum controller_id dest)
{
	slurm_addr *addr;
	slurm_fd rc;

	if (slurm_api_set_default_config() < 0) {
		debug3("Error: Unable to set default config");
		return SLURM_ERROR;
	}
		
	addr = (dest == PRIMARY_CONTROLLER) ? 
		  &proto_conf->primary_controller : 
		  &proto_conf->secondary_controller;

	if (!addr) return SLURM_ERROR;

	rc = slurm_open_msg_conn(addr);
	if (rc == -1)
		_remap_slurmctld_errno();
	return rc;
}

/* In the bsd implmentation maps directly to a accept call 
 * IN open_fd		- file descriptor to accept connection on
 * OUT slurm_address	- slurm_addr of the accepted connection
 * RET slurm_fd		- file descriptor of the connection created
 */
slurm_fd slurm_accept_msg_conn(slurm_fd open_fd,
			       slurm_addr * slurm_address)
{
	return _slurm_accept_msg_conn(open_fd, slurm_address);
}

/* In the bsd implmentation maps directly to a close call, to close 
 *	the socket that was accepted
 * IN open_fd		- an open file descriptor to close
 * RET int		- the return code
 */
int slurm_close_accepted_conn(slurm_fd open_fd)
{
	return _slurm_close_accepted_conn(open_fd);
}

/**********************************************************************\
 * receive message functions
\**********************************************************************/

/*
 * NOTE: memory is allocated for the returned msg and must be freed at 
 *	some point using the slurm_free_functions
 * IN open_fd	- file descriptor to receive msg on
 * OUT msg	- a slurm_msg struct to be filled in by the function
 * RET int	- size of msg received in bytes before being unpacked
 */
List slurm_receive_msg(slurm_fd fd, slurm_msg_t *msg, int timeout)
{
	char *buf = NULL;
	size_t buflen = 0;
	header_t header;
	int rc;
	void *auth_cred = NULL;
	Buf buffer;
	int count = 0;
	forward_struct_t *forward_struct = NULL;
	ret_types_t *ret_type = NULL;
	ListIterator itr;
	int i=0;
	int16_t fwd_cnt = 0;
	char addrbuf[INET_ADDRSTRLEN];
	int steps = 0;

	List ret_list = list_create(destroy_ret_types);
	
	xassert(fd >= 0);
	
	if ((timeout*=1000) == 0)
		timeout = SLURM_MESSAGE_TIMEOUT_MSEC_STATIC;
	/*
	 * Receive a msg. slurm_msg_recvfrom() will read the message
	 *  length and allocate space on the heap for a buffer containing
	 *  the message. 
	 */
	if (_slurm_msg_recvfrom_timeout(fd, &buf, &buflen, 0, timeout) < 0) {
		forward_init(&header.forward, NULL);
		rc = errno;		
		goto total_return;
	}
	
#if	_DEBUG
	_print_data (buftemp, rc);
#endif
	buffer = create_buf(buf, buflen);

	if(unpack_header(&header, buffer) == SLURM_ERROR) {
		free_buf(buffer);
		rc = SLURM_COMMUNICATIONS_RECEIVE_ERROR;
		goto total_return;
	}
	
	if (check_header_version(&header) < 0) {
		free_buf(buffer);
		rc = SLURM_PROTOCOL_VERSION_ERROR;
		goto total_return;
	}
	//info("ret_cnt = %d",header.ret_cnt);
	if(header.ret_cnt > 0) {
		while((ret_type = list_pop(header.ret_list)) != NULL)
			list_push(ret_list, ret_type);
		header.ret_cnt = 0;
		list_destroy(header.ret_list);
		header.ret_list = NULL;
	}
	
	if(header.orig_addr.sin_addr.s_addr != 0) {
		memcpy(&msg->orig_addr, &header.orig_addr, sizeof(slurm_addr));
	} else {
		memcpy(&header.orig_addr, &msg->orig_addr, sizeof(slurm_addr));
	}
	fwd_cnt = header.forward.cnt;
	/* Forward message to other nodes */
	if(fwd_cnt > 0) {
		forward_struct = xmalloc(sizeof(forward_struct_t));
		forward_struct->buf_len = remaining_buf(buffer);
		forward_struct->buf = 
			xmalloc(sizeof(char) * forward_struct->buf_len);
		memcpy(forward_struct->buf, 
		       &buffer->head[buffer->processed], 
		       forward_struct->buf_len);
		forward_struct->ret_list = ret_list;

		forward_struct->timeout = timeout - header.forward.timeout;
		
		debug3("forwarding messages to %d nodes!!!!", fwd_cnt);
		
		if(forward_msg(forward_struct, &header) == SLURM_ERROR) {
			error("problem with forward msg");
		}
	}
	
	if ((auth_cred = g_slurm_auth_unpack(buffer)) == NULL) {
		error( "authentication: %s ",
			g_slurm_auth_errstr(g_slurm_auth_errno(NULL)));
		free_buf(buffer);
		rc = ESLURM_PROTOCOL_INCOMPLETE_PACKET;
		goto total_return;
	}
	rc = g_slurm_auth_verify( auth_cred, NULL, 2 );
	
	if (rc != SLURM_SUCCESS) {
		error( "authentication: %s ",
		       g_slurm_auth_errstr(g_slurm_auth_errno(auth_cred)));
		(void) slurm_free_cred(auth_cred);
		free_buf(buffer);
		rc = SLURM_PROTOCOL_AUTHENTICATION_ERROR;
		goto total_return;
	}	

	/*
	 * Unpack message body 
	 */
	msg->msg_type = header.msg_type;
	msg->srun_node_id = header.srun_node_id;
	
	if ( (header.body_length > remaining_buf(buffer)) ||
	     (unpack_msg(msg, buffer) != SLURM_SUCCESS) ) {
		(void) slurm_free_cred(auth_cred);
		free_buf(buffer);
		rc = ESLURM_PROTOCOL_INCOMPLETE_PACKET;
		goto total_return;
	}
	
	msg->cred = (void *) auth_cred;

	free_buf(buffer);
	rc = SLURM_SUCCESS;
	if(forward_struct) {
		slurm_mutex_lock(&forward_struct->forward_mutex);
		count = 0;
		itr = list_iterator_create(ret_list);
		while((ret_type = (ret_types_t *) list_next(itr)) 
		      != NULL) {
			count += list_count(ret_type->ret_data_list);
		}
		list_iterator_destroy(itr);
		debug3("Got back %d",count);
		while((count < fwd_cnt)) {
			pthread_cond_wait(&forward_struct->notify, 
					  &forward_struct->forward_mutex);
			count = 0;
			itr = list_iterator_create(ret_list);
			while((ret_type = (ret_types_t *) list_next(itr)) 
			      != NULL) {
				count += list_count(ret_type->ret_data_list);
			}
			list_iterator_destroy(itr);
			debug3("Got back %d",count);
				
		}
		debug2("Got them all");
		slurm_mutex_unlock(&forward_struct->forward_mutex);
	}
	
total_return:
	destroy_forward(&header.forward);
	destroy_forward_struct(forward_struct);
	
	if(rc != SLURM_SUCCESS) {
		error("slurm_receive_msg: %s", slurm_strerror(rc));
	}
	//info("rc= %d count of ret is %d",rc, list_count(ret_list));
	errno = rc;
	return ret_list;
		
}

/**********************************************************************\
 * send message functions
\**********************************************************************/

/*
 *  Do the wonderful stuff that needs be done to pack msg
 *  and hdr into buffer
 */
static void
_pack_msg(slurm_msg_t *msg, header_t *hdr, Buf buffer)
{
	unsigned int tmplen, msglen;

	tmplen = get_buf_offset(buffer);
	pack_msg(msg, buffer);
	msglen = get_buf_offset(buffer) - tmplen;

	/* update header with correct cred and msg lengths */
	update_header(hdr, msglen);
	
	/* repack updated header */
	tmplen = get_buf_offset(buffer);
	set_buf_offset(buffer, 0);
	pack_header(hdr, buffer);
	set_buf_offset(buffer, tmplen);
}

/* 
 *  Adds header to buffer and sends message over an open file descriptor `fd'
 *    Returns the size of the message sent in bytes, or -1 on failure.
 */
int slurm_add_header_and_send(slurm_fd fd, slurm_msg_t *msg)
{
	Buf send_buf = NULL;
	header_t header;
	int buf_len = get_buf_offset(msg->buffer);
	int rc = SLURM_SUCCESS;
	void *auth_cred = NULL;

	init_header(&header, msg, SLURM_PROTOCOL_NO_FLAGS);
	header.body_length = buf_len;
	send_buf = init_buf(BUF_SIZE+buf_len);
	
	/*
	 * Pack header into buffer for transmission
	 */
	pack_header(&header, send_buf);
	/* 
	 * Initialize header with Auth credential and message type.
	 */
	auth_cred = g_slurm_auth_create(NULL, 2);
	if (auth_cred == NULL) {
		error("authentication: %s",
		       g_slurm_auth_errstr(g_slurm_auth_errno(NULL)) );
		errno = SLURM_PROTOCOL_AUTHENTICATION_ERROR;
		goto ret_error;
	}
	/* 
	 * Pack auth credential
	 */
	rc = g_slurm_auth_pack(auth_cred, send_buf);
	(void) slurm_free_cred(auth_cred);
	if (rc) {
		error("authentication: %s",
		       g_slurm_auth_errstr(g_slurm_auth_errno(auth_cred)));
		free_buf(send_buf);
		send_buf = NULL;
		errno = SLURM_PROTOCOL_AUTHENTICATION_ERROR;
		goto ret_error;
	}

	/* add to send_buf */
	if (remaining_buf(send_buf) < buf_len) {
		send_buf->size += (buf_len + BUF_SIZE);
		xrealloc(send_buf->head, send_buf->size);
	}
	if (buf_len) {
		memcpy(&send_buf->head[send_buf->processed],
		       msg->buffer->head, buf_len);
		send_buf->processed += buf_len;
	}
	
	/*
	 * send message
	 */
	rc = _slurm_msg_sendto(fd, 
			       get_buf_data(send_buf), 
			       get_buf_offset(send_buf),
			       SLURM_PROTOCOL_NO_SEND_RECV_FLAGS);
ret_error:
	if (rc < 0) 
		error("slurm_add_header_and_send: %m");
	free_buf(send_buf);
	return rc;

}


/* 
 *  Send a slurm message over an open file descriptor `fd'
 *    Returns the size of the message sent in bytes, or -1 on failure.
 */
int slurm_send_node_msg(slurm_fd fd, slurm_msg_t * msg)
{
	header_t header;
	Buf      buffer;
	int      rc;
	void *   auth_cred;
	
	/* 
	 * Initialize header with Auth credential and message type.
	 */
	auth_cred = g_slurm_auth_create(NULL, 2);
	if (auth_cred == NULL) {
		error("authentication: %s",
		       g_slurm_auth_errstr(g_slurm_auth_errno(NULL)) );
		slurm_seterrno_ret(SLURM_PROTOCOL_AUTHENTICATION_ERROR);
	}

	if(msg->forward.init != FORWARD_INIT) {
		forward_init(&msg->forward, NULL);
		msg->ret_list = NULL;
	}
	init_header(&header, msg, SLURM_PROTOCOL_NO_FLAGS);

	/*
	 * Pack header into buffer for transmission
	 */
	buffer = init_buf(0);
	pack_header(&header, buffer);
	
	/* 
	 * Pack auth credential
	 */
	rc = g_slurm_auth_pack(auth_cred, buffer);
	(void) slurm_free_cred(auth_cred);
	if (rc) {
		error("authentication: %s",
		       g_slurm_auth_errstr(g_slurm_auth_errno(auth_cred)));
		free_buf(buffer);
		slurm_seterrno_ret(SLURM_PROTOCOL_AUTHENTICATION_ERROR);
	}
	
	/*
	 * Pack message into buffer
	 */
	_pack_msg(msg, &header, buffer);

#if	_DEBUG
	_print_data (get_buf_data(buffer),get_buf_offset(buffer));
#endif
	/*
	 * Send message
	 */
	rc = _slurm_msg_sendto( fd, get_buf_data(buffer), 
				get_buf_offset(buffer),
				SLURM_PROTOCOL_NO_SEND_RECV_FLAGS );

	if (rc < 0) 
		error("slurm_msg_sendto: %m");
		
	free_buf(buffer);
	return rc;
}

/**********************************************************************\
 * stream functions
\**********************************************************************/

/* slurm_listen_stream
 * opens a stream server and listens on it
 * IN slurm_address	- slurm_addr to bind the server stream to
 * RET slurm_fd		- file descriptor of the stream created
 */
slurm_fd slurm_listen_stream(slurm_addr * slurm_address)
{
	return _slurm_listen_stream(slurm_address);
}

/* slurm_accept_stream
 * accepts a incomming stream connection on a stream server slurm_fd 
 * IN open_fd		- file descriptor to accept connection on
 * OUT slurm_address	- slurm_addr of the accepted connection
 * RET slurm_fd		- file descriptor of the accepted connection 
 */
slurm_fd slurm_accept_stream(slurm_fd open_fd, slurm_addr * slurm_address)
{
	return _slurm_accept_stream(open_fd, slurm_address);
}

/* slurm_open_stream
 * opens a client connection to stream server
 * IN slurm_address     - slurm_addr of the connection destination
 * RET slurm_fd         - file descriptor of the connection created
 * NOTE: Retry with various ports as needed if connection is refused
 */
slurm_fd slurm_open_stream(slurm_addr * slurm_address)
{
	return _slurm_open_stream(slurm_address, true);
}

/* slurm_write_stream
 * writes a buffer out a stream file descriptor
 * IN open_fd		- file descriptor to write on
 * IN buffer		- buffer to send
 * IN size		- size of buffer send
 * IN timeout		- how long to wait in milliseconds
 * RET size_t		- bytes sent , or -1 on errror
 */
size_t slurm_write_stream(slurm_fd open_fd, char *buffer, size_t size)
{
	return _slurm_send_timeout(open_fd, buffer, size,
				   SLURM_PROTOCOL_NO_SEND_RECV_FLAGS,
				   SLURM_MESSAGE_TIMEOUT_MSEC_STATIC);
}
size_t slurm_write_stream_timeout(slurm_fd open_fd, char *buffer,
				  size_t size, int timeout)
{
	return _slurm_send_timeout(open_fd, buffer, size,
				   SLURM_PROTOCOL_NO_SEND_RECV_FLAGS,
				   timeout);
}

/* slurm_read_stream
 * read into buffer grom a stream file descriptor
 * IN open_fd	- file descriptor to read from
 * OUT buffer   - buffer to receive into
 * IN size	- size of buffer
 * IN timeout	- how long to wait in milliseconds
 * RET size_t	- bytes read , or -1 on errror
 */
size_t slurm_read_stream(slurm_fd open_fd, char *buffer, size_t size)
{
	return _slurm_recv_timeout(open_fd, buffer, size,
				   SLURM_PROTOCOL_NO_SEND_RECV_FLAGS,
				   SLURM_MESSAGE_TIMEOUT_MSEC_STATIC);
}
size_t slurm_read_stream_timeout(slurm_fd open_fd, char *buffer,
				 size_t size, int timeout)
{
	return _slurm_recv_timeout(open_fd, buffer, size,
				   SLURM_PROTOCOL_NO_SEND_RECV_FLAGS,
				   timeout);
}

/* slurm_get_stream_addr
 * esentially a encapsilated get_sockname  
 * IN open_fd		- file descriptor to retreive slurm_addr for
 * OUT address		- address that open_fd to bound to
 */
int slurm_get_stream_addr(slurm_fd open_fd, slurm_addr * address)
{
	return _slurm_get_stream_addr(open_fd, address);
}

/* slurm_close_stream
 * closes either a server or client stream file_descriptor
 * IN open_fd	- an open file descriptor to close
 * RET int	- the return code
 */
int slurm_close_stream(slurm_fd open_fd)
{
	return _slurm_close_stream(open_fd);
}

/* make an open slurm connection blocking or non-blocking
 *	(i.e. wait or do not wait for i/o completion )
 * IN open_fd	- an open file descriptor to change the effect
 * RET int	- the return code
 */
int slurm_set_stream_non_blocking(slurm_fd open_fd)
{
	return _slurm_set_stream_non_blocking(open_fd);
}
int slurm_set_stream_blocking(slurm_fd open_fd)
{
	return _slurm_set_stream_blocking(open_fd);
}

/**********************************************************************\
 * address conversion and management functions
\**********************************************************************/

/* slurm_set_addr_uint
 * initializes the slurm_address with the supplied port and ip_address
 * OUT slurm_address	- slurm_addr to be filled in
 * IN port		- port in host order
 * IN ip_address	- ipv4 address in uint32 host order form
 */
void slurm_set_addr_uint(slurm_addr * slurm_address, uint16_t port,
			 uint32_t ip_address)
{
	_slurm_set_addr_uint(slurm_address, port, ip_address);
}

/* slurm_set_addr_any
 * initialized the slurm_address with the supplied port on INADDR_ANY
 * OUT slurm_address	- slurm_addr to be filled in
 * IN port		- port in host order
 */
void slurm_set_addr_any(slurm_addr * slurm_address, uint16_t port)
{
	_slurm_set_addr_uint(slurm_address, port, SLURM_INADDR_ANY);
}

/* slurm_set_addr
 * initializes the slurm_address with the supplied port and host name
 * OUT slurm_address	- slurm_addr to be filled in
 * IN port		- port in host order
 * IN host		- hostname or dns name 
 */
void slurm_set_addr(slurm_addr * slurm_address, uint16_t port, char *host)
{
	_slurm_set_addr_char(slurm_address, port, host);
}

/* reset_slurm_addr
 * resets the address field of a slurm_addr, port and family unchanged
 * OUT slurm_address	- slurm_addr to be reset in
 * IN new_address	- source of address to write into slurm_address
 */
void reset_slurm_addr(slurm_addr * slurm_address, slurm_addr new_address)
{
	_reset_slurm_addr(slurm_address, new_address);
}

/* slurm_set_addr_char
 * initializes the slurm_address with the supplied port and host
 * OUT slurm_address	- slurm_addr to be filled in
 * IN port		- port in host order
 * IN host		- hostname or dns name 
 */
void slurm_set_addr_char(slurm_addr * slurm_address, uint16_t port,
			 char *host)
{
	_slurm_set_addr_char(slurm_address, port, host);
}

/* slurm_get_addr 
 * given a slurm_address it returns its port and hostname
 * IN slurm_address	- slurm_addr to be queried
 * OUT port		- port number
 * OUT host		- hostname
 * IN buf_len		- length of hostname buffer
 */
void slurm_get_addr(slurm_addr * slurm_address, uint16_t * port,
		    char *host, unsigned int buf_len)
{
	_slurm_get_addr(slurm_address, port, host, buf_len);
}

/* slurm_get_ip_str 
 * given a slurm_address it returns its port and ip address string
 * IN slurm_address	- slurm_addr to be queried
 * OUT port		- port number
 * OUT ip		- ip address in dotted-quad string form
 * IN buf_len		- length of ip buffer
 */
void slurm_get_ip_str(slurm_addr * slurm_address, uint16_t * port,
		      char *ip, unsigned int buf_len)
{
	unsigned char *uc = (unsigned char *)&slurm_address->sin_addr.s_addr;
	*port = slurm_address->sin_port;
	snprintf(ip, buf_len, "%u.%u.%u.%u", uc[0], uc[1], uc[2], uc[3]);
}

/* slurm_get_peer_addr
 * get the slurm address of the peer connection, similar to getpeeraddr
 * IN fd		- an open connection
 * OUT slurm_address	- place to park the peer's slurm_addr
 */
int slurm_get_peer_addr(slurm_fd fd, slurm_addr * slurm_address)
{
	struct sockaddr name;
	socklen_t namelen = (socklen_t) sizeof(struct sockaddr);
	int rc;

	if ((rc = _slurm_getpeername((int) fd, &name, &namelen)))
		return rc;
	memcpy(slurm_address, &name, sizeof(slurm_addr));
	return 0;
}

/* slurm_print_slurm_addr
 * prints a slurm_addr into a buf
 * IN address		- slurm_addr to print
 * IN buf		- space for string representation of slurm_addr
 * IN n			- max number of bytes to write (including NUL)
 */
void slurm_print_slurm_addr(slurm_addr * address, char *buf, size_t n)
{
	_slurm_print_slurm_addr(address, buf, n);
}

/**********************************************************************\
 * slurm_addr pack routines
\**********************************************************************/

/* 
 *  Pack just the message with no header and send back the buffer.
 */
Buf slurm_pack_msg_no_header(slurm_msg_t * msg)
{
	Buf      buffer = NULL;
	void *   auth_cred;
	int rc;

	buffer = init_buf(0);
	
	/*
	 * Pack message into buffer
	 */
	pack_msg(msg, buffer);
	
ret_error:

	return buffer;
}

/* slurm_pack_slurm_addr
 * packs a slurm_addr into a buffer to serialization transport
 * IN slurm_address	- slurm_addr to pack
 * IN/OUT buffer	- buffer to pack the slurm_addr into
 */
void slurm_pack_slurm_addr(slurm_addr * slurm_address, Buf buffer)
{
	_slurm_pack_slurm_addr(slurm_address, buffer);
}

/* slurm_pack_slurm_addr
 * unpacks a buffer into a slurm_addr after serialization transport
 * OUT slurm_address	- slurm_addr to unpack to
 * IN/OUT buffer	- buffer to upack the slurm_addr from
 * returns		- SLURM error code
 */
int slurm_unpack_slurm_addr_no_alloc(slurm_addr * slurm_address,
				     Buf buffer)
{
	return _slurm_unpack_slurm_addr_no_alloc(slurm_address, buffer);
}

/**********************************************************************\
 * simplified communication routines 
 * They open a connection do work then close the connection all within 
 * the function
\**********************************************************************/

/* slurm_send_rc_msg
 * given the original request message this function sends a 
 *	slurm_return_code message back to the client that made the request
 * IN request_msg	- slurm_msg the request msg
 * IN rc		- the return_code to send back to the client
 */
int slurm_send_rc_msg(slurm_msg_t *msg, int rc)
{
	slurm_msg_t resp_msg;
	return_code_msg_t rc_msg;
	
	if (msg->conn_fd < 0)
		return (ENOTCONN);
	rc_msg.return_code = rc;

	resp_msg.address  = msg->address;
	resp_msg.msg_type = RESPONSE_SLURM_RC;
	resp_msg.data     = &rc_msg;
	resp_msg.forward = msg->forward;
	resp_msg.ret_list = msg->ret_list;
	resp_msg.orig_addr = msg->orig_addr;
	resp_msg.srun_node_id = msg->srun_node_id;
	
	/* send message */
	return slurm_send_node_msg(msg->conn_fd, &resp_msg);
}

/*
 * Send and recv a slurm request and response on the open slurm descriptor
 */
static List 
_send_and_recv_msg(slurm_fd fd, slurm_msg_t *req, 
		   slurm_msg_t *resp, int timeout)
{
	int err = SLURM_SUCCESS;
	int retry = 0;
	List ret_list = NULL;
	int steps = 0;

	if (slurm_send_node_msg(fd, req) < 0)
		err = errno;
			
	if(err == SLURM_SUCCESS) {
		if ((timeout*=1000) == 0)
			timeout = SLURM_MESSAGE_TIMEOUT_MSEC_STATIC;
		
		if(req->forward.cnt>0) {
			steps = req->forward.cnt/slurmctld_conf.tree_width;
			steps += 1;
			timeout += (req->forward.timeout*steps);
		}
		ret_list = slurm_receive_msg(fd, resp, timeout);
	}
	
	if(!ret_list || list_count(ret_list) == 0) {
		no_resp_forwards(&req->forward, &ret_list, errno);
	}

	/* 
	 *  Attempt to close an open connection
	 */
	while ((slurm_shutdown_msg_conn(fd) < 0) && (errno == EINTR) ) {
		if (retry++ > MAX_SHUTDOWN_RETRY) {
			err = errno;
			break;
		}
	}
	if (err) 
		errno = err; 

	return ret_list;
}

/*
 * slurm_send_recv_controller_msg
 * opens a connection to the controller, sends the controller a message, 
 * listens for the response, then closes the connection
 * IN request_msg	- slurm_msg request
 * OUT response_msg     - slurm_msg response
 * RET int              - return code
 */
int slurm_send_recv_controller_msg(slurm_msg_t *req, slurm_msg_t *resp)
{
	slurm_fd fd = -1;
	int rc = SLURM_SUCCESS;
	time_t start_time = time(NULL);
	List ret_list = NULL;
	int retry = 1;

	if ((fd = slurm_open_controller_conn()) < 0) {
		rc = SLURM_SOCKET_ERROR;
		goto cleanup;
	}
	forward_init(&req->forward, NULL);
	req->ret_list = NULL;
	req->orig_addr.sin_addr.s_addr = 0; 
	//info("here 2");
	
	while(retry) {
		retry = 0;
		/* If the backup controller is in the process of assuming 
		 * control, we sleep and retry later */
		
		ret_list = _send_and_recv_msg(fd, req, resp, 0);
		
		if(errno == SLURM_SUCCESS)
			slurm_free_cred(resp->cred);

		rc = errno;
		if(ret_list) {
			if(list_count(ret_list)>0) {
				error("We didn't do things correctly "
				      "missed %d responses",
				      list_count(ret_list));
			}
			list_destroy(ret_list);
		}
		
		if((rc == SLURM_SUCCESS) &&
		   (resp->msg_type == RESPONSE_SLURM_RC) &&
		   ((((return_code_msg_t *) resp->data)->return_code) 
		    == ESLURM_IN_STANDBY_MODE) &&
		   (req->msg_type 
		    != MESSAGE_NODE_REGISTRATION_STATUS) && 
		   (slurmctld_conf.backup_controller) &&
		   (difftime(time(NULL), start_time) < 
		    (slurmctld_conf.slurmctld_timeout +
		     slurmctld_conf.heartbeat_interval))) {
			debug("Neither primary nor backup controller "
			      "responding, sleep and retry");
			slurm_free_return_code_msg(resp->data);
			sleep(30);
			if ((fd = slurm_open_controller_conn()) < 0) 
				rc = SLURM_SOCKET_ERROR;
			else 
				retry = 1;
		}
		if(rc != SLURM_SUCCESS)
			break;
	}
			
      cleanup:
	if (rc != SLURM_SUCCESS) 
 		_remap_slurmctld_errno(); 
	return rc;
}

/* slurm_send_recv_node_msg
 * opens a connection to node, sends the node a message, listens 
 * for the response, then closes the connection
 * IN request_msg	- slurm_msg request
 * OUT response_msg	- slurm_msg response
 * RET List		- return list from multiple nodes
 */
List slurm_send_recv_node_msg(slurm_msg_t *req, slurm_msg_t *resp, int timeout)
{
	slurm_fd fd = -1;

	if ((fd = slurm_open_msg_conn(&req->address)) < 0)
		return NULL; //SLURM_SOCKET_ERROR;
	//info("here 3");
	
	return _send_and_recv_msg(fd, req, resp, timeout);

}

/* slurm_send_only_controller_msg
 * opens a connection to the controller, sends the controller a 
 * message then, closes the connection
 * IN request_msg	- slurm_msg request
 * RET int		- return code
 */
int slurm_send_only_controller_msg(slurm_msg_t *req)
{
	int      rc = SLURM_SUCCESS;
	int      retry = 0;
	slurm_fd fd = -1;

	/*
	 *  Open connection to SLURM controller:
	 */
	if ((fd = slurm_open_controller_conn()) < 0) {
		rc = SLURM_SOCKET_ERROR;
		goto cleanup;
	}

	rc = slurm_send_node_msg(fd, req);

	/* 
	 *  Attempt to close an open connection
	 */
	while ( (slurm_shutdown_msg_conn(fd) < 0) && (errno == EINTR) ) {
		if (retry++ > MAX_SHUTDOWN_RETRY) {
			rc = SLURM_SOCKET_ERROR;
			goto cleanup;
		}
	}

      cleanup:
	if (rc != SLURM_SUCCESS)
		_remap_slurmctld_errno();
	return rc;
}

/* 
 *  Open a connection to the "address" specified in the slurm msg `req'
 *   Then, immediately close the connection w/out waiting for a reply.
 *
 *   Returns SLURM_SUCCESS on success SLURM_FAILURE (< 0) for failure.
 */
int slurm_send_only_node_msg(slurm_msg_t *req)
{
	int      rc = SLURM_SUCCESS;
	int      retry = 0;
	slurm_fd fd = -1;
	
	if ((fd = slurm_open_msg_conn(&req->address)) < 0)
		return SLURM_SOCKET_ERROR;

	rc = slurm_send_node_msg(fd, req);

	/* 
	 *  Attempt to close an open connection
	 */
	while ( (slurm_shutdown_msg_conn(fd) < 0) && (errno == EINTR) ) {
		if (retry++ > MAX_SHUTDOWN_RETRY)
			return SLURM_SOCKET_ERROR;
	}

	return rc;
}


/*
 *  Send message and recv "return code" message on an already open
 *    slurm file descriptor
 */
static List _send_recv_rc_msg(slurm_fd fd, slurm_msg_t *req, int timeout)
{
	int		retval = SLURM_SUCCESS;
	slurm_msg_t	msg;
	List ret_list = NULL;
	ListIterator itr = NULL;
	ret_types_t *ret_type = NULL;
	ret_data_info_t *ret_data_info = NULL;
	int msg_rc;
	int set = 0;
	int err;
	//info("here 1");
	
	ret_list = _send_and_recv_msg(fd, req, &msg, timeout);
	if(!ret_list) {
		return ret_list;
	} 
	
	err = errno;
	if(errno != SLURM_SUCCESS) 
		msg_rc = SLURM_ERROR;	
	else {
		msg_rc = ((return_code_msg_t *)msg.data)->return_code;
		slurm_free_return_code_msg(msg.data);
		slurm_free_cred(msg.cred);
	}
	ret_data_info = xmalloc(sizeof(ret_data_info_t));
	ret_data_info->node_name = xstrdup("localhost");
	ret_data_info->data = NULL;

	itr = list_iterator_create(ret_list);		
	while((ret_type = list_next(itr)) != NULL) {
		if(ret_type->err == err) {
			list_push(ret_type->ret_data_list, ret_data_info);
			set = 1;
			break;
		}
	}
	list_iterator_destroy(itr);
	if(!set) {
		ret_type = xmalloc(sizeof(ret_types_t));
		list_push(ret_list, ret_type);
		ret_type->type = msg.msg_type;
		ret_type->msg_rc = msg_rc;
		ret_type->err = err;
		ret_type->ret_data_list = list_create(destroy_data_info);
		list_push(ret_type->ret_data_list, ret_data_info);
	}
	errno = err;
	return ret_list;
}

List slurm_send_recv_rc_packed_msg(slurm_msg_t *msg, int timeout)
{
	slurm_fd fd = -1;
	int err = SLURM_SUCCESS;
	int steps = 0;
	slurm_msg_t resp;
	List ret_list = NULL;
	ListIterator itr;
	ret_types_t *ret_type = NULL;
	ret_data_info_t *ret_data_info = NULL;
	int msg_rc;
	int set = 0;
	int retry = 0;

	if(!msg->buffer) {
		err = SLURMCTLD_COMMUNICATIONS_SEND_ERROR;
		goto failed;	
	}

	if ((fd = slurm_open_msg_conn(&msg->address)) < 0) {
		err = SLURM_SOCKET_ERROR;
		goto failed;
	}

	if(slurm_add_header_and_send(fd, msg) >= 0) {
		if ((timeout*=1000) == 0)
			timeout = SLURM_MESSAGE_TIMEOUT_MSEC_STATIC;
		
		if(msg->forward.cnt>0) {
			steps = msg->forward.cnt/slurmctld_conf.tree_width;
			steps += 1;
			timeout += (msg->forward.timeout*steps);
		}

		ret_list = slurm_receive_msg(fd, &resp, timeout);
	}
	err = errno;
			
	/* 
	 *  Attempt to close an open connection
	 */
	while ((slurm_shutdown_msg_conn(fd) < 0) && (errno == EINTR) ) {
		if (retry++ > MAX_RETRIES) {
			err = errno;
			break;
		}
	}
	
failed:
	if(!ret_list || list_count(ret_list) == 0) {
		no_resp_forwards(&msg->forward, &ret_list, err);
	}
	
	if(err != SLURM_SUCCESS) 
		msg_rc = SLURM_ERROR;	
	else {
		msg_rc = ((return_code_msg_t *)resp.data)->return_code;
		slurm_free_return_code_msg(resp.data);
		slurm_free_cred(resp.cred);
	}
	ret_data_info = xmalloc(sizeof(ret_data_info_t));
	ret_data_info->node_name = xstrdup("localhost");
	ret_data_info->data = NULL;

	itr = list_iterator_create(ret_list);		
	while((ret_type = list_next(itr)) != NULL) {
		if(ret_type->msg_rc == msg_rc) {
			list_push(ret_type->ret_data_list, ret_data_info);
			set = 1;
			break;
		}
	}
	list_iterator_destroy(itr);
	if(!set) {
		ret_type = xmalloc(sizeof(ret_types_t));
		list_push(ret_list, ret_type);
		ret_type->type = resp.msg_type;
		ret_type->msg_rc = msg_rc;
		ret_type->err = err;
		ret_type->ret_data_list = list_create(destroy_data_info);
		list_push(ret_type->ret_data_list, ret_data_info);
	}
	
	return ret_list; 
}

/*
 *  Open a connection to the "address" specified in the the slurm msg "req"
 *    Then read back an "rc" message returning the "return_code" specified
 *    in the response in the "rc" parameter.
 */
List slurm_send_recv_rc_msg(slurm_msg_t *req, int timeout)
{
	slurm_fd fd = -1;
	
	if ((fd = slurm_open_msg_conn(&req->address)) < 0) {
		errno = SLURM_SOCKET_ERROR;
		return NULL;
	}
	return _send_recv_rc_msg(fd, req, timeout);
}
/*
 *  Open a connection to the "address" specified in the the slurm msg "req"
 *    Then read back an "rc" message returning the "return_code" specified
 *    in the response in the "rc" parameter.
 */
int slurm_send_recv_rc_msg_only_one(slurm_msg_t *req, int *rc, int timeout)
{
	slurm_fd fd = -1;
	List ret_list = NULL;
	ret_types_t *ret_type = NULL;
	int ret_c = SLURM_SUCCESS;

	if ((fd = slurm_open_msg_conn(&req->address)) < 0) {
		return SLURM_SOCKET_ERROR;
	}

	forward_init(&req->forward, NULL);
	req->ret_list = NULL;
	req->orig_addr.sin_addr.s_addr = 0;
		
	ret_list = _send_recv_rc_msg(fd, req, timeout);
	if(ret_list) {
		if(list_count(ret_list)>1) 
			error("Got %d, expecting 1 from message receiving",
			      list_count(ret_list));

		ret_type = list_pop(ret_list);
	
		if(ret_type) {
			*rc = ret_type->msg_rc;
			ret_c = ret_type->err;
			destroy_ret_types(ret_type);
		}
		list_destroy(ret_list);
	} else 
		ret_c = SLURM_ERROR;
	return ret_c;
}

/*
 *  Same as above, but send message to controller
 */
int slurm_send_recv_controller_rc_msg(slurm_msg_t *req, int *rc)
{
	slurm_fd fd = -1;
	List ret_list = NULL;
	ret_types_t *ret_type = NULL;
	int ret_val = SLURM_SUCCESS;

	forward_init(&req->forward, NULL);
	req->ret_list = NULL;
	req->orig_addr.sin_addr.s_addr = 0; 
	
	if ((fd = slurm_open_controller_conn()) < 0)
		return SLURM_SOCKET_ERROR;
	ret_list = _send_recv_rc_msg(fd, req, 0);
	
	if(ret_list) {
		if(list_count(ret_list)>1)
			error("controller_rc_msg: Got %d instead of 1 back",
			      list_count(ret_list));
		ret_type = list_pop(ret_list);
		
		if(ret_type) {
			*rc = ret_type->msg_rc;
			ret_val = ret_type->err;
			destroy_ret_types(ret_type);
		}
		list_destroy(ret_list);
	} else 
		ret_val = SLURM_ERROR;
	return ret_val;
}

extern int *set_span(int total)
{
	int *span = xmalloc(sizeof(int)*slurmctld_conf.tree_width);
	int left = total;
	int i = 0;
	//info("span count = %d",slurmctld_conf.tree_width);
	memset(span,0,slurmctld_conf.tree_width);
	if(total <= slurmctld_conf.tree_width) {
		return span;
	} 
	
	while(left>0) {
		for(i=0; i<slurmctld_conf.tree_width; i++) {
			if((slurmctld_conf.tree_width-i)>=left) {
				if(span[i] == 0) {
					left = 0;
					break;
				} else {
					span[i] += left;
					left = 0;
					break;
				}
			} else if(left<=slurmctld_conf.tree_width) {
				span[i]+=left;
				left = 0;
				break;
			}
			span[i] += slurmctld_conf.tree_width;
			left -= slurmctld_conf.tree_width;
		}
	}
	return span;
}

/*
 * Free a slurm message
 */
void slurm_free_msg(slurm_msg_t * msg)
{
	(void) slurm_free_cred(msg->cred);
	if(msg->ret_list) {
		list_destroy(msg->ret_list);
		msg->ret_list = NULL;
	}
	xfree(msg);
}

/* 
 * Free just the credential of a message
 */
void slurm_free_cred(void *cred)
{
	(void) g_slurm_auth_destroy(cred);
}

int convert_to_kilo(int number, char *tmp)
{
	int i;
	if(number >= 1024) {
		i = number % 1024;
		if(i > 0) {
			i *= 10;
			i /= 1024;
			sprintf(tmp, "%d.%dk\0", number/1024, i);
		} else 
			sprintf(tmp, "%dk\0", number/1024);
	} else
		sprintf(tmp, "%d\0", number);

	return SLURM_SUCCESS;
}

#if _DEBUG

static void _print_data(char *data, int len)
{
	int i;
	for (i = 0; i < len; i++) {
		if ((i % 10 == 0) && (i != 0))
			printf("\n");
		printf("%2.2x ", ((int) data[i] & 0xff));
		if (i >= 200)
			break;
	}
	printf("\n\n");
}

#endif


/*
 * vi: shiftwidth=8 tabstop=8 expandtab
 */
