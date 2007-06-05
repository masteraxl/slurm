/****************************************************************************\
 *  config_info.c - get/print the system configuration information of slurm
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> and Kevin Tew <tew1@llnl.gov>.
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

#include <errno.h>
#include <stdio.h>

#include <slurm/slurm.h>

#include "src/api/job_info.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"

/*
 * slurm_api_version - Return a single number reflecting the SLURM API's 
 *      version number. Use the macros SLURM_VERSION_NUM, SLURM_VERSION_MAJOR, 
 *      SLURM_VERSION_MINOR, and SLURM_VERSION_MICRO to work with this value
 * RET API's version number
 */
extern long slurm_api_version (void)
{
	return (long) SLURM_API_VERSION;
}

static char *
_select_info(uint16_t select_type_param)
{
	switch (select_type_param) {
		case SELECT_TYPE_INFO_NONE:
			return "NONE";
		case CR_CPU:
			return "CR_CPU";
		case CR_SOCKET:
			return "CR_SOCKET";
		case CR_CORE:
			return "CR_CORE";
		case CR_MEMORY:
			return "CR_MEMORY";
		case CR_SOCKET_MEMORY:
			return "CR_SOCKET_MEMORY";
		case CR_CORE_MEMORY:
			return "CR_CORE_MEMORY";
		case CR_CPU_MEMORY:
			return "CR_CPU_MEMORY";
		default:
			return "unknown";
	}
}

static char *_task_plugin_param(uint16_t task_plugin_param)
{
	switch(task_plugin_param) {
		case TASK_PARAM_NONE:
			return "none";
		case TASK_PARAM_CPUSETS:
			return "cpusets";
		case TASK_PARAM_SCHED:
			return "sched";
		default:
			return "unknown";
	}
}

/*
 * slurm_print_ctl_conf - output the contents of slurm control configuration 
 *	message as loaded using slurm_load_ctl_conf
 * IN out - file to write to
 * IN slurm_ctl_conf_ptr - slurm control configuration pointer
 */
void slurm_print_ctl_conf ( FILE* out, 
                            slurm_ctl_conf_info_msg_t * slurm_ctl_conf_ptr )
{
	char time_str[32];

	if ( slurm_ctl_conf_ptr == NULL )
		return ;

	slurm_make_time_str ((time_t *)&slurm_ctl_conf_ptr->last_update, 
			     time_str, sizeof(time_str));
	fprintf(out, "Configuration data as of %s\n", time_str);
	fprintf(out, "AuthType          = %s\n", 
		slurm_ctl_conf_ptr->authtype);
	fprintf(out, "BackupAddr        = %s\n", 
		slurm_ctl_conf_ptr->backup_addr);
	fprintf(out, "BackupController  = %s\n", 
		slurm_ctl_conf_ptr->backup_controller);
	fprintf(out, "CacheGroups       = %u\n", 
		slurm_ctl_conf_ptr->cache_groups);
	fprintf(out, "CheckpointType    = %s\n",
		slurm_ctl_conf_ptr->checkpoint_type);
	fprintf(out, "ControlAddr       = %s\n", 
		slurm_ctl_conf_ptr->control_addr);
	fprintf(out, "ControlMachine    = %s\n", 
		slurm_ctl_conf_ptr->control_machine);
	fprintf(out, "Epilog            = %s\n", 
		slurm_ctl_conf_ptr->epilog);
	fprintf(out, "FastSchedule      = %u\n", 
		slurm_ctl_conf_ptr->fast_schedule);
	fprintf(out, "FirstJobId        = %u\n", 
		slurm_ctl_conf_ptr->first_job_id);
#ifdef HAVE_XCPU
	fprintf(out, "HAVE_XCPU         = %d\n", HAVE_XCPU);
#endif
	fprintf(out, "InactiveLimit     = %u\n", 
		slurm_ctl_conf_ptr->inactive_limit);
	fprintf(out, "JobAcctLogFile    = %s\n", 
		slurm_ctl_conf_ptr->job_acct_logfile);
	fprintf(out, "JobAcctFrequency  = %u\n",
		slurm_ctl_conf_ptr->job_acct_freq);
	fprintf(out, "JobAcctType       = %s\n", 
		slurm_ctl_conf_ptr->job_acct_type);
	fprintf(out, "JobCompLoc        = %s\n", 
		slurm_ctl_conf_ptr->job_comp_loc);
	fprintf(out, "JobCompType       = %s\n", 
		slurm_ctl_conf_ptr->job_comp_type);
	fprintf(out, "JobCredentialPrivateKey = %s\n", 
		slurm_ctl_conf_ptr->job_credential_private_key);
	fprintf(out, "JobCredentialPublicCertificate = %s\n", 
		slurm_ctl_conf_ptr->job_credential_public_certificate);
	fprintf(out, "KillWait          = %u\n", 
		slurm_ctl_conf_ptr->kill_wait);
	fprintf(out, "MailProg          = %s\n",
		slurm_ctl_conf_ptr->mail_prog);
	fprintf(out, "MaxJobCount       = %u\n", 
		slurm_ctl_conf_ptr->max_job_cnt);
	fprintf(out, "MessageTimeout    = %u\n",
		slurm_ctl_conf_ptr->msg_timeout);
	fprintf(out, "MinJobAge         = %u\n", 
		slurm_ctl_conf_ptr->min_job_age);
	fprintf(out, "MpiDefault        = %s\n",
		slurm_ctl_conf_ptr->mpi_default);
#ifdef MULTIPLE_SLURMD
	fprintf(out, "MULTIPLE_SLURMD   = %d\n", MULTIPLE_SLURMD);
#endif
	fprintf(out, "NEXT_JOB_ID       = %u\n",
		slurm_ctl_conf_ptr->next_job_id);
	fprintf(out, "PluginDir         = %s\n", 
		slurm_ctl_conf_ptr->plugindir);
	fprintf(out, "PlugStackConfig   = %s\n",
		slurm_ctl_conf_ptr->plugstack);
	fprintf(out, "PrivateData       = %u\n",
		slurm_ctl_conf_ptr->private_data);
	fprintf(out, "ProctrackType     = %s\n",
		slurm_ctl_conf_ptr->proctrack_type);
	fprintf(out, "Prolog            = %s\n", 
		slurm_ctl_conf_ptr->prolog);
	fprintf(out, "PropagatePrioProcess = %u\n",
		slurm_ctl_conf_ptr->propagate_prio_process);
        fprintf(out, "PropagateResourceLimits = %s\n",
                slurm_ctl_conf_ptr->propagate_rlimits);
        fprintf(out, "PropagateResourceLimitsExcept = %s\n", 
                slurm_ctl_conf_ptr->propagate_rlimits_except);
	fprintf(out, "ReturnToService   = %u\n", 
		slurm_ctl_conf_ptr->ret2service);
	fprintf(out, "SchedulerPort     = %u\n",
		slurm_ctl_conf_ptr->schedport);
	fprintf(out, "SchedulerRootFilter = %u\n",
		slurm_ctl_conf_ptr->schedrootfltr);
	fprintf(out, "SchedulerType     = %s\n",
		slurm_ctl_conf_ptr->schedtype);
	fprintf(out, "SelectType        = %s\n",
		slurm_ctl_conf_ptr->select_type);
	if (slurm_ctl_conf_ptr->select_type_param) {
		fprintf(out, "SelectTypeParameters = %s\n",
			_select_info(slurm_ctl_conf_ptr->
			select_type_param));
	}
	fprintf(out, "SlurmUser         = %s(%u)\n", 
		slurm_ctl_conf_ptr->slurm_user_name,
		slurm_ctl_conf_ptr->slurm_user_id);
	fprintf(out, "SlurmctldDebug    = %u\n", 
		slurm_ctl_conf_ptr->slurmctld_debug);
	fprintf(out, "SlurmctldLogFile  = %s\n", 
		slurm_ctl_conf_ptr->slurmctld_logfile);
	fprintf(out, "SlurmctldPidFile  = %s\n", 
		slurm_ctl_conf_ptr->slurmctld_pidfile);
	fprintf(out, "SlurmctldPort     = %u\n", 
		slurm_ctl_conf_ptr->slurmctld_port);
	fprintf(out, "SlurmctldTimeout  = %u\n", 
		slurm_ctl_conf_ptr->slurmctld_timeout);
	fprintf(out, "SlurmdDebug       = %u\n", 
		slurm_ctl_conf_ptr->slurmd_debug);
	fprintf(out, "SlurmdLogFile     = %s\n", 
		slurm_ctl_conf_ptr->slurmd_logfile);
	fprintf(out, "SlurmdPidFile     = %s\n", 
		slurm_ctl_conf_ptr->slurmd_pidfile);
#ifndef MULTIPLE_SLURMD
	fprintf(out, "SlurmdPort        = %u\n", 
		slurm_ctl_conf_ptr->slurmd_port);
#endif
	fprintf(out, "SlurmdSpoolDir    = %s\n", 
		slurm_ctl_conf_ptr->slurmd_spooldir);
	fprintf(out, "SlurmdTimeout     = %u\n", 
		slurm_ctl_conf_ptr->slurmd_timeout);
	fprintf(out, "SLURM_CONFIG_FILE = %s\n", 
		slurm_ctl_conf_ptr->slurm_conf);
	fprintf(out, "SLURM_VERSION     = %s\n", SLURM_VERSION);
	fprintf(out, "SrunProlog        = %s\n",
		slurm_ctl_conf_ptr->srun_prolog);
	fprintf(out, "SrunEpilog        = %s\n",
		slurm_ctl_conf_ptr->srun_epilog);
	fprintf(out, "StateSaveLocation = %s\n", 
		slurm_ctl_conf_ptr->state_save_location);
	fprintf(out, "SwitchType        = %s\n",
		slurm_ctl_conf_ptr->switch_type);
	fprintf(out, "TaskEpilog        = %s\n",
		slurm_ctl_conf_ptr->task_epilog);
	fprintf(out, "TaskPlugin        = %s\n",
		 slurm_ctl_conf_ptr->task_plugin);
	fprintf(out, "TaskPluginParam   = %s\n",
		_task_plugin_param(slurm_ctl_conf_ptr->task_plugin_param));
	fprintf(out, "TaskProlog        = %s\n",
		slurm_ctl_conf_ptr->task_prolog);
	fprintf(out, "TmpFS             = %s\n", 
		slurm_ctl_conf_ptr->tmp_fs);
	fprintf(out, "TreeWidth         = %u\n",
		slurm_ctl_conf_ptr->tree_width);
	fprintf(out, "UsePam            = %u\n",
		slurm_ctl_conf_ptr->use_pam);
	fprintf(out, "WaitTime          = %u\n", 
		slurm_ctl_conf_ptr->wait_time);
}

/*
 * slurm_load_ctl_conf - issue RPC to get slurm control configuration  
 *	information if changed since update_time 
 * IN update_time - time of current configuration data
 * IN slurm_ctl_conf_ptr - place to store slurm control configuration 
 *	pointer
 * RET 0 on success, otherwise return -1 and set errno to indicate the error
 * NOTE: free the response using slurm_free_ctl_conf
 */
int
slurm_load_ctl_conf (time_t update_time, slurm_ctl_conf_t **confp)
{
	int rc;
	slurm_msg_t req_msg;
	slurm_msg_t resp_msg;
        last_update_msg_t req; 
	
	slurm_msg_t_init(&req_msg);
	slurm_msg_t_init(&resp_msg);

	req.last_update  = update_time;
	req_msg.msg_type = REQUEST_BUILD_INFO;
	req_msg.data     = &req;

	if (slurm_send_recv_controller_msg(&req_msg, &resp_msg) < 0) 
		return SLURM_ERROR;

	switch (resp_msg.msg_type) {
	case RESPONSE_BUILD_INFO:
		*confp = (slurm_ctl_conf_info_msg_t *) resp_msg.data;
		break;
	case RESPONSE_SLURM_RC:
		rc = ((return_code_msg_t *) resp_msg.data)->return_code;
		slurm_free_return_code_msg(resp_msg.data);	
		if (rc) 
			slurm_seterrno_ret(rc);
		break;
	default:
		slurm_seterrno_ret(SLURM_UNEXPECTED_MSG_ERROR);
		break;
	}
        return SLURM_PROTOCOL_SUCCESS;
}

