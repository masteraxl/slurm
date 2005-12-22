/****************************************************************************\
 *  config_info.c - get/print the system configuration information of slurm
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette1@llnl.gov> and Kevin Tew <tew1@llnl.gov>.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <stdio.h>

#include <slurm/slurm.h>

#include "src/api/job_info.h"
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

/*
 * slurm_print_ctl_conf - output the contents of slurm control configuration 
 *	message as loaded using slurm_load_ctl_conf
 * IN out - file to write to
 * IN slurm_ctl_conf_ptr - slurm control configuration pointer
 */
void slurm_print_ctl_conf ( FILE* out, 
                            slurm_ctl_conf_info_msg_t * slurm_ctl_conf_ptr )
{
	char time_str[16];

	if ( slurm_ctl_conf_ptr == NULL )
		return ;

	slurm_make_time_str ((time_t *)&slurm_ctl_conf_ptr->last_update, 
			     time_str);
	fprintf(out, "Configuration data as of %s\n", time_str);
	fprintf(out, "AuthType          = %s\n", 
		slurm_ctl_conf_ptr->authtype);
	fprintf(out, "BackupAddr        = %s\n", 
		slurm_ctl_conf_ptr->backup_addr);
	fprintf(out, "BackupController  = %s\n", 
		slurm_ctl_conf_ptr->backup_controller);
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
	fprintf(out, "HeartbeatInterval = %u\n", 
		slurm_ctl_conf_ptr->heartbeat_interval);
	fprintf(out, "InactiveLimit     = %u\n", 
		slurm_ctl_conf_ptr->inactive_limit);
	fprintf(out, "JobAcctLoc        = %s\n", 
		slurm_ctl_conf_ptr->job_acct_loc);
	fprintf(out, "JobAcctParameters = %s\n",
		slurm_ctl_conf_ptr->job_acct_parameters);
	fprintf(out, "JobAcctType       = %s\n", 
		slurm_ctl_conf_ptr->job_acct_type);
	fprintf(out, "JobCompLoc        = %s\n", 
		slurm_ctl_conf_ptr->job_comp_loc);
	fprintf(out, "JobCompType       = %s\n", 
		slurm_ctl_conf_ptr->job_comp_type);
	fprintf(out, "JobCredPrivateKey = %s\n", 
		slurm_ctl_conf_ptr->job_credential_private_key);
	fprintf(out, "JobCredPublicKey  = %s\n", 
		slurm_ctl_conf_ptr->job_credential_public_certificate);
	fprintf(out, "KillWait          = %u\n", 
		slurm_ctl_conf_ptr->kill_wait);
	fprintf(out, "MaxJobCnt         = %u\n", 
		slurm_ctl_conf_ptr->max_job_cnt);
	fprintf(out, "MinJobAge         = %u\n", 
		slurm_ctl_conf_ptr->min_job_age);
	fprintf(out, "MpiDefault        = %s\n",
		slurm_ctl_conf_ptr->mpi_default);
	fprintf(out, "PluginDir         = %s\n", 
		slurm_ctl_conf_ptr->plugindir);
	fprintf(out, "ProctrackType     = %s\n",
		slurm_ctl_conf_ptr->proctrack_type);
	fprintf(out, "Prolog            = %s\n", 
		slurm_ctl_conf_ptr->prolog);
        fprintf(out, "PropagateResourceLimits = %s\n",
                slurm_ctl_conf_ptr->propagate_rlimits);
        fprintf(out, "PropagateResourceLimitsExcept = %s\n", 
                slurm_ctl_conf_ptr->propagate_rlimits_except);
	fprintf(out, "ReturnToService   = %u\n", 
		slurm_ctl_conf_ptr->ret2service);
	fprintf(out, "SchedulerAuth     = %s\n",
		slurm_ctl_conf_ptr->schedauth);
	fprintf(out, "SchedulerPort     = %u\n",
		slurm_ctl_conf_ptr->schedport);
	fprintf(out, "SchedulerRootFilter = %u\n",
		slurm_ctl_conf_ptr->schedrootfltr);
	fprintf(out, "SchedulerType     = %s\n",
		slurm_ctl_conf_ptr->schedtype);
	fprintf(out, "SelectType        = %s\n",
		slurm_ctl_conf_ptr->select_type);
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
	fprintf(out, "SlurmdPort        = %u\n", 
		slurm_ctl_conf_ptr->slurmd_port);
	fprintf(out, "SlurmdSpoolDir    = %s\n", 
		slurm_ctl_conf_ptr->slurmd_spooldir);
	fprintf(out, "SlurmdTimeout     = %u\n", 
		slurm_ctl_conf_ptr->slurmd_timeout);
	fprintf(out, "SLURM_CONFIG_FILE = %s\n", 
		slurm_ctl_conf_ptr->slurm_conf);
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
	fprintf(out, "TaskProlog        = %s\n",
		slurm_ctl_conf_ptr->task_prolog);
	fprintf(out, "TmpFS             = %s\n", 
		slurm_ctl_conf_ptr->tmp_fs);
	fprintf(out, "WaitTime          = %u\n", 
		slurm_ctl_conf_ptr->wait_time);
}

/*
 * slurm_load_ctl_conf - issue RPC to get slurm control configuration  
 *	information if changed since update_time 
 * IN update_time - time of current configuration data
 * IN slurm_ctl_conf_ptr - place to store slurm control configuration 
 *	pointer
 * RET 0 on success or slurm error code
 * NOTE: free the response using slurm_free_ctl_conf
 */
int
slurm_load_ctl_conf (time_t update_time, slurm_ctl_conf_t **confp)
{
	int rc;
	slurm_msg_t req_msg;
	slurm_msg_t resp_msg;
        last_update_msg_t req; 
	
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

