/*****************************************************************************\
 *  account_storage_none.c - account interface to none.
 *
 *  $Id: account_storage_none.c 13061 2008-01-22 21:23:56Z da $
 *****************************************************************************
 *  Copyright (C) 2004-2008 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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

#include "src/common/slurm_account_storage.h"

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
const char plugin_name[] = "Account storage NOT INVOKED plugin";
const char plugin_type[] = "account_storage/none";
const uint32_t plugin_version = 100;

/*
 * init() is called when the plugin is loaded, before any other functions
 * are called.  Put global initialization here.
 */
extern int init ( void )
{
	verbose("%s loaded", plugin_name);
	return SLURM_SUCCESS;
}

extern int fini ( void )
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_add_users(List user_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_add_coord(char *account, List user_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_add_accounts(List account_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_add_clusters(List cluster_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_add_associations(List association_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_modify_users(List user_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_modify_user_admin_level(
	account_admin_level_t level, List user_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_modify_accounts(List account_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_modify_clusters(List cluster_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_modify_associations(List association_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_remove_users(List user_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_remove_coord(char *account, List user_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_remove_accounts(List account_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_remove_clusters(List cluster_list)
{
	return SLURM_SUCCESS;
}

extern int account_storage_p_remove_associations(List association_list)
{
	return SLURM_SUCCESS;
}

extern List account_storage_p_get_users(List selected_users,
					void *params)
{
	return NULL;
}

extern List account_storage_p_get_accounts(List selected_accounts,
					   void *params)
{
	return NULL;
}

extern List account_storage_p_get_clusters(List selected_clusters,
					   void *params)
{
	return NULL;
}

extern List account_storage_p_get_associations(List selected_users,
					       List selected_accounts,
					       List selected_parts,
					       char *cluster,
					       void *params)
{
	return NULL;
}

extern int account_storage_p_get_hourly_usage(
	account_association_rec_t *acct_rec,
	time_t start, time_t end, void *params)
{
	int rc = SLURM_SUCCESS;

	return rc;
}

extern int account_storage_p_get_daily_usage(
	account_association_rec_t *acct_rec,
	time_t start, time_t end, void *params)
{
	int rc = SLURM_SUCCESS;

	return rc;
}

extern int account_storage_p_get_monthly_usage(
	account_association_rec_t *acct_rec,
	time_t start, time_t end, void *params)
{
	int rc = SLURM_SUCCESS;
	return rc;
}
