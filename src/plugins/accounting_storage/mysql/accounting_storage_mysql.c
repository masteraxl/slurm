/*****************************************************************************\
 *  accounting_storage_mysql.c - accounting interface to mysql.
 *
 *  $Id: accounting_storage_mysql.c 13061 2008-01-22 21:23:56Z da $
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Copyright (C) 2008 Lawrence Livermore National Security.
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
 *****************************************************************************
 * Notes on mysql configuration 
 *	Assumes mysql is installed as user root
 *	Assumes SlurmUser is configured as user slurm
 * # mysqladmin create <db_name>
 *	The <db_name> goes into slurmdbd.conf as StorageLoc
 * # mysql --user=root -p
 * mysql> GRANT ALL ON *.* TO 'slurm'@'localhost' IDENTIFIED BY PASSWORD 'pw';
 * mysql> GRANT SELECT, INSERT ON *.* TO 'slurm'@'localhost';
\*****************************************************************************/

#include <strings.h>
#include "mysql_jobacct_process.h"
#include "mysql_rollup.h"
#include "src/common/slurmdbd_defs.h"
#include "src/common/slurm_auth.h"
#include "src/common/uid.h"

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
 * the plugin (e.g., "accounting_storage" for SLURM job completion
 * logging) and <method>
 * is a description of how this plugin satisfies that application.  SLURM will
 * only load job completion logging plugins if the plugin_type string has a 
 * prefix of "accounting_storage/".
 *
 * plugin_version - an unsigned 32-bit integer giving the version number
 * of the plugin.  If major and minor revisions are desired, the major
 * version number may be multiplied by a suitable magnitude constant such
 * as 100 or 1000.  Various SLURM versions will likely require a certain
 * minimum versions for their plugins as the job accounting API 
 * matures.
 */
const char plugin_name[] = "Accounting storage MYSQL plugin";
const char plugin_type[] = "accounting_storage/mysql";
const uint32_t plugin_version = 100;

#ifdef HAVE_MYSQL

static mysql_db_info_t *mysql_db_info = NULL;
static char *mysql_db_name = NULL;

#define DEFAULT_ACCT_DB "slurm_acct_db"
#define DELETE_SEC_BACK 86400

char *acct_coord_table = "acct_coord_table";
char *acct_table = "acct_table";
char *assoc_day_table = "assoc_day_usage_table";
char *assoc_hour_table = "assoc_hour_usage_table";
char *assoc_month_table = "assoc_month_usage_table";
char *assoc_table = "assoc_table";
char *cluster_day_table = "cluster_day_usage_table";
char *cluster_hour_table = "cluster_hour_usage_table";
char *cluster_month_table = "cluster_month_usage_table";
char *cluster_table = "cluster_table";
char *event_table = "cluster_event_table";
char *job_table = "job_table";
char *qos_table = "qos_table";
char *step_table = "step_table";
char *txn_table = "txn_table";
char *user_table = "user_table";
char *last_ran_table = "last_ran_table";
char *suspend_table = "suspend_table";

static int normal_qos_id = NO_VAL;

extern int acct_storage_p_commit(mysql_conn_t *mysql_conn, bool commit);

extern int acct_storage_p_add_associations(mysql_conn_t *mysql_conn,
					   uint32_t uid, 
					   List association_list);

extern List acct_storage_p_get_associations(
	mysql_conn_t *mysql_conn, uid_t uid, 
	acct_association_cond_t *assoc_cond);

extern int acct_storage_p_get_usage(mysql_conn_t *mysql_conn, uid_t uid,
				    acct_association_rec_t *acct_assoc,
				    time_t start, time_t end);

extern int clusteracct_storage_p_get_usage(
	mysql_conn_t *mysql_conn, uid_t uid,
	acct_cluster_rec_t *cluster_rec, time_t start, time_t end);

extern List acct_storage_p_remove_coord(mysql_conn_t *mysql_conn, uint32_t uid, 
					List acct_list,
					acct_user_cond_t *user_cond);

/* This should be added to the beginning of each function to make sure
 * we have a connection to the database before we try to use it.
 */
static int _check_connection(mysql_conn_t *mysql_conn)
{
	if(!mysql_conn) {
		error("We need a connection to run this");
		return SLURM_ERROR;
	} else if(!mysql_conn->db_conn
		  || mysql_db_ping(mysql_conn->db_conn) != 0) {
		if(mysql_get_db_connection(&mysql_conn->db_conn,
					   mysql_db_name, mysql_db_info)
		   != SLURM_SUCCESS) {
			error("unable to re-connect to mysql database");
			return SLURM_ERROR;
		}
	}
	return SLURM_SUCCESS;
}

static int _setup_association_limits(acct_association_rec_t *assoc,
				     char **in_cols, char **in_vals,
				     char **in_extra, bool get_qos)
{
	char *cols = (*in_cols), *vals = (*in_vals), *extra = (*in_extra);
	
	if(!assoc)
		return SLURM_ERROR;

	if((int)assoc->fairshare >= 0) {
		xstrcat(cols, ", fairshare");
		xstrfmtcat(vals, ", %u", assoc->fairshare);
		xstrfmtcat(extra, ", fairshare=%u", assoc->fairshare);
	} else if ((int)assoc->fairshare == INFINITE) {
		xstrcat(cols, ", fairshare");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", fairshare=NULL");		
	}
	
	if((int)assoc->grp_cpu_hours >= 0) {
		xstrcat(cols, ", grp_cpu_hours");
		xstrfmtcat(vals, ", %llu", assoc->grp_cpu_hours);
		xstrfmtcat(extra, ", grp_cpu_mins=%llu", assoc->grp_cpu_hours);
	} else if((int)assoc->grp_cpu_hours == INFINITE) {
		xstrcat(cols, ", grp_cpu_hours");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", grp_cpu_mins=NULL");
	}
		
	if((int)assoc->grp_cpus >= 0) {
		xstrcat(cols, ", grp_cpus");
		xstrfmtcat(vals, ", %u", assoc->grp_cpus);
		xstrfmtcat(extra, ", grp_cpus=%u", assoc->grp_cpus);
	} else if((int)assoc->grp_cpus == INFINITE) {
		xstrcat(cols, ", grp_cpus");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", grp_cpus=NULL");
	}

	if((int)assoc->grp_jobs >= 0) {
		xstrcat(cols, ", grp_jobs");
		xstrfmtcat(vals, ", %u", assoc->grp_jobs);
		xstrfmtcat(extra, ", grp_jobs=%u", assoc->grp_jobs);
	} else if((int)assoc->grp_jobs == INFINITE) {
		xstrcat(cols, ", grp_jobs");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", grp_jobs=NULL");
	}

	if((int)assoc->grp_nodes >= 0) {
		xstrcat(cols, ", grp_nodes");
		xstrfmtcat(vals, ", %u", assoc->grp_nodes);
		xstrfmtcat(extra, ", grp_nodes=%u", assoc->grp_nodes);
	} else if((int)assoc->grp_nodes == INFINITE) {
		xstrcat(cols, ", grp_nodes");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", grp_nodes=NULL");
	}

	if((int)assoc->grp_submit_jobs >= 0) {
		xstrcat(cols, ", grp_submit_jobs");
		xstrfmtcat(vals, ", %u",
			   assoc->grp_submit_jobs);
		xstrfmtcat(extra, ", grp_submit_jobs=%u",
			   assoc->grp_submit_jobs);
	} else if((int)assoc->grp_submit_jobs == INFINITE) {
		xstrcat(cols, ", grp_submit_jobs");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", grp_submit_jobs=NULL");
	}

	if((int)assoc->grp_wall >= 0) {
		xstrcat(cols, ", grp_wall");
		xstrfmtcat(vals, ", %u", assoc->grp_wall);
		xstrfmtcat(extra, ", grp_wall=%u",
			   assoc->grp_wall);
	} else if((int)assoc->grp_wall == INFINITE) {
		xstrcat(cols, ", grp_wall");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", grp_wall=NULL");
	}

	if((int)assoc->max_cpu_mins_pj >= 0) {
		xstrcat(cols, ", max_cpu_mins_per_job");
		xstrfmtcat(vals, ", %llu", assoc->max_cpu_mins_pj);
		xstrfmtcat(extra, ", max_cpu_mins_per_job=%u",
			   assoc->max_cpu_mins_pj);
	} else if((int)assoc->max_cpu_mins_pj == INFINITE) {
		xstrcat(cols, ", max_cpu_mins_per_job");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", max_cpu_mins_per_job=NULL");
	}

	if((int)assoc->max_cpus_pj >= 0) {
		xstrcat(cols, ", max_cpus_per_job");
		xstrfmtcat(vals, ", %u", assoc->max_cpus_pj);
		xstrfmtcat(extra, ", max_cpus_per_job=%u",
			   assoc->max_cpus_pj);
	} else if((int)assoc->max_cpus_pj == INFINITE) {
		xstrcat(cols, ", max_cpus_per_job");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", max_cpus_per_job=NULL");
	}
		
	if((int)assoc->max_jobs >= 0) {
		xstrcat(cols, ", max_jobs");
		xstrfmtcat(vals, ", %u", assoc->max_jobs);
		xstrfmtcat(extra, ", max_jobs=%u",
			   assoc->max_jobs);
	} else if((int)assoc->max_jobs == INFINITE) {
		xstrcat(cols, ", max_jobs");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", max_jobs=NULL");		
	}

	if((int)assoc->max_nodes_pj >= 0) {
		xstrcat(cols, ", max_nodes_per_job");
		xstrfmtcat(vals, ", %u", assoc->max_nodes_pj);
		xstrfmtcat(extra, ", max_nodes_per_job=%u",
			   assoc->max_nodes_pj);
	} else if((int)assoc->max_nodes_pj == INFINITE) {
		xstrcat(cols, ", max_nodes_per_job");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", max_nodes_per_job=NULL");
	}

	if((int)assoc->max_submit_jobs >= 0) {
		xstrcat(cols, ", max_submit_jobs");
		xstrfmtcat(vals, ", %u", assoc->max_submit_jobs);
		xstrfmtcat(extra, ", max_submit_jobs=%u",
			   assoc->max_submit_jobs);
	} else if((int)assoc->max_submit_jobs == INFINITE) {
		xstrcat(cols, ", max_submit_jobs");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", max_submit_jobs=NULL");
	}

	if((int)assoc->max_wall_pj >= 0) {
		xstrcat(cols, ", max_wall_duration_per_job");
		xstrfmtcat(vals, ", %u", assoc->max_wall_pj);
		xstrfmtcat(extra, ", max_wall_duration_per_job=%u",
			   assoc->max_wall_pj);
	} else if((int)assoc->max_wall_pj == INFINITE) {
		xstrcat(cols, ", max_wall_duration_per_job");
		xstrcat(vals, ", NULL");
		xstrcat(extra, ", max_wall_duration_per_job=NULL");
	}

	if(assoc->qos_list && list_count(assoc->qos_list)) {
		char *qos_val = NULL;
		char *tmp_char = NULL;
		ListIterator qos_itr = list_iterator_create(assoc->qos_list);
			
		xstrcat(cols, ", qos");
			
		while((tmp_char = list_next(qos_itr))) 
			xstrfmtcat(qos_val, ",%s", tmp_char);
			
		list_iterator_destroy(qos_itr);
			
		xstrfmtcat(vals, ", '%s'", qos_val); 		
		xstrfmtcat(extra, ", qos='%s'", qos_val); 
		xfree(qos_val);
	} else if(get_qos && normal_qos_id != NO_VAL) { 
		/* Add normal qos to the account */
		xstrcat(cols, ", qos");
		xstrfmtcat(vals, ", ',%d'", normal_qos_id);
		xstrfmtcat(extra, ", qos=',%d'", normal_qos_id);
	}

	return SLURM_SUCCESS;

}

static int _setup_association_cond_limits(acct_association_cond_t *assoc_cond,
					  char **in_extra)
{
	char *extra = (*in_extra);
	int set = 0;
	ListIterator itr = NULL;
	char *object = NULL;

	if(!assoc_cond)
		return 0;

	if(assoc_cond->acct_list && list_count(assoc_cond->acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "acct='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->cluster_list && list_count(assoc_cond->cluster_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->cluster_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "cluster='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->fairshare_list
	   && list_count(assoc_cond->fairshare_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->fairshare_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "fairshare='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->grp_cpu_hours_list
	   && list_count(assoc_cond->grp_cpu_hours_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->grp_cpu_hours_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "grp_cpu_hours='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->grp_cpus_list
	   && list_count(assoc_cond->grp_cpus_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->grp_cpus_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "grp_cpus='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->grp_jobs_list
	   && list_count(assoc_cond->grp_jobs_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->grp_jobs_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "grp_jobs='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->grp_nodes_list
	   && list_count(assoc_cond->grp_nodes_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->grp_nodes_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "grp_nodes='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->grp_submit_jobs_list
	   && list_count(assoc_cond->grp_submit_jobs_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->grp_submit_jobs_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "grp_submit_jobs='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->grp_wall_list
	   && list_count(assoc_cond->grp_wall_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->grp_wall_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "grp_wall='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->max_cpu_mins_pj_list
	   && list_count(assoc_cond->max_cpu_mins_pj_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->max_cpu_mins_pj_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "max_cpu_mins_pj='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->max_cpus_pj_list
	   && list_count(assoc_cond->max_cpus_pj_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->max_cpus_pj_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "max_cpus_pj='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->max_jobs_list
	   && list_count(assoc_cond->max_jobs_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->max_jobs_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "max_jobs='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->max_nodes_pj_list
	   && list_count(assoc_cond->max_nodes_pj_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->max_nodes_pj_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "max_nodes_pj='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->max_submit_jobs_list
	   && list_count(assoc_cond->max_submit_jobs_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->max_submit_jobs_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "max_submit_jobs='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->max_wall_pj_list
	   && list_count(assoc_cond->max_wall_pj_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->max_wall_pj_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "max_wall_pj='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->user_list && list_count(assoc_cond->user_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->user_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "user='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	} else if (!assoc_cond->user_list) {
		debug4("no user specified looking at accounts");
		xstrcat(extra, " && user = '' ");
	} else {
		debug4("no user specified looking at users");
		xstrcat(extra, " && user != '' ");
	}

	if(assoc_cond->partition_list 
	   && list_count(assoc_cond->partition_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->partition_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "partition='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(assoc_cond->id_list && list_count(assoc_cond->id_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->id_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "id=%s", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(assoc_cond->qos_list && list_count(assoc_cond->qos_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->qos_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, 
				   "(qos like '%%,%s' || qos like '%%,%s,%%')",
				   object, object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(assoc_cond->parent_acct_list
	   && list_count(assoc_cond->parent_acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(assoc_cond->parent_acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "parent_acct=%s", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	return set;
}
/* This function will take the object given and free it later so it
 * needed to be removed from a list if in one before 
 */
static int _addto_update_list(List update_list, acct_update_type_t type,
			      void *object)
{
	acct_update_object_t *update_object = NULL;
	ListIterator itr = NULL;
	if(!update_list) {
		error("no update list given");
		return SLURM_ERROR;
	}

	itr = list_iterator_create(update_list);
	while((update_object = list_next(itr))) {
		if(update_object->type == type)
			break;
	}
	list_iterator_destroy(itr);

	if(update_object) {
		list_append(update_object->objects, object);
		return SLURM_SUCCESS;
	} 
	update_object = xmalloc(sizeof(acct_update_object_t));

	list_append(update_list, update_object);

	update_object->type = type;
	
	switch(type) {
	case ACCT_MODIFY_USER:
	case ACCT_ADD_USER:
	case ACCT_REMOVE_USER:
	case ACCT_ADD_COORD:
	case ACCT_REMOVE_COORD:
		update_object->objects = list_create(destroy_acct_user_rec);
		break;
	case ACCT_ADD_ASSOC:
	case ACCT_MODIFY_ASSOC:
	case ACCT_REMOVE_ASSOC:
		update_object->objects = list_create(
			destroy_acct_association_rec);
		break;
	case ACCT_ADD_QOS:
	case ACCT_REMOVE_QOS:
		update_object->objects = list_create(
			destroy_acct_qos_rec);
		break;
	case ACCT_UPDATE_NOTSET:
	default:
		error("unknown type set in update_object: %d", type);
		return SLURM_ERROR;
	}
	list_append(update_object->objects, object);
	return SLURM_SUCCESS;
}

/* This should take care of all the lft and rgts when you move an
 * account.  This handles deleted associations also.
 */
static int _move_account(mysql_conn_t *mysql_conn, uint32_t lft, uint32_t rgt,
			 char *cluster,
			 char *id, char *parent)
{
	int rc = SLURM_SUCCESS;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	int par_left = 0;
	int diff = 0;
	int width = 0;
	char *query = xstrdup_printf(
		"SELECT lft from %s " 
		"where cluster='%s' && acct='%s' && user='';",
		assoc_table,
		cluster, parent);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);
	if(!(row = mysql_fetch_row(result))) {
		error("no row");
		mysql_free_result(result);
		return SLURM_ERROR;
	}
	par_left = atoi(row[0]);
	mysql_free_result(result);

	diff = ((par_left + 1) - lft);

	if(diff == 0) {
		debug3("Trying to move association to the same position?  "
		       "Nothing to do.");
		return rc;
	}
	
	width = (rgt - lft + 1);

	/* every thing below needs to be a %d not a %u because we are
	   looking for -1 */
	xstrfmtcat(query,
		   "update %s set deleted = deleted + 2, "
		   "lft = lft + %d, rgt = rgt + %d "
		   "WHERE lft BETWEEN %d AND %d;",
		   assoc_table, diff, diff, lft, rgt);

	xstrfmtcat(query,
		   "UPDATE %s SET rgt = rgt + %d WHERE "
		   "rgt > %d && deleted < 2;"
		   "UPDATE %s SET lft = lft + %d WHERE "
		   "lft > %d && deleted < 2;",
		   assoc_table, width,
		   par_left,
		   assoc_table, width,
		   par_left);

	xstrfmtcat(query,
		   "UPDATE %s SET rgt = rgt - %d WHERE "
		   "(%d < 0 && rgt > %d && deleted < 2) "
		   "|| (%d > 0 && rgt > %d);"
		   "UPDATE %s SET lft = lft - %d WHERE "
		   "(%d < 0 && lft > %d && deleted < 2) "
		   "|| (%d > 0 && lft > %d);",
		   assoc_table, width,
		   diff, rgt,
		   diff, lft,
		   assoc_table, width,
		   diff, rgt,
		   diff, lft);

	xstrfmtcat(query,
		   "update %s set deleted = deleted - 2 WHERE deleted > 1;",
		   assoc_table);
	xstrfmtcat(query,
		   "update %s set parent_acct='%s' where id = %s;",
		   assoc_table, parent, id);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);

	return rc;
}


/* This code will move an account from one parent to another.  This
 * should work either way in the tree.  (i.e. move child to be parent
 * of current parent, and parent to be child of child.)
 */
static int _move_parent(mysql_conn_t *mysql_conn, uid_t uid,
			uint32_t lft, uint32_t rgt,
			char *cluster,
			char *id, char *old_parent, char *new_parent)
{
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *query = NULL;
	int rc = SLURM_SUCCESS;
	List assoc_list = NULL;
	ListIterator itr = NULL;
	acct_association_rec_t *assoc = NULL;
		
	/* first we need to see if we are going to make a child of this
	 * account the new parent.  If so we need to move that child to this
	 * accounts parent and then do the move.
	 */
	query = xstrdup_printf(
		"select id, lft, rgt from %s where lft between %d and %d "
		"&& acct='%s' && user='' order by lft;",
		assoc_table, lft, rgt,
		new_parent);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = 
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	if((row = mysql_fetch_row(result))) {
		debug4("%s(%s) %s,%s is a child of %s",
		       new_parent, row[0], row[1], row[2], id);
		rc = _move_account(mysql_conn, atoi(row[1]), atoi(row[2]),
				   cluster, row[0], old_parent);
	}

	mysql_free_result(result);

	if(rc == SLURM_ERROR) 
		return rc;
	
	/* now move the one we wanted to move in the first place 
	 * We need to get the new lft and rgts though since they may
	 * have changed.
	 */
	query = xstrdup_printf(
		"select lft, rgt from %s where id=%s;",
		assoc_table, id);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = 
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	if((row = mysql_fetch_row(result))) {
		rc = _move_account(mysql_conn, atoi(row[0]), atoi(row[1]),
				   cluster, id, new_parent);
	} else {
		error("can't find parent? we were able to a second ago.");
		rc = SLURM_ERROR;
	}
	mysql_free_result(result);

	if(rc == SLURM_ERROR) 
		return rc;
	
	/* now we need to send the update of the new parents and
	 * limits, so just to be safe, send the whole tree
	 */
	assoc_list = acct_storage_p_get_associations(mysql_conn, uid, NULL);
	/* NOTE: you can not use list_pop, or list_push
	   anywhere either, since mysql is
	   exporting something of the same type as a macro,
	   which messes everything up (my_list.h is the bad boy).
	   So we are just going to delete each item as it
	   comes out since we are moving it to the update_list.
	*/
	itr = list_iterator_create(assoc_list);
	while((assoc = list_next(itr))) {
		if(_addto_update_list(mysql_conn->update_list, 
				      ACCT_MODIFY_ASSOC,
				      assoc) == SLURM_SUCCESS) 
			list_remove(itr);
	}
	list_iterator_destroy(itr);
	list_destroy(assoc_list);
	return rc;
}

/* Let me know if the last statement had rows that were affected.
 */
static int _last_affected_rows(MYSQL *mysql_db)
{
	int status=0, rows=0;
	MYSQL_RES *result = NULL;

	do {
		result = mysql_store_result(mysql_db);
		if (result) 
			mysql_free_result(result);
		else 
			if (mysql_field_count(mysql_db) == 0) {
				status = mysql_affected_rows(mysql_db);
				if(status > 0)
					rows = status;
			}
		if ((status = mysql_next_result(mysql_db)) > 0)
			debug3("Could not execute statement\n");
	} while (status == 0);
	
	return rows;
}

/* This is called by most modify functions to alter the table and
 * insert a new line in the transaction table.
 */
static int _modify_common(mysql_conn_t *mysql_conn,
			  uint16_t type,
			  time_t now,
			  char *user_name,
			  char *table,
			  char *cond_char,
			  char *vals) 
{
	char *query = NULL;
	int rc = SLURM_SUCCESS;

	xstrfmtcat(query, 
		   "update %s set mod_time=%d%s "
		   "where deleted=0 && %s;",
		   table, now, vals,
		   cond_char);
	xstrfmtcat(query, 	
		   "insert into %s "
		   "(timestamp, action, name, actor, info) "
		   "values (%d, %d, \"%s\", '%s', \"%s\");",
		   txn_table,
		   now, type, cond_char, user_name, vals);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);		
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);

	if(rc != SLURM_SUCCESS) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		
		return SLURM_ERROR;
	}

	return SLURM_SUCCESS;
}

/* Used to get all the users inside a lft and rgt set.  This is just
 * to send the user all the associations that are being modified from
 * a previous change to it's parent.    
 */
static int _modify_unset_users(mysql_conn_t *mysql_conn,
			       acct_association_rec_t *assoc,
			       char *acct,
			       uint32_t lft, uint32_t rgt,
			       List ret_list)
{
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *query = NULL, *object = NULL;
	int i;

	char *assoc_req_inx[] = {
		"id",
		"user",
		"acct",
		"cluster",
		"partition",
		"max_jobs",
		"max_nodes_per_job",
		"max_wall_duration_per_job",
		"max_cpu_mins_per_job",
		"lft",
		"rgt"
	};
	
	enum {
		ASSOC_ID,
		ASSOC_USER,
		ASSOC_ACCT,
		ASSOC_CLUSTER,
		ASSOC_PART,
		ASSOC_MJ,
		ASSOC_MNPJ,
		ASSOC_MWPJ,
		ASSOC_MCPJ,
		ASSOC_LFT,
		ASSOC_RGT,
		ASSOC_COUNT
	};

	if(!ret_list || !acct)
		return SLURM_ERROR;

	for(i=0; i<ASSOC_COUNT; i++) {
		if(i) 
			xstrcat(object, ", ");
		xstrcat(object, assoc_req_inx[i]);
	}

	/* We want all the sub accounts and user accounts */
	query = xstrdup_printf("select distinct %s from %s where deleted=0 "
			       "&& lft between %d and %d && "
			       "((user = '' && parent_acct = '%s') || "
			       "(user != '' && acct = '%s')) "
			       "order by lft;",
			       object, assoc_table, lft, rgt, acct, acct);
	xfree(object);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result =
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	while((row = mysql_fetch_row(result))) {
		acct_association_rec_t *mod_assoc = NULL;
		int modified = 0;

		mod_assoc = xmalloc(sizeof(acct_association_rec_t));
		mod_assoc->id = atoi(row[ASSOC_ID]);

		if(!row[ASSOC_MJ] && assoc->max_jobs != NO_VAL) {
			mod_assoc->max_jobs = assoc->max_jobs;
			modified = 1;
		} else
			mod_assoc->max_jobs = NO_VAL;
		
		if(!row[ASSOC_MNPJ] &&
		   assoc->max_nodes_pj != NO_VAL) {
			mod_assoc->max_nodes_pj =
				assoc->max_nodes_pj;
			modified = 1;
		} else 
			mod_assoc->max_nodes_pj = NO_VAL;

		
		if(!row[ASSOC_MWPJ] && 
		   assoc->max_wall_pj != NO_VAL) {
			mod_assoc->max_wall_pj =
				assoc->max_wall_pj;
			modified = 1;
		} else 
			mod_assoc->max_wall_pj = NO_VAL;
					
		if(!row[ASSOC_MCPJ] && 
		   assoc->max_cpu_mins_pj != NO_VAL) {
			mod_assoc->max_cpu_mins_pj = 
				assoc->max_cpu_mins_pj;
			modified = 1;
		} else
			mod_assoc->max_cpu_mins_pj = NO_VAL;
		
		/* We only want to add those that are modified here */
		if(modified) {
			/* Since we aren't really changing this non
			 * user association we don't want to send it.
			 */
			if(!row[ASSOC_USER][0]) {
				/* This is a sub account so run it
				 * through as if it is a parent.
				 */
				_modify_unset_users(mysql_conn,
						    mod_assoc,
						    row[ASSOC_ACCT],
						    atoi(row[ASSOC_LFT]),
						    atoi(row[ASSOC_RGT]),
						    ret_list);
				destroy_acct_association_rec(mod_assoc);
				continue;
			}
			/* We do want to send all user accounts though */
			mod_assoc->fairshare = NO_VAL;
			if(row[ASSOC_PART][0]) { 
				// see if there is a partition name
				object = xstrdup_printf(
					"C = %-10s A = %-20s U = %-9s P = %s",
					row[ASSOC_CLUSTER], row[ASSOC_ACCT],
					row[ASSOC_USER], row[ASSOC_PART]);
			} else {
				object = xstrdup_printf(
					"C = %-10s A = %-20s U = %-9s",
					row[ASSOC_CLUSTER], 
					row[ASSOC_ACCT], 
					row[ASSOC_USER]);
			}
			
			list_append(ret_list, object);
			
			if(_addto_update_list(mysql_conn->update_list, 
					      ACCT_MODIFY_ASSOC,
					      mod_assoc) != SLURM_SUCCESS) 
				error("couldn't add to the update list");
		} else {
			xfree(mod_assoc);
		}
	}
	mysql_free_result(result);

	return SLURM_SUCCESS;
}

/* this function is here to see if any of what we are trying to remove
 * has jobs that are or were once running.  So if we have jobs and the
 * object is less than a day old we don't want to delete it only set
 * the deleted flag.
 */
static bool _check_jobs_before_remove(mysql_conn_t *mysql_conn,
				      char *assoc_char)
{
	char *query = NULL;
	bool rc = 0;
	MYSQL_RES *result = NULL;

	query = xstrdup_printf("select t0.associd from %s as t0, %s as t1, "
			       "%s as t2 where t1.lft between "
			       "t2.lft and t2.rgt && (%s)"
			       "and t0.associd=t1.id limit 1;",
			       job_table, assoc_table, assoc_table, 
			       assoc_char);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return rc;
	}
	xfree(query);

	if(mysql_num_rows(result)) {
		debug4("We have jobs for this combo");
		rc = true;
	}

	mysql_free_result(result);
	return rc;
}

/* Same as above but for associations instead of other tables */
static bool _check_jobs_before_remove_assoc(mysql_conn_t *mysql_conn,
					    char *assoc_char)
{
	char *query = NULL;
	bool rc = 0;
	MYSQL_RES *result = NULL;

	query = xstrdup_printf("select t1.associd from %s as t1, "
			       "%s as t2 where (%s)"
			       "and t1.associd=t2.id limit 1;",
			       job_table, assoc_table, 
			       assoc_char);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return rc;
	}
	xfree(query);

	if(mysql_num_rows(result)) {
		debug4("We have jobs for this combo");
		rc = true;
	}

	mysql_free_result(result);
	return rc;
}

/* Every option in assoc_char should have a 't1.' infront of it. */
static int _remove_common(mysql_conn_t *mysql_conn,
			  uint16_t type,
			  time_t now,
			  char *user_name,
			  char *table,
			  char *name_char,
			  char *assoc_char) 
{
	int rc = SLURM_SUCCESS;
	char *query = NULL;
	char *loc_assoc_char = NULL;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	time_t day_old = now - DELETE_SEC_BACK;
	bool has_jobs = false;

	/* If we have jobs associated with this we do not want to
	 * really delete it for accounting purposes.  This is for
	 * corner cases most of the time this won't matter.
	 */
	if(table == acct_coord_table || table == qos_table) {
		/* This doesn't apply for these tables since we are
		 * only looking for association type tables.
		 */
	} else if(table != assoc_table) {
		has_jobs = _check_jobs_before_remove(mysql_conn, assoc_char);
	} else {
		has_jobs = _check_jobs_before_remove_assoc(mysql_conn,
							   name_char);	
	}
	/* we want to remove completely all that is less than a day old */
	if(!has_jobs && table != assoc_table) {
		query = xstrdup_printf("delete from %s where creation_time>%d "
				       "&& (%s);",
				       table, day_old, name_char);
	}

	if(table != assoc_table)
		xstrfmtcat(query,
			   "update %s set mod_time=%d, deleted=1 "
			   "where deleted=0 && (%s);",
			   table, now, name_char);
	
	xstrfmtcat(query, 	
		   "insert into %s (timestamp, action, name, actor) "
		   "values (%d, %d, \"%s\", '%s');",
		   txn_table,
		   now, type, name_char, user_name);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	if(rc != SLURM_SUCCESS) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		
		return SLURM_ERROR;
	}
	
	if(table == qos_table) {
		/* remove this qos from all the users/accts that have it
		 */
		xstrfmtcat(query,
			   "update %s set mod_time=%d, %s "
			   "where deleted=0;",
			   assoc_table, now, assoc_char);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			if(mysql_conn->rollback) {
				mysql_db_rollback(mysql_conn->db_conn);
			}
			list_flush(mysql_conn->update_list);
			
			return SLURM_ERROR;
		}
		/* now get what we changed and set the update */
		xstrfmtcat(query,
			   "select id, qos from %s where "
			   "mod_time=%d and deleted=0;",
			   assoc_table, now);
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			if(mysql_conn->rollback) {
				mysql_db_rollback(mysql_conn->db_conn);
			}
			list_flush(mysql_conn->update_list);
			
			return SLURM_ERROR;
		}
		
		rc = 0;
		while((row = mysql_fetch_row(result))) {
			acct_association_rec_t *assoc_rec = 
				xmalloc(sizeof(acct_association_rec_t));
			assoc_rec->id = atoi(row[0]);
			assoc_rec->qos_list = list_create(slurm_destroy_char);
			slurm_addto_char_list(assoc_rec->qos_list, row[1]);
			_addto_update_list(mysql_conn->update_list,
					   ACCT_MODIFY_ASSOC,
					   assoc_rec);
		}
		mysql_free_result(result);
		
		return SLURM_SUCCESS;
	} else if(table == acct_coord_table)
		return SLURM_SUCCESS;

	/* mark deleted=1 or remove completely the
	   accounting tables
	*/
	if(table != assoc_table) {
		if(!assoc_char) {
			error("no assoc_char");
			if(mysql_conn->rollback) {
				mysql_db_rollback(mysql_conn->db_conn);
			}
			list_flush(mysql_conn->update_list);
			return SLURM_ERROR;
		}

		/* If we are doing this on an assoc_table we have
		   already done this, so don't */
/* 		query = xstrdup_printf("select lft, rgt " */
/* 				       "from %s as t2 where %s order by lft;", */
/* 				       assoc_table, assoc_char); */
		query = xstrdup_printf("select distinct t1.id "
				       "from %s as t1, %s as t2 "
				       "where (%s) && t1.lft between "
				       "t2.lft and t2.rgt && t1.deleted=0 "
				       " && t2.deleted=0;",
				       assoc_table, assoc_table, assoc_char);
		
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			if(mysql_conn->rollback) {
				mysql_db_rollback(mysql_conn->db_conn);
			}
			list_flush(mysql_conn->update_list);
			return SLURM_ERROR;
		}
		xfree(query);

		rc = 0;
		loc_assoc_char = NULL;
		while((row = mysql_fetch_row(result))) {
			acct_association_rec_t *rem_assoc = NULL;
			if(!rc) {
				xstrfmtcat(loc_assoc_char, "id=%s", row[0]);
				rc = 1;
			} else {
				xstrfmtcat(loc_assoc_char,
					   " || id=%s", row[0]);
			}
			rem_assoc = xmalloc(sizeof(acct_association_rec_t));
			rem_assoc->id = atoi(row[0]);
			if(_addto_update_list(mysql_conn->update_list, 
					      ACCT_REMOVE_ASSOC,
					      rem_assoc) != SLURM_SUCCESS) 
				error("couldn't add to the update list");
		}
		mysql_free_result(result);
	} else 
		loc_assoc_char = assoc_char;

	if(!loc_assoc_char) {
		debug2("No associations with object being deleted\n");
		return rc;
	}

	/* We should not have to delete from usage table, only flag since we
	 * only delete things that are typos.
	 */ 
	xstrfmtcat(query,
		   "update %s set mod_time=%d, deleted=1 where (%s);"
		   "update %s set mod_time=%d, deleted=1 where (%s);"
		   "update %s set mod_time=%d, deleted=1 where (%s);",
		   assoc_day_table, now, loc_assoc_char,
		   assoc_hour_table, now, loc_assoc_char,
		   assoc_month_table, now, loc_assoc_char);

	debug3("%d(%d) query\n%s %d",
	       mysql_conn->conn, __LINE__, query, strlen(query));
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	if(rc != SLURM_SUCCESS) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		return SLURM_ERROR;
	}

	/* If we have jobs that have ran don't go through the logic of
	 * removing the associations. Since we may want them for
	 * reports in the future since jobs had ran.
	 */
	if(has_jobs)
		goto just_update;

	/* remove completely all the associations for this added in the last
	 * day, since they are most likely nothing we really wanted in
	 * the first place.
	 */
	query = xstrdup_printf("select id from %s as t1 where "
			       "creation_time>%d && (%s);",
			       assoc_table, day_old, loc_assoc_char);
	
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		return SLURM_ERROR;
	}
	xfree(query);

	while((row = mysql_fetch_row(result))) {
		MYSQL_RES *result2 = NULL;
		MYSQL_ROW row2;
		
		/* we have to do this one at a time since the lft's and rgt's
		   change. If you think you need to remove this make
		   sure your new way can handle changing lft and rgt's
		   in the association. */
		xstrfmtcat(query,
			   "SELECT lft, rgt, (rgt - lft + 1) "
			   "FROM %s WHERE id = %s;",
			   assoc_table, row[0]);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		if(!(result2 = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			rc = SLURM_ERROR;
			break;
		}
		xfree(query);
		if(!(row2 = mysql_fetch_row(result2))) {
			mysql_free_result(result2);
			continue;
		}

		xstrfmtcat(query,
			   "delete quick from %s where lft between "
			   "%s AND %s;",
			   assoc_table,
			   row2[0], row2[1]);
		
		xstrfmtcat(query,
			   "UPDATE %s SET rgt = rgt - %s WHERE "
			   "rgt > %s;"
			   "UPDATE %s SET lft = lft - %s WHERE "
			   "lft > %s;",
			   assoc_table, row2[2],
			   row2[1],
			   assoc_table, row2[2],
			   row2[1]);
		
		mysql_free_result(result2);

		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("couldn't remove assoc");
			break;
		}
	}
	mysql_free_result(result);
	if(rc == SLURM_ERROR) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		return rc;
	}

just_update:
	/* now update the associations themselves that are still
	 * around clearing all the limits since if we add them back
	 * we don't want any residue from past associations lingering
	 * around.
	 */
	query = xstrdup_printf("update %s as t1 set mod_time=%d, deleted=1, "
			       "fairshare=1, max_jobs=NULL, "
			       "max_nodes_per_job=NULL, "
			       "max_wall_duration_per_job=NULL, "
			       "max_cpu_mins_per_job=NULL "
			       "where (%s);",
			       assoc_table, now,
			       loc_assoc_char);

	if(table != assoc_table)
		xfree(loc_assoc_char);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	if(rc != SLURM_SUCCESS) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
	}
	
	return rc;
}

/* Fill in all the users that are coordinator for this account.  This
 * will fill in if there are coordinators from a parent account also.
 */
static int _get_account_coords(mysql_conn_t *mysql_conn, 
			       acct_account_rec_t *acct)
{
	char *query = NULL;
	acct_coord_rec_t *coord = NULL;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	
	if(!acct) {
		error("We need a account to fill in.");
		return SLURM_ERROR;
	}

	if(!acct->coordinators)
		acct->coordinators = list_create(destroy_acct_coord_rec);
			
	query = xstrdup_printf(
		"select user from %s where acct='%s' && deleted=0",
		acct_coord_table, acct->name);
			
	if(!(result =
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);
	while((row = mysql_fetch_row(result))) {
		coord = xmalloc(sizeof(acct_coord_rec_t));
		list_append(acct->coordinators, coord);
		coord->name = xstrdup(row[0]);
		coord->direct = 1;
	}
	mysql_free_result(result);
		
	query = xstrdup_printf("select distinct t0.user from %s as t0, "
			       "%s as t1, %s as t2 where t0.acct=t1.acct && "
			       "t1.lft<t2.lft && t1.rgt>t2.lft && "
			       "t1.user='' && t2.acct='%s' && t1.acct!='%s' && "
			       "!t0.deleted;",
			       acct_coord_table, assoc_table, assoc_table,
			       acct->name, acct->name);
	if(!(result =
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);
	while((row = mysql_fetch_row(result))) {
		coord = xmalloc(sizeof(acct_coord_rec_t));
		list_append(acct->coordinators, coord);
		coord->name = xstrdup(row[0]);
		coord->direct = 0;
	}
	return SLURM_SUCCESS;
}

/* Fill in all the accounts this user is coordinator over.  This
 * will fill in all the sub accounts they are coordinator over also.
 */
static int _get_user_coords(mysql_conn_t *mysql_conn, acct_user_rec_t *user)
{
	char *query = NULL;
	acct_coord_rec_t *coord = NULL;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	ListIterator itr = NULL;

	if(!user) {
		error("We need a user to fill in.");
		return SLURM_ERROR;
	}

	if(!user->coord_accts)
		user->coord_accts = list_create(destroy_acct_coord_rec);
			
	query = xstrdup_printf(
		"select acct from %s where user='%s' && deleted=0",
		acct_coord_table, user->name);
			
	if(!(result =
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);
	while((row = mysql_fetch_row(result))) {
		coord = xmalloc(sizeof(acct_coord_rec_t));
		list_append(user->coord_accts, coord);
		coord->name = xstrdup(row[0]);
		coord->direct = 1;
		if(query) 
			xstrcat(query, " || ");
		else 
			query = xstrdup_printf(
				"select distinct t1.acct from "
				"%s as t1, %s as t2 where t1.deleted=0 && ",
				assoc_table, assoc_table);
		/* Make sure we don't get the same
		 * account back since we want to keep
		 * track of the sub-accounts.
		 */
		xstrfmtcat(query, "(t2.acct='%s' "
			   "&& t1.lft between t2.lft "
			   "and t2.rgt && t1.user='' "
			   "&& t1.acct!='%s')",
			   coord->name, coord->name);
	}
	mysql_free_result(result);

	if(query) {
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			return SLURM_ERROR;
		}
		xfree(query);

		itr = list_iterator_create(user->coord_accts);
		while((row = mysql_fetch_row(result))) {

			while((coord = list_next(itr))) {
				if(!strcmp(coord->name, row[0]))
					break;
			}
			list_iterator_reset(itr);
			if(coord) 
				continue;
					
			coord = xmalloc(sizeof(acct_coord_rec_t));
			list_append(user->coord_accts, coord);
			coord->name = xstrdup(row[0]);
			coord->direct = 0;
		}
		list_iterator_destroy(itr);
		mysql_free_result(result);
	}
	return SLURM_SUCCESS;
}

/* Used in job functions for getting the database index based off the
 * submit time, job and assoc id.  0 is returned if none is found
 */
static int _get_db_index(MYSQL *db_conn, 
			 time_t submit, uint32_t jobid, uint32_t associd)
{
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	int db_index = 0;
	char *query = xstrdup_printf("select id from %s where "
				     "submit=%d and jobid=%u and associd=%u",
				     job_table, (int)submit, jobid, associd);

	if(!(result = mysql_db_query_ret(db_conn, query, 0))) {
		xfree(query);
		return 0;
	}
	xfree(query);

	row = mysql_fetch_row(result);
	if(!row) {
		mysql_free_result(result);
		error("We can't get a db_index for this combo, "
		      "submit=%d and jobid=%u and associd=%u.",
		      (int)submit, jobid, associd);
		return 0;
	}
	db_index = atoi(row[0]);
	mysql_free_result(result);
	
	return db_index;
}

static mysql_db_info_t *_mysql_acct_create_db_info()
{
	mysql_db_info_t *db_info = xmalloc(sizeof(mysql_db_info_t));
	db_info->port = slurm_get_accounting_storage_port();
	if(!db_info->port) 
		db_info->port = 3306;
	db_info->host = slurm_get_accounting_storage_host();	
	db_info->user = slurm_get_accounting_storage_user();	
	db_info->pass = slurm_get_accounting_storage_pass();	
	return db_info;
}

/* Any time a new table is added set it up here */
static int _mysql_acct_check_tables(MYSQL *db_conn)
{
	int rc = SLURM_SUCCESS;
	storage_field_t acct_coord_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "acct", "tinytext not null" },
		{ "user", "tinytext not null" },
		{ NULL, NULL}		
	};

	storage_field_t acct_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "name", "tinytext not null" },
		{ "description", "text not null" },
		{ "organization", "text not null" },
		{ NULL, NULL}		
	};

	storage_field_t assoc_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "id", "int not null auto_increment" },
		{ "user", "tinytext not null default ''" },
		{ "acct", "tinytext not null" },
		{ "cluster", "tinytext not null" },
		{ "partition", "tinytext not null default ''" },
		{ "parent_acct", "tinytext not null default ''" },
		{ "lft", "int not null" },
		{ "rgt", "int not null" },
		{ "fairshare", "int default 1 not null" },
		{ "max_jobs", "int default NULL" },
		{ "max_submit_jobs", "int default NULL" },
		{ "max_cpus_per_job", "int default NULL" },
		{ "max_nodes_per_job", "int default NULL" },
		{ "max_wall_duration_per_job", "int default NULL" },
		{ "max_cpu_mins_per_job", "bigint default NULL" },
		{ "grp_jobs", "int default NULL" },
		{ "grp_submit_jobs", "int default NULL" },
		{ "grp_cpus", "int default NULL" },
		{ "grp_nodes", "int default NULL" },
		{ "grp_wall", "int default NULL" },
		{ "grp_cpu_hours", "bigint default NULL" },
		{ "qos", "blob not null default ''" },
		{ NULL, NULL}		
	};

	storage_field_t assoc_usage_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "id", "int not null" },
		{ "period_start", "int unsigned not null" },
		{ "alloc_cpu_secs", "bigint default 0" },
		{ NULL, NULL}		
	};

	storage_field_t cluster_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "name", "tinytext not null" },
		{ "control_host", "tinytext not null default ''" },
		{ "control_port", "mediumint not null default 0" },
		{ "rpc_version", "mediumint not null default 0" },
		{ NULL, NULL}		
	};

	storage_field_t cluster_usage_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "cluster", "tinytext not null" },
		{ "period_start", "int unsigned not null" },
		{ "cpu_count", "int default 0" },
		{ "alloc_cpu_secs", "bigint default 0" },
		{ "down_cpu_secs", "bigint default 0" },
		{ "idle_cpu_secs", "bigint default 0" },
		{ "resv_cpu_secs", "bigint default 0" },
		{ "over_cpu_secs", "bigint default 0" },
		{ NULL, NULL}		
	};

	storage_field_t event_table_fields[] = {
		{ "node_name", "tinytext default '' not null" },
		{ "cluster", "tinytext not null" },
		{ "cpu_count", "int not null" },
		{ "period_start", "int unsigned not null" },
		{ "period_end", "int unsigned default 0 not null" },
		{ "reason", "tinytext not null" },
		{ NULL, NULL}		
	};

	storage_field_t job_table_fields[] = {
		{ "id", "int not null auto_increment" },
		{ "jobid", "mediumint unsigned not null" },
		{ "associd", "mediumint unsigned not null" },
		{ "uid", "smallint unsigned not null" },
		{ "gid", "smallint unsigned not null" },
		{ "partition", "tinytext not null" },
		{ "blockid", "tinytext" },
		{ "account", "tinytext" },
		{ "eligible", "int unsigned default 0 not null" },
		{ "submit", "int unsigned default 0 not null" },
		{ "start", "int unsigned default 0 not null" },
		{ "end", "int unsigned default 0 not null" },
		{ "suspended", "int unsigned default 0 not null" },
		{ "name", "tinytext not null" }, 
		{ "track_steps", "tinyint not null" },
		{ "state", "smallint not null" }, 
		{ "comp_code", "int default 0 not null" },
		{ "priority", "int unsigned not null" },
		{ "req_cpus", "mediumint unsigned not null" }, 
		{ "alloc_cpus", "mediumint unsigned not null" }, 
		{ "nodelist", "text" },
		{ "kill_requid", "smallint default -1 not null" },
		{ "qos", "smallint default 0" },
		{ NULL, NULL}
	};

	storage_field_t last_ran_table_fields[] = {
		{ "hourly_rollup", "int unsigned default 0 not null" },
		{ "daily_rollup", "int unsigned default 0 not null" },
		{ "monthly_rollup", "int unsigned default 0 not null" },
		{ NULL, NULL}		
	};

	storage_field_t qos_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "id", "int not null auto_increment" },
		{ "name", "tinytext not null" }, 
		{ "description", "text" }, 
		{ NULL, NULL}		
	};

	storage_field_t step_table_fields[] = {
		{ "id", "int not null" },
		{ "stepid", "smallint not null" },
		{ "start", "int unsigned default 0 not null" },
		{ "end", "int unsigned default 0 not null" },
		{ "suspended", "int unsigned default 0 not null" },
		{ "name", "text not null" },
		{ "nodelist", "text not null" },
		{ "state", "smallint not null" },
		{ "kill_requid", "smallint default -1 not null" },
		{ "comp_code", "int default 0 not null" },
		{ "cpus", "mediumint unsigned not null" },
		{ "user_sec", "int unsigned default 0 not null" },
		{ "user_usec", "int unsigned default 0 not null" },
		{ "sys_sec", "int unsigned default 0 not null" },
		{ "sys_usec", "int unsigned default 0 not null" },
		{ "max_vsize", "int unsigned default 0 not null" },
		{ "max_vsize_task", "smallint unsigned default 0 not null" },
		{ "max_vsize_node", "mediumint unsigned default 0 not null" },
		{ "ave_vsize", "float default 0.0 not null" },
		{ "max_rss", "int unsigned default 0 not null" },
		{ "max_rss_task", "smallint unsigned default 0 not null" },
		{ "max_rss_node", "mediumint unsigned default 0 not null" },
		{ "ave_rss", "float default 0.0 not null" },
		{ "max_pages", "mediumint unsigned default 0 not null" },
		{ "max_pages_task", "smallint unsigned default 0 not null" },
		{ "max_pages_node", "mediumint unsigned default 0 not null" },
		{ "ave_pages", "float default 0.0 not null" },
		{ "min_cpu", "mediumint unsigned default 0 not null" },
		{ "min_cpu_task", "smallint unsigned default 0 not null" },
		{ "min_cpu_node", "mediumint unsigned default 0 not null" },
		{ "ave_cpu", "float default 0.0 not null" },
		{ NULL, NULL}
	};

	storage_field_t suspend_table_fields[] = {
		{ "id", "int not null" },
		{ "associd", "mediumint not null" },
		{ "start", "int unsigned default 0 not null" },
		{ "end", "int unsigned default 0 not null" },
		{ NULL, NULL}		
	};

	storage_field_t txn_table_fields[] = {
		{ "id", "int not null auto_increment" },
		{ "timestamp", "int unsigned default 0 not null" },
		{ "action", "smallint not null" },
		{ "name", "tinytext not null" },
		{ "actor", "tinytext not null" },
		{ "info", "text" },
		{ NULL, NULL}		
	};

	storage_field_t user_table_fields[] = {
		{ "creation_time", "int unsigned not null" },
		{ "mod_time", "int unsigned default 0 not null" },
		{ "deleted", "tinyint default 0" },
		{ "name", "tinytext not null" },
		{ "default_acct", "tinytext not null" },
		{ "admin_level", "smallint default 1 not null" },
		{ NULL, NULL}		
	};

	char *get_parent_proc = 
		"drop procedure if exists get_parent_limits; "
		"create procedure get_parent_limits("
		"my_table text, acct text, cluster text, without_limits int) "
		"begin "
		"set @par_id = NULL; "
		"set @mj = NULL; "
		"set @msj = NULL; "
		"set @mcpj = NULL; "
		"set @mnpj = NULL; "
		"set @mwpj = NULL; "
		"set @mcmpj = NULL; "
		"set @qos = NULL; "
		"set @my_acct = acct; "
		"if without_limits then "
		"set @mj = 0; " 
		"set @msj = 0; " 
		"set @mcpj = 0; "
		"set @mnpj = 0; "
		"set @mwpj = 0; "
		"set @mcmpj = 0; "
		"set @qos = 0; "
		"end if; "
		"REPEAT "
		"set @s = 'select '; "
		"if @par_id is NULL then set @s = CONCAT("
		"@s, '@par_id := id, '); "
		"end if; "
		"if @mj is NULL then set @s = CONCAT("
		"@s, '@mj := max_jobs, '); "
		"end if; "
		"if @msj is NULL then set @s = CONCAT("
		"@s, '@msj := max_submit_jobs, '); "
		"end if; "
		"if @mcpj is NULL then set @s = CONCAT("
		"@s, '@mcpj := max_cpus_per_job, ') ;"
		"end if; "
		"if @mnpj is NULL then set @s = CONCAT("
		"@s, '@mnpj := max_nodes_per_job, ') ;"
		"end if; "
		"if @mwpj is NULL then set @s = CONCAT("
		"@s, '@mwpj := max_wall_duration_per_job, '); "
		"end if; "
		"if @mcmpj is NULL then set @s = CONCAT("
		"@s, '@mcmpj := max_cpu_mins_per_job, '); "
		"end if; "
		"if @qos is NULL then set @s = CONCAT("
		"@s, '@qos := qos, '); "
		"end if; "
		"set @s = concat(@s, ' @my_acct := parent_acct from ', "
		"my_table, ' where acct = \"', @my_acct, '\" && "
		"cluster = \"', cluster, '\" && user=\"\"'); "
		"prepare query from @s; "
		"execute query; "
		"deallocate prepare query; "
		"UNTIL (@mj != -1 && @msj != -1 && @mcpj != -1 "
		"&& @mnpj != -1 && @mwpj != -1 "
		"&& @mcmpj != -1 && @qos != '') || @my_acct = '' END REPEAT; "
		"END;";
	char *query = NULL;
	time_t now = time(NULL);

	if(mysql_db_create_table(db_conn, acct_coord_table,
				 acct_coord_table_fields,
				 ", primary key (acct(20), user(20)))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, acct_table, acct_table_fields,
				 ", primary key (name(20)))") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, assoc_day_table,
				 assoc_usage_table_fields,
				 ", primary key (id, period_start))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, assoc_hour_table,
				 assoc_usage_table_fields,
				 ", primary key (id, period_start))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, assoc_month_table,
				 assoc_usage_table_fields,
				 ", primary key (id, period_start))") 
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, assoc_table, assoc_table_fields,
				 ", primary key (id), "
				 " unique index (user(20), acct(20), "
				 "cluster(20), partition(20)))"
/* 				 " unique index (lft), " */
				 /* 				 " unique index (rgt))" */)
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, cluster_day_table,
				 cluster_usage_table_fields,
				 ", primary key (cluster(20), period_start))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, cluster_hour_table,
				 cluster_usage_table_fields,
				 ", primary key (cluster(20), period_start))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, cluster_month_table,
				 cluster_usage_table_fields,
				 ", primary key (cluster(20), period_start))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, cluster_table,
				 cluster_table_fields,
				 ", primary key (name(20)))") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, event_table,
				 event_table_fields,
				 ", primary key (node_name(20), cluster(20), "
				 "period_start))") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, job_table, job_table_fields,
				 ", primary key (id), "
				 "unique index (jobid, associd, submit))")
	   == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, last_ran_table,
				 last_ran_table_fields, 
				 ")") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, qos_table,
				 qos_table_fields, 
				 ", primary key (id), "
				 "unique index (name(20)))")
	   == SLURM_ERROR)
		return SLURM_ERROR;
	else {
		query = xstrdup_printf(
			"insert into %s "
			"(creation_time, mod_time, name, description) "
			"values (%d, %d, 'normal', 'Normal QOS default') "
			"on duplicate key update id=LAST_INSERT_ID(id), "
			"deleted=0;",
			qos_table, now, now);
		//debug3("%s", query);
		normal_qos_id = mysql_insert_ret_id(db_conn, query);
		xfree(query);		
	}
	if(mysql_db_create_table(db_conn, step_table,
				 step_table_fields, 
				 ", primary key (id, stepid))") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, suspend_table,
				 suspend_table_fields, 
				 ")") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, txn_table, txn_table_fields,
				 ", primary key (id))") == SLURM_ERROR)
		return SLURM_ERROR;

	if(mysql_db_create_table(db_conn, user_table, user_table_fields,
				 ", primary key (name(20)))") == SLURM_ERROR)
		return SLURM_ERROR;

	rc = mysql_db_query(db_conn, get_parent_proc);

	/* Add user root to be a user by default and have this default
	 * account be root.  If already there just update
	 * name='root'.  That way if the admins delete it it will
	 * remained deleted. Creation time will be 0 so it will never
	 * really be deleted.
	 */
	query = xstrdup_printf(
		"insert into %s (creation_time, mod_time, name, default_acct, "
		"admin_level) values (0, %d, 'root', 'root', %u) "
		"on duplicate key update name='root';",
		user_table, now, ACCT_ADMIN_SUPER_USER, now);
	xstrfmtcat(query, 
		   "insert into %s (creation_time, mod_time, name, "
		   "description, organization) values (0, %d, 'root', "
		   "'default root account', 'root') on duplicate key "
		   "update name='root';",
		   acct_table, now); 

	//debug3("%s", query);
	mysql_db_query(db_conn, query);
	xfree(query);		

	return rc;
}
#endif

/*
 * init() is called when the plugin is loaded, before any other functions
 * are called.  Put global initialization here.
 */
extern int init ( void )
{
	static int first = 1;
	int rc = SLURM_SUCCESS;
#ifdef HAVE_MYSQL
	MYSQL *db_conn = NULL;
	char *location = NULL;
#else
	fatal("No MySQL database was found on the machine. "
	      "Please check the config.log from the run of configure "
	      "and run again.");
#endif

	/* since this can be loaded from many different places
	   only tell us once. */
	if(!first)
		return SLURM_SUCCESS;

	first = 0;

#ifdef HAVE_MYSQL
	mysql_db_info = _mysql_acct_create_db_info();

	location = slurm_get_accounting_storage_loc();
	if(!location)
		mysql_db_name = xstrdup(DEFAULT_ACCT_DB);
	else {
		int i = 0;
		while(location[i]) {
			if(location[i] == '.' || location[i] == '/') {
				debug("%s doesn't look like a database "
				      "name using %s",
				      location, DEFAULT_ACCT_DB);
				break;
			}
			i++;
		}
		if(location[i]) {
			mysql_db_name = xstrdup(DEFAULT_ACCT_DB);
			xfree(location);
		} else
			mysql_db_name = location;
	}

	debug2("mysql_connect() called for db %s", mysql_db_name);
	
	mysql_get_db_connection(&db_conn, mysql_db_name, mysql_db_info);
		
	rc = _mysql_acct_check_tables(db_conn);

	mysql_close_db_connection(&db_conn);
	
#endif		

	if(rc == SLURM_SUCCESS)
		verbose("%s loaded", plugin_name);
	else 
		verbose("%s failed", plugin_name);
	
	return rc;
}

extern int fini ( void )
{
#ifdef HAVE_MYSQL
	destroy_mysql_db_info(mysql_db_info);		
	xfree(mysql_db_name);
	mysql_cleanup();
	return SLURM_SUCCESS;
#else
	return SLURM_ERROR;
#endif
}

extern void *acct_storage_p_get_connection(bool make_agent, bool rollback)
{
#ifdef HAVE_MYSQL
	mysql_conn_t *mysql_conn = xmalloc(sizeof(mysql_conn_t));
	static int conn = 0;
	if(!mysql_db_info)
		init();

	debug2("acct_storage_p_get_connection: request new connection");
	
	mysql_get_db_connection(&mysql_conn->db_conn,
				mysql_db_name, mysql_db_info);
	mysql_conn->rollback = rollback;
	if(rollback) {
		mysql_autocommit(mysql_conn->db_conn, 0);
	}
	mysql_conn->conn = conn++;
	mysql_conn->update_list = list_create(destroy_acct_update_object);
	return (void *)mysql_conn;
#else
	return NULL;
#endif
}

extern int acct_storage_p_close_connection(mysql_conn_t **mysql_conn)
{
#ifdef HAVE_MYSQL

	if(!mysql_conn || !(*mysql_conn))
		return SLURM_SUCCESS;

	acct_storage_p_commit((*mysql_conn), 0);
	mysql_close_db_connection(&(*mysql_conn)->db_conn);
	list_destroy((*mysql_conn)->update_list);
	xfree((*mysql_conn));

	return SLURM_SUCCESS;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_commit(mysql_conn_t *mysql_conn, bool commit)
{
#ifdef HAVE_MYSQL
	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	debug4("got %d commits", list_count(mysql_conn->update_list));

	if(mysql_conn->rollback) {
		if(!commit) {
			if(mysql_db_rollback(mysql_conn->db_conn))
				error("rollback failed");
		} else {
			if(mysql_db_commit(mysql_conn->db_conn))
				error("commit failed");
		}
	}
	
	if(commit && list_count(mysql_conn->update_list)) {
		int rc;
		char *query = NULL;
		MYSQL_RES *result = NULL;
		MYSQL_ROW row;
		accounting_update_msg_t msg;
		slurm_msg_t req;
		slurm_msg_t resp;
		ListIterator itr = NULL;
		acct_update_object_t *object = NULL;
		
		memset(&msg, 0, sizeof(accounting_update_msg_t));
		msg.update_list = mysql_conn->update_list;
		
		xstrfmtcat(query, "select control_host, control_port, "
			   "name, rpc_version "
			   "from %s where deleted=0 && control_port != 0",
			   cluster_table);
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			goto skip;
		}
		xfree(query);
		while((row = mysql_fetch_row(result))) {
			info("sending to %s at %s(%s) ver %s",
			     row[2], row[0], row[1], row[3]);
			msg.rpc_version = atoi(row[3]);
			slurm_msg_t_init(&req);
			slurm_set_addr_char(&req.address, atoi(row[1]), row[0]);
			req.msg_type = ACCOUNTING_UPDATE_MSG;
			req.flags = SLURM_GLOBAL_AUTH_KEY;
			req.data = &msg;			
			slurm_msg_t_init(&resp);
			
			rc = slurm_send_recv_node_msg(&req, &resp, 0);
			if ((rc != 0) || !resp.auth_cred) {
				error("update cluster: %m to %s at %s(%s)",
				      row[2], row[0], row[1]);
				if (resp.auth_cred)
					g_slurm_auth_destroy(resp.auth_cred);
				rc = SLURM_ERROR;
			}
			if (resp.auth_cred)
				g_slurm_auth_destroy(resp.auth_cred);
			
			switch (resp.msg_type) {
			case RESPONSE_SLURM_RC:
				rc = ((return_code_msg_t *)resp.data)->
					return_code;
				slurm_free_return_code_msg(resp.data);	
				break;
			default:
				break;
			}	
			//info("got rc of %d", rc);
		}
		mysql_free_result(result);
	skip:
		/* NOTE: you can not use list_pop, or list_push
		   anywhere either, since mysql is
		   exporting something of the same type as a macro,
		   which messes everything up (my_list.h is the bad boy).
		   So we are just going to delete each item as it
		   comes out.
		*/
		itr = list_iterator_create(mysql_conn->update_list);
		while((object = list_next(itr))) {
			if(!object->objects || !list_count(object->objects)) {
				list_delete_item(itr);
				continue;
			}
			switch(object->type) {
			case ACCT_MODIFY_USER:
			case ACCT_ADD_USER:
			case ACCT_REMOVE_USER:
			case ACCT_ADD_COORD:
			case ACCT_REMOVE_COORD:
				rc = assoc_mgr_update_local_users(object);
				break;
			case ACCT_ADD_ASSOC:
			case ACCT_MODIFY_ASSOC:
			case ACCT_REMOVE_ASSOC:
				rc = assoc_mgr_update_local_assocs(object);
				break;
			case ACCT_ADD_QOS:
			case ACCT_REMOVE_QOS:
				rc = assoc_mgr_update_local_qos(object);
				break;
			case ACCT_UPDATE_NOTSET:
			default:
				error("unknown type set in "
				      "update_object: %d",
				      object->type);
				break;
			}
			list_delete_item(itr);
		}
		list_iterator_destroy(itr);
	}
	list_flush(mysql_conn->update_list);

	return SLURM_SUCCESS;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_add_users(mysql_conn_t *mysql_conn, uint32_t uid,
				    List user_list)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	int rc = SLURM_SUCCESS;
	acct_user_rec_t *object = NULL;
	char *cols = NULL, *vals = NULL, *query = NULL, *txn_query = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	char *extra = NULL;
	int affect_rows = 0;
	List assoc_list = list_create(destroy_acct_association_rec);

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	user_name = uid_to_string((uid_t) uid);
	itr = list_iterator_create(user_list);
	while((object = list_next(itr))) {
		if(!object->name || !object->default_acct) {
			error("We need a user name and "
			      "default acct to add.");
			rc = SLURM_ERROR;
			continue;
		}
		xstrcat(cols, "creation_time, mod_time, name, default_acct");
		xstrfmtcat(vals, "%d, %d, '%s', '%s'", 
			   now, now, object->name, object->default_acct); 
		xstrfmtcat(extra, ", default_acct='%s'", object->default_acct);
		
		if(object->admin_level != ACCT_ADMIN_NOTSET) {
			xstrcat(cols, ", admin_level");
			xstrfmtcat(vals, ", %u", object->admin_level);
			xstrfmtcat(extra, ", admin_level=%u", 
				   object->admin_level); 		
		}

		query = xstrdup_printf(
			"insert into %s (%s) values (%s) "
			"on duplicate key update deleted=0, mod_time=%d %s;",
			user_table, cols, vals,
			now, extra);

		xfree(cols);
		xfree(vals);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add user %s", object->name);
			xfree(extra);
			continue;
		}

		affect_rows = _last_affected_rows(mysql_conn->db_conn);
		if(!affect_rows) {
			debug("nothing changed");
			xfree(extra);
			continue;
		}

		if(_addto_update_list(mysql_conn->update_list, ACCT_ADD_USER,
				      object) == SLURM_SUCCESS) 
			list_remove(itr);
			

		if(txn_query)
			xstrfmtcat(txn_query, 	
				   ", (%d, %u, '%s', '%s', \"%s\")",
				   now, DBD_ADD_USERS, object->name,
				   user_name, extra);
		else
			xstrfmtcat(txn_query, 	
				   "insert into %s "
				   "(timestamp, action, name, actor, info) "
				   "values (%d, %u, '%s', '%s', \"%s\")",
				   txn_table,
				   now, DBD_ADD_USERS, object->name,
				   user_name, extra);
		xfree(extra);
		
		if(!object->assoc_list)
			continue;

		list_transfer(assoc_list, object->assoc_list);
	}
	list_iterator_destroy(itr);
	xfree(user_name);

	if(rc != SLURM_ERROR) {
		if(txn_query) {
			xstrcat(txn_query, ";");
			rc = mysql_db_query(mysql_conn->db_conn,
					    txn_query);
			xfree(txn_query);
			if(rc != SLURM_SUCCESS) {
				error("Couldn't add txn");
				rc = SLURM_SUCCESS;
			}
		}
	} else
		xfree(txn_query);

	if(list_count(assoc_list)) {
		if(acct_storage_p_add_associations(mysql_conn, uid, assoc_list)
		   == SLURM_ERROR) {
			error("Problem adding user associations");
			rc = SLURM_ERROR;
		}
	}
	list_destroy(assoc_list);

	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_add_coord(mysql_conn_t *mysql_conn, uint32_t uid, 
				    List acct_list, acct_user_cond_t *user_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL, *user = NULL, *acct = NULL;
	char *user_name = NULL, *txn_query = NULL;
	ListIterator itr, itr2;
	time_t now = time(NULL);
	int rc = SLURM_SUCCESS;
	acct_user_rec_t *user_rec = NULL;
	
	if(!user_cond || !user_cond->assoc_cond 
	   || !user_cond->assoc_cond->user_list 
	   || !list_count(user_cond->assoc_cond->user_list) 
	   || !acct_list || !list_count(acct_list)) {
		error("we need something to add");
		return SLURM_ERROR;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	user_name = uid_to_string((uid_t) uid);
	itr = list_iterator_create(user_cond->assoc_cond->user_list);
	itr2 = list_iterator_create(acct_list);
	while((user = list_next(itr))) {
		while((acct = list_next(itr2))) {
			if(query) 
				xstrfmtcat(query, ", (%d, %d, '%s', '%s')",
					   now, now, acct, user);
			else
				query = xstrdup_printf(
					"insert into %s (creation_time, "
					"mod_time, acct, user) values "
					"(%d, %d, '%s', '%s')",
					acct_coord_table, 
					now, now, acct, user); 

			if(txn_query)
				xstrfmtcat(txn_query, 	
					   ", (%d, %u, '%s', '%s', '%s')",
					   now, DBD_ADD_ACCOUNT_COORDS, user,
					   user_name, acct);
			else
				xstrfmtcat(txn_query, 	
					   "insert into %s "
					   "(timestamp, action, name, "
					   "actor, info) "
					   "values (%d, %u, '%s', "
					   "'%s', \"%s\")",
					   txn_table,
					   now, DBD_ADD_ACCOUNT_COORDS, user,
					   user_name, acct);
		}
		list_iterator_reset(itr2);
	}
	xfree(user_name);
	list_iterator_destroy(itr);
	list_iterator_destroy(itr2);

	if(query) {
		xstrfmtcat(query, 
			   " on duplicate key update mod_time=%d, deleted=0;%s",
			   now, txn_query);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		xfree(txn_query);
		
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add cluster hour rollup");
			return rc;
		}
		/* get the update list set */
		itr = list_iterator_create(user_cond->assoc_cond->user_list);
		while((user = list_next(itr))) {
			user_rec = xmalloc(sizeof(acct_user_rec_t));
			user_rec->name = xstrdup(user);
			_get_user_coords(mysql_conn, user_rec);
			_addto_update_list(mysql_conn->update_list, 
					   ACCT_ADD_COORD, user_rec);
		}
		list_iterator_destroy(itr);
	}
	
	return SLURM_SUCCESS;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_add_accts(mysql_conn_t *mysql_conn, uint32_t uid, 
				    List acct_list)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	int rc = SLURM_SUCCESS;
	acct_account_rec_t *object = NULL;
	char *cols = NULL, *vals = NULL, *query = NULL, *txn_query = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	char *extra = NULL;
	int affect_rows = 0;
	List assoc_list = list_create(destroy_acct_association_rec);

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	user_name = uid_to_string((uid_t) uid);
	itr = list_iterator_create(acct_list);
	while((object = list_next(itr))) {
		if(!object->name || !object->description
		   || !object->organization) {
			error("We need an account name, description, and "
			      "organization to add. %s %s %s", 
			      object->name, object->description,
			      object->organization);
			rc = SLURM_ERROR;
			continue;
		}
		xstrcat(cols, "creation_time, mod_time, name, "
			"description, organization");
		xstrfmtcat(vals, "%d, %d, '%s', '%s', '%s'", 
			   now, now, object->name, 
			   object->description, object->organization); 
		xstrfmtcat(extra, ", description='%s', organization='%s'",
			   object->description, object->organization); 		
		
		query = xstrdup_printf(
			"insert into %s (%s) values (%s) "
			"on duplicate key update deleted=0, mod_time=%d %s;",
			acct_table, cols, vals,
			now, extra);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(cols);
		xfree(vals);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add acct");
			xfree(extra);
			continue;
		}
		affect_rows = _last_affected_rows(mysql_conn->db_conn);
/* 		debug3("affected %d", affect_rows); */

		if(!affect_rows) {
			debug3("nothing changed");
			xfree(extra);
			continue;
		}

		if(txn_query)
			xstrfmtcat(txn_query, 	
				   ", (%d, %u, '%s', '%s', \"%s\")",
				   now, DBD_ADD_ACCOUNTS, object->name,
				   user_name, extra);
		else
			xstrfmtcat(txn_query, 	
				   "insert into %s "
				   "(timestamp, action, name, actor, info) "
				   "values (%d, %u, '%s', '%s', \"%s\")",
				   txn_table,
				   now, DBD_ADD_ACCOUNTS, object->name,
				   user_name, extra);
		xfree(extra);
		
		if(!object->assoc_list)
			continue;

		list_transfer(assoc_list, object->assoc_list);
	}
	list_iterator_destroy(itr);
	xfree(user_name);
	
	if(rc != SLURM_ERROR) {
		if(txn_query) {
			xstrcat(txn_query, ";");
			rc = mysql_db_query(mysql_conn->db_conn,
					    txn_query);
			xfree(txn_query);
			if(rc != SLURM_SUCCESS) {
				error("Couldn't add txn");
				rc = SLURM_SUCCESS;
			}
		}
	} else
		xfree(txn_query);

	if(list_count(assoc_list)) {
		if(acct_storage_p_add_associations(mysql_conn, uid, assoc_list)
		   == SLURM_ERROR) {
			error("Problem adding user associations");
			rc = SLURM_ERROR;
		}
	}
	list_destroy(assoc_list);

	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_add_clusters(mysql_conn_t *mysql_conn, uint32_t uid, 
				       List cluster_list)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	int rc = SLURM_SUCCESS;
	acct_cluster_rec_t *object = NULL;
	char *cols = NULL, *vals = NULL, *extra = NULL, *query = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int affect_rows = 0;
	int added = 0;
	List assoc_list = NULL;
	acct_association_rec_t *assoc = NULL;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	assoc_list = list_create(destroy_acct_association_rec);

	user_name = uid_to_string((uid_t) uid);
	itr = list_iterator_create(cluster_list);
	while((object = list_next(itr))) {
		if(!object->name) {
			error("We need a cluster name to add.");
			rc = SLURM_ERROR;
			continue;
		}

		xstrcat(cols, "creation_time, mod_time, acct, cluster");
		xstrfmtcat(vals, "%d, %d, 'root', '%s'",
			   now, now, object->name);
		xstrfmtcat(extra, ", mod_time=%d", now);
		if(object->root_assoc)
			_setup_association_limits(object->root_assoc, &cols, 
						  &vals, &extra, 1);
		xstrfmtcat(query, 
			   "insert into %s (creation_time, mod_time, name) "
			   "values (%d, %d, '%s') "
			   "on duplicate key update deleted=0, mod_time=%d, "
			   "control_host='', control_port=0;",
			   cluster_table, 
			   now, now, object->name,
			   now);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add cluster %s", object->name);
			xfree(extra);
			xfree(cols);
			xfree(vals);
			added=0;
			break;
		}

		affect_rows = _last_affected_rows(mysql_conn->db_conn);

		if(!affect_rows) {
			debug2("nothing changed %d", affect_rows);
			xfree(extra);
			xfree(cols);
			xfree(vals);
			continue;
		}

		xstrfmtcat(query,
			   "SELECT @MyMax := coalesce(max(rgt), 0) FROM %s "
			   "FOR UPDATE;",
			   assoc_table);
		xstrfmtcat(query,
			   "insert into %s (%s, lft, rgt) "
			   "values (%s, @MyMax+1, @MyMax+2) "
			   "on duplicate key update deleted=0, "
			   "id=LAST_INSERT_ID(id)%s;",
			   assoc_table, cols,
			   vals,
			   extra);
		
		xfree(cols);
		xfree(vals);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);

		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);

		if(rc != SLURM_SUCCESS) {
			error("Couldn't add cluster root assoc");
			xfree(extra);
			added=0;
			break;
		}
		xstrfmtcat(query,
			   "insert into %s "
			   "(timestamp, action, name, actor, info) "
			   "values (%d, %u, '%s', '%s', \"%s\");",
			   txn_table, now, DBD_ADD_CLUSTERS, 
			   object->name, user_name, extra);
		xfree(extra);			
		debug4("query\n%s",query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add txn");
		} else
			added++;

		/* Add user root by default to run from the root
		 * association.  This gets popped off so we need to
		 * readd it every time here. 
		 */
		assoc = xmalloc(sizeof(acct_association_rec_t));
		list_append(assoc_list, assoc);
		
		assoc->cluster = xstrdup(object->name);
		assoc->user = xstrdup("root");
		assoc->acct = xstrdup("root");
		assoc->fairshare = NO_VAL;
		assoc->max_cpu_mins_pj = NO_VAL;
		assoc->max_jobs = NO_VAL;
		assoc->max_nodes_pj = NO_VAL;
		assoc->max_wall_pj = NO_VAL;

		if(acct_storage_p_add_associations(mysql_conn, uid, assoc_list)
		   == SLURM_ERROR) {
			error("Problem adding root user association");
			rc = SLURM_ERROR;
		}

	}
	list_iterator_destroy(itr);
	xfree(user_name);

	list_destroy(assoc_list);

	if(!added) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
	}

	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_add_associations(mysql_conn_t *mysql_conn,
					   uint32_t uid, 
					   List association_list)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	int rc = SLURM_SUCCESS;
	int i=0;
	acct_association_rec_t *object = NULL;
	char *cols = NULL, *vals = NULL, *txn_query = NULL,
		*extra = NULL, *query = NULL, *update = NULL;
	char *parent = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	char *tmp_char = NULL;
	int assoc_id = 0;
	int incr = 0, my_left = 0;
	int affect_rows = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *old_parent = NULL, *old_cluster = NULL;
	char *massoc_req_inx[] = {
		"id",
		"parent_acct",
		"lft",
		"rgt",
		"deleted"
	};
	
	enum {
		MASSOC_ID,
		MASSOC_PACCT,
		MASSOC_LFT,
		MASSOC_RGT,
		MASSOC_DELETED,
		MASSOC_COUNT
	};

	if(!association_list) {
		error("No association list given");
		return SLURM_ERROR;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	user_name = uid_to_string((uid_t) uid);
	itr = list_iterator_create(association_list);
	while((object = list_next(itr))) {
		if(!object->cluster || !object->acct) {
			error("We need a association cluster and "
			      "acct to add one.");
			rc = SLURM_ERROR;
			continue;
		}

		if(object->parent_acct) {
			parent = object->parent_acct;
		} else if(object->user) {
			parent = object->acct;
		} else {
			parent = "root";
		}

		xstrcat(cols, "creation_time, mod_time, cluster, acct");
		xstrfmtcat(vals, "%d, %d, '%s', '%s'", 
			   now, now, object->cluster, object->acct); 
		xstrfmtcat(update, "where id>=0 && cluster='%s' && acct='%s'",
			   object->cluster, object->acct); 

		xstrfmtcat(extra, ", mod_time=%d", now);
		if(!object->user) {
			xstrcat(cols, ", parent_acct");
			xstrfmtcat(vals, ", '%s'", parent);
			xstrfmtcat(extra, ", parent_acct='%s'", parent);
			xstrfmtcat(update, " && user=''"); 
		} else {
			char *part = object->partition;
			xstrcat(cols, ", user");
			xstrfmtcat(vals, ", '%s'", object->user); 		
			xstrfmtcat(update, " && user='%s'",
				   object->user); 

			/* We need to give a partition wiether it be
			 * '' or the actual partition name given
			 */
			if(!part)
				part = "";
			xstrcat(cols, ", partition");
			xstrfmtcat(vals, ", '%s'", part);
			xstrfmtcat(update, " && partition='%s'", part);
		}

		_setup_association_limits(object, &cols, &vals, &extra, 1);

		for(i=0; i<MASSOC_COUNT; i++) {
			if(i) 
				xstrcat(tmp_char, ", ");
			xstrcat(tmp_char, massoc_req_inx[i]);
		}
		
		xstrfmtcat(query, 
			   "select distinct %s from %s %s order by lft "
			   "FOR UPDATE;",
			   tmp_char, assoc_table, update);
		xfree(tmp_char);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			xfree(cols);
			xfree(vals);
			xfree(extra);
			xfree(update);
			error("couldn't query the database");
			rc = SLURM_ERROR;
			break;
		}
		xfree(query);

		assoc_id = 0;
		if(!(row = mysql_fetch_row(result))) {
			/* This code speeds up the add process quite a bit
			 * here we are only doing an update when we are done
			 * adding to a specific group (cluster/account) other
			 * than that we are adding right behind what we were
			 * so just total them up and then do one update
			 * instead of the slow ones that require an update
			 * every time.  There is a incr check outside of the
			 * loop to catch everything on the last spin of the
			 * while. 
			 */ 
			if(!old_parent || !old_cluster
			   || strcasecmp(parent, old_parent) 
			   || strcasecmp(object->cluster, old_cluster)) {
				char *sel_query = xstrdup_printf(
					"SELECT lft FROM %s WHERE "
					"acct = '%s' and cluster = '%s' "
					"and user = '' order by lft;",
					assoc_table,
					parent, object->cluster);
				MYSQL_RES *sel_result = NULL;
				
				if(incr) {
					char *up_query = xstrdup_printf(
						"UPDATE %s SET rgt = rgt+%d "
						"WHERE rgt > %d && deleted < 2;"
						"UPDATE %s SET lft = lft+%d "
						"WHERE lft > %d "
						"&& deleted < 2;"
						"UPDATE %s SET deleted = 0 "
						"WHERE deleted = 2;",
						assoc_table, incr,
						my_left,
						assoc_table, incr,
						my_left,
						assoc_table);
					debug3("%d query\n%s", mysql_conn->conn,
					       up_query);
					rc = mysql_db_query(
						mysql_conn->db_conn,
						up_query);
					xfree(up_query);
					if(rc != SLURM_SUCCESS) {
						error("Couldn't do update");
						xfree(cols);
						xfree(vals);
						xfree(update);
						xfree(extra);
						xfree(sel_query);
						break;
					}
				}

				debug3("%d query\n%s", mysql_conn->conn,
				       sel_query);
				if(!(sel_result = mysql_db_query_ret(
					     mysql_conn->db_conn,
					     sel_query, 0))) {
					xfree(cols);
					xfree(vals);
					xfree(update);
					xfree(extra);
					xfree(sel_query);
					rc = SLURM_ERROR;
					break;
				}
				
				if(!(row = mysql_fetch_row(sel_result))) {
					error("Couldn't get left from query\n",
					      sel_query);
					mysql_free_result(sel_result);
					xfree(cols);
					xfree(vals);
					xfree(update);
					xfree(extra);
					xfree(sel_query);
					rc = SLURM_ERROR;
					break;
				}
				xfree(sel_query);

				my_left = atoi(row[0]);
				mysql_free_result(sel_result);
				//info("left is %d", my_left);
				xfree(old_parent);
				xfree(old_cluster);
				old_parent = xstrdup(parent);
				old_cluster = xstrdup(object->cluster);
				incr = 0;
			}
			incr += 2;
			xstrfmtcat(query,
				   "insert into %s (%s, lft, rgt, deleted) "
				   "values (%s, %d, %d, 2);",
				   assoc_table, cols,
				   vals, my_left+(incr-1), my_left+incr);
			
			/* definantly works but slow */
/* 			xstrfmtcat(query, */
/* 				   "SELECT @myLeft := lft FROM %s WHERE " */
/* 				   "acct = '%s' " */
/* 				   "and cluster = '%s' and user = '';", */
/* 				   assoc_table, */
/* 				   parent, */
/* 				   object->cluster); */
/* 			xstrfmtcat(query, */
/* 				   "UPDATE %s SET rgt = rgt+2 " */
/* 				   "WHERE rgt > @myLeft;" */
/* 				   "UPDATE %s SET lft = lft+2 " */
/* 				   "WHERE lft > @myLeft;", */
/* 				   assoc_table, */
/* 				   assoc_table); */
/* 			xstrfmtcat(query, */
/* 				   "insert into %s (%s, lft, rgt) " */
/* 				   "values (%s, @myLeft+1, @myLeft+2);", */
/* 				   assoc_table, cols, */
/* 				   vals); */
		} else if(!atoi(row[MASSOC_DELETED])) {
			/* We don't need to do anything here */
			debug("This account was added already");
			xfree(cols);
			xfree(vals);
			xfree(update);
			mysql_free_result(result);
			xfree(extra);
			continue;
		} else {
			/* If it was once deleted we have kept the lft
			 * and rgt's consant while it was deleted and
			 * so we can just unset the deleted flag,
			 * check for the parent and move if needed.
			 */
			assoc_id = atoi(row[MASSOC_ID]);
			if(object->parent_acct 
			   && strcasecmp(object->parent_acct,
					 row[MASSOC_PACCT])) {
				
				/* We need to move the parent! */
				if(_move_parent(mysql_conn, uid,
						atoi(row[MASSOC_LFT]),
						atoi(row[MASSOC_RGT]),
						object->cluster,
						row[MASSOC_ID],
						row[MASSOC_PACCT],
						object->parent_acct)
				   == SLURM_ERROR)
					continue;
			}


			affect_rows = 2;
			xstrfmtcat(query,
				   "update %s set deleted=0, "
				   "id=LAST_INSERT_ID(id)%s %s;",
				   assoc_table, 
				   extra, update);
		}
		mysql_free_result(result);

		xfree(cols);
		xfree(vals);
		xfree(update);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add assoc");
			xfree(extra);
			break;
		}
		/* see if this was an insert or update.  On an update
		 * the assoc_id will already be set
		 */
		if(!assoc_id) {
			affect_rows = _last_affected_rows(
				mysql_conn->db_conn);
			assoc_id = mysql_insert_id(mysql_conn->db_conn);
			//info("last id was %d", assoc_id);
		}

		object->id = assoc_id;

		if(_addto_update_list(mysql_conn->update_list, ACCT_ADD_ASSOC,
				      object) == SLURM_SUCCESS) {
			list_remove(itr);
		}

		if(txn_query)
			xstrfmtcat(txn_query, 	
				   ", (%d, %d, '%d', '%s', \"%s\")",
				   now, DBD_ADD_ASSOCS, assoc_id, user_name,
				   extra);
		else
			xstrfmtcat(txn_query, 	
				   "insert into %s "
				   "(timestamp, action, name, actor, info) "
				   "values (%d, %d, '%d', '%s', \"%s\")",
				   txn_table,
				   now, DBD_ADD_ASSOCS, assoc_id, user_name, 
				   extra);
		xfree(extra);
	}
	list_iterator_destroy(itr);
	xfree(user_name);

	if(rc != SLURM_SUCCESS)
		goto end_it;

	if(incr) {
		char *up_query = xstrdup_printf(
			"UPDATE %s SET rgt = rgt+%d "
			"WHERE rgt > %d && deleted < 2;"
			"UPDATE %s SET lft = lft+%d "
			"WHERE lft > %d "
			"&& deleted < 2;"
			"UPDATE %s SET deleted = 0 "
			"WHERE deleted = 2;",
			assoc_table, incr,
			my_left,
			assoc_table, incr,
			my_left,
			assoc_table);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, up_query);
		rc = mysql_db_query(mysql_conn->db_conn, up_query);
		xfree(up_query);
		if(rc != SLURM_SUCCESS)
			error("Couldn't do update 2");
		
	}

end_it:
	if(rc != SLURM_ERROR) {
		if(txn_query) {
			xstrcat(txn_query, ";");
			rc = mysql_db_query(mysql_conn->db_conn,
					    txn_query);
			xfree(txn_query);
			if(rc != SLURM_SUCCESS) {
				error("Couldn't add txn");
				rc = SLURM_SUCCESS;
			}
		}
	} else {
		xfree(txn_query);
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
	}

	xfree(old_parent);
	xfree(old_cluster);
					
	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_add_qos(mysql_conn_t *mysql_conn, uint32_t uid, 
				  List qos_list)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	int rc = SLURM_SUCCESS;
	acct_qos_rec_t *object = NULL;
	char *query = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int affect_rows = 0;
	int added = 0;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	user_name = uid_to_string((uid_t) uid);
	itr = list_iterator_create(qos_list);
	while((object = list_next(itr))) {
		if(!object->name) {
			error("We need a qos name to add.");
			rc = SLURM_ERROR;
			continue;
		}

		xstrfmtcat(query, 
			   "insert into %s (creation_time, mod_time, "
			   "name, description) "
			   "values (%d, %d, '%s', '%s') "
			   "on duplicate key update deleted=0, mod_time=%d;",
			   qos_table, 
			   now, now, object->name, object->description,
			   now);
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add qos %s", object->name);
			added=0;
			break;
		}

		affect_rows = _last_affected_rows(mysql_conn->db_conn);

		if(!affect_rows) {
			debug2("nothing changed %d", affect_rows);
			continue;
		}
		xstrfmtcat(query,
			   "insert into %s "
			   "(timestamp, action, name, actor, info) "
			   "values (%d, %u, '%s', '%s', \"%s\");",
			   txn_table,
			   now, DBD_ADD_QOS, object->name, user_name,
			   object->description);

		debug4("query\n%s",query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
		if(rc != SLURM_SUCCESS) {
			error("Couldn't add txn");
		} else {
			if(_addto_update_list(mysql_conn->update_list, 
					      ACCT_ADD_QOS,
					      object) == SLURM_SUCCESS) 
				list_remove(itr);
			added++;
		}
		
	}
	list_iterator_destroy(itr);
	xfree(user_name);

	if(!added) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
	}

	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern List acct_storage_p_modify_users(mysql_conn_t *mysql_conn, uint32_t uid, 
					acct_user_cond_t *user_cond,
					acct_user_rec_t *user)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *vals = NULL, *extra = NULL, *query = NULL, *name_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	if(!user_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");
	if(user_cond->assoc_cond && user_cond->assoc_cond->user_list
	   && list_count(user_cond->assoc_cond->user_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(user_cond->assoc_cond->user_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(user_cond->def_acct_list && list_count(user_cond->def_acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(user_cond->def_acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "default_acct='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(user_cond->admin_level != ACCT_ADMIN_NOTSET) {
		xstrfmtcat(extra, " && admin_level=%u", user_cond->admin_level);
	}

	if(user->default_acct)
		xstrfmtcat(vals, ", default_acct='%s'", user->default_acct);

	if(user->admin_level != ACCT_ADMIN_NOTSET)
		xstrfmtcat(vals, ", admin_level=%u", user->admin_level);

	if(!extra || !vals) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		error("Nothing to change");
		return NULL;
	}
	query = xstrdup_printf("select name, qos from %s %s;",
			       user_table, extra);
	xfree(extra);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}

	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		acct_user_rec_t *user_rec = NULL;
		
		object = xstrdup(row[0]);
		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(name_char, "(name='%s'", object);
			rc = 1;
		} else  {
			xstrfmtcat(name_char, " || name='%s'", object);
		}
		user_rec = xmalloc(sizeof(acct_user_rec_t));
		user_rec->name = xstrdup(object);
		user_rec->default_acct = xstrdup(user->default_acct);
		user_rec->admin_level = user->admin_level;
		_addto_update_list(mysql_conn->update_list, ACCT_MODIFY_USER,
				   user_rec);
	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(vals);
		xfree(query);
		return ret_list;
	}
	xfree(query);
	xstrcat(name_char, ")");

	user_name = uid_to_string((uid_t) uid);
	rc = _modify_common(mysql_conn, DBD_MODIFY_USERS, now,
			    user_name, user_table, name_char, vals);
	xfree(user_name);
	xfree(name_char);
	xfree(vals);
	if (rc == SLURM_ERROR) {
		error("Couldn't modify users");
		list_destroy(ret_list);
		ret_list = NULL;
	}
				
	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_modify_accounts(
	mysql_conn_t *mysql_conn, uint32_t uid, 
	acct_account_cond_t *acct_cond,
	acct_account_rec_t *acct)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *vals = NULL, *extra = NULL, *query = NULL, *name_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	if(!acct_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");
	if(acct_cond->assoc_cond 
	   && acct_cond->assoc_cond->acct_list 
	   && list_count(acct_cond->assoc_cond->acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->assoc_cond->acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(acct_cond->description_list 
	   && list_count(acct_cond->description_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->description_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "description='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(acct_cond->organization_list 
	   && list_count(acct_cond->organization_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->organization_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "organization='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(acct->description)
		xstrfmtcat(vals, ", description='%s'", acct->description);
	if(acct->organization)
		xstrfmtcat(vals, ", organization='%s'", acct->organization);

	if(!extra || !vals) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		error("Nothing to change");
		return NULL;
	}

	query = xstrdup_printf("select name from %s %s;", acct_table, extra);
	xfree(extra);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		xfree(vals);
		return NULL;
	}

	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		object = xstrdup(row[0]);
		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(name_char, "(name='%s'", object);
			rc = 1;
		} else  {
			xstrfmtcat(name_char, " || name='%s'", object);
		}

	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(query);
		xfree(vals);
		return ret_list;
	}
	xfree(query);
	xstrcat(name_char, ")");

	user_name = uid_to_string((uid_t) uid);
	rc = _modify_common(mysql_conn, DBD_MODIFY_ACCOUNTS, now,
			    user_name, acct_table, name_char, vals);
	xfree(user_name);
	if (rc == SLURM_ERROR) {
		error("Couldn't modify accounts");
		list_destroy(ret_list);
		errno = SLURM_ERROR;
		ret_list = NULL;
	}
		
	xfree(name_char);
	xfree(vals);

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_modify_clusters(mysql_conn_t *mysql_conn, 
					   uint32_t uid, 
					   acct_cluster_cond_t *cluster_cond,
					   acct_cluster_rec_t *cluster)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *vals = NULL, *extra = NULL, *query = NULL,
		*name_char = NULL, *assoc_char= NULL, *send_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	/* If you need to alter the default values of the cluster use
	 * modify_associations since this is used only for registering
	 * the controller when it loads 
	 */

	if(!cluster_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");
	if(cluster_cond->cluster_list
	   && list_count(cluster_cond->cluster_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(cluster_cond->cluster_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	set = 0;
	if(cluster->control_host) {
		xstrfmtcat(vals, ", control_host='%s'", cluster->control_host);
		set++;
	}

	if(cluster->control_port) {
		xstrfmtcat(vals, ", control_port=%u", cluster->control_port);
		set++;
	}

	if(cluster->rpc_version) {
		xstrfmtcat(vals, ", rpc_version=%u", cluster->rpc_version);
		set++;
	}

	if(!vals) {
		xfree(extra);
		errno = SLURM_NO_CHANGE_IN_DATA;
		error("Nothing to change");
		return NULL;
	} else if(set != 3) {
		xfree(vals);
		xfree(extra);
		errno = EFAULT;
		error("Need control host, port and rpc version "
		      "to register a cluster");
		return NULL;
	}


	xstrfmtcat(query, "select name, control_port from %s %s;",
		   cluster_table, extra);

	xfree(extra);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		xfree(vals);
		error("no result given for %s", extra);
		return NULL;
	}

	/* Set here is used to ask for jobs and nodes in anything
	 * other than up state, so it you reset it later make sure
	 * this is accounted for before you do
	 */
	set = 1;
	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		object = xstrdup(row[0]);

		/* check to see if this is the first time to register */
		if(row[1][0] == '0')
			set = 0;

		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(name_char, "name='%s'", object);
			rc = 1;
		} else  {
			xstrfmtcat(name_char, " || name='%s'", object);
		}
	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(vals);
		xfree(query);
		return ret_list;
	}
	xfree(query);

	if(vals) {
		send_char = xstrdup_printf("(%s)", name_char);
		user_name = uid_to_string((uid_t) uid);
		rc = _modify_common(mysql_conn, DBD_MODIFY_CLUSTERS, now,
				    user_name, cluster_table, send_char, vals);
		xfree(user_name);
		if (rc == SLURM_ERROR) {
			error("Couldn't modify cluster 1");
			list_destroy(ret_list);
			ret_list = NULL;
			goto end_it;
		}
	}

	/* Get all nodes in a down state and jobs pending or running.
	 * This is for the first time a cluster registers
	 */

	if(!set && slurmdbd_conf) {
		/* This only happens here with the slurmdbd.  If
		 * calling this plugin directly we do this in
		 * clusteracct_storage_p_cluster_procs.
		 */
		slurm_addr ctld_address;
		slurm_fd fd;

		info("First time to register cluster requesting "
		     "running jobs and system information.");

		slurm_set_addr_char(&ctld_address, cluster->control_port,
				    cluster->control_host);
		fd =  slurm_open_msg_conn(&ctld_address);
		if (fd < 0) {
			error("can not open socket back to slurmctld");
		} else {
			slurm_msg_t out_msg;
			slurm_msg_t_init(&out_msg);
			out_msg.msg_type = ACCOUNTING_FIRST_REG;
			out_msg.flags = SLURM_GLOBAL_AUTH_KEY;
			slurm_send_node_msg(fd, &out_msg);
			/* We probably need to add matching recv_msg function
			 * for an arbitray fd or should these be fire
			 * and forget?  For this, that we can probably
			 * forget about it */
			slurm_close_stream(fd);
		}
	}

end_it:
	xfree(name_char);
	xfree(assoc_char);
	xfree(vals);
	xfree(send_char);

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_modify_associations(
	mysql_conn_t *mysql_conn, uint32_t uid, 
	acct_association_cond_t *assoc_cond,
	acct_association_rec_t *assoc)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *vals = NULL, *extra = NULL, *query = NULL, *name_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0, i = 0, is_admin=0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	acct_user_rec_t user;
	int replace_qos = 0;

	char *massoc_req_inx[] = {
		"id",
		"acct",
		"parent_acct",
		"cluster",
		"user",
		"partition",
		"lft",
		"rgt"
	};
	
	enum {
		MASSOC_ID,
		MASSOC_ACCT,
		MASSOC_PACCT,
		MASSOC_CLUSTER,
		MASSOC_USER,
		MASSOC_PART,
		MASSOC_LFT,
		MASSOC_RGT,
		MASSOC_COUNT
	};

	if(!assoc_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	/* This only works when running though the slurmdbd.
	 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
	 * SLURMDBD!
	 */
	if(slurmdbd_conf) {
		/* we have to check the authentication here in the
		 * plugin since we don't know what accounts are being
		 * referenced until after the query.  Here we will
		 * set if they are an operator or greater and then
		 * check it below after the query.
		 */
		if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
		   || assoc_mgr_get_admin_level(mysql_conn, uid) 
		   >= ACCT_ADMIN_OPERATOR) 
			is_admin = 1;	
		else {
			if(assoc_mgr_fill_in_user(mysql_conn, &user, 1)
			   != SLURM_SUCCESS) {
				error("couldn't get information for this user");
				errno = SLURM_ERROR;
				return NULL;
			}
			if(!user.coord_accts || !list_count(user.coord_accts)) {
				error("This user doesn't have any "
				      "coordinator abilities");
				errno = ESLURM_ACCESS_DENIED;
				return NULL;
			}
		}
	} else {
		/* Setting this here just makes it easier down below
		 * since user will not be filled in.
		 */
		is_admin = 1;
	}

	set = _setup_association_cond_limits(assoc_cond, &extra);
	
	if((int)assoc->fairshare >= 0) 
		xstrfmtcat(vals, ", fairshare=%u", assoc->fairshare);
	else if((int)assoc->fairshare == INFINITE) {
		xstrfmtcat(vals, ", fairshare=1");
		assoc->fairshare = 1;
	}
	if((int)assoc->max_cpu_mins_pj >= 0) 
		xstrfmtcat(vals, ", max_cpu_mins_per_job=%u",
			   assoc->max_cpu_mins_pj);
	else if((int)assoc->max_cpu_mins_pj == INFINITE) {
		xstrfmtcat(vals, ", max_cpu_mins_per_job=NULL");
	}
	if((int)assoc->max_cpus_pj >= 0) 
		xstrfmtcat(vals, ", max_cpus_per_job=%u",
			   assoc->max_cpus_pj);
	else if((int)assoc->max_cpus_pj == INFINITE) {
		xstrfmtcat(vals, ", max_cpus_per_job=NULL");
	}
	if((int)assoc->max_jobs >= 0) 
		xstrfmtcat(vals, ", max_jobs=%u", assoc->max_jobs);
	else if((int)assoc->max_jobs == INFINITE) {
		xstrfmtcat(vals, ", max_jobs=NULL");
	}
	if((int)assoc->max_nodes_pj >= 0) 
		xstrfmtcat(vals, ", max_nodes_per_job=%u",
			   assoc->max_nodes_pj);
	else if((int)assoc->max_nodes_pj == INFINITE) {
		xstrfmtcat(vals, ", max_nodes_per_job=NULL");
	}
	if((int)assoc->max_submit_jobs >= 0) 
		xstrfmtcat(vals, ", max_submit_jobs=%u",
			   assoc->max_submit_jobs);
	else if((int)assoc->max_submit_jobs == INFINITE) {
		xstrfmtcat(vals, ", max_submit_jobs=NULL");
	}
	if((int)assoc->max_wall_pj >= 0) 
		xstrfmtcat(vals, ", max_wall_duration_per_job=%u",
			   assoc->max_wall_pj);
	else if((int)assoc->max_wall_pj == INFINITE) {
		xstrfmtcat(vals, ", max_wall_duration_per_job=NULL");
	}

	if(assoc->qos_list && list_count(assoc->qos_list)) {
		char *tmp_qos = NULL;
		set = 0;
		itr = list_iterator_create(assoc->qos_list);
		while((object = list_next(itr))) {
			/* when adding we need to make sure we don't
			 * already have it so we remove it and then add
			 * it.
			 */
			if(object[0] == '-') {
				xstrfmtcat(vals,
					   ", qos=replace(qos, ',%s', '')",
					   object+1);
			} else if(object[0] == '+') {
				xstrfmtcat(vals,
					   ", qos=concat_ws(',', "
					   "replace(qos, ',%s', ''), '%s')",
					   object+1, object+1);
			} else {
				xstrfmtcat(tmp_qos, ",%s", object);
			}
		}
		list_iterator_destroy(itr);
		if(tmp_qos) {
			xstrfmtcat(vals, ", qos='%s'", tmp_qos);
			xfree(tmp_qos);
			replace_qos = 1;
		}
	}


	if(!extra || (!vals && !assoc->parent_acct)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		error("Nothing to change");
		return NULL;
	}

	for(i=0; i<MASSOC_COUNT; i++) {
		if(i) 
			xstrcat(object, ", ");
		xstrcat(object, massoc_req_inx[i]);
	}

	query = xstrdup_printf("select distinct %s from %s where deleted=0%s "
			       "order by lft FOR UPDATE;",
			       object, assoc_table, extra);
	xfree(object);
	xfree(extra);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	rc = SLURM_SUCCESS;
	set = 0;
	extra = NULL;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		acct_association_rec_t *mod_assoc = NULL;
		int account_type=0;
/* 		MYSQL_RES *result2 = NULL; */
/* 		MYSQL_ROW row2; */

		if(!is_admin) {
			acct_coord_rec_t *coord = NULL;
			char *account = row[MASSOC_ACCT];

			/* Here we want to see if the person
			 * is a coord of the parent account
			 * since we don't want him to be able
			 * to alter the limits of the account
			 * he is directly coord of.  They
			 * should be able to alter the
			 * sub-accounts though. If no parent account
			 * that means we are talking about a user
			 * association so account is really the parent
			 * of the user a coord can change that all day long.
			 */
			if(row[MASSOC_PACCT][0])
				account = row[MASSOC_PACCT];

			if(!user.coord_accts) { // This should never
						// happen
				error("We are here with no coord accts.");
				if(mysql_conn->rollback) {
					mysql_db_rollback(
						mysql_conn->db_conn);
				}
				errno = ESLURM_ACCESS_DENIED;
				mysql_free_result(result);
				xfree(vals);
				list_destroy(ret_list);
				return NULL;
			}
			itr = list_iterator_create(user.coord_accts);
			while((coord = list_next(itr))) {
				if(!strcasecmp(coord->name, account))
					break;
			}
			list_iterator_destroy(itr);

			if(!coord) {
				if(row[MASSOC_PACCT][0])
					error("User %s(%d) can not modify "
					      "account (%s) because they "
					      "are not coordinators of "
					      "parent account '%s'.",
					      user.name, user.uid,
					      row[MASSOC_ACCT], 
					      row[MASSOC_PACCT]);
				else
					error("User %s(%d) does not have the "
					      "ability to modify the account "
					      "(%s).",
					      user.name, user.uid, 
					      row[MASSOC_ACCT]);
					
				if(mysql_conn->rollback) {
					mysql_db_rollback(
						mysql_conn->db_conn);
				}
				errno = ESLURM_ACCESS_DENIED;
				mysql_free_result(result);
				xfree(vals);
				list_destroy(ret_list);
				return NULL;
			}
		}

		if(row[MASSOC_PART][0]) { 
			// see if there is a partition name
			object = xstrdup_printf(
				"C = %-10s A = %-20s U = %-9s P = %s",
				row[MASSOC_CLUSTER], row[MASSOC_ACCT],
				row[MASSOC_USER], row[MASSOC_PART]);
		} else if(row[MASSOC_USER][0]){
			object = xstrdup_printf(
				"C = %-10s A = %-20s U = %-9s",
				row[MASSOC_CLUSTER], row[MASSOC_ACCT], 
				row[MASSOC_USER]);
		} else {
			if(row[MASSOC_PACCT][0]) {
				object = xstrdup_printf(
					"C = %-10s A = %s of %s",
					row[MASSOC_CLUSTER], row[MASSOC_ACCT],
					row[MASSOC_PACCT]);
			} else {
				object = xstrdup_printf(
					"C = %-10s A = %s",
					row[MASSOC_CLUSTER], row[MASSOC_ACCT]);
			}
			if(assoc->parent_acct) {
				if(!strcasecmp(row[MASSOC_ACCT],
					       assoc->parent_acct)) {
					error("You can't make an account be a "
					      "child of it's self");
					xfree(object);
					continue;
				}

				if(_move_parent(mysql_conn, uid,
						atoi(row[MASSOC_LFT]),
						atoi(row[MASSOC_RGT]),
						row[MASSOC_CLUSTER],
						row[MASSOC_ID],
						row[MASSOC_PACCT],
						assoc->parent_acct)
				   == SLURM_ERROR)
					break;
			}
			account_type = 1;
		}
		list_append(ret_list, object);

		if(!set) {
			xstrfmtcat(name_char, "(id=%s", row[MASSOC_ID]);
			set = 1;
		} else {
			xstrfmtcat(name_char, " || id=%s", row[MASSOC_ID]);
		}
		
		mod_assoc = xmalloc(sizeof(acct_association_rec_t));
		mod_assoc->id = atoi(row[MASSOC_ID]);

		mod_assoc->fairshare = assoc->fairshare;

		mod_assoc->grp_cpus = assoc->grp_cpus;
		mod_assoc->grp_cpu_hours = assoc->grp_cpu_hours;
		mod_assoc->grp_jobs = assoc->grp_jobs;
		mod_assoc->grp_nodes = assoc->grp_nodes;
		mod_assoc->grp_submit_jobs = assoc->grp_submit_jobs;
		mod_assoc->grp_wall = assoc->grp_wall;

		mod_assoc->max_cpus_pj = assoc->max_cpus_pj;
		mod_assoc->max_cpu_mins_pj = assoc->max_cpu_mins_pj;
		mod_assoc->max_jobs = assoc->max_jobs;
		mod_assoc->max_nodes_pj = assoc->max_nodes_pj;
		mod_assoc->max_submit_jobs = assoc->max_submit_jobs;
		mod_assoc->max_wall_pj = assoc->max_wall_pj;
		if(!row[MASSOC_USER][0])
			mod_assoc->parent_acct = xstrdup(assoc->parent_acct);
		if(assoc->qos_list) {
			ListIterator new_qos_itr = 
				list_iterator_create(assoc->qos_list);
			ListIterator curr_qos_itr = NULL;
			char *new_qos = NULL, *curr_qos = NULL;

			mod_assoc->qos_list = list_create(slurm_destroy_char);
			if(!replace_qos)
				slurm_addto_char_list(mod_assoc->qos_list,
						      row[1]);
			curr_qos_itr = 
				list_iterator_create(mod_assoc->qos_list);
			
			while((new_qos = list_next(new_qos_itr))) {
				char *tmp_char = NULL;
				if(new_qos[0] == '-') {
					tmp_char = xstrdup(object+1);
					while((curr_qos =
					       list_next(curr_qos_itr))) {
						if(!strcmp(curr_qos,
							   tmp_char)) {
							list_delete_item(
								curr_qos_itr);
							break;
						}
					}
					xfree(tmp_char);
					list_iterator_reset(curr_qos_itr);
				} else if(new_qos[0] == '+') {
					tmp_char = xstrdup(object+1);
					while((curr_qos =
					       list_next(curr_qos_itr))) {
						if(!strcmp(curr_qos,
							   tmp_char)) {
							break;
						}
					}
					if(!curr_qos)
						list_append(mod_assoc->qos_list,
							    tmp_char);
					else
						xfree(tmp_char);
					list_iterator_reset(curr_qos_itr);
				} else {
					list_append(mod_assoc->qos_list,
						    xstrdup(object));
				}
			}
			list_iterator_destroy(curr_qos_itr);
			list_iterator_destroy(new_qos_itr);			
		}

		if(_addto_update_list(mysql_conn->update_list, 
				      ACCT_MODIFY_ASSOC,
				      mod_assoc) != SLURM_SUCCESS) 
			error("couldn't add to the update list");
		if(account_type) {
			_modify_unset_users(mysql_conn,
					    mod_assoc,
					    row[MASSOC_ACCT],
					    atoi(row[MASSOC_LFT]),
					    atoi(row[MASSOC_RGT]),
					    ret_list);
		}
	}
	mysql_free_result(result);

	if(assoc->parent_acct) {
		if(rc != SLURM_SUCCESS) {
			if(mysql_conn->rollback) {
				mysql_db_rollback(mysql_conn->db_conn);
			}
			list_flush(mysql_conn->update_list);
			list_destroy(ret_list);
			xfree(vals);
			errno = rc;
			return NULL;
		}
	}


	if(!list_count(ret_list)) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything");
		xfree(vals);
		return ret_list;
	}
	xstrcat(name_char, ")");

	if(vals) {
		user_name = uid_to_string((uid_t) uid);
		rc = _modify_common(mysql_conn, DBD_MODIFY_ASSOCS, now,
				    user_name, assoc_table, name_char, vals);
		xfree(user_name);
		if (rc == SLURM_ERROR) {
			if(mysql_conn->rollback) {
				mysql_db_rollback(mysql_conn->db_conn);
			}
			list_flush(mysql_conn->update_list);
			error("Couldn't modify associations");
			list_destroy(ret_list);
			ret_list = NULL;
			goto end_it;
		}
	}

end_it:
	xfree(name_char);
	xfree(vals);

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_remove_users(mysql_conn_t *mysql_conn, uint32_t uid, 
					acct_user_cond_t *user_cond)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	List coord_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *extra = NULL, *query = NULL,
		*name_char = NULL, *assoc_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	acct_user_cond_t user_coord_cond;
	acct_association_cond_t assoc_cond;

	if(!user_cond) {
		error("we need something to remove");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");

	if(user_cond->assoc_cond && user_cond->assoc_cond->user_list
	   && list_count(user_cond->assoc_cond->user_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(user_cond->assoc_cond->user_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(user_cond->def_acct_list && list_count(user_cond->def_acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(user_cond->def_acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "default_acct='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(user_cond->admin_level != ACCT_ADMIN_NOTSET) {
		xstrfmtcat(extra, " && admin_level=%u", user_cond->admin_level);
	}

	if(!extra) {
		error("Nothing to remove");
		return NULL;
	}

	query = xstrdup_printf("select name from %s %s;", user_table, extra);
	xfree(extra);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}

	memset(&user_coord_cond, 0, sizeof(acct_user_cond_t));
	memset(&assoc_cond, 0, sizeof(acct_association_cond_t));
	/* we do not need to free the objects we put in here since
	   they are also placed in a list that will be freed
	*/
	assoc_cond.user_list = list_create(NULL);
	user_coord_cond.assoc_cond = &assoc_cond;

	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		char *object = xstrdup(row[0]);
		acct_user_rec_t *user_rec = NULL;
		
		list_append(ret_list, object);
		list_append(assoc_cond.user_list, object);

		if(!rc) {
			xstrfmtcat(name_char, "name='%s'", object);
			xstrfmtcat(assoc_char, "t2.user='%s'", object);
			rc = 1;
		} else {
			xstrfmtcat(name_char, " || name='%s'", object);
			xstrfmtcat(assoc_char, " || t2.user='%s'", object);
		}
		user_rec = xmalloc(sizeof(acct_user_rec_t));
		user_rec->name = xstrdup(object);
		_addto_update_list(mysql_conn->update_list, ACCT_REMOVE_USER,
				   user_rec);

	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(query);
		list_destroy(assoc_cond.user_list);
		return ret_list;
	}
	xfree(query);

	/* We need to remove these accounts from the coord's that have it */
	coord_list = acct_storage_p_remove_coord(
		mysql_conn, uid, NULL, &user_coord_cond);
	if(coord_list)
		list_destroy(coord_list);
	list_destroy(assoc_cond.user_list);

	user_name = uid_to_string((uid_t) uid);
	rc = _remove_common(mysql_conn, DBD_REMOVE_USERS, now,
			    user_name, user_table, name_char, assoc_char);
	xfree(user_name);
	xfree(name_char);
	if (rc == SLURM_ERROR) {
		list_destroy(ret_list);
		xfree(assoc_char);
		return NULL;
	}

	query = xstrdup_printf(
		"update %s as t2 set deleted=1, mod_time=%d where %s",
		acct_coord_table, now, assoc_char);
	xfree(assoc_char);

	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	if(rc != SLURM_SUCCESS) {
		error("Couldn't remove user coordinators");
		list_destroy(ret_list);
		return NULL;
	}		

	return ret_list;

#else
	return NULL;
#endif
}

extern List acct_storage_p_remove_coord(mysql_conn_t *mysql_conn, uint32_t uid, 
					List acct_list,
					acct_user_cond_t *user_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL, *object = NULL, *extra = NULL, *last_user = NULL;
	char *user_name = NULL;
	time_t now = time(NULL);
	int set = 0, is_admin=0, rc;
	ListIterator itr = NULL;
	acct_user_rec_t *user_rec = NULL;
	List ret_list = NULL;
	List user_list = NULL;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	acct_user_rec_t user;

	if(!user_cond && !acct_list) {
		error("we need something to remove");
		return NULL;
	} else if(user_cond && user_cond->assoc_cond)
		user_list = user_cond->assoc_cond->user_list;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	/* This only works when running though the slurmdbd.
	 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
	 * SLURMDBD!
	 */
	if(slurmdbd_conf) {
		/* we have to check the authentication here in the
		 * plugin since we don't know what accounts are being
		 * referenced until after the query.  Here we will
		 * set if they are an operator or greater and then
		 * check it below after the query.
		 */
		if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
		   || assoc_mgr_get_admin_level(mysql_conn, uid) 
		   >= ACCT_ADMIN_OPERATOR) 
			is_admin = 1;	
		else {
			if(assoc_mgr_fill_in_user(mysql_conn, &user, 1)
			   != SLURM_SUCCESS) {
				error("couldn't get information for this user");
				errno = SLURM_ERROR;
				return NULL;
			}
			if(!user.coord_accts || !list_count(user.coord_accts)) {
				error("This user doesn't have any "
				      "coordinator abilities");
				errno = ESLURM_ACCESS_DENIED;
				return NULL;
			}
		}
	} else {
		/* Setting this here just makes it easier down below
		 * since user will not be filled in.
		 */
		is_admin = 1;
	}

	/* Leave it this way since we are using extra below */

	if(user_list && list_count(user_list)) {
		set = 0;
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, "(");
			
		itr = list_iterator_create(user_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "user='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(acct_list && list_count(acct_list)) {
		set = 0;
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, "(");

		itr = list_iterator_create(acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "acct='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(!extra) {
		errno = SLURM_ERROR;
		debug3("No conditions given");
		return NULL;
	}

	query = xstrdup_printf(
		"select user, acct from %s where deleted=0 && %s order by user",
		acct_coord_table, extra);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result =
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		xfree(extra);
		errno = SLURM_ERROR;
		return NULL;
	}
	xfree(query);
	ret_list = list_create(slurm_destroy_char);
	user_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		if(!is_admin) {
			acct_coord_rec_t *coord = NULL;
			if(!user.coord_accts) { // This should never
						// happen
				error("We are here with no coord accts");
				errno = ESLURM_ACCESS_DENIED;
				list_destroy(ret_list);
				list_destroy(user_list);
				xfree(extra);
				mysql_free_result(result);
				return NULL;
			}
			itr = list_iterator_create(user.coord_accts);
			while((coord = list_next(itr))) {
				if(!strcasecmp(coord->name, row[1]))
					break;
			}
			list_iterator_destroy(itr);

			if(!coord) {
				error("User %s(%d) does not have the "
				      "ability to change this account (%s)",
				      user.name, user.uid, row[1]);
				errno = ESLURM_ACCESS_DENIED;
				list_destroy(ret_list);
				list_destroy(user_list);
				xfree(extra);
				mysql_free_result(result);
				return NULL;
			}
		}
		if(!last_user || strcasecmp(last_user, row[0])) {
			list_append(user_list, xstrdup(row[0]));
			last_user = row[0];
		}
		list_append(ret_list, xstrdup_printf("U = %-9s A = %-10s", 
						     row[0], row[1]));
	}
	mysql_free_result(result);
	
	user_name = uid_to_string((uid_t) uid);
	rc = _remove_common(mysql_conn, DBD_REMOVE_ACCOUNT_COORDS, now,
			    user_name, acct_coord_table, extra, NULL);
	xfree(user_name);
	xfree(extra);
	if (rc == SLURM_ERROR) {
		list_destroy(ret_list);
		list_destroy(user_list);
		errno = SLURM_ERROR;
		return NULL;
	}

	/* get the update list set */
	itr = list_iterator_create(user_list);
	while((last_user = list_next(itr))) {
		user_rec = xmalloc(sizeof(acct_user_rec_t));
		user_rec->name = xstrdup(last_user);
		_get_user_coords(mysql_conn, user_rec);
		_addto_update_list(mysql_conn->update_list, 
				   ACCT_REMOVE_COORD, user_rec);
	}
	list_iterator_destroy(itr);
	list_destroy(user_list);

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_remove_accts(mysql_conn_t *mysql_conn, uint32_t uid, 
					acct_account_cond_t *acct_cond)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	List coord_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *extra = NULL, *query = NULL,
		*name_char = NULL, *assoc_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	if(!acct_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");
	if(acct_cond->assoc_cond 
	   && acct_cond->assoc_cond->acct_list 
	   && list_count(acct_cond->assoc_cond->acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->assoc_cond->acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(acct_cond->description_list 
	   && list_count(acct_cond->description_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->description_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "description='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(acct_cond->organization_list
	   && list_count(acct_cond->organization_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->organization_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "organization='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(!extra) {
		error("Nothing to remove");
		return NULL;
	}

	query = xstrdup_printf("select name from %s %s;", acct_table, extra);
	xfree(extra);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}

	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		char *object = xstrdup(row[0]);
		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(name_char, "name='%s'", object);
			xstrfmtcat(assoc_char, "t2.acct='%s'", object);
			rc = 1;
		} else  {
			xstrfmtcat(name_char, " || name='%s'", object);
			xstrfmtcat(assoc_char, " || t2.acct='%s'", object);
		}
	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(query);
		return ret_list;
	}
	xfree(query);

	/* We need to remove these accounts from the coord's that have it */
	coord_list = acct_storage_p_remove_coord(
		mysql_conn, uid, ret_list, NULL);
	if(coord_list)
		list_destroy(coord_list);

	user_name = uid_to_string((uid_t) uid);
	rc = _remove_common(mysql_conn, DBD_REMOVE_ACCOUNTS, now,
			    user_name, acct_table, name_char, assoc_char);
	xfree(user_name);
	xfree(name_char);
	xfree(assoc_char);
	if (rc == SLURM_ERROR) {
		list_destroy(ret_list);
		return NULL;
	}

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_remove_clusters(mysql_conn_t *mysql_conn,
					   uint32_t uid, 
					   acct_cluster_cond_t *cluster_cond)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *extra = NULL, *query = NULL,
		*name_char = NULL, *assoc_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	if(!cluster_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");
	if(cluster_cond->cluster_list
	   && list_count(cluster_cond->cluster_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(cluster_cond->cluster_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(!extra) {
		error("Nothing to remove");
		return NULL;
	}

	query = xstrdup_printf("select name from %s %s;", cluster_table, extra);
	xfree(extra);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		char *object = xstrdup(row[0]);
		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(name_char, "name='%s'", object);
			xstrfmtcat(extra, "t2.cluster='%s'", object);
			xstrfmtcat(assoc_char, "cluster='%s'", object);
			rc = 1;
		} else  {
			xstrfmtcat(name_char, " || name='%s'", object);
			xstrfmtcat(extra, " || t2.cluster='%s'", object);
			xstrfmtcat(assoc_char, " || cluster='%s'", object);
		}
	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(query);
		return ret_list;
	}
	xfree(query);

	/* We should not need to delete any cluster usage just set it
	 * to deleted */
	xstrfmtcat(query,
		   "update %s set mod_time=%d, deleted=1 where (%s);"
		   "update %s set mod_time=%d, deleted=1 where (%s);"
		   "update %s set mod_time=%d, deleted=1 where (%s);",
		   cluster_day_table, now, assoc_char,
		   cluster_hour_table, now, assoc_char,
		   cluster_month_table, now, assoc_char);
	xfree(assoc_char);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	if(rc != SLURM_SUCCESS) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		list_destroy(ret_list);
		xfree(name_char);
		xfree(extra);
		return NULL;
	}

	assoc_char = xstrdup_printf("t2.acct='root' && (%s)", extra);
	xfree(extra);

	user_name = uid_to_string((uid_t) uid);
	rc = _remove_common(mysql_conn, DBD_REMOVE_CLUSTERS, now,
			    user_name, cluster_table, name_char, assoc_char);
	xfree(user_name);
	xfree(name_char);
	xfree(assoc_char);
	if (rc  == SLURM_ERROR) {
		list_destroy(ret_list);
		return NULL;
	}

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_remove_associations(
	mysql_conn_t *mysql_conn, uint32_t uid, 
	acct_association_cond_t *assoc_cond)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *extra = NULL, *query = NULL,
		*name_char = NULL, *assoc_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0, i = 0, is_admin=0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	acct_user_rec_t user;

	/* if this changes you will need to edit the corresponding 
	 * enum below also t1 is step_table */
	char *rassoc_req_inx[] = {
		"id",
		"acct",
		"parent_acct",
		"cluster",
		"user",
		"partition"
	};
	
	enum {
		RASSOC_ID,
		RASSOC_ACCT,
		RASSOC_PACCT,
		RASSOC_CLUSTER,
		RASSOC_USER,
		RASSOC_PART,
		RASSOC_COUNT
	};

	if(!assoc_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	/* This only works when running though the slurmdbd.
	 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
	 * SLURMDBD!
	 */
	if(slurmdbd_conf) {
		/* we have to check the authentication here in the
		 * plugin since we don't know what accounts are being
		 * referenced until after the query.  Here we will
		 * set if they are an operator or greater and then
		 * check it below after the query.
		 */
		if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
		   || assoc_mgr_get_admin_level(mysql_conn, uid) 
		   >= ACCT_ADMIN_OPERATOR) 
			is_admin = 1;	
		else {
			if(assoc_mgr_fill_in_user(mysql_conn, &user, 1)
			   != SLURM_SUCCESS) {
				error("couldn't get information for this user");
				errno = SLURM_ERROR;
				return NULL;
			}
			if(!user.coord_accts || !list_count(user.coord_accts)) {
				error("This user doesn't have any "
				      "coordinator abilities");
				errno = ESLURM_ACCESS_DENIED;
				return NULL;
			}
		}
	} else {
		/* Setting this here just makes it easier down below
		 * since user will not be filled in.
		 */
		is_admin = 1;
	}

	xstrcat(extra, "where id>0 && deleted=0");

	set = _setup_association_cond_limits(assoc_cond, &extra);

	for(i=0; i<RASSOC_COUNT; i++) {
		if(i) 
			xstrcat(object, ", ");
		xstrcat(object, rassoc_req_inx[i]);
	}

	query = xstrdup_printf("select lft, rgt from %s %s order by lft "
			       "FOR UPDATE;",
			       assoc_table, extra);
	xfree(extra);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
		
	rc = 0;
	while((row = mysql_fetch_row(result))) {
		if(!rc) {
			xstrfmtcat(name_char, "lft between %s and %s",
				   row[0], row[1]);
			rc = 1;
		} else {
			xstrfmtcat(name_char, " || lft between %s and %s",
				   row[0], row[1]);
		}
	}
	mysql_free_result(result);

	if(!name_char) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(query);
		return ret_list;
	}

	xfree(query);
	query = xstrdup_printf("select distinct %s "
			       "from %s where (%s) order by lft;",
			       object,
			       assoc_table, name_char);
	xfree(extra);
	xfree(object);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		if(mysql_conn->rollback) {
			mysql_db_rollback(mysql_conn->db_conn);
		}
		list_flush(mysql_conn->update_list);
		xfree(query);
		xfree(name_char);
		return NULL;
	}

	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		acct_association_rec_t *rem_assoc = NULL;
		if(!is_admin) {
			acct_coord_rec_t *coord = NULL;
			if(!user.coord_accts) { // This should never
						// happen
				error("We are here with no coord accts");
				errno = ESLURM_ACCESS_DENIED;
				goto end_it;
			}
			itr = list_iterator_create(user.coord_accts);
			while((coord = list_next(itr))) {
				if(!strcasecmp(coord->name,
					       row[RASSOC_ACCT]))
					break;
			}
			list_iterator_destroy(itr);

			if(!coord) {
				error("User %s(%d) does not have the "
				      "ability to change this account (%s)",
				      user.name, user.uid, row[RASSOC_ACCT]);
				errno = ESLURM_ACCESS_DENIED;
				goto end_it;
			}
		}
		if(row[RASSOC_PART][0]) { 
			// see if there is a partition name
			object = xstrdup_printf(
				"C = %-10s A = %-10s U = %-9s P = %s",
				row[RASSOC_CLUSTER], row[RASSOC_ACCT],
				row[RASSOC_USER], row[RASSOC_PART]);
		} else if(row[RASSOC_USER][0]){
			object = xstrdup_printf(
				"C = %-10s A = %-10s U = %-9s",
				row[RASSOC_CLUSTER], row[RASSOC_ACCT], 
				row[RASSOC_USER]);
		} else {
			if(row[RASSOC_PACCT][0]) {
				object = xstrdup_printf(
					"C = %-10s A = %s of %s",
					row[RASSOC_CLUSTER], row[RASSOC_ACCT],
					row[RASSOC_PACCT]);
			} else {
				object = xstrdup_printf(
					"C = %-10s A = %s",
					row[RASSOC_CLUSTER], row[RASSOC_ACCT]);
			}
		}
		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(assoc_char, "id=%s", row[RASSOC_ID]);
			rc = 1;
		} else {
			xstrfmtcat(assoc_char, " || id=%s", row[RASSOC_ID]);
		}

		rem_assoc = xmalloc(sizeof(acct_association_rec_t));
		rem_assoc->id = atoi(row[RASSOC_ID]);
		if(_addto_update_list(mysql_conn->update_list, 
				      ACCT_REMOVE_ASSOC,
				      rem_assoc) != SLURM_SUCCESS) 
			error("couldn't add to the update list");

	}
	mysql_free_result(result);

	user_name = uid_to_string((uid_t) uid);
	rc = _remove_common(mysql_conn, DBD_REMOVE_ASSOCS, now,
			    user_name, assoc_table, name_char, assoc_char);
	xfree(user_name);
	xfree(name_char);
	xfree(assoc_char);
	if (rc  == SLURM_ERROR)
		goto end_it;

	return ret_list;
end_it:
	if(mysql_conn->rollback) {
		mysql_db_rollback(mysql_conn->db_conn);
	}
	list_flush(mysql_conn->update_list);
	
	if(ret_list) {
		list_destroy(ret_list);
		ret_list = NULL;
	}
	mysql_free_result(result);

	return NULL;
#else
	return NULL;
#endif
}

extern List acct_storage_p_remove_qos(mysql_conn_t *mysql_conn, uint32_t uid, 
				      acct_qos_cond_t *qos_cond)
{
#ifdef HAVE_MYSQL
	ListIterator itr = NULL;
	List ret_list = NULL;
	int rc = SLURM_SUCCESS;
	char *object = NULL;
	char *extra = NULL, *query = NULL,
		*name_char = NULL, *assoc_char = NULL;
	time_t now = time(NULL);
	char *user_name = NULL;
	int set = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	if(!qos_cond) {
		error("we need something to change");
		return NULL;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	xstrcat(extra, "where deleted=0");
	if(qos_cond->description_list 
	   && list_count(qos_cond->description_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(qos_cond->description_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "description='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(qos_cond->id_list 
	   && list_count(qos_cond->id_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(qos_cond->id_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "id='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(qos_cond->name_list
	   && list_count(qos_cond->name_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(qos_cond->name_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(!extra) {
		error("Nothing to remove");
		return NULL;
	}

	query = xstrdup_printf("select id from %s %s;", qos_table, extra);
	xfree(extra);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}

	rc = 0;
	ret_list = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		char *object = xstrdup(row[0]);
		acct_qos_rec_t *qos_rec = NULL;

		list_append(ret_list, object);
		if(!rc) {
			xstrfmtcat(name_char, "id='%s'", object);
			xstrfmtcat(assoc_char, "qos=replace(qos, ',%s', '')",
				   object);
			rc = 1;
		} else  {
			xstrfmtcat(name_char, " || id='%s'", object); 
			xstrfmtcat(assoc_char, ", qos=replace(qos, ',%s', '')",
				   object);
		}
		qos_rec = xmalloc(sizeof(acct_qos_rec_t));
		qos_rec->name = xstrdup(object);
		_addto_update_list(mysql_conn->update_list, ACCT_REMOVE_QOS,
				   qos_rec);
	}
	mysql_free_result(result);

	if(!list_count(ret_list)) {
		errno = SLURM_NO_CHANGE_IN_DATA;
		debug3("didn't effect anything\n%s", query);
		xfree(query);
		return ret_list;
	}
	xfree(query);

	user_name = uid_to_string((uid_t) uid);
	rc = _remove_common(mysql_conn, DBD_REMOVE_ACCOUNTS, now,
			    user_name, qos_table, name_char, assoc_char);
	xfree(name_char);
	xfree(user_name);
	if (rc == SLURM_ERROR) {
		list_destroy(ret_list);
		return NULL;
	}

	return ret_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_get_users(mysql_conn_t *mysql_conn, uid_t uid, 
				     acct_user_cond_t *user_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL;	
	char *extra = NULL;	
	char *tmp = NULL;	
	List user_list = NULL;
	ListIterator itr = NULL;
	char *object = NULL;
	int set = 0;
	int i=0, is_admin=1;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	uint16_t private_data = 0;
	acct_user_rec_t user;

	/* if this changes you will need to edit the corresponding enum */
	char *user_req_inx[] = {
		"name",
		"default_acct",
		"admin_level"
	};
	enum {
		USER_REQ_NAME,
		USER_REQ_DA,
		USER_REQ_AL,
		USER_REQ_COUNT
	};

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	private_data = slurm_get_private_data();
	if (private_data & PRIVATE_DATA_USERS) {
		/* This only works when running though the slurmdbd.
		 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
		 * SLURMDBD!
		 */
		if(slurmdbd_conf) {
			is_admin = 0;
			/* we have to check the authentication here in the
			 * plugin since we don't know what accounts are being
			 * referenced until after the query.  Here we will
			 * set if they are an operator or greater and then
			 * check it below after the query.
			 */
			if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
			   || assoc_mgr_get_admin_level(mysql_conn, uid) 
			   >= ACCT_ADMIN_OPERATOR) 
				is_admin = 1;	
			else {
				assoc_mgr_fill_in_user(mysql_conn, &user, 1);
			}
		}
	}
	
	if(!user_cond) {
		xstrcat(extra, "where deleted=0");
		goto empty;
	} 
	
	if(user_cond->with_deleted) 
		xstrcat(extra, "where (deleted=0 || deleted=1)");
	else
		xstrcat(extra, "where deleted=0");
		

	if(user_cond->assoc_cond && 
	   user_cond->assoc_cond->user_list
	   && list_count(user_cond->assoc_cond->user_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(user_cond->assoc_cond->user_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(user_cond->def_acct_list && list_count(user_cond->def_acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(user_cond->def_acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "default_acct='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(user_cond->admin_level != ACCT_ADMIN_NOTSET) {
		xstrfmtcat(extra, " && admin_level=%u",
			   user_cond->admin_level);
	}
empty:
	/* This is here to make sure we are looking at only this user
	 * if this flag is set. 
	 */
	if(!is_admin && (private_data & PRIVATE_DATA_USERS)) {
		xstrfmtcat(extra, " && name='%s'", user.name);
	}

	xfree(tmp);
	xstrfmtcat(tmp, "%s", user_req_inx[i]);
	for(i=1; i<USER_REQ_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", user_req_inx[i]);
	}

	query = xstrdup_printf("select %s from %s %s", tmp, user_table, extra);
	xfree(tmp);
	xfree(extra);
	
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	user_list = list_create(destroy_acct_user_rec);

	if(user_cond && user_cond->with_assocs) {
		/* We are going to be freeing the inners of
		   this list in the user->name so we don't
		   free it here
		*/
		if(user_cond->assoc_cond->user_list)
			list_destroy(user_cond->assoc_cond->user_list);
		user_cond->assoc_cond->user_list = list_create(NULL);
	}

	while((row = mysql_fetch_row(result))) {
		acct_user_rec_t *user = xmalloc(sizeof(acct_user_rec_t));
/* 		uid_t pw_uid; */
		list_append(user_list, user);

		user->name =  xstrdup(row[USER_REQ_NAME]);
		user->default_acct = xstrdup(row[USER_REQ_DA]);
		user->admin_level = atoi(row[USER_REQ_AL]);
		
		/* user id will be set on the client since this could be on a
		 * different machine where this user may not exist or
		 * may have a different uid
		 */
/* 		pw_uid = uid_from_string(user->name); */
/* 		if(pw_uid == (uid_t) -1)  */
/* 			user->uid = (uint32_t)NO_VAL; */
/* 		else */
/* 			user->uid = passwd_ptr->pw_uid; */

		if(user_cond && user_cond->with_coords) {
			_get_user_coords(mysql_conn, user);
		}

		if(user_cond && user_cond->with_assocs) {
			if(!user_cond->assoc_cond) {
				user_cond->assoc_cond = xmalloc(
					sizeof(acct_association_cond_t));
			}

			list_append(user_cond->assoc_cond->user_list,
				    user->name);
		}
	}
	mysql_free_result(result);

	if(user_cond && user_cond->with_assocs) {
		ListIterator assoc_itr = NULL;
		acct_user_rec_t *user = NULL;
		acct_association_rec_t *assoc = NULL;
		List assoc_list = acct_storage_p_get_associations(
			mysql_conn, uid, user_cond->assoc_cond);

		if(!assoc_list) {
			error("no associations");
			return user_list;
		}

		itr = list_iterator_create(user_list);
		assoc_itr = list_iterator_create(assoc_list);
		while((user = list_next(itr))) {
			while((assoc = list_next(assoc_itr))) {
				if(strcmp(assoc->user, user->name)) 
					continue;
				
				if(!user->assoc_list)
					user->assoc_list = list_create(
						destroy_acct_association_rec);
				list_append(user->assoc_list, assoc);
				list_remove(assoc_itr);
			}
			list_iterator_reset(assoc_itr);
		}
		list_iterator_destroy(itr);
		list_iterator_destroy(assoc_itr);

		list_destroy(assoc_list);
	}

	return user_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_get_accts(mysql_conn_t *mysql_conn, uid_t uid,
				     acct_account_cond_t *acct_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL;	
	char *extra = NULL;	
	char *tmp = NULL;	
	List acct_list = NULL;
	ListIterator itr = NULL;
	char *object = NULL;
	int set = 0;
	int i=0, is_admin=1;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	uint16_t private_data = 0;
	acct_user_rec_t user;

	/* if this changes you will need to edit the corresponding enum */
	char *acct_req_inx[] = {
		"name",
		"description",
		"organization"
	};
	enum {
		ACCT_REQ_NAME,
		ACCT_REQ_DESC,
		ACCT_REQ_ORG,
		ACCT_REQ_COUNT
	};

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	private_data = slurm_get_private_data();

	if (private_data & PRIVATE_DATA_ACCOUNTS) {
		/* This only works when running though the slurmdbd.
		 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
		 * SLURMDBD!
		 */
		if(slurmdbd_conf) {
			is_admin = 0;
			/* we have to check the authentication here in the
			 * plugin since we don't know what accounts are being
			 * referenced until after the query.  Here we will
			 * set if they are an operator or greater and then
			 * check it below after the query.
			 */
			if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
			   || assoc_mgr_get_admin_level(mysql_conn, uid) 
			   >= ACCT_ADMIN_OPERATOR) 
				is_admin = 1;	
			else {
				assoc_mgr_fill_in_user(mysql_conn, &user, 1);
			}

			if(!is_admin && (!user.coord_accts 
					 || !list_count(user.coord_accts))) {
				errno = ESLURM_ACCESS_DENIED;
				return NULL;
			}
		}
	}
	
	if(!acct_cond) {
		xstrcat(extra, "where deleted=0");
		goto empty;
	} 

	if(acct_cond->with_deleted) 
		xstrcat(extra, "where (deleted=0 || deleted=1)");
	else
		xstrcat(extra, "where deleted=0");

	if(acct_cond->assoc_cond 
	   && acct_cond->assoc_cond->acct_list 
	   && list_count(acct_cond->assoc_cond->acct_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->assoc_cond->acct_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(acct_cond->description_list
	   && list_count(acct_cond->description_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->description_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "description='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(acct_cond->organization_list 
	   && list_count(acct_cond->organization_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(acct_cond->organization_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "organization='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
empty:

	xfree(tmp);
	xstrfmtcat(tmp, "%s", acct_req_inx[i]);
	for(i=1; i<ACCT_REQ_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", acct_req_inx[i]);
	}

	/* This is here to make sure we are looking at only this user
	 * if this flag is set.  We also include any accounts they may be
	 * coordinator of.
	 */
	if(!is_admin && (private_data & PRIVATE_DATA_ACCOUNTS)) {
		acct_coord_rec_t *coord = NULL;
		set = 0;
		itr = list_iterator_create(user.coord_accts);
		while((coord = list_next(itr))) {
			if(set) {
				xstrfmtcat(extra, " || name='%s'", coord->name);
			} else {
				set = 1;
				xstrfmtcat(extra, " && (name='%s'",coord->name);
			}
		}		
		list_iterator_destroy(itr);
		if(set)
			xstrcat(extra,")");
	}

	query = xstrdup_printf("select %s from %s %s", tmp, acct_table, extra);
	xfree(tmp);
	xfree(extra);
	
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	acct_list = list_create(destroy_acct_account_rec);
	
	if(acct_cond && acct_cond->with_assocs) {
		/* We are going to be freeing the inners of
			   this list in the acct->name so we don't
			   free it here
			*/
		if(acct_cond->assoc_cond->acct_list) 
			list_destroy(acct_cond->assoc_cond->acct_list);
		acct_cond->assoc_cond->acct_list = list_create(NULL);
	}

	while((row = mysql_fetch_row(result))) {
		acct_account_rec_t *acct = xmalloc(sizeof(acct_account_rec_t));
		list_append(acct_list, acct);

		acct->name =  xstrdup(row[ACCT_REQ_NAME]);
		acct->description = xstrdup(row[ACCT_REQ_DESC]);
		acct->organization = xstrdup(row[ACCT_REQ_ORG]);

		if(acct_cond && acct_cond->with_coords) {
			_get_account_coords(mysql_conn, acct);
		}

		if(acct_cond && acct_cond->with_assocs) {
			if(!acct_cond->assoc_cond) {
				acct_cond->assoc_cond = xmalloc(
					sizeof(acct_association_cond_t));
			}

			list_append(acct_cond->assoc_cond->acct_list,
				    acct->name);
		}
	}
	mysql_free_result(result);

	if(acct_cond && acct_cond->with_assocs) {
		ListIterator assoc_itr = NULL;
		acct_account_rec_t *acct = NULL;
		acct_association_rec_t *assoc = NULL;
		List assoc_list = acct_storage_p_get_associations(
			mysql_conn, uid, acct_cond->assoc_cond);

		if(!assoc_list) {
			error("no associations");
			return acct_list;
		}

		itr = list_iterator_create(acct_list);
		assoc_itr = list_iterator_create(assoc_list);
		while((acct = list_next(itr))) {
			while((assoc = list_next(assoc_itr))) {
				if(strcmp(assoc->acct, acct->name)) 
					continue;
				
				if(!acct->assoc_list)
					acct->assoc_list = list_create(
						destroy_acct_association_rec);
				list_append(acct->assoc_list, assoc);
				list_remove(assoc_itr);
			}
			list_iterator_reset(assoc_itr);
		}
		list_iterator_destroy(itr);
		list_iterator_destroy(assoc_itr);

		list_destroy(assoc_list);
	}

	return acct_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_get_clusters(mysql_conn_t *mysql_conn, uid_t uid, 
					acct_cluster_cond_t *cluster_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL;	
	char *extra = NULL;	
	char *tmp = NULL;	
	List cluster_list = NULL;
	ListIterator itr = NULL;
	char *object = NULL;
	int set = 0;
	int i=0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	acct_association_cond_t assoc_cond;
	ListIterator assoc_itr = NULL;
	acct_cluster_rec_t *cluster = NULL;
	acct_association_rec_t *assoc = NULL;
	List assoc_list = NULL;

	/* if this changes you will need to edit the corresponding enum */
	char *cluster_req_inx[] = {
		"name",
		"control_host",
		"control_port",
		"rpc_version"
	};
	enum {
		CLUSTER_REQ_NAME,
		CLUSTER_REQ_CH,
		CLUSTER_REQ_CP,
		CLUSTER_REQ_VERSION,
		CLUSTER_REQ_COUNT
	};

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

		
	if(!cluster_cond) {
		xstrcat(extra, "where deleted=0");
		goto empty;
	}

	if(cluster_cond->with_deleted) 
		xstrcat(extra, "where (deleted=0 || deleted=1)");
	else
		xstrcat(extra, "where deleted=0");

	if(cluster_cond->cluster_list 
	   && list_count(cluster_cond->cluster_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(cluster_cond->cluster_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

empty:

	xfree(tmp);
	i=0;
	xstrfmtcat(tmp, "%s", cluster_req_inx[i]);
	for(i=1; i<CLUSTER_REQ_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", cluster_req_inx[i]);
	}

	query = xstrdup_printf("select %s from %s %s", 
			       tmp, cluster_table, extra);
	xfree(tmp);
	xfree(extra);
	
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	cluster_list = list_create(destroy_acct_cluster_rec);

	memset(&assoc_cond, 0, sizeof(acct_association_cond_t));

	assoc_cond.cluster_list = list_create(NULL);
	assoc_cond.acct_list = list_create(NULL);
	list_append(assoc_cond.acct_list, "root");

	assoc_cond.user_list = list_create(NULL);
	list_append(assoc_cond.user_list, "");

	while((row = mysql_fetch_row(result))) {
		cluster = xmalloc(sizeof(acct_cluster_rec_t));
		list_append(cluster_list, cluster);

		cluster->name = xstrdup(row[CLUSTER_REQ_NAME]);

		list_append(assoc_cond.cluster_list, cluster->name);

		/* get the usage if requested */
		if(cluster_cond->with_usage) {
			clusteracct_storage_p_get_usage(
				mysql_conn, uid, cluster,
				cluster_cond->usage_start,
				cluster_cond->usage_end);
		}

		cluster->control_host = xstrdup(row[CLUSTER_REQ_CH]);
		cluster->control_port = atoi(row[CLUSTER_REQ_CP]);
		cluster->rpc_version = atoi(row[CLUSTER_REQ_VERSION]);
	}
	mysql_free_result(result);
	
	assoc_list = acct_storage_p_get_associations(mysql_conn,
						     uid, &assoc_cond);
	list_destroy(assoc_cond.cluster_list);
	list_destroy(assoc_cond.acct_list);
	list_destroy(assoc_cond.user_list);

	if(!assoc_list) 
		return cluster_list;
	
	
	itr = list_iterator_create(cluster_list);
	assoc_itr = list_iterator_create(assoc_list);
	while((cluster = list_next(itr))) {
		while((assoc = list_next(assoc_itr))) {
			if(strcmp(assoc->cluster, cluster->name)) 
				continue;
			
			if(cluster->root_assoc) {
				debug("This cluster %s already has "
				      "an association.");
				continue;
			}
			cluster->root_assoc = assoc;
			list_remove(assoc_itr);
		}
		list_iterator_reset(assoc_itr);
	}
	list_iterator_destroy(itr);
	list_iterator_destroy(assoc_itr);
	if(list_count(assoc_list))
		info("I have %d left over associations", 
		     list_count(assoc_list));
	list_destroy(assoc_list);

	return cluster_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_get_associations(mysql_conn_t *mysql_conn,
					    uid_t uid, 
					    acct_association_cond_t *assoc_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL;	
	char *extra = NULL;	
	char *tmp = NULL;	
	List assoc_list = NULL;
	ListIterator itr = NULL;
	int set = 0;
	int i=0, is_admin=1;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	uint32_t parent_mj = INFINITE;
        uint32_t parent_msj = INFINITE;
	uint32_t parent_mcpj = INFINITE;
	uint32_t parent_mnpj = INFINITE;
	uint32_t parent_mwpj = INFINITE;
	uint64_t parent_mcmpj = INFINITE;
	char *parent_qos = NULL;
	char *last_acct = NULL;
	char *last_acct_parent = NULL;
	char *last_cluster = NULL;
	char *last_cluster2 = NULL;
	uint32_t user_parent_id = 0;
	uint32_t acct_parent_id = 0;
	uint16_t private_data = 0;
	acct_user_rec_t user;

	/* needed if we don't have an assoc_cond */
	uint16_t without_parent_info = 0;
	uint16_t without_parent_limits = 0;
	uint16_t with_usage = 0;

	/* if this changes you will need to edit the corresponding enum */
	char *assoc_req_inx[] = {
		"id",
		"lft",
		"rgt",
		"user",
		"acct",
		"cluster",
		"partition",
		"parent_acct",
		"fairshare",
		"grp_jobs",
		"grp_submit_jobs",
		"grp_cpus",
		"grp_nodes",
		"grp_wall",
		"grp_cpu_hours",
		"max_jobs",
		"max_submit_jobs",
		"max_cpus_per_job",
		"max_nodes_per_job",
		"max_wall_duration_per_job",
		"max_cpu_mins_per_job",
		"qos",
	};
	enum {
		ASSOC_REQ_ID,
		ASSOC_REQ_LFT,
		ASSOC_REQ_RGT,
		ASSOC_REQ_USER,
		ASSOC_REQ_ACCT,
		ASSOC_REQ_CLUSTER,
		ASSOC_REQ_PART,
		ASSOC_REQ_PARENT,
		ASSOC_REQ_FS,
		ASSOC_REQ_GJ,
		ASSOC_REQ_GSJ,
		ASSOC_REQ_GC,
		ASSOC_REQ_GN,
		ASSOC_REQ_GW,
		ASSOC_REQ_GCH,
		ASSOC_REQ_MJ,
		ASSOC_REQ_MSJ,
		ASSOC_REQ_MCPJ,
		ASSOC_REQ_MNPJ,
		ASSOC_REQ_MWPJ,
		ASSOC_REQ_MCMPJ,
		ASSOC_REQ_QOS,
		ASSOC_REQ_COUNT
	};

	enum {
		ASSOC2_REQ_PARENT_ID,
		ASSOC2_REQ_MJ,
		ASSOC2_REQ_MSJ,
		ASSOC2_REQ_MCPJ,
		ASSOC2_REQ_MNPJ,
		ASSOC2_REQ_MWPJ,
		ASSOC2_REQ_MCMPJ,
		ASSOC2_REQ_QOS,
	};

	if(!assoc_cond) {
		xstrcat(extra, "where deleted=0");
		goto empty;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	private_data = slurm_get_private_data();
	if (private_data & PRIVATE_DATA_USERS) {
		/* This only works when running though the slurmdbd.
		 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
		 * SLURMDBD!
		 */
		if(slurmdbd_conf) {
			is_admin = 0;
			/* we have to check the authentication here in the
			 * plugin since we don't know what accounts are being
			 * referenced until after the query.  Here we will
			 * set if they are an operator or greater and then
			 * check it below after the query.
			 */
			if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
			   || assoc_mgr_get_admin_level(mysql_conn, uid) 
			   >= ACCT_ADMIN_OPERATOR) 
				is_admin = 1;	
			else {
				assoc_mgr_fill_in_user(mysql_conn, &user, 1);
			}
		}
	}

	if(assoc_cond->with_deleted) 
		xstrcat(extra, "where (deleted=0 || deleted=1)");
	else
		xstrcat(extra, "where deleted=0");

	set = _setup_association_cond_limits(assoc_cond, &extra);

	with_usage = assoc_cond->with_usage;
	without_parent_limits = assoc_cond->without_parent_limits;
	without_parent_info = assoc_cond->without_parent_info;
empty:
	xfree(tmp);
	xstrfmtcat(tmp, "%s", assoc_req_inx[i]);
	for(i=1; i<ASSOC_REQ_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", assoc_req_inx[i]);
	}
	
	/* this is here to make sure we are looking at only this user
	 * if this flag is set.  We also include any accounts they may be
	 * coordinator of.
	 */
	if(!is_admin && (private_data & PRIVATE_DATA_USERS)) {
		query = xstrdup_printf("select lft from %s where user='%s'", 
				       assoc_table, user.name);
		if(user.coord_accts) {
			acct_coord_rec_t *coord = NULL;
			itr = list_iterator_create(user.coord_accts);
			while((coord = list_next(itr))) {
				xstrfmtcat(query, " || acct='%s'",
					   coord->name);
			}
			list_iterator_destroy(itr);
		}
		debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(extra);
			xfree(query);
			return NULL;
		}
		xfree(query);
		set = 0;
		while((row = mysql_fetch_row(result))) {
			if(set) {
				xstrfmtcat(extra,
					   " || (%s between lft and rgt)",
					   row[0]);
			} else {
				set = 1;
				xstrfmtcat(extra,
					" && ((%s between lft and rgt)",
					row[0]);
			}
		}		
		if(set)
			xstrcat(extra,")");
		mysql_free_result(result);
	}
	
	query = xstrdup_printf("select %s from %s %s order by lft;", 
			       tmp, assoc_table, extra);
	xfree(tmp);
	xfree(extra);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	assoc_list = list_create(destroy_acct_association_rec);

	while((row = mysql_fetch_row(result))) {
		acct_association_rec_t *assoc =
			xmalloc(sizeof(acct_association_rec_t));
		MYSQL_RES *result2 = NULL;
		MYSQL_ROW row2;
		
		list_append(assoc_list, assoc);
		
		assoc->id = atoi(row[ASSOC_REQ_ID]);
		assoc->lft = atoi(row[ASSOC_REQ_LFT]);
		assoc->rgt = atoi(row[ASSOC_REQ_RGT]);

		if(row[ASSOC_REQ_USER][0])
			assoc->user = xstrdup(row[ASSOC_REQ_USER]);
		assoc->acct = xstrdup(row[ASSOC_REQ_ACCT]);
		assoc->cluster = xstrdup(row[ASSOC_REQ_CLUSTER]);

		if(row[ASSOC_REQ_GJ])
			assoc->grp_jobs = atoi(row[ASSOC_REQ_GJ]);
		else
			assoc->grp_jobs = INFINITE;

		if(row[ASSOC_REQ_GSJ])
			assoc->grp_submit_jobs = atoi(row[ASSOC_REQ_GSJ]);
		else
			assoc->grp_submit_jobs = INFINITE;

		if(row[ASSOC_REQ_GC])
			assoc->grp_cpus = atoi(row[ASSOC_REQ_GC]);
		else
			assoc->grp_cpus = INFINITE;

		if(row[ASSOC_REQ_GN])
			assoc->grp_nodes = atoi(row[ASSOC_REQ_GN]);
		else
			assoc->grp_nodes = INFINITE;
		if(row[ASSOC_REQ_GW])
			assoc->grp_wall = atoi(row[ASSOC_REQ_GW]);
		else
			assoc->grp_wall = INFINITE;

		if(row[ASSOC_REQ_GCH])
			assoc->grp_cpu_hours = atoll(row[ASSOC_REQ_GCH]);
		else
			assoc->grp_cpu_hours = INFINITE;

		/* get the usage if requested */
		if(with_usage) {
			acct_storage_p_get_usage(mysql_conn, uid, assoc,
						 assoc_cond->usage_start,
						 assoc_cond->usage_end);
		}

		if(!without_parent_info 
		   && row[ASSOC_REQ_PARENT][0]) {
/* 			info("got %s?=%s and %s?=%s", */
/* 			     row[ASSOC_REQ_PARENT], last_acct_parent, */
/* 			     row[ASSOC_REQ_CLUSTER], last_cluster); */
			if(!last_acct_parent || !last_cluster
			   || strcmp(row[ASSOC_REQ_PARENT], last_acct_parent)
			   || strcmp(row[ASSOC_REQ_CLUSTER], last_cluster)) {
				query = xstrdup_printf(
					"select id from %s where user='' "
					"and deleted = 0 and acct='%s' "
					"and cluster='%s';", 
					assoc_table, row[ASSOC_REQ_PARENT],
					row[ASSOC_REQ_CLUSTER]);
				debug4("%d(%d) query\n%s",
				       mysql_conn->conn, __LINE__, query);

				if(!(result2 = mysql_db_query_ret(
					     mysql_conn->db_conn,
					     query, 1))) {
					xfree(query);
					break;
				}
				xfree(query);
				row2 = mysql_fetch_row(result2);
				last_acct_parent = row[ASSOC_REQ_PARENT];
				last_cluster = row[ASSOC_REQ_CLUSTER];
				acct_parent_id = atoi(row2[0]);	
				mysql_free_result(result2);
			}
			assoc->parent_acct = xstrdup(row[ASSOC_REQ_PARENT]);
			assoc->parent_id = acct_parent_id;
		} 

		if(row[ASSOC_REQ_PART][0])
			assoc->partition = xstrdup(row[ASSOC_REQ_PART]);
		if(row[ASSOC_REQ_FS])
			assoc->fairshare = atoi(row[ASSOC_REQ_FS]);
		else
			assoc->fairshare = 1;

		if((!last_acct || !last_cluster2 
		    || strcmp(row[ASSOC_REQ_ACCT], last_acct)
		    || strcmp(row[ASSOC_REQ_CLUSTER], last_cluster2))) {
			query = xstrdup_printf(
				"call get_parent_limits('%s', '%s', '%s', %u);"
				"select @par_id, @mj, @msj, @mcpj, "
				"@mnpj, @mwpj, @mcmpj, @qos;", 
				assoc_table, row[ASSOC_REQ_ACCT],
				row[ASSOC_REQ_CLUSTER],
				without_parent_limits);
			debug4("%d(%d) query\n%s",
			       mysql_conn->conn, __LINE__, query);
			if(!(result2 = mysql_db_query_ret(
				     mysql_conn->db_conn, query, 1))) {
				xfree(query);
				break;
			}
			xfree(query);
			
			row2 = mysql_fetch_row(result2);
			user_parent_id = atoi(row2[ASSOC2_REQ_PARENT_ID]);
			if(!without_parent_limits) {
				if(row2[ASSOC2_REQ_MCMPJ])
					parent_mcmpj =
						atoi(row2[ASSOC2_REQ_MCMPJ]);
				else
					parent_mcmpj = INFINITE;
				
				if(row2[ASSOC2_REQ_MCPJ])
					parent_mcpj =
						atoi(row2[ASSOC2_REQ_MCPJ]);
				else
					parent_mcpj = INFINITE;
				
				if(row2[ASSOC2_REQ_MJ])
					parent_mj = atoi(row2[ASSOC2_REQ_MJ]);
				else
					parent_mj = INFINITE;
				
				if(row2[ASSOC2_REQ_MNPJ])
					parent_mnpj =
						atoi(row2[ASSOC2_REQ_MNPJ]);
				else
					parent_mnpj = INFINITE;
				
				if(row2[ASSOC2_REQ_MWPJ])
					parent_mwpj =
						atoi(row2[ASSOC2_REQ_MWPJ]);
				else
					parent_mwpj = INFINITE;
				
				if(row2[ASSOC2_REQ_MCMPJ])
					parent_mcmpj =
						atoll(row2[ASSOC2_REQ_MCMPJ]);
				else 
					parent_mcmpj = INFINITE;

				xfree(parent_qos);
				if(row2[ASSOC2_REQ_QOS][0])
					parent_qos =
						xstrdup(row2[ASSOC2_REQ_QOS]);
				else 
					parent_qos = NULL;
			
				if(row2[ASSOC2_REQ_MSJ])
					parent_msj = atoi(row2[ASSOC2_REQ_MSJ]);
				else
					parent_msj = INFINITE;
			}
			last_acct = row[ASSOC_REQ_ACCT];
			last_cluster2 = row[ASSOC_REQ_CLUSTER];
			mysql_free_result(result2);
		}
		if(row[ASSOC_REQ_MJ])
			assoc->max_jobs = atoi(row[ASSOC_REQ_MJ]);
		else
			assoc->max_jobs = parent_mj;

		if(row[ASSOC_REQ_MSJ])
			assoc->max_submit_jobs = atoi(row[ASSOC_REQ_MSJ]);
		else
			assoc->max_submit_jobs = parent_msj;

		if(row[ASSOC_REQ_MCPJ])
			assoc->max_cpus_pj = 
				atoi(row[ASSOC_REQ_MCPJ]);
		else
			assoc->max_cpus_pj = parent_mcpj;

		if(row[ASSOC_REQ_MNPJ])
			assoc->max_nodes_pj = 
				atoi(row[ASSOC_REQ_MNPJ]);
		else
			assoc->max_nodes_pj = parent_mnpj;

		if(row[ASSOC_REQ_MWPJ])
			assoc->max_wall_pj = 
				atoi(row[ASSOC_REQ_MWPJ]);
		else
			assoc->max_wall_pj = parent_mwpj;

		if(row[ASSOC_REQ_MCMPJ])
			assoc->max_cpu_mins_pj = 
				atoi(row[ASSOC_REQ_MCMPJ]);
		else
			assoc->max_cpu_mins_pj = parent_mcmpj;

		assoc->qos_list = list_create(slurm_destroy_char);
		if(row[ASSOC_REQ_QOS][0]) 
			slurm_addto_char_list(assoc->qos_list,
					      row[ASSOC_REQ_QOS]);
		else if(parent_qos) 			
			slurm_addto_char_list(assoc->qos_list, parent_qos);
		
		/* don't do this unless this is an user association */
		if(assoc->user && assoc->parent_id != acct_parent_id) 
			assoc->parent_id = user_parent_id;

		//info("parent id is %d", assoc->parent_id);
		//log_assoc_rec(assoc);
	}
	mysql_free_result(result);

	xfree(parent_qos);

	return assoc_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_get_qos(mysql_conn_t *mysql_conn, uid_t uid,
				   acct_qos_cond_t *qos_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL;	
	char *extra = NULL;	
	char *tmp = NULL;	
	List qos_list = NULL;
	ListIterator itr = NULL;
	char *object = NULL;
	int set = 0;
	int i=0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	/* if this changes you will need to edit the corresponding enum */
	char *qos_req_inx[] = {
		"name",
		"description",
		"id"
	};
	enum {
		QOS_REQ_NAME,
		QOS_REQ_DESC,
		QOS_REQ_ID,
		QOS_REQ_COUNT
	};

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;


	
	if(!qos_cond) {
		xstrcat(extra, "where deleted=0");
		goto empty;
	} 
	
	if(qos_cond->with_deleted) 
		xstrcat(extra, "where (deleted=0 || deleted=1)");
	else
		xstrcat(extra, "where deleted=0");
		

	if(qos_cond->description_list 
	   && list_count(qos_cond->description_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(qos_cond->description_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "description='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(qos_cond->id_list 
	   && list_count(qos_cond->id_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(qos_cond->id_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "id='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}
	
	if(qos_cond->name_list
	   && list_count(qos_cond->name_list)) {
		set = 0;
		xstrcat(extra, " && (");
		itr = list_iterator_create(qos_cond->name_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "name='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

empty:

	xfree(tmp);
	xstrfmtcat(tmp, "%s", qos_req_inx[i]);
	for(i=1; i<QOS_REQ_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", qos_req_inx[i]);
	}

	query = xstrdup_printf("select %s from %s %s", tmp, qos_table, extra);
	xfree(tmp);
	xfree(extra);
	
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	qos_list = list_create(destroy_acct_qos_rec);

	while((row = mysql_fetch_row(result))) {
		acct_qos_rec_t *qos = xmalloc(sizeof(acct_qos_rec_t));
		list_append(qos_list, qos);

		qos->description = xstrdup(row[QOS_REQ_DESC]);
		qos->id = atoi(row[QOS_REQ_ID]);
		qos->name =  xstrdup(row[QOS_REQ_NAME]);
	}
	mysql_free_result(result);

	return qos_list;
#else
	return NULL;
#endif
}

extern List acct_storage_p_get_txn(mysql_conn_t *mysql_conn, uid_t uid,
				   acct_txn_cond_t *txn_cond)
{
#ifdef HAVE_MYSQL
	char *query = NULL;	
	char *extra = NULL;	
	char *tmp = NULL;	
	List txn_list = NULL;
	ListIterator itr = NULL;
	char *object = NULL;
	int set = 0;
	int i=0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

	/* if this changes you will need to edit the corresponding enum */
	char *txn_req_inx[] = {
		"id",
		"timestamp",
		"action",
		"name",
		"actor",
		"info"
	};
	enum {
		TXN_REQ_ID,
		TXN_REQ_TS,
		TXN_REQ_ACTION,
		TXN_REQ_NAME,
		TXN_REQ_ACTOR,
		TXN_REQ_INFO,
		TXN_REQ_COUNT
	};

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;

	if(!txn_cond) 
		goto empty;

	if(txn_cond->action_list && list_count(txn_cond->action_list)) {
		set = 0;
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, " where (");
		itr = list_iterator_create(txn_cond->action_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "action='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(txn_cond->actor_list && list_count(txn_cond->actor_list)) {
		set = 0;
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, " where (");
		itr = list_iterator_create(txn_cond->actor_list);
		while((object = list_next(itr))) {
			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "actor='%s'", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(txn_cond->id_list && list_count(txn_cond->id_list)) {
		set = 0;
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, " where (");
		itr = list_iterator_create(txn_cond->id_list);
		while((object = list_next(itr))) {
			char *ptr = NULL;
			long num = strtol(object, &ptr, 10);
			if ((num == 0) && ptr && ptr[0]) {
				error("Invalid value for txn id (%s)",
				      object);
				xfree(extra);
				list_iterator_destroy(itr);
				return NULL;
			}

			if(set) 
				xstrcat(extra, " || ");
			xstrfmtcat(extra, "id=%s", object);
			set = 1;
		}
		list_iterator_destroy(itr);
		xstrcat(extra, ")");
	}

	if(txn_cond->time_start && txn_cond->time_end) {
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, " where (");
		xstrfmtcat(extra, "timestamp < %d && timestamp >= %d)", 
			   txn_cond->time_end, txn_cond->time_start);
	} else if(txn_cond->time_start) {
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, " where (");
		xstrfmtcat(extra, "timestamp >= %d)", txn_cond->time_start);
		
	} else if(txn_cond->time_end) {
		if(extra)
			xstrcat(extra, " && (");
		else
			xstrcat(extra, " where (");
		xstrfmtcat(extra, "timestamp < %d)", txn_cond->time_end);
	}
empty:
	xfree(tmp);
	xstrfmtcat(tmp, "%s", txn_req_inx[i]);
	for(i=1; i<TXN_REQ_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", txn_req_inx[i]);
	}

	query = xstrdup_printf("select %s from %s", tmp, txn_table);

	if(extra) {
		xstrfmtcat(query, "%s", extra);
		xfree(extra);
	}
	xstrcat(query, " order by timestamp;");

	xfree(tmp);

	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return NULL;
	}
	xfree(query);

	txn_list = list_create(destroy_acct_txn_rec);
	
	while((row = mysql_fetch_row(result))) {
		acct_txn_rec_t *txn = xmalloc(sizeof(acct_txn_rec_t));

		list_append(txn_list, txn);
		
		txn->action = atoi(row[TXN_REQ_ACTION]);
		txn->actor_name = xstrdup(row[TXN_REQ_ACTOR]);
		txn->id = atoi(row[TXN_REQ_ID]);
		txn->set_info = xstrdup(row[TXN_REQ_INFO]);
		txn->timestamp = atoi(row[TXN_REQ_TS]);
		txn->where_query = xstrdup(row[TXN_REQ_NAME]);
	}
	mysql_free_result(result);

	return txn_list;
#else
	return NULL;
#endif
}

extern int acct_storage_p_get_usage(mysql_conn_t *mysql_conn, uid_t uid,
				    acct_association_rec_t *acct_assoc,
				    time_t start, time_t end)
{
#ifdef HAVE_MYSQL
	int rc = SLURM_SUCCESS;
	int i=0, is_admin=1;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *tmp = NULL;
	char *my_usage_table = assoc_day_table;
	time_t my_time = time(NULL);
	struct tm start_tm;
	struct tm end_tm;
	char *query = NULL;
	uint16_t private_data = 0;
	acct_user_rec_t user;

	char *assoc_req_inx[] = {
		"t1.id",
		"t1.period_start",
		"t1.alloc_cpu_secs"
	};
	
	enum {
		ASSOC_ID,
		ASSOC_START,
		ASSOC_ACPU,
		ASSOC_COUNT
	};

	if(!acct_assoc->id) {
		error("We need a assoc id to set data for");
		return SLURM_ERROR;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	memset(&user, 0, sizeof(acct_user_rec_t));
	user.uid = uid;

	private_data = slurm_get_private_data();
	if (private_data & PRIVATE_DATA_USAGE) {
		/* This only works when running though the slurmdbd.
		 * THERE IS NO AUTHENTICATION WHEN RUNNNING OUT OF THE
		 * SLURMDBD!
		 */
		if(slurmdbd_conf) {
			is_admin = 0;
			/* we have to check the authentication here in the
			 * plugin since we don't know what accounts are being
			 * referenced until after the query.  Here we will
			 * set if they are an operator or greater and then
			 * check it below after the query.
			 */
			if((uid == slurmdbd_conf->slurm_user_id || uid == 0)
			   || assoc_mgr_get_admin_level(mysql_conn, uid) 
			   >= ACCT_ADMIN_OPERATOR) 
				is_admin = 1;	
			else {
				assoc_mgr_fill_in_user(mysql_conn, &user, 1);
			}
			
			if(!is_admin) {
				ListIterator itr = NULL;
				acct_coord_rec_t *coord = NULL;

				if(acct_assoc->user && 
				   !strcmp(acct_assoc->user, user.name)) 
					goto is_user;
				
				if(!user.coord_accts) {
					debug4("This user isn't a coord.");
					goto bad_user;
				}

				if(!acct_assoc->acct) {
					debug("No account name given "
					      "in association.");
					goto bad_user;				
				}
				
				itr = list_iterator_create(user.coord_accts);
				while((coord = list_next(itr))) {
					if(!strcasecmp(coord->name, 
						       acct_assoc->acct))
						break;
				}
				list_iterator_destroy(itr);
				
				if(coord) 
					goto is_user;
				
			bad_user:
				errno = ESLURM_ACCESS_DENIED;
				return SLURM_ERROR;
			}
		}
	}
is_user:

	/* Default is going to be the last day */
	if(!end) {
		if(!localtime_r(&my_time, &end_tm)) {
			error("Couldn't get localtime from end %d",
			      my_time);
			return SLURM_ERROR;
		}
		end_tm.tm_hour = 0;
		end = mktime(&end_tm);		
	} else {
		if(!localtime_r(&end, &end_tm)) {
			error("Couldn't get localtime from user end %d",
			      my_time);
			return SLURM_ERROR;
		}
	}
	end_tm.tm_sec = 0;
	end_tm.tm_min = 0;
	end_tm.tm_isdst = -1;
	end = mktime(&end_tm);		

	if(!start) {
		if(!localtime_r(&my_time, &start_tm)) {
			error("Couldn't get localtime from start %d",
			      my_time);
			return SLURM_ERROR;
		}
		start_tm.tm_hour = 0;
		start_tm.tm_mday--;
		start = mktime(&start_tm);		
	} else {
		if(!localtime_r(&start, &start_tm)) {
			error("Couldn't get localtime from user start %d",
			      my_time);
			return SLURM_ERROR;
		}
	}
	start_tm.tm_sec = 0;
	start_tm.tm_min = 0;
	start_tm.tm_isdst = -1;
	start = mktime(&start_tm);		

	if(end-start < 3600) {
		end = start + 3600;
		if(!localtime_r(&end, &end_tm)) {
			error("2 Couldn't get localtime from user end %d",
			      my_time);
			return SLURM_ERROR;
		}
	}
	/* check to see if we are off day boundaries or on month
	 * boundaries other wise use the day table.
	 */
	if(start_tm.tm_hour || end_tm.tm_hour || (end-start < 86400)) 
		my_usage_table = assoc_hour_table;
	else if(start_tm.tm_mday == 0 && end_tm.tm_mday == 0 
		&& (end-start > 86400))
		my_usage_table = assoc_month_table;
		
	xfree(tmp);
	i=0;
	xstrfmtcat(tmp, "%s", assoc_req_inx[i]);
	for(i=1; i<ASSOC_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", assoc_req_inx[i]);
	}

	query = xstrdup_printf(
		"select %s from %s as t1, %s as t2, %s as t3 "
		"where (t1.period_start < %d && t1.period_start >= %d) "
		"&& t1.id=t2.id && t3.id=%u && "
		"t2.lft between t3.lft and t3.rgt "
		"order by t1.id, period_start;",
		tmp, my_usage_table, assoc_table, assoc_table, end, start,
		acct_assoc->id);
	xfree(tmp);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	if(!acct_assoc->accounting_list)
		acct_assoc->accounting_list =
			list_create(destroy_acct_accounting_rec);

	while((row = mysql_fetch_row(result))) {
		acct_accounting_rec_t *accounting_rec =
			xmalloc(sizeof(acct_accounting_rec_t));
		accounting_rec->assoc_id = atoi(row[ASSOC_ID]);
		accounting_rec->period_start = atoi(row[ASSOC_START]);
		accounting_rec->alloc_secs = atoll(row[ASSOC_ACPU]);
		list_append(acct_assoc->accounting_list, accounting_rec);
	}
	mysql_free_result(result);
	
	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int acct_storage_p_roll_usage(mysql_conn_t *mysql_conn, 
				     time_t sent_start)
{
#ifdef HAVE_MYSQL
	int rc = SLURM_SUCCESS;
	int i = 0;
	time_t my_time = time(NULL);
	struct tm start_tm;
	struct tm end_tm;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *query = NULL;
	char *tmp = NULL;
	time_t last_hour = sent_start;
	time_t last_day = sent_start;
	time_t last_month = sent_start;
	time_t start_time = 0;
  	time_t end_time = 0;
	DEF_TIMERS;

	char *update_req_inx[] = {
		"hourly_rollup",
		"daily_rollup",
		"monthly_rollup"
	};
	
	enum {
		UPDATE_HOUR,
		UPDATE_DAY,
		UPDATE_MONTH,
		UPDATE_COUNT
	};

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	if(!sent_start) {
		i=0;
		xstrfmtcat(tmp, "%s", update_req_inx[i]);
		for(i=1; i<UPDATE_COUNT; i++) {
			xstrfmtcat(tmp, ", %s", update_req_inx[i]);
		}
		query = xstrdup_printf("select %s from %s",
				       tmp, last_ran_table);
		xfree(tmp);
		
		if(!(result = mysql_db_query_ret(
			     mysql_conn->db_conn, query, 0))) {
			xfree(query);
			return SLURM_ERROR;
		}
		
		xfree(query);
		row = mysql_fetch_row(result);
		if(row) {
			last_hour = atoi(row[UPDATE_HOUR]);
			last_day = atoi(row[UPDATE_DAY]);
			last_month = atoi(row[UPDATE_MONTH]);		
			mysql_free_result(result);
		} else {
			time_t now = time(NULL);
			/* If we don't have any events like adding a
			 * cluster this will not work correctly, so we
			 * will insert now as a starting point.
			 */
			query = xstrdup_printf(
				"set @PS = %d;"
				"select @PS := period_start from %s limit 1;"
				"insert into %s "
				"(hourly_rollup, daily_rollup, monthly_rollup) "
				"values (@PS, @PS, @PS);",
				now, event_table, last_ran_table);
			
			debug3("%d(%d) query\n%s", mysql_conn->conn, 
			       __LINE__, query);
			mysql_free_result(result);
			if(!(result = mysql_db_query_ret(
				     mysql_conn->db_conn, query, 0))) {
				xfree(query);
				return SLURM_ERROR;
			}
			xfree(query);
			row = mysql_fetch_row(result);
			if(!row) {
				debug("No clusters have been added "
				      "not doing rollup");
				mysql_free_result(result);
				return SLURM_SUCCESS;
			}
			
			last_hour = last_day = last_month = atoi(row[0]);
			mysql_free_result(result);
		}
	}
	
	/* test month gap */
/* 	last_hour = 1212299999; */
/* 	last_day = 1212217200; */
/* 	last_month = 1212217200; */
/* 	my_time = 1212307200; */

/* 	last_hour = 1211475599; */
/* 	last_day = 1211475599; */
/* 	last_month = 1211475599; */

//	last_hour = 1211403599;
	//	last_hour = 1206946800;
//	last_day = 1207033199;
//	last_day = 1197033199;
//	last_month = 1204358399;

	if(!localtime_r(&last_hour, &start_tm)) {
		error("Couldn't get localtime from hour start %d", last_hour);
		return SLURM_ERROR;
	}
	
	if(!localtime_r(&my_time, &end_tm)) {
		error("Couldn't get localtime from hour end %d", my_time);
		return SLURM_ERROR;
	}

	/* below and anywhere in a rollup plugin when dealing with
	 * epoch times we need to set the tm_isdst = -1 so we don't
	 * have to worry about the time changes.  Not setting it to -1
	 * will cause problems in the day and month with the date change.
	 */

	start_tm.tm_sec = 0;
	start_tm.tm_min = 0;
	start_tm.tm_isdst = -1;
	start_time = mktime(&start_tm);
	end_tm.tm_sec = 0;
	end_tm.tm_min = 0;
	end_tm.tm_isdst = -1;
	end_time = mktime(&end_tm);

/* 	info("hour start %s", ctime(&start_time)); */
/* 	info("hour end %s", ctime(&end_time)); */
/* 	info("diff is %d", end_time-start_time); */

	if(end_time-start_time > 0) {
		START_TIMER;
		if((rc = mysql_hourly_rollup(mysql_conn, start_time, end_time)) 
		   != SLURM_SUCCESS)
			return rc;
		END_TIMER2("hourly_rollup");
		query = xstrdup_printf("update %s set hourly_rollup=%d",
				       last_ran_table, end_time);
	} else {
		debug2("no need to run this hour %d <= %d", 
		       end_time, start_time);
	}

	if(!localtime_r(&last_day, &start_tm)) {
		error("Couldn't get localtime from day %d", last_day);
		return SLURM_ERROR;
	}
	start_tm.tm_sec = 0;
	start_tm.tm_min = 0;
	start_tm.tm_hour = 0;
	start_tm.tm_isdst = -1;
	start_time = mktime(&start_tm);
	end_tm.tm_hour = 0;
	end_tm.tm_isdst = -1;
	end_time = mktime(&end_tm);

/* 	info("day start %s", ctime(&start_time)); */
/* 	info("day end %s", ctime(&end_time)); */
/* 	info("diff is %d", end_time-start_time); */

	if(end_time-start_time > 0) {
		START_TIMER;
		if((rc = mysql_daily_rollup(mysql_conn, start_time, end_time)) 
		   != SLURM_SUCCESS)
			return rc;
		END_TIMER2("daily_rollup");
		if(query) 
			xstrfmtcat(query, ", daily_rollup=%d", end_time);
		else 
			query = xstrdup_printf("update %s set daily_rollup=%d",
					       last_ran_table, end_time);
	} else {
		debug2("no need to run this day %d <= %d",
		       end_time, start_time);
	}

	if(!localtime_r(&last_month, &start_tm)) {
		error("Couldn't get localtime from month %d", last_month);
		return SLURM_ERROR;
	}

	start_tm.tm_sec = 0;
	start_tm.tm_min = 0;
	start_tm.tm_hour = 0;
	start_tm.tm_mday = 1;
	start_tm.tm_isdst = -1;
	start_time = mktime(&start_tm);
	end_time = mktime(&end_tm);

	end_tm.tm_sec = 0;
	end_tm.tm_min = 0;
	end_tm.tm_hour = 0;
	end_tm.tm_mday = 1;
	end_tm.tm_isdst = -1;
	end_time = mktime(&end_tm);

/* 	info("month start %s", ctime(&start_time)); */
/* 	info("month end %s", ctime(&end_time)); */
/* 	info("diff is %d", end_time-start_time); */

	if(end_time-start_time > 0) {
		START_TIMER;
		if((rc = mysql_monthly_rollup(
			    mysql_conn, start_time, end_time)) != SLURM_SUCCESS)
			return rc;
		END_TIMER2("monthly_rollup");

		if(query) 
			xstrfmtcat(query, ", monthly_rollup=%d", end_time);
		else 
			query = xstrdup_printf(
				"update %s set monthly_rollup=%d",
				last_ran_table, end_time);
	} else {
		debug2("no need to run this month %d <= %d",
		       end_time, start_time);
	}
	
	if(query) {
		debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
	}
	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int clusteracct_storage_p_node_down(mysql_conn_t *mysql_conn, 
					   char *cluster,
					   struct node_record *node_ptr,
					   time_t event_time, char *reason)
{
#ifdef HAVE_MYSQL
	uint16_t cpus;
	int rc = SLURM_SUCCESS;
	char *query = NULL;
	char *my_reason;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	if(!node_ptr) {
		error("No node_ptr given!");
		return SLURM_ERROR;
	}

	if (slurmctld_conf.fast_schedule && !slurmdbd_conf)
		cpus = node_ptr->config_ptr->cpus;
	else
		cpus = node_ptr->cpus;

	if (reason)
		my_reason = reason;
	else
		my_reason = node_ptr->reason;
	
	debug2("inserting %s(%s) with %u cpus", node_ptr->name, cluster, cpus);

	query = xstrdup_printf(
		"update %s set period_end=%d where cluster='%s' "
		"and period_end=0 and node_name='%s';",
		event_table, event_time, cluster, node_ptr->name);
	/* If you are clean-restarting the controller over and over again you
	 * could get records that are duplicates in the database.  If
	 * this is the case we will zero out the period_end we are
	 * just filled in.  This will cause the last time to be erased
	 * from the last restart, but if you are restarting things
	 * this often the pervious one didn't mean anything anyway.
	 * This way we only get one for the last time we let it run.
	 */
	xstrfmtcat(query,
		   "insert into %s "
		   "(node_name, cluster, cpu_count, period_start, reason) "
		   "values ('%s', '%s', %u, %d, \"%s\") on duplicate key "
		   "update period_end=0;",
		   event_table, node_ptr->name, cluster, 
		   cpus, event_time, my_reason);
	debug4("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);

	return rc;
#else
	return SLURM_ERROR;
#endif
}
extern int clusteracct_storage_p_node_up(mysql_conn_t *mysql_conn, 
					 char *cluster,
					 struct node_record *node_ptr,
					 time_t event_time)
{
#ifdef HAVE_MYSQL
	char* query;
	int rc = SLURM_SUCCESS;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	query = xstrdup_printf(
		"update %s set period_end=%d where cluster='%s' "
		"and period_end=0 and node_name='%s';",
		event_table, event_time, cluster, node_ptr->name);
	debug4("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int clusteracct_storage_p_register_ctld(char *cluster,
					       uint16_t port)
{
	return SLURM_SUCCESS;
}

extern int clusteracct_storage_p_cluster_procs(mysql_conn_t *mysql_conn, 
					       char *cluster,
					       uint32_t procs,
					       time_t event_time)
{
#ifdef HAVE_MYSQL
	char* query;
	int rc = SLURM_SUCCESS;
	int first = 0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;

 	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	/* Record the processor count */
	query = xstrdup_printf(
		"select cpu_count from %s where cluster='%s' "
		"and period_end=0 and node_name='' limit 1",
		event_table, cluster);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	/* we only are checking the first one here */
	if(!(row = mysql_fetch_row(result))) {
		debug("We don't have an entry for this machine %s "
		      "most likely a first time running.", cluster);

		/* Get all nodes in a down state and jobs pending or running.
		 * This is for the first time a cluster registers
		 *
		 * This only happens here when calling the plugin directly.  If
		 * calling this plugin throught the slurmdbd we do this in
		 * acct_storage_p_modify_clusters.
		 */
		if(!slurmdbd_conf) {
			/* We will return ACCOUNTING_FIRST_REG so this
			   is taken care of since the message thread
			   may not be up when we run this in the controller.
			*/
			first = 1;
		}
		goto add_it;
	}

	if(atoi(row[0]) == procs) {
		debug3("we have the same procs as before no need to "
		       "update the database.");
		goto end_it;
	}
	debug("%s has changed from %s cpus to %u", cluster, row[0], procs);   

	query = xstrdup_printf(
		"update %s set period_end=%d where cluster='%s' "
		"and period_end=0 and node_name=''",
		event_table, event_time, cluster);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	if(rc != SLURM_SUCCESS)
		goto end_it;
add_it:
	query = xstrdup_printf(
		"insert into %s (cluster, cpu_count, period_start, reason) "
		"values ('%s', %u, %d, 'Cluster processor count')",
		event_table, cluster, procs, event_time);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
end_it:
	mysql_free_result(result);
	if(first && rc == SLURM_SUCCESS)
		rc = ACCOUNTING_FIRST_REG;

	return rc;
#else
	return SLURM_ERROR;
#endif
}

extern int clusteracct_storage_p_get_usage(
	mysql_conn_t *mysql_conn, uid_t uid,
	acct_cluster_rec_t *cluster_rec, time_t start, time_t end)
{
#ifdef HAVE_MYSQL
	int rc = SLURM_SUCCESS;
	int i=0;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *tmp = NULL;
	char *my_usage_table = cluster_day_table;
	time_t my_time = time(NULL);
	struct tm start_tm;
	struct tm end_tm;
	char *query = NULL;
	char *cluster_req_inx[] = {
		"alloc_cpu_secs",
		"down_cpu_secs",
		"idle_cpu_secs",
		"resv_cpu_secs",
		"over_cpu_secs",
		"cpu_count",
		"period_start"
	};
	
	enum {
		CLUSTER_ACPU,
		CLUSTER_DCPU,
		CLUSTER_ICPU,
		CLUSTER_RCPU,
		CLUSTER_OCPU,
		CLUSTER_CPU_COUNT,
		CLUSTER_START,
		CLUSTER_COUNT
	};

	if(!cluster_rec->name) {
		error("We need a cluster name to set data for");
		return SLURM_ERROR;
	}

	/* Default is going to be the last day */
	if(!end) {
		if(!localtime_r(&my_time, &end_tm)) {
			error("Couldn't get localtime from end %d",
			      my_time);
			return SLURM_ERROR;
		}
		end_tm.tm_hour = 0;
		end = mktime(&end_tm);		
	} else {
		if(!localtime_r(&end, &end_tm)) {
			error("Couldn't get localtime from user end %d",
			      my_time);
			return SLURM_ERROR;
		}
	}
	end_tm.tm_sec = 0;
	end_tm.tm_min = 0;
	end_tm.tm_isdst = -1;
	end = mktime(&end_tm);		

	if(!start) {
		if(!localtime_r(&my_time, &start_tm)) {
			error("Couldn't get localtime from start %d",
			      my_time);
			return SLURM_ERROR;
		}
		start_tm.tm_hour = 0;
		start_tm.tm_mday--;
		start = mktime(&start_tm);		
	} else {
		if(!localtime_r(&start, &start_tm)) {
			error("Couldn't get localtime from user start %d",
			      my_time);
			return SLURM_ERROR;
		}
	}
	start_tm.tm_sec = 0;
	start_tm.tm_min = 0;
	start_tm.tm_isdst = -1;
	start = mktime(&start_tm);		

	if(end-start < 3600) {
		end = start + 3600;
		if(!localtime_r(&end, &end_tm)) {
			error("2 Couldn't get localtime from user end %d",
			      my_time);
			return SLURM_ERROR;
		}
	}
	/* check to see if we are off day boundaries or on month
	 * boundaries other wise use the day table.
	 */
	if(start_tm.tm_hour || end_tm.tm_hour || (end-start < 86400)) 
		my_usage_table = cluster_hour_table;
	else if(start_tm.tm_mday == 0 && end_tm.tm_mday == 0 
		&& (end-start > 86400))
		my_usage_table = cluster_month_table;

	xfree(tmp);
	i=0;
	xstrfmtcat(tmp, "%s", cluster_req_inx[i]);
	for(i=1; i<CLUSTER_COUNT; i++) {
		xstrfmtcat(tmp, ", %s", cluster_req_inx[i]);
	}

	query = xstrdup_printf(
		"select %s from %s where (period_start < %d "
		"&& period_start >= %d) and cluster='%s'",
		tmp, my_usage_table, end, start, cluster_rec->name);

	xfree(tmp);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	if(!(result = mysql_db_query_ret(
		     mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	if(!cluster_rec->accounting_list)
		cluster_rec->accounting_list =
			list_create(destroy_cluster_accounting_rec);
	
	while((row = mysql_fetch_row(result))) {
		cluster_accounting_rec_t *accounting_rec =
			xmalloc(sizeof(cluster_accounting_rec_t));
		accounting_rec->alloc_secs = atoll(row[CLUSTER_ACPU]);
		accounting_rec->down_secs = atoll(row[CLUSTER_DCPU]);
		accounting_rec->idle_secs = atoll(row[CLUSTER_ICPU]);
		accounting_rec->over_secs = atoll(row[CLUSTER_OCPU]);
		accounting_rec->resv_secs = atoll(row[CLUSTER_RCPU]);
		accounting_rec->cpu_count = atoi(row[CLUSTER_CPU_COUNT]);
		accounting_rec->period_start = atoi(row[CLUSTER_START]);
		list_append(cluster_rec->accounting_list, accounting_rec);
	}
	mysql_free_result(result);

	return rc;
#else
	return SLURM_ERROR;
#endif
}

/* 
 * load into the storage the start of a job
 */
extern int jobacct_storage_p_job_start(mysql_conn_t *mysql_conn, 
				       struct job_record *job_ptr)
{
#ifdef HAVE_MYSQL
	int	rc=SLURM_SUCCESS;
	char	*jname = NULL, *nodes = NULL;
	long	priority;
	int track_steps = 0;
	char *block_id = NULL;
	char *query = NULL;
	int reinit = 0;

	if (!job_ptr->details || !job_ptr->details->submit_time) {
		error("jobacct_storage_p_job_start: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;
	
	debug2("mysql_jobacct_job_start() called");
	priority = (job_ptr->priority == NO_VAL) ?
		-1L : (long) job_ptr->priority;

	if (job_ptr->name && job_ptr->name[0]) {
		int i;
		jname = xmalloc(strlen(job_ptr->name) + 1);
		for (i=0; job_ptr->name[i]; i++) {
			if (isalnum(job_ptr->name[i]))
				jname[i] = job_ptr->name[i];
			else
				jname[i] = '_';
		}
	} else {
		jname = xstrdup("allocation");
		track_steps = 1;
	}

	if (job_ptr->nodes && job_ptr->nodes[0])
		nodes = job_ptr->nodes;
	else
		nodes = "None assigned";

	if(job_ptr->batch_flag)
		track_steps = 1;

	if(slurmdbd_conf) {
		block_id = xstrdup(job_ptr->comment);
	} else {
		select_g_get_jobinfo(job_ptr->select_jobinfo, 
				     SELECT_DATA_BLOCK_ID, 
				     &block_id);
	}

	job_ptr->requid = -1; /* force to -1 for sacct to know this
			       * hasn't been set yet */
	
	/* We need to put a 0 for 'end' incase of funky job state
	 * files from a hot start of the controllers we call
	 * job_start on jobs we may still know about after
	 * job_flush has been called so we need to restart
	 * them by zeroing out the end.
	 */
	if(!job_ptr->db_index) {
		query = xstrdup_printf(
			"insert into %s "
			"(jobid, associd, uid, gid, nodelist, ",
			job_table);

		if(job_ptr->account) 
			xstrcat(query, "account, ");
		if(job_ptr->partition) 
			xstrcat(query, "partition, ");
		if(block_id) 
			xstrcat(query, "blockid, ");
		
		xstrfmtcat(query, 
			   "eligible, submit, start, name, track_steps, "
			   "state, priority, req_cpus, alloc_cpus) "
			   "values (%u, %u, %u, %u, '%s', ",
			   job_ptr->job_id, job_ptr->assoc_id,
			   job_ptr->user_id, job_ptr->group_id, nodes);
		
		if(job_ptr->account) 
			xstrfmtcat(query, "'%s', ", job_ptr->account);
		if(job_ptr->partition) 
			xstrfmtcat(query, "'%s', ", job_ptr->partition);
		if(block_id) 
			xstrfmtcat(query, "'%s', ", block_id);
		
		xstrfmtcat(query, 
			   "%d, %d, %d, '%s', %u, %u, %u, %u, %u) "
			   "on duplicate key update "
			   "id=LAST_INSERT_ID(id), state=%u, associd=%u",
			   (int)job_ptr->details->begin_time,
			   (int)job_ptr->details->submit_time,
			   (int)job_ptr->start_time,
			   jname, track_steps,
			   job_ptr->job_state & (~JOB_COMPLETING),
			   priority, job_ptr->num_procs,
			   job_ptr->total_procs, 
			   job_ptr->job_state & (~JOB_COMPLETING),
			   job_ptr->assoc_id);

		if(job_ptr->account) 
			xstrfmtcat(query, ", account='%s'", job_ptr->account);
		if(job_ptr->partition) 
			xstrfmtcat(query, ", partition='%s'",
				   job_ptr->partition);
		if(block_id)
			xstrfmtcat(query, ", blockid='%s'", block_id);
		
		debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	try_again:
		if(!(job_ptr->db_index = mysql_insert_ret_id(
			     mysql_conn->db_conn, query))) {
			if(!reinit) {
				error("It looks like the storage has gone "
				      "away trying to reconnect");
				mysql_close_db_connection(
					&mysql_conn->db_conn);
				mysql_get_db_connection(
					&mysql_conn->db_conn,
					mysql_db_name, mysql_db_info);
				reinit = 1;
				goto try_again;
			} else
				rc = SLURM_ERROR;
		}
	} else {
		query = xstrdup_printf("update %s set nodelist='%s', ", 
				       job_table, nodes);

		if(job_ptr->account) 
			xstrfmtcat(query, "account='%s', ",
				   job_ptr->account);
		if(job_ptr->partition) 
			xstrfmtcat(query, "partition='%s', ",
				   job_ptr->partition);
		if(block_id)
			xstrfmtcat(query, "blockid='%s', ", block_id);

		xstrfmtcat(query, "start=%d, name='%s', state=%u, "
			   "alloc_cpus=%u, associd=%d where id=%d",
			   (int)job_ptr->start_time,
			   jname, job_ptr->job_state & (~JOB_COMPLETING),
			   job_ptr->total_procs, job_ptr->assoc_id,
			   job_ptr->db_index);
		debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
		rc = mysql_db_query(mysql_conn->db_conn, query);
	}

	xfree(block_id);
	xfree(jname);

	xfree(query);

	return rc;
#else
	return SLURM_ERROR;
#endif
}

/* 
 * load into the storage the end of a job
 */
extern int jobacct_storage_p_job_complete(mysql_conn_t *mysql_conn, 
					  struct job_record *job_ptr)
{
#ifdef HAVE_MYSQL
	char *query = NULL, *nodes = NULL;
	int rc=SLURM_SUCCESS;

	if (!job_ptr->db_index 
	    && (!job_ptr->details || !job_ptr->details->submit_time)) {
		error("jobacct_storage_p_job_complete: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;
	debug2("mysql_jobacct_job_complete() called");
	
	/* If we get an error with this just fall through to avoid an
	 * infinite loop
	 */
	if (job_ptr->end_time == 0) {
		debug("mysql_jobacct: job %u never started", job_ptr->job_id);
		return SLURM_SUCCESS;
	}	
	
	if (job_ptr->nodes && job_ptr->nodes[0])
		nodes = job_ptr->nodes;
	else
		nodes = "None assigned";

	if(!job_ptr->db_index) {
		if(!(job_ptr->db_index =
		     _get_db_index(mysql_conn->db_conn,
				   job_ptr->details->submit_time,
				   job_ptr->job_id,
				   job_ptr->assoc_id))) {
			/* If we get an error with this just fall
			 * through to avoid an infinite loop
			 */
			if(jobacct_storage_p_job_start(mysql_conn, job_ptr)
			   == SLURM_ERROR) {
				error("couldn't add job %u at job completion",
				      job_ptr->job_id);
				return SLURM_SUCCESS;
			}
		}
	}

	query = xstrdup_printf("update %s set start=%u, end=%u, state=%d, "
			       "nodelist='%s', comp_code=%u, "
			       "kill_requid=%u where id=%u",
			       job_table, (int)job_ptr->start_time,
			       (int)job_ptr->end_time, 
			       job_ptr->job_state & (~JOB_COMPLETING),
			       nodes, job_ptr->exit_code,
			       job_ptr->requid, job_ptr->db_index);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	
	return rc;
#else
	return SLURM_ERROR;
#endif
}

/* 
 * load into the storage the start of a job step
 */
extern int jobacct_storage_p_step_start(mysql_conn_t *mysql_conn, 
					struct step_record *step_ptr)
{
#ifdef HAVE_MYSQL
	int cpus = 0;
	int rc=SLURM_SUCCESS;
	char node_list[BUFFER_SIZE];
#ifdef HAVE_BG
	char *ionodes = NULL;
#endif
	char *query = NULL;
	
	if (!step_ptr->job_ptr->db_index 
	    && (!step_ptr->job_ptr->details
		|| !step_ptr->job_ptr->details->submit_time)) {
		error("jobacct_storage_p_step_start: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;
	if(slurmdbd_conf) {
		cpus = step_ptr->job_ptr->total_procs;
		snprintf(node_list, BUFFER_SIZE, "%s",
			 step_ptr->job_ptr->nodes);
	} else {
#ifdef HAVE_BG
		cpus = step_ptr->job_ptr->num_procs;
		select_g_get_jobinfo(step_ptr->job_ptr->select_jobinfo, 
				     SELECT_DATA_IONODES, 
				     &ionodes);
		if(ionodes) {
			snprintf(node_list, BUFFER_SIZE, 
				 "%s[%s]", step_ptr->job_ptr->nodes, ionodes);
			xfree(ionodes);
		} else
			snprintf(node_list, BUFFER_SIZE, "%s",
				 step_ptr->job_ptr->nodes);
		
#else
		if(!step_ptr->step_layout || !step_ptr->step_layout->task_cnt) {
			cpus = step_ptr->job_ptr->total_procs;
			snprintf(node_list, BUFFER_SIZE, "%s",
				 step_ptr->job_ptr->nodes);
		} else {
			cpus = step_ptr->step_layout->task_cnt;
			snprintf(node_list, BUFFER_SIZE, "%s", 
				 step_ptr->step_layout->node_list);
		}
#endif
	}

	step_ptr->job_ptr->requid = -1; /* force to -1 for sacct to know this
					 * hasn't been set yet  */

	if(!step_ptr->job_ptr->db_index) {
		if(!(step_ptr->job_ptr->db_index = 
		     _get_db_index(mysql_conn->db_conn,
				   step_ptr->job_ptr->details->submit_time,
				   step_ptr->job_ptr->job_id,
				   step_ptr->job_ptr->assoc_id))) {
			/* If we get an error with this just fall
			 * through to avoid an infinite loop
			 */
			if(jobacct_storage_p_job_start(mysql_conn,
						       step_ptr->job_ptr)
			   == SLURM_ERROR) {
				error("couldn't add job %u at step start",
				      step_ptr->job_ptr->job_id);
				return SLURM_SUCCESS;
			}
		}
	}
	/* we want to print a -1 for the requid so leave it a
	   %d */
	query = xstrdup_printf(
		"insert into %s (id, stepid, start, name, state, "
		"cpus, nodelist) "
		"values (%d, %u, %d, '%s', %d, %u, '%s') "
		"on duplicate key update cpus=%u, end=0, state=%u",
		step_table, step_ptr->job_ptr->db_index,
		step_ptr->step_id, 
		(int)step_ptr->start_time, step_ptr->name,
		JOB_RUNNING, cpus, node_list, cpus, JOB_RUNNING);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);

	return rc;
#else
	return SLURM_ERROR;
#endif
}

/* 
 * load into the storage the end of a job step
 */
extern int jobacct_storage_p_step_complete(mysql_conn_t *mysql_conn, 
					   struct step_record *step_ptr)
{
#ifdef HAVE_MYSQL
	time_t now;
	int elapsed;
	int comp_status;
	int cpus = 0;
	struct jobacctinfo *jobacct = (struct jobacctinfo *)step_ptr->jobacct;
	struct jobacctinfo dummy_jobacct;
	float ave_vsize = 0, ave_rss = 0, ave_pages = 0;
	float ave_cpu = 0, ave_cpu2 = 0;
	char *query = NULL;
	int rc =SLURM_SUCCESS;
	uint32_t exit_code = 0;

	if (!step_ptr->job_ptr->db_index 
	    && (!step_ptr->job_ptr->details
		|| !step_ptr->job_ptr->details->submit_time)) {
		error("jobacct_storage_p_step_complete: "
		      "Not inputing this job, it has no submit time.");
		return SLURM_ERROR;
	}

	if (jobacct == NULL) {
		/* JobAcctGather=jobacct_gather/none, no data to process */
		bzero(&dummy_jobacct, sizeof(dummy_jobacct));
		jobacct = &dummy_jobacct;
	}

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	if(slurmdbd_conf) {
		now = step_ptr->job_ptr->end_time;
		cpus = step_ptr->job_ptr->total_procs;

	} else {
		now = time(NULL);
#ifdef HAVE_BG
		cpus = step_ptr->job_ptr->num_procs;
		
#else
		if(!step_ptr->step_layout || !step_ptr->step_layout->task_cnt)
			cpus = step_ptr->job_ptr->total_procs;
		else 
			cpus = step_ptr->step_layout->task_cnt;
#endif
	}
	
	if ((elapsed=now-step_ptr->start_time)<0)
		elapsed=0;	/* For *very* short jobs, if clock is wrong */
	
	exit_code = step_ptr->exit_code;
	if (exit_code == NO_VAL) {
		comp_status = JOB_CANCELLED;
		exit_code = 0;
	} else if (exit_code)
		comp_status = JOB_FAILED;
	else
		comp_status = JOB_COMPLETE;
       
	/* figure out the ave of the totals sent */
	if(cpus > 0) {
		ave_vsize = jobacct->tot_vsize;
		ave_vsize /= cpus;
		ave_rss = jobacct->tot_rss;
		ave_rss /= cpus;
		ave_pages = jobacct->tot_pages;
		ave_pages /= cpus;
		ave_cpu = jobacct->tot_cpu;
		ave_cpu /= cpus;	
		ave_cpu /= 100;
	}
 
	if(jobacct->min_cpu != NO_VAL) {
		ave_cpu2 = jobacct->min_cpu;
		ave_cpu2 /= 100;
	}

	if(!step_ptr->job_ptr->db_index) {
		if(!(step_ptr->job_ptr->db_index = 
		     _get_db_index(mysql_conn->db_conn,
				   step_ptr->job_ptr->details->submit_time,
				   step_ptr->job_ptr->job_id,
				   step_ptr->job_ptr->assoc_id))) {
			/* If we get an error with this just fall
			 * through to avoid an infinite loop
			 */
			if(jobacct_storage_p_job_start(mysql_conn,
						       step_ptr->job_ptr)
			   == SLURM_ERROR) {
				error("couldn't add job %u "
				      "at step completion",
				      step_ptr->job_ptr->job_id);
				return SLURM_SUCCESS;
			}
		}
	}

	query = xstrdup_printf(
		"update %s set end=%d, state=%d, "
		"kill_requid=%u, comp_code=%u, "
		"user_sec=%ld, user_usec=%ld, "
		"sys_sec=%ld, sys_usec=%ld, "
		"max_vsize=%u, max_vsize_task=%u, "
		"max_vsize_node=%u, ave_vsize=%.2f, "
		"max_rss=%u, max_rss_task=%u, "
		"max_rss_node=%u, ave_rss=%.2f, "
		"max_pages=%u, max_pages_task=%u, "
		"max_pages_node=%u, ave_pages=%.2f, "
		"min_cpu=%.2f, min_cpu_task=%u, "
		"min_cpu_node=%u, ave_cpu=%.2f "
		"where id=%u and stepid=%u",
		step_table, (int)now,
		comp_status,
		step_ptr->job_ptr->requid, 
		exit_code,
		/* user seconds */
		jobacct->user_cpu_sec,	
		/* user microseconds */
		jobacct->user_cpu_usec,
		/* system seconds */
		jobacct->sys_cpu_sec,
		/* system microsecs */
		jobacct->sys_cpu_usec,
		jobacct->max_vsize,	/* max vsize */
		jobacct->max_vsize_id.taskid,	/* max vsize task */
		jobacct->max_vsize_id.nodeid,	/* max vsize node */
		ave_vsize,	/* ave vsize */
		jobacct->max_rss,	/* max vsize */
		jobacct->max_rss_id.taskid,	/* max rss task */
		jobacct->max_rss_id.nodeid,	/* max rss node */
		ave_rss,	/* ave rss */
		jobacct->max_pages,	/* max pages */
		jobacct->max_pages_id.taskid,	/* max pages task */
		jobacct->max_pages_id.nodeid,	/* max pages node */
		ave_pages,	/* ave pages */
		ave_cpu2,	/* min cpu */
		jobacct->min_cpu_id.taskid,	/* min cpu task */
		jobacct->min_cpu_id.nodeid,	/* min cpu node */
		ave_cpu,	/* ave cpu */
		step_ptr->job_ptr->db_index, step_ptr->step_id);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
	rc = mysql_db_query(mysql_conn->db_conn, query);
	xfree(query);
	 
	return rc;
#else
	return SLURM_ERROR;
#endif
}

/* 
 * load into the storage a suspention of a job
 */
extern int jobacct_storage_p_suspend(mysql_conn_t *mysql_conn, 
				     struct job_record *job_ptr)
{
#ifdef HAVE_MYSQL
	char *query = NULL;
	int rc = SLURM_SUCCESS;
	bool suspended = false;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;
	if(!job_ptr->db_index) {
		if(!(job_ptr->db_index =
		     _get_db_index(mysql_conn->db_conn,
				   job_ptr->details->submit_time,
				   job_ptr->job_id,
				   job_ptr->assoc_id))) {
			/* If we get an error with this just fall
			 * through to avoid an infinite loop
			 */
			if(jobacct_storage_p_job_start(mysql_conn, job_ptr)
			   == SLURM_ERROR) {
				error("couldn't suspend job %u",
				      job_ptr->job_id);
				return SLURM_SUCCESS;
			}
		}
	}

	if (job_ptr->job_state == JOB_SUSPENDED)
		suspended = true;

	xstrfmtcat(query,
		   "update %s set suspended=%d-suspended, state=%d "
		   "where id=%u;",
		   job_table, (int)job_ptr->suspend_time, 
		   job_ptr->job_state & (~JOB_COMPLETING),
		   job_ptr->db_index);
	if(suspended)
		xstrfmtcat(query,
			   "insert into %s (id, associd, start, end) "
			   "values (%u, %u, %d, 0);",
			   suspend_table, job_ptr->db_index, job_ptr->assoc_id,
			   (int)job_ptr->suspend_time);
	else
		xstrfmtcat(query,
			   "update %s set end=%d where id=%u && end=0;",
			   suspend_table, (int)job_ptr->suspend_time, 
			   job_ptr->db_index);
	debug3("%d(%d) query\n%s", mysql_conn->conn, __LINE__, query);
				
	rc = mysql_db_query(mysql_conn->db_conn, query);

	xfree(query);
	if(rc != SLURM_ERROR) {
		xstrfmtcat(query,
			   "update %s set suspended=%u-suspended, "
			   "state=%d where id=%u and end=0",
			   step_table, (int)job_ptr->suspend_time, 
			   job_ptr->job_state, job_ptr->db_index);
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
	}
	
	return rc;
#else
	return SLURM_ERROR;
#endif
}

/* 
 * get info from the storage 
 * returns List of job_rec_t *
 * note List needs to be freed when called
 */
extern List jobacct_storage_p_get_jobs(mysql_conn_t *mysql_conn, uid_t uid, 
				       List selected_steps,
				       List selected_parts,
				       sacct_parameters_t *params)
{
	List job_list = NULL;
#ifdef HAVE_MYSQL
	acct_job_cond_t job_cond;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;
	memset(&job_cond, 0, sizeof(acct_job_cond_t));

	job_cond.acct_list = selected_steps;
	job_cond.step_list = selected_steps;
	job_cond.partition_list = selected_parts;
	job_cond.cluster_list = params->opt_cluster_list;

	if (params->opt_uid >=0) {
		char *temp = xstrdup_printf("%u", params->opt_uid);
		job_cond.userid_list = list_create(NULL);
		list_append(job_cond.userid_list, temp);
	}	

	if (params->opt_gid >=0) {
		char *temp = xstrdup_printf("%u", params->opt_gid);
		job_cond.groupid_list = list_create(NULL);
		list_append(job_cond.groupid_list, temp);
	}	

	job_list = mysql_jobacct_process_get_jobs(mysql_conn, uid, &job_cond);

	if(job_cond.userid_list)
		list_destroy(job_cond.userid_list);
	if(job_cond.groupid_list)
		list_destroy(job_cond.groupid_list);
		
#endif
	return job_list;
}

/* 
 * get info from the storage 
 * returns List of job_rec_t *
 * note List needs to be freed when called
 */
extern List jobacct_storage_p_get_jobs_cond(mysql_conn_t *mysql_conn, 
					    uid_t uid, 
					    acct_job_cond_t *job_cond)
{
	List job_list = NULL;
#ifdef HAVE_MYSQL
	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return NULL;
	job_list = mysql_jobacct_process_get_jobs(mysql_conn, uid, job_cond);	
#endif
	return job_list;
}

/* 
 * expire old info from the storage 
 */
extern void jobacct_storage_p_archive(mysql_conn_t *mysql_conn, 
				      List selected_parts,
				      void *params)
{
#ifdef HAVE_MYSQL
	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return;
	mysql_jobacct_process_archive(mysql_conn,
				      selected_parts, params);
#endif
	return;
}

extern int acct_storage_p_update_shares_used(mysql_conn_t *mysql_conn, 
					     List shares_used)
{
	/* This definitely needs to be fleshed out.
	 * Go through the list of shares_used_object_t objects and store them */
	return SLURM_SUCCESS;
}

extern int acct_storage_p_flush_jobs_on_cluster(
	mysql_conn_t *mysql_conn, char *cluster, time_t event_time)
{
	int rc = SLURM_SUCCESS;
#ifdef HAVE_MYSQL
	/* put end times for a clean start */
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	char *query = NULL;
	char *id_char = NULL;
	char *suspended_char = NULL;

	if(_check_connection(mysql_conn) != SLURM_SUCCESS)
		return SLURM_ERROR;

	/* First we need to get the id's and states so we can clean up
	 * the suspend table and the step table 
	 */
	query = xstrdup_printf("select t1.id, t1.state from %s as t1, %s as t2 "
			       "where ((t2.id=t1.associd and t2.cluster='%s') "
			       "|| !t1.associd) && t1.end=0;",
			       job_table, assoc_table, cluster);
	if(!(result =
	     mysql_db_query_ret(mysql_conn->db_conn, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	while((row = mysql_fetch_row(result))) {
		int state = atoi(row[1]);
		if(state == JOB_SUSPENDED) {
			if(suspended_char) 
				xstrfmtcat(suspended_char, " || id=%s", row[0]);
			else
				xstrfmtcat(suspended_char, "id=%s", row[0]);
		}
		
		if(id_char) 
			xstrfmtcat(id_char, " || id=%s", row[0]);
		else
			xstrfmtcat(id_char, "id=%s", row[0]);
	}
	mysql_free_result(result);
	
	if(suspended_char) {
		xstrfmtcat(query,
			   "update %s set suspended=%d-suspended where %s;",
			   job_table, event_time, suspended_char);
		xstrfmtcat(query,
			   "update %s set suspended=%d-suspended where %s;",
			   step_table, event_time, suspended_char);
		xstrfmtcat(query,
			   "update %s set end=%d where (%s) && end=0;",
			   suspend_table, event_time, suspended_char);
		xfree(suspended_char);
	}
	if(id_char) {
		xstrfmtcat(query,
			   "update %s set state=%d, end=%u where %s;",
			   job_table, JOB_CANCELLED, event_time, id_char);
		xstrfmtcat(query,
			   "update %s set state=%d, end=%u where %s;",
			   step_table, JOB_CANCELLED, event_time, id_char);
		xfree(id_char);
	}
/* 	query = xstrdup_printf("update %s as t1, %s as t2 set " */
/* 			       "t1.state=%u, t1.end=%u where " */
/* 			       "t2.id=t1.associd and t2.cluster='%s' " */
/* 			       "&& t1.end=0;", */
/* 			       job_table, assoc_table, JOB_CANCELLED,  */
/* 			       event_time, cluster); */
	if(query) {
		debug3("%d(%d) query\n%s",
		       mysql_conn->conn, __LINE__, query);
		
		rc = mysql_db_query(mysql_conn->db_conn, query);
		xfree(query);
	}
#endif

	return rc;
}
