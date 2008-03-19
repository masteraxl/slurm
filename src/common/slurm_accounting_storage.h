/*****************************************************************************\
 *  slurm_accounting_storage.h - Define accounting storage plugin functions.
 *****************************************************************************
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Copyright (C) 2008 Lawrence Livermore National Security.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
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

#ifndef _SLURM_ACCOUNTING_STORAGE_H 
#define _SLURM_ACCOUNTING_STORAGE_H

#include "src/common/list.h"
#include "src/slurmctld/slurmctld.h"
#include <slurm/slurm.h>
#include <slurm/slurm_errno.h>
#include <sys/types.h>
#include <pwd.h>

typedef enum {
	ACCT_USAGE_NOTSET,
	ACCT_USAGE_HOUR,
	ACCT_USAGE_DAY,
	ACCT_USAGE_MONTH
} acct_usage_type_t;

typedef enum {
	ACCT_ADMIN_NOTSET,
	ACCT_ADMIN_NONE,
	ACCT_ADMIN_OPERATOR,
	ACCT_ADMIN_SUPER_USER
} acct_admin_level_t;

typedef enum {
	ACCT_QOS_NOTSET,
	ACCT_QOS_NORMAL,
	ACCT_QOS_EXPEDITE,
	ACCT_QOS_STANDBY,
	ACCT_QOS_EXEMPT	
} acct_qos_level_t;

typedef struct {
	List acct_list; /* list of char * */
	List description_list; /* list of char * */
	acct_qos_level_t qos;	
	List organization_list; /* list of char * */
	uint16_t with_assocs; 
} acct_account_cond_t;

typedef struct {
	List assoc_list; /* list of acct_association_rec_t *'s */
	List coordinators; /* list of char *'s */
	char *description;
	acct_qos_level_t qos;
	char *name;
	char *organization;
} acct_account_rec_t;

typedef struct {
	uint32_t alloc_secs; /* number of cpu seconds allocated */
	time_t period_start; 
} acct_accounting_rec_t;

typedef struct {
	List acct_list; /* list of char * */
	List cluster_list; /* list of char * */
	List id_list; /* list of char */
	List partition_list; /* list of char * */
	char *parent_acct; /* name of parent account */
	List user_list; /* list of char * */
} acct_association_cond_t;

typedef struct {
	List accounting_list; 	/* list of acct_accounting_rec_t *'s */
	char *acct;		/* account/project associated to association */
	char *cluster;		/* cluster associated to association */
	uint32_t fairshare;	/* fairshare number */
	uint32_t id;		/* id identifing a combination of
				 * user-account-cluster(-partition) */
	uint32_t max_cpu_secs_per_job; /* max number of cpu seconds this 
					   * association can have per job */
	uint32_t max_jobs;	/* max number of jobs this association can run
				 * at one time */
	uint32_t max_nodes_per_job; /* max number of nodes this
				     * association can allocate per job */
	uint32_t max_wall_duration_per_job; /* longest time this
					     * association can run a job */
	char *parent_acct;	/* name of parent account */
	char *partition;	/* optional partition in a cluster 
				 * associated to association */
	uint32_t uid;		/* user ID */
	char *user;		/* user associated to association */
} acct_association_rec_t;

typedef struct {
	List cluster_list; /* list of char * */
} acct_cluster_cond_t;

typedef struct {
	List accounting_list; /* list of cluster_accounting_rec_t *'s */
	char *control_host;
	uint32_t control_port;
	uint32_t default_fairshare;	/* fairshare number */
	uint32_t default_max_cpu_secs_per_job; /* max number of cpu seconds this 
					* association can have per job */
	uint32_t default_max_jobs;	/* max number of jobs this association can run
				 * at one time */
	uint32_t default_max_nodes_per_job; /* max number of nodes this
				     * association can allocate per job */
	uint32_t default_max_wall_duration_per_job; /* longest time this
					     * association can run a job */
	char *name;

} acct_cluster_rec_t;

typedef struct {
	char *acct_name;
	uint16_t sub_acct;
} acct_coord_rec_t;

typedef struct {
	acct_admin_level_t admin_level;
	List def_acct_list; /* list of char * */
	acct_qos_level_t qos;	
	List user_list; /* list of char * */
	uint16_t with_assocs; 
} acct_user_cond_t;

typedef struct {
	acct_admin_level_t admin_level;
	List assoc_list; /* list of acct_association_rec_t *'s */
	List coord_accts; /* list of acct_coord_rec_t *'s */
	char *default_acct;
	acct_qos_level_t qos;
	char *name;
	uint32_t uid;
} acct_user_rec_t;

typedef struct {
	uint32_t alloc_secs; /* number of cpu seconds allocated */
	uint32_t cpu_count; /* number of cpus during time period */
	uint32_t down_secs; /* number of cpu seconds down */
	uint32_t idle_secs; /* number of cpu seconds idle */
	time_t period_start; /* when this record was started */
	uint32_t resv_secs; /* number of cpu seconds reserved */	
} cluster_accounting_rec_t;

extern void destroy_acct_user_rec(void *object);
extern void destroy_acct_account_rec(void *object);
extern void destroy_acct_coord_rec(void *object);
extern void destroy_cluster_accounting_rec(void *object);
extern void destroy_acct_cluster_rec(void *object);
extern void destroy_acct_accounting_rec(void *object);
extern void destroy_acct_association_rec(void *object);

extern void destroy_acct_user_cond(void *object);
extern void destroy_acct_account_cond(void *object);
extern void destroy_acct_cluster_cond(void *object);
extern void destroy_acct_association_cond(void *object);

/* pack functions */
extern void pack_acct_user_rec(void *object, Buf buffer);
extern int unpack_acct_user_rec(void **object, Buf buffer);
extern void pack_acct_account_rec(void *object, Buf buffer);
extern int unpack_acct_account_rec(void **object, Buf buffer);
extern void pack_acct_coord_rec(void *object, Buf buffer);
extern int unpack_acct_coord_rec(void **object, Buf buffer);
extern void pack_cluster_accounting_rec(void *object, Buf buffer);
extern int unpack_cluster_accounting_rec(void **object, Buf buffer);
extern void pack_acct_cluster_rec(void *object, Buf buffer);
extern int unpack_acct_cluster_rec(void **object, Buf buffer);
extern void pack_acct_accounting_rec(void *object, Buf buffer);
extern int unpack_acct_accounting_rec(void **object, Buf buffer);
extern void pack_acct_association_rec(void *object, Buf buffer);
extern int unpack_acct_association_rec(void **object, Buf buffer);

extern void pack_acct_user_cond(void *object, Buf buffer);
extern int unpack_acct_user_cond(void **object, Buf buffer);
extern void pack_acct_account_cond(void *object, Buf buffer);
extern int unpack_acct_account_cond(void **object, Buf buffer);
extern void pack_acct_cluster_cond(void *object, Buf buffer);
extern int unpack_acct_cluster_cond(void **object, Buf buffer);
extern void pack_acct_association_cond(void *object, Buf buffer);
extern int unpack_acct_association_cond(void **object, Buf buffer);

extern char *acct_qos_str(acct_qos_level_t level);
extern acct_qos_level_t str_2_acct_qos(char *level);
extern char *acct_admin_level_str(acct_admin_level_t level);
extern acct_admin_level_t str_2_acct_admin_level(char *level);

extern int slurm_acct_storage_init(char *loc); /* load the plugin */
extern int slurm_acct_storage_fini(void); /* unload the plugin */

/*
 * get a new connection to the storage unit
 * RET: pointer used to access db 
 */
extern void *acct_storage_g_get_connection(bool rollback);

/*
 * release connection to the storage unit
 * IN: void * pointer returned from acct_storage_g_get_connection()
 * RET: SLURM_SUCCESS on success SLURM_ERROR else 
 */
extern int acct_storage_g_close_connection(void **db_conn, bool commit);


/* 
 * add users to accounting system 
 * IN:  user_list List of acct_user_rec_t *
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_add_users(void *db_conn, uint32_t uid, 
				    List user_list);

/* 
 * add users as account coordinators 
 * IN:  acct name of account
 * IN:  acct_user_cond_t *user_q
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_add_coord(void *db_conn, uint32_t uid,
				    char *acct, acct_user_cond_t *user_q);


/* 
 * add accounts to accounting system 
 * IN:  account_list List of acct_account_rec_t *
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_add_accounts(void *db_conn, uint32_t uid,
				       List acct_list);

/* 
 * add clusters to accounting system 
 * IN:  cluster_list List of acct_cluster_rec_t *
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_add_clusters(void *db_conn, uint32_t uid,
				       List cluster_list);

/* 
 * add accts to accounting system 
 * IN:  association_list List of acct_association_rec_t *
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_add_associations(void *db_conn, uint32_t uid, 
					   List association_list);

/* 
 * modify existing users in the accounting system 
 * IN:  acct_user_cond_t *user_q
 * IN:  acct_user_rec_t *user
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_modify_users(void *db_conn, uint32_t uid, 
				       acct_user_cond_t *user_q,
				       acct_user_rec_t *user);

/* 
 * modify existing accounts in the accounting system 
 * IN:  acct_acct_cond_t *acct_q
 * IN:  acct_account_rec_t *acct
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_modify_accounts(void *db_conn, uint32_t uid, 
					  acct_account_cond_t *acct_q,
					  acct_account_rec_t *acct);

/* 
 * modify existing clusters in the accounting system 
 * IN:  acct_cluster_cond_t *cluster_q
 * IN:  acct_cluster_rec_t *cluster
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_modify_clusters(void *db_conn, uint32_t uid, 
					  acct_cluster_cond_t *cluster_q,
					  acct_cluster_rec_t *cluster);

/* 
 * modify existing associations in the accounting system 
 * IN:  acct_association_cond_t *assoc_q
 * IN:  acct_association_rec_t *assoc
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_modify_associations(void *db_conn, uint32_t uid, 
					      acct_association_cond_t *assoc_q,
					      acct_association_rec_t *assoc);

/* 
 * remove users from accounting system 
 * IN:  acct_user_cond_t *user_q
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_remove_users(void *db_conn, uint32_t uid, 
				       acct_user_cond_t *user_q);

/* 
 * remove users from being a coordinator of an account
 * IN: acct name of acct
 * IN: acct_user_cond_t *user_q
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_remove_coord(void *db_conn, uint32_t uid, 
				       char *acct, acct_user_cond_t *user_q);

/* 
 * remove accounts from accounting system 
 * IN:  acct_account_cond_t *acct_q
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_remove_accounts(void *db_conn, uint32_t uid, 
					  acct_account_cond_t *acct_q);

/* 
 * remove clusters from accounting system 
 * IN:  acct_cluster_cond_t *cluster_q
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_remove_clusters(void *db_conn, uint32_t uid, 
					  acct_cluster_cond_t *cluster_q);

/* 
 * remove associations from accounting system 
 * IN:  acct_association_cond_t *assoc_q
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern List acct_storage_g_remove_associations(void *db_conn, uint32_t uid, 
					      acct_association_cond_t *assoc_q);

/* 
 * get info from the storage 
 * IN:  acct_user_cond_t *
 * IN:  params void *
 * returns List of acct_user_rec_t *
 * note List needs to be freed when called
 */
extern List acct_storage_g_get_users(void *db_conn, 
				     acct_user_cond_t *user_q);

/* 
 * get info from the storage 
 * IN:  acct_account_cond_t *
 * IN:  params void *
 * returns List of acct_account_rec_t *
 * note List needs to be freed when called
 */
extern List acct_storage_g_get_accounts(void *db_conn, 
					acct_account_cond_t *acct_q);

/* 
 * get info from the storage 
 * IN:  acct_cluster_cond_t *
 * IN:  params void *
 * returns List of acct_cluster_rec_t *
 * note List needs to be freed when called
 */
extern List acct_storage_g_get_clusters(void *db_conn, 
					acct_cluster_cond_t *cluster_q);

/* 
 * get info from the storage 
 * IN:  acct_association_cond_t *
 * RET: List of acct_association_rec_t *
 * note List needs to be freed when called
 */
extern List acct_storage_g_get_associations(void *db_conn, 
					    acct_association_cond_t *assoc_q);

/* 
 * get info from the storage 
 * IN:  type period specifier
 * IN/OUT:  assoc void * (acct_association_rec_t *) with the id set
 * IN:  start time stamp for records >=
 * IN:  end time stamp for records <
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_get_usage(
	void *db_conn, acct_usage_type_t type, void *assoc,
	time_t start, time_t end);
/* 
 * roll up data in the storage 
 * IN:  type period specifier
 * IN:  start time stamp for records >=
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int acct_storage_g_roll_usage(
	void *db_conn, acct_usage_type_t type, time_t start);

/*********************** CLUSTER ACCOUNTING STORAGE **************************/

extern int clusteracct_storage_g_node_down(void *db_conn, 
					   char *cluster,
					   struct node_record *node_ptr,
					   time_t event_time,
					   char *reason);

extern int clusteracct_storage_g_node_up(void *db_conn, 
					 char *cluster,
					 struct node_record *node_ptr,
					 time_t event_time);

extern int clusteracct_storage_g_cluster_procs(void *db_conn, 
					       char *cluster,
					       uint32_t procs,
					       time_t event_time);

extern int clusteracct_storage_g_register_ctld(char *cluster, uint16_t port);

/* 
 * get info from the storage 
 * IN:  type period specifier
 * IN/OUT:  cluster_rec void * (acct_cluster_rec_t *) with the name set
 * IN:  start time stamp for records >=
 * IN:  end time stamp for records <
 * IN:  params void *
 * RET: SLURM_SUCCESS on success SLURM_ERROR else
 */
extern int clusteracct_storage_g_get_usage(
	void *db_conn, acct_usage_type_t type, void *cluster_rec,
	time_t start, time_t end);

/* 
 * load into the storage the start of a job
 */
extern int jobacct_storage_g_job_start (void *db_conn, 
					struct job_record *job_ptr);

/* 
 * load into the storage the end of a job
 */
extern int jobacct_storage_g_job_complete (void *db_conn, 
					   struct job_record *job_ptr);

/* 
 * load into the storage the start of a job step
 */
extern int jobacct_storage_g_step_start (void *db_conn, 
					 struct step_record *step_ptr);

/* 
 * load into the storage the end of a job step
 */
extern int jobacct_storage_g_step_complete (void *db_conn, 
					    struct step_record *step_ptr);

/* 
 * load into the storage a suspention of a job
 */
extern int jobacct_storage_g_job_suspend (void *db_conn, 
					  struct job_record *job_ptr);

/* 
 * get info from the storage 
 * returns List of jobacct_job_rec_t *
 * note List needs to be freed when called
 */
extern List jobacct_storage_g_get_jobs(void *db_conn, 
				       List selected_steps,
				       List selected_parts,
				       void *params);

/* 
 * expire old info from the storage 
 */
extern void jobacct_storage_g_archive(void *db_conn, 
				      List selected_parts,
				      void *params);

#endif /*_SLURM_ACCOUNTING_STORAGE_H*/
