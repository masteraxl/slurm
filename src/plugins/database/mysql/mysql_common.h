/*****************************************************************************\
 *  mysql_common.h - common functions for the the mysql database plugin.
 *****************************************************************************
 *
 *  Copyright (C) 2004-2007 The Regents of the University of California.
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
 *
 *  This file is patterned after jobcomp_linux.c, written by Morris Jette and
 *  Copyright (C) 2002 The Regents of the University of California.
\*****************************************************************************/
#ifndef _HAVE_MYSQL_COMMON_H
#define _HAVE_MYSQL_COMMON_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#if HAVE_STDINT_H
#  include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#include <stdio.h>
#include <slurm/slurm_errno.h>
#include "src/slurmctld/slurmctld.h"

#ifdef HAVE_MYSQL
#include <mysql/mysql.h>
#include <mysql/mysqld_error.h>

typedef struct {
	uint32_t port;	
	char *host;	
	char *user;	
	char *pass;	
} mysql_db_info_t;
#endif

extern bool thread_safe;
extern pthread_mutex_t mysql_lock;

#ifdef HAVE_MYSQL
extern mysql_db_info_t *create_mysql_db_info();

extern int *destroy_mysql_db_info(mysql_db_info_t *db_info);

extern int mysql_create_db(MYSQL *mysql_db, char *db_name,
			   mysql_db_info_t *db_info);

extern int mysql_get_db_connection(MYSQL **mysql_db, char *db_name,
				   mysql_db_info_t *db_info,
				   int *database_init);
extern int mysql_db_query(MYSQL *mysql_db, int database_init, char *query);

extern MYSQL_RES *mysql_db_query_ret(MYSQL *mysql_db, int database_init,
				     char *query);

extern int mysql_insert_ret_id(MYSQL *mysql_db, int database_init, 
			       char *query);
#endif

#endif
