/*****************************************************************************\
 *  pgsql_common.c - common functions for the the pgsql storage plugin.
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

#include "pgsql_common.h"
#include <stdlib.h>

pthread_mutex_t pgsql_lock = PTHREAD_MUTEX_INITIALIZER;

#ifdef HAVE_PGSQL

static int rollback_started = 0;

extern int *destroy_pgsql_db_info(pgsql_db_info_t *db_info)
{
	if(db_info) {
		xfree(db_info->host);
		xfree(db_info->user);
		xfree(db_info->pass);
		xfree(db_info);
	}
	return SLURM_SUCCESS;
}

extern int pgsql_create_db(PGconn *pgsql_db, char *db_name,
			   pgsql_db_info_t *db_info)
{
	char create_line[50];
	char *connect_line = xstrdup_printf("dbname = 'postgres'"
					    " host = '%s'"
					    " port = '%u'"
					    " user = '%s'"
					    " password = '%s'",
					    db_info->host,
					    db_info->port,
					    db_info->user,
					    db_info->pass);

	pgsql_db = PQconnectdb(connect_line);

	if (PQstatus(pgsql_db) == CONNECTION_OK) {
		PGresult *result = NULL;
		snprintf(create_line, sizeof(create_line),
			 "create database %s", db_name);
		result = PQexec(pgsql_db, create_line);
		if (PQresultStatus(result) != PGRES_COMMAND_OK) {
			fatal("PQexec failed: %d %s\n%s",
			     PQresultStatus(result), PQerrorMessage(pgsql_db), create_line);
		}
		PQclear(result);
	} else {
		info("Connection failed to %s", connect_line);
		fatal("Status was: %d %s",
		      PQstatus(pgsql_db), PQerrorMessage(pgsql_db));
	}
	xfree(connect_line);
	return SLURM_SUCCESS;
}

extern int pgsql_get_db_connection(PGconn **pgsql_db, char *db_name,
				   pgsql_db_info_t *db_info, bool rollback)
{
	int rc = SLURM_SUCCESS;
	bool storage_init = false;
	char *connect_line = xstrdup_printf("dbname = '%s'"
					    " host = '%s'"
					    " port = '%u'"
					    " user = '%s'"
					    " password = '%s'",
					    db_name,
					    db_info->host,
					    db_info->port,
					    db_info->user,
					    db_info->pass);

	while(!storage_init) {
		*pgsql_db = PQconnectdb(connect_line);
		
		if(PQstatus(*pgsql_db) != CONNECTION_OK) {
			if(!strcmp(PQerrorMessage(*pgsql_db),
				   "no password supplied")) {
				PQfinish(*pgsql_db);
				fatal("This Postgres connection needs "
				      "a password.  It doesn't appear to "
				      "like blank ones");
			} 
			
			info("Database %s not created. Creating", db_name);
			PQfinish(*pgsql_db);
			pgsql_create_db(*pgsql_db, db_name, db_info);		
		} else {
			storage_init = true;
			debug2("connected to %s", db_name);
			if(rollback || rollback_started) {
				rollback_started = 1;
				PQexec(*pgsql_db, "BEGIN WORK");
			}
		} 
	}
	xfree(connect_line);
	return rc;
}

extern int pgsql_close_db_connection(PGconn *pgsql_db, bool commit)
{
	if(rollback_started) {
		if(commit) 
			PQexec(pgsql_db, "COMMIT WORK");
		else 
			PQexec(pgsql_db, "ROLLBACK WORK");
	}

	PQfinish(pgsql_db);	      
	return SLURM_SUCCESS;
}

extern int pgsql_db_query(PGconn *pgsql_db, char *query)
{
	PGresult *result = NULL;
	
	if(!pgsql_db)
		fatal("You haven't inited this storage yet.");
	
	if(!(result = pgsql_db_query_ret(pgsql_db, query))) 
		return SLURM_ERROR;
	
	PQclear(result);
	return SLURM_SUCCESS;
}

extern PGresult *pgsql_db_query_ret(PGconn *pgsql_db, char *query)
{
	PGresult *result = NULL;
	
	if(!pgsql_db)
		fatal("You haven't inited this storage yet.");

	result = PQexec(pgsql_db, query);

	if(PQresultStatus(result) != PGRES_COMMAND_OK
	   && PQresultStatus(result) != PGRES_TUPLES_OK) {
		error("PQexec failed: %d %s", PQresultStatus(result), 
		      PQerrorMessage(pgsql_db));
		info("query was %s", query);
		PQclear(result);
		return NULL;
	}
	return result;
}

extern int pgsql_insert_ret_id(PGconn *pgsql_db, char *sequence_name,
			       char *query)
{
	int new_id = 0;
	PGresult *result = NULL;

	slurm_mutex_lock(&pgsql_lock);
	if(pgsql_db_query(pgsql_db, query) != SLURM_ERROR)  {
		char *new_query = xstrdup_printf(
			"select last_value from %s", sequence_name);
		
		if((result = pgsql_db_query_ret(pgsql_db, new_query))) {
			new_id = atoi(PQgetvalue(result, 0, 0));
			PQclear(result);		
		}
		xfree(new_query);
		if(!new_id) {
			/* should have new id */
			error("We should have gotten a new id: %s", 
			      PQerrorMessage(pgsql_db));
		}
	}
	slurm_mutex_unlock(&pgsql_lock);
	
	return new_id;
	
}

extern int pgsql_db_create_table(PGconn *pgsql_db,  
				 char *table_name, storage_field_t *fields,
				 char *ending)
{
	char *query = NULL;
	char *tmp = NULL;
	char *next = NULL;
	int i = 0;

	query = xstrdup_printf("create table %s (", table_name);
	i=0;
	while(fields && fields->name) {
		next = xstrdup_printf(" %s %s",
				      fields->name, 
				      fields->options);
		if(i) 
			xstrcat(tmp, ",");
		xstrcat(tmp, next);
		xfree(next);
		fields++;
		i++;
	}
	xstrcat(query, tmp);
	xfree(tmp);
	xstrcat(query, ending);

	if(pgsql_db_query(pgsql_db, query) == SLURM_ERROR) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);

	return SLURM_SUCCESS;
}

extern int pgsql_db_make_table_current(PGconn *pgsql_db, char *table_name,
				       storage_field_t *fields)
{
	char *query = NULL, *opt_part = NULL, *temp_char = NULL;
	char *type = NULL;
	int not_null = 0;
	char *default_str = NULL, *default_query = NULL, *null_query = NULL;
	char* original_ptr = NULL;
	int i = 0;

	while(fields[i].name) {
		if(!strcmp("serial", fields[i].options)) {
			i++;
			continue;
		} 
		opt_part = xstrdup(fields[i].options);
		original_ptr = opt_part;
		opt_part = strtok_r(opt_part, " ", &temp_char);
		if(opt_part) {
			type = xstrdup(opt_part);
			opt_part = temp_char;
			opt_part = strtok_r(opt_part, " ", &temp_char);
			while(opt_part) {
				if(!strcmp("not null", opt_part)) {
					not_null = 1;
					opt_part = temp_char;
					opt_part = strtok_r(opt_part,
							    " ", &temp_char);
				} else if(!strcmp("default", opt_part)){
					opt_part = temp_char;
					opt_part = strtok_r(opt_part,
							    " ", &temp_char);
					default_str = xstrdup(opt_part);
				}
				if(opt_part) {
					opt_part = temp_char;
					opt_part = strtok_r(opt_part,
							    " ", &temp_char);
				}
			}
		} else {
			type = xstrdup(fields[i].options);
		}
		xfree(original_ptr);

		if(default_str) 
			default_query = xstrdup_printf(
				", alter column %s set default %s",
				fields[i].name, default_str);
		else 
			default_query = xstrdup_printf(
				", alter column %s drop default",
				fields[i].name);

		if(not_null) 
			null_query = xstrdup_printf(
				", alter column %s set not null",
				fields[i].name);
		else 
			null_query = xstrdup_printf(
				", alter column %s drop not null",
				fields[i].name);
		
		
		query = xstrdup_printf("alter table %s alter column "
				       "%s type %s%s%s",
				       table_name, fields[i].name, type,
				       default_query, null_query);
		xfree(default_query);
		xfree(null_query);
			
		if(pgsql_db_query(pgsql_db, query)) {
			info("adding column %s", fields[i].name);
			if(default_str) 
				default_query = xstrdup_printf(
					" default %s", default_str);
						
			if(not_null) 
				null_query = xstrdup_printf(" not null");
			
			xfree(query);
			query = xstrdup_printf(
				"alter table %s add %s %s",
				table_name, fields[i].name,
				type);
			if(default_query) {
				xstrcat(query, default_query);
				xfree(default_query);
			}

			if(null_query) {
				xstrcat(query, null_query);
				xfree(null_query);
			}

			if(pgsql_db_query(pgsql_db, query)) {
				xfree(default_str);
				xfree(query);
				xfree(type);
				return SLURM_ERROR;
			}
			

		}
		xfree(default_str);
		xfree(query);
		xfree(type);
		i++;
	}
	
	return SLURM_SUCCESS;
}


#endif

