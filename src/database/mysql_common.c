/*****************************************************************************\
 *  mysql_common.c - common functions for the the mysql storage plugin.
 *****************************************************************************
 *
 *  Copyright (C) 2004-2007 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
 *  
 *  This file is part of SLURM, a resource management program.
 *  For details, see <https://computing.llnl.gov/linux/slurm/>.
 *  Please also read the included file: DISCLAIMER.
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

#include "mysql_common.h"
#include "src/common/xmalloc.h"
#include "src/common/timers.h"
#include "src/common/slurm_protocol_api.h"

#ifdef MYSQL_NOT_THREAD_SAFE
pthread_mutex_t mysql_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

static char *table_defs_table = "table_defs_table";

static MYSQL_RES *_get_first_result(MYSQL *mysql_db)
{
	MYSQL_RES *result = NULL;
	int rc = 0;
	do {
		/* did current statement return data? */
		if((result = mysql_store_result(mysql_db)))
			return result;
		
		/* more results? -1 = no, >0 = error, 0 = yes (keep looping) */
		if ((rc = mysql_next_result(mysql_db)) > 0)
			debug3("error: Could not execute statement %d\n", rc);
		
	} while (rc == 0);
	
	return NULL;
}

static MYSQL_RES *_get_last_result(MYSQL *mysql_db)
{
	MYSQL_RES *result = NULL;
	MYSQL_RES *last_result = NULL;
	int rc = 0;
	do {
		/* did current statement return data? */
		if((result = mysql_store_result(mysql_db))) {
			if(last_result)
				mysql_free_result(last_result);
			last_result = result;
		}
		/* more results? -1 = no, >0 = error, 0 = yes (keep looping) */
		if ((rc = mysql_next_result(mysql_db)) > 0)
			debug3("error: Could not execute statement %d\n", rc);
	} while (rc == 0);
	
	return last_result;
}

static int _mysql_make_table_current(MYSQL *mysql_db, char *table_name,
				     storage_field_t *fields, char *ending)
{
	char *query = NULL;
	char *correct_query = NULL;
	MYSQL_RES *result = NULL;
	MYSQL_ROW row;
	int i = 0;
	List columns = NULL;
	ListIterator itr = NULL;
	char *col = NULL;
	int adding = 0;
	int run_update = 0;
	char *primary_key = NULL;
	char *unique_index = NULL;
	int old_primary = 0;
	char *old_index = NULL;
	char *temp = NULL;

	DEF_TIMERS;
	
	/* figure out the keys in the table */
	query = xstrdup_printf("show index from %s", table_name);
	if(!(result = mysql_db_query_ret(mysql_db, query, 0))) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);
	while((row = mysql_fetch_row(result))) {
		// row[2] is the key name
		if(!strcasecmp(row[2], "PRIMARY"))
			old_primary = 1;
		else if(!old_index)
			old_index = xstrdup(row[2]);
	}
	mysql_free_result(result);

	/* figure out the existing columns in the table */
	query = xstrdup_printf("show columns from %s", table_name);
	if(!(result = mysql_db_query_ret(mysql_db, query, 0))) {
		xfree(query);
		xfree(old_index);
		return SLURM_ERROR;
	}
	xfree(query);
	columns = list_create(slurm_destroy_char);
	while((row = mysql_fetch_row(result))) {
		col = xstrdup(row[0]); //Field
		list_append(columns, col);
	}
	mysql_free_result(result);


	itr = list_iterator_create(columns);
	query = xstrdup_printf("alter table %s", table_name);
	correct_query = xstrdup_printf("alter table %s", table_name);
	START_TIMER;
	while(fields[i].name) {
		int found = 0;
		list_iterator_reset(itr);
		while((col = list_next(itr))) {
			if(!strcmp(col, fields[i].name)) {
				xstrfmtcat(query, " modify %s %s,",
					   fields[i].name,
					   fields[i].options);
				xstrfmtcat(correct_query, " modify %s %s,",
					   fields[i].name,
					   fields[i].options);
				list_delete_item(itr);
				found = 1;
				break;
			}
		}
		if(!found) {
			if(i) {
				info("adding column %s after %s in table %s",
				     fields[i].name,
				     fields[i-1].name,
				     table_name);
				xstrfmtcat(query, " add %s %s after %s,",
					   fields[i].name,
					   fields[i].options,
					   fields[i-1].name);
				xstrfmtcat(correct_query, " modify %s %s,",
					   fields[i].name,
					   fields[i].options);
			} else {
				info("adding column %s at the beginning "
				     "of table %s",
				     fields[i].name,
				     fields[i-1].name,
				     table_name);
				xstrfmtcat(query, " add %s %s first,",
					   fields[i].name,
					   fields[i].options);
				xstrfmtcat(correct_query, " modify %s %s,",
					   fields[i].name,
					   fields[i].options);
			}
			adding = 1;
		}

		i++;
	}

	list_iterator_reset(itr);
	while((col = list_next(itr))) {
		adding = 1;
		info("dropping column %s from table %s", col, table_name);
		xstrfmtcat(query, " drop %s,", col);
	}
	
	list_iterator_destroy(itr);
	list_destroy(columns);
	
	if((temp = strstr(ending, "primary key ("))) {
		int open = 0, close =0;
		int end = 0;
		while(temp[end++]) {
			if(temp[end] == '(')
				open++;
			else if(temp[end] == ')')
				close++;
			else
				continue;
			if(open == close)
				break;
		}
		if(temp[end]) {
			end++;
			primary_key = xstrndup(temp, end);
			if(old_primary) {
				xstrcat(query, " drop primary key,");
				xstrcat(correct_query, " drop primary key,");
			}
			xstrfmtcat(query, " add %s,",  primary_key);
			xstrfmtcat(correct_query, " add %s,",  primary_key);
			
			xfree(primary_key);
		}
	}

	if((temp = strstr(ending, "unique index ("))) {
		int open = 0, close =0;
		int end = 0;
		while(temp[end++]) {
			if(temp[end] == '(')
				open++;
			else if(temp[end] == ')')
				close++;
			else
				continue;
			if(open == close)
				break;
		}
		if(temp[end]) {
			end++;
			unique_index = xstrndup(temp, end);
			if(old_index) {
				xstrfmtcat(query, " drop index %s,",
					   old_index);
				xstrfmtcat(correct_query, " drop index %s,",
					   old_index);
			}
			xstrfmtcat(query, " add %s,", unique_index);
			xstrfmtcat(correct_query, " add %s,", unique_index);
			xfree(unique_index);
		}
	}
	xfree(old_index);

	query[strlen(query)-1] = ';';
	correct_query[strlen(correct_query)-1] = ';';
	//info("%d query\n%s", __LINE__, query);

	/* see if we have already done this definition */
	if(!adding) {
		char *query2 = xstrdup_printf("select table_name from "
					      "%s where definition=\"%s\"",
					      table_defs_table, query);
		MYSQL_RES *result = NULL;
		MYSQL_ROW row;
		
		run_update = 1;
		
		if((result = mysql_db_query_ret(mysql_db, query2, 0))) {
			if((row = mysql_fetch_row(result)))
				run_update = 0;
			mysql_free_result(result);
		}
		xfree(query2);
	}

	/* if something has changed run the alter line */
	if(run_update || adding) {
		time_t now = time(NULL);
		char *query2 = NULL;
	
		debug("Table %s has changed.  Updating...", table_name);

		if(mysql_db_query(mysql_db, query)) {
			xfree(query);
			return SLURM_ERROR;
		}
		
		query2 = xstrdup_printf("insert into %s (creation_time, "
					"mod_time, table_name, definition) "
					"values (%d, %d, \"%s\", \"%s\") "
					"on duplicate key update "
					"definition=\"%s\", mod_time=%d;",
					table_defs_table, now, now,
					table_name, correct_query,
					correct_query, now);
		if(mysql_db_query(mysql_db, query2)) {
			xfree(query2);
			return SLURM_ERROR;
		}
		xfree(query2);
	}

	xfree(query);
	xfree(correct_query);
	query = xstrdup_printf("make table current %s", table_name);
	END_TIMER2(query);
	xfree(query);
	return SLURM_SUCCESS;
}

static int _create_db(char *db_name, mysql_db_info_t *db_info)
{
	char create_line[50];
	MYSQL *mysql_db = NULL;
	int rc = SLURM_ERROR;
	
	MYSQL *db_ptr = NULL;
	char *db_host = NULL;

	while(rc == SLURM_ERROR) {
		rc = SLURM_SUCCESS;
#ifdef MYSQL_NOT_THREAD_SAFE
		slurm_mutex_lock(&mysql_lock);
#endif
		if(!(mysql_db = mysql_init(mysql_db)))
			fatal("mysql_init failed: %s", mysql_error(mysql_db));
		
		db_host = db_info->host;
		db_ptr = mysql_real_connect(mysql_db,
					    db_host, db_info->user,
					    db_info->pass, NULL,
					    db_info->port, NULL, 0);

		if (!db_ptr && db_info->backup) {
			info("Connection failed to host = %s "
			     "user = %s port = %u",
			     db_host, db_info->user,
			     db_info->port);  
			db_host = db_info->backup;
			db_ptr = mysql_real_connect(mysql_db, db_host,
						    db_info->user,
						    db_info->pass, NULL,
						    db_info->port, NULL, 0);
		}

		if (db_ptr) {
			snprintf(create_line, sizeof(create_line),
				 "create database %s", db_name);
			if(mysql_query(mysql_db, create_line)) {
				fatal("mysql_real_query failed: %d %s\n%s",
				      mysql_errno(mysql_db),
				      mysql_error(mysql_db), create_line);
			}
			mysql_close_db_connection(&mysql_db);
		} else {
			info("Connection failed to host = %s "
			     "user = %s port = %u",
			     db_host, db_info->user,
			     db_info->port);
#ifdef MYSQL_NOT_THREAD_SAFE
			slurm_mutex_unlock(&mysql_lock);
#endif
			error("mysql_real_connect failed: %d %s\n",
			      mysql_errno(mysql_db),
			      mysql_error(mysql_db));
			rc = SLURM_ERROR;
		}
#ifdef MYSQL_NOT_THREAD_SAFE
		slurm_mutex_unlock(&mysql_lock);
#endif
		if(rc == SLURM_ERROR)
			sleep(3);
	}
	return rc;
}

extern int *destroy_mysql_db_info(mysql_db_info_t *db_info)
{
	if(db_info) {
		xfree(db_info->backup);
		xfree(db_info->host);
		xfree(db_info->user);
		xfree(db_info->pass);
		xfree(db_info);
	}
	return SLURM_SUCCESS;
}

extern int mysql_get_db_connection(MYSQL **mysql_db, char *db_name,
				   mysql_db_info_t *db_info)
{
	int rc = SLURM_SUCCESS;
	bool storage_init = false;
	
	char *db_host = db_info->host;

	if(!(*mysql_db = mysql_init(*mysql_db)))
		fatal("mysql_init failed: %s", mysql_error(*mysql_db));
	else {
		unsigned int my_timeout = 30;
#ifdef MYSQL_OPT_RECONNECT
		my_bool reconnect = 1;
		/* make sure reconnect is on */
		mysql_options(*mysql_db, MYSQL_OPT_RECONNECT, &reconnect);
#endif
		mysql_options(*mysql_db, MYSQL_OPT_CONNECT_TIMEOUT,
			      (char *)&my_timeout);
		while(!storage_init) {
			if(!mysql_real_connect(*mysql_db, db_host,
					       db_info->user, db_info->pass,
					       db_name, db_info->port,
					       NULL, CLIENT_MULTI_STATEMENTS)) {
				if(mysql_errno(*mysql_db) == ER_BAD_DB_ERROR) {
					debug("Database %s not created.  "
					      "Creating", db_name);
					rc = _create_db(db_name, db_info);
				} else {
					error("mysql_real_connect failed: "
					      "%d %s",
					      mysql_errno(*mysql_db),
					      mysql_error(*mysql_db));
					if ((db_host == db_info->host)
					    && db_info->backup) {
						db_host = db_info->backup;
						continue;
					}
					rc = SLURM_ERROR;
					break;
				}
			} else {
				storage_init = true;
			}
		}
	}
	return rc;
}

extern int mysql_close_db_connection(MYSQL **mysql_db)
{
	if(mysql_db && *mysql_db) {
		if(mysql_thread_safe())
			mysql_thread_end();
		mysql_close(*mysql_db);
		*mysql_db = NULL;
	}

	return SLURM_SUCCESS;
}

extern int mysql_cleanup()
{
	debug3("starting mysql cleaning up");

#ifdef mysql_library_end
	mysql_library_end();
#else
	mysql_server_end();
#endif

	debug3("finished mysql cleaning up");
	return SLURM_SUCCESS;
}

extern int mysql_clear_results(MYSQL *mysql_db)
{
	MYSQL_RES *result = NULL;
	int rc = 0;
	do {
		/* did current statement return data? */
		if((result = mysql_store_result(mysql_db)))
			mysql_free_result(result);
		
		/* more results? -1 = no, >0 = error, 0 = yes (keep looping) */
		if ((rc = mysql_next_result(mysql_db)) > 0)
			error("Could not execute statement %d %s\n",
			      mysql_errno(mysql_db),
			      mysql_error(mysql_db));
	} while (rc == 0);

	if(rc > 0) {
		errno = rc;
		return SLURM_ERROR;
	} 

	return SLURM_SUCCESS;
}

extern int mysql_db_query(MYSQL *mysql_db, char *query)
{
	if(!mysql_db)
		fatal("You haven't inited this storage yet.");
#ifdef MYSQL_NOT_THREAD_SAFE
	slurm_mutex_lock(&mysql_lock);
#endif
	/* clear out the old results so we don't get a 2014 error */
	mysql_clear_results(mysql_db);		
//try_again:
	if(mysql_query(mysql_db, query)) {
		error("mysql_query failed: %d %s\n%s",
		      mysql_errno(mysql_db),
		      mysql_error(mysql_db), query);
		errno = mysql_errno(mysql_db);
#ifdef MYSQL_NOT_THREAD_SAFE
		slurm_mutex_unlock(&mysql_lock);
#endif
		/* FIXME: If we get ER_LOCK_WAIT_TIMEOUT here we need
		to restart the connections, but it appears restarting
		the calling program is the only way to handle this.
		If anyone in the future figures out a way to handle
		this, super.  Until then we will need to restart the
		calling program if you ever get this error. 
		*/

		return SLURM_ERROR;
	}

#ifdef MYSQL_NOT_THREAD_SAFE
	slurm_mutex_unlock(&mysql_lock);
#endif
	return SLURM_SUCCESS;
}

extern int mysql_db_ping(MYSQL *mysql_db)
{
	/* clear out the old results so we don't get a 2014 error */
	mysql_clear_results(mysql_db);		
	return mysql_ping(mysql_db);
}

extern int mysql_db_commit(MYSQL *mysql_db)
{
#ifdef MYSQL_NOT_THREAD_SAFE
	slurm_mutex_lock(&mysql_lock);
#endif
	/* clear out the old results so we don't get a 2014 error */
	mysql_clear_results(mysql_db);		
	if(mysql_commit(mysql_db)) {
		error("mysql_commit failed: %d %s",
		      mysql_errno(mysql_db),
		      mysql_error(mysql_db));
		errno = mysql_errno(mysql_db);
#ifdef MYSQL_NOT_THREAD_SAFE
		slurm_mutex_unlock(&mysql_lock);
#endif
		return SLURM_ERROR;
	}
#ifdef MYSQL_NOT_THREAD_SAFE
	slurm_mutex_unlock(&mysql_lock);
#endif
	return SLURM_SUCCESS;
}

extern int mysql_db_rollback(MYSQL *mysql_db)
{
#ifdef MYSQL_NOT_THREAD_SAFE
	slurm_mutex_lock(&mysql_lock);
#endif
	/* clear out the old results so we don't get a 2014 error */
	mysql_clear_results(mysql_db);		
	if(mysql_rollback(mysql_db)) {
		error("mysql_commit failed: %d %s",
		      mysql_errno(mysql_db),
		      mysql_error(mysql_db));
		errno = mysql_errno(mysql_db);
#ifdef MYSQL_NOT_THREAD_SAFE
		slurm_mutex_unlock(&mysql_lock);
#endif
		return SLURM_ERROR;
	}
	//mysql_db_query(mysql_db, "unlock tables;");
#ifdef MYSQL_NOT_THREAD_SAFE
	slurm_mutex_unlock(&mysql_lock);
#endif
	return SLURM_SUCCESS;

}

extern MYSQL_RES *mysql_db_query_ret(MYSQL *mysql_db, char *query, bool last)
{
	MYSQL_RES *result = NULL;
	
	if(mysql_db_query(mysql_db, query) != SLURM_ERROR)  {
		if(last)
			result = _get_last_result(mysql_db);
		else
			result = _get_first_result(mysql_db);
		if(!result && mysql_field_count(mysql_db)) {
			/* should have returned data */
			error("We should have gotten a result: %s", 
			      mysql_error(mysql_db));
		}
	}

	return result;
}

extern int mysql_db_query_check_after(MYSQL *mysql_db, char *query)
{
	int rc = SLURM_SUCCESS;
		
	if((rc = mysql_db_query(mysql_db, query)) != SLURM_ERROR)  
		rc = mysql_clear_results(mysql_db);
	
	return rc;
}

extern int mysql_insert_ret_id(MYSQL *mysql_db, char *query)
{
	int new_id = 0;
	
	if(mysql_db_query(mysql_db, query) != SLURM_ERROR)  {
		new_id = mysql_insert_id(mysql_db);
		if(!new_id) {
			/* should have new id */
			error("We should have gotten a new id: %s", 
			      mysql_error(mysql_db));
		}
	}

	return new_id;
	
}

extern int mysql_db_create_table(MYSQL *mysql_db, char *table_name,
				 storage_field_t *fields, char *ending)
{
	char *query = NULL;
	int i = 0;
	storage_field_t *first_field = fields;

	if(!fields || !fields->name) {
		error("Not creating an empty table");
		return SLURM_ERROR;
	}

	/* We have an internal table called table_defs_table which
	 * contains the definition of each table in the database.  To
	 * speed things up we just check against that to see if
	 * anything has changed.
	 */
	query = xstrdup_printf("create table if not exists %s "
			       "(creation_time int unsigned not null, "
			       "mod_time int unsigned default 0 not null, "
			       "table_name text not null, "
			       "definition text not null, "
			       "primary key (table_name(50))) engine='innodb'",
			       table_defs_table);

	if(mysql_db_query(mysql_db, query) == SLURM_ERROR) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);	

	query = xstrdup_printf("create table if not exists %s (%s %s",
			       table_name, fields->name, fields->options);
	i=1;
	fields++;
		
	while(fields && fields->name) {
		xstrfmtcat(query, ", %s %s", fields->name, fields->options);
		fields++;
		i++;
	}
	xstrcat(query, ending);

	/* make sure we can do a rollback */
	xstrcat(query, " engine='innodb'");

	if(mysql_db_query(mysql_db, query) == SLURM_ERROR) {
		xfree(query);
		return SLURM_ERROR;
	}
	xfree(query);	
	
	return _mysql_make_table_current(mysql_db, table_name,
					 first_field, ending);
}
