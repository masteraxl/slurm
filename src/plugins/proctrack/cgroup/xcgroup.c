/*****************************************************************************\
 *  xcgroup.c - cgroup related primitives
 *****************************************************************************
 *  Copyright (C) 2009 CEA/DAM/DIF
 *  Written by Matthieu Hautreux <matthieu.hautreux@cea.fr>
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
\*****************************************************************************/

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDINT_H
#  include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include <slurm/slurm.h>
#include <slurm/slurm_errno.h>
#include "src/common/log.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/slurmd/slurmstepd/slurmstepd_job.h"

#include "xcgroup.h"

#ifndef PATH_MAX
#define PATH_MAX 256
#endif

/* internal functions */
size_t _file_getsize(int fd);
int _file_read_uint32s(char* file_path,uint32_t** pvalues,int* pnb);
int _file_write_uint32s(char* file_path,uint32_t* values,int nb);
int _file_read_uint64s(char* file_path,uint64_t** pvalues,int* pnb);
int _file_write_uint64s(char* file_path,uint64_t* values,int nb);
int _file_read_content(char* file_path,char** content,size_t *csize);
int _file_write_content(char* file_path,char* content,size_t csize);
int _xcgroup_cpuset_init(char* file_path);

/* xcgroup primitives */
int xcgroup_is_available(void)
{
	char* value;
	size_t s;

	if ( xcgroup_get_param(CGROUP_BASEDIR,"release_agent",
			       &value,&s) != XCGROUP_SUCCESS )
		return 0;
	else {
		xfree(value);
		return 1;
	}

}

int xcgroup_set_release_agent(char* agent)
{
	int fstatus;
	char* rag;
	char* value;
	size_t s;

	if ( agent == NULL )
		return XCGROUP_ERROR;

	rag = (char*) xstrdup("release_agent=");
	fstatus = xcgroup_get_param(CGROUP_BASEDIR,"release_agent",
				    &value,&s);
	if (  fstatus == XCGROUP_SUCCESS ) {
		if ( strcmp(value,agent) != 0 ) {
			xstrcat(rag,agent);
			fstatus = xcgroup_set_params(CGROUP_BASEDIR,rag);
		}
		xfree(value);
	}

	xfree(rag);
	return fstatus;
}

int xcgroup_mount(char* mount_opts)
{
	char* mount_cmd_fmt;
	char mount_cmd[1024];

	mode_t cmask;
	mode_t omask;

	cmask = S_IWGRP | S_IWOTH;
	omask = umask(cmask);

	if ( mkdir(CGROUP_BASEDIR,0755) && errno != EEXIST) {
		debug("unable to create cgroup directory '%s'"
		      " : %m",CGROUP_BASEDIR);
		umask(omask);
		return XCGROUP_ERROR;
	}
	umask(omask);

	if ( mount_opts == NULL ||
	     strlen(mount_opts) == 0 ) {
		mount_cmd_fmt="/bin/mount -t cgroup none " CGROUP_BASEDIR;
	}
	else
		mount_cmd_fmt="/bin/mount -o %s -t cgroup none " CGROUP_BASEDIR;

	if ( snprintf(mount_cmd,1024,mount_cmd_fmt,
		      mount_opts) >= 1024 ) {
		debug2("unable to build mount cmd line");
		return XCGROUP_ERROR;
	}
	else
		debug3("cgroup mount cmd line is '%s'",mount_cmd);

	if ( system(mount_cmd) )
		return XCGROUP_ERROR;
	else
		return XCGROUP_SUCCESS;

}

int xcgroup_create(char* file_path,xcgroup_opts_t* opts)
{
	int fstatus;
	uid_t uid;
	gid_t gid;
	int create_only;
	int notify;
	
	mode_t cmask;
	mode_t omask;

	uid=opts->uid;
	gid=opts->gid;
	create_only=opts->create_only;
	notify=opts->notify;
	
	fstatus = XCGROUP_ERROR;

	/* save current mask and apply working one */
	cmask = S_IWGRP | S_IWOTH;
	omask = umask(cmask);

	/* build cgroup */
	if ( mkdir(file_path,0755) ) {
		if ( create_only || errno != EEXIST ) {
			debug2("unable to create cgroup '%s' : %m",
			       file_path);
			umask(omask);
			return fstatus;
		}
	}
	umask(omask);

	/* initialize cpuset support (if enabled in cgroup ) */
	if ( _xcgroup_cpuset_init(file_path) != XCGROUP_SUCCESS ) {
		debug2("unable to initialize cpuset cgroup component");
		rmdir(file_path);
		return fstatus;
	}
	
	/* change cgroup ownership as requested */
	if ( chown(file_path,uid,gid) ) {
		debug2("unable to chown %d:%d cgroup '%s' : %m",
		       uid,gid,file_path);
		return fstatus;
	}

	/* following operations failure might not result in a general
	 * failure so set output status to success */
	fstatus = XCGROUP_SUCCESS;

	/* set notify on release flag */
	if ( notify == 1 )
		xcgroup_set_params(file_path,"notify_on_release=1");
	else if ( notify == 0 )
		xcgroup_set_params(file_path,"notify_on_release=0");

	return fstatus;
}

int xcgroup_destroy(char* file_path)
{

	/* 
	 * nothing to be done here, notify_on_release was set
	 * so hope that all will works perfectly...
	 *
	 * with memory cgroup some pages can still be accounted
	 * to the cgroup but no more processes are present, this results
	 * in a directory not being removed until the pages are accounted
	 * to an other cgroup...
	 * echoing 1 into memory.force_empty can purge this memory but 
	 * as slurmstepd is still present in the cgroup and use pages, 
	 * this is not sufficient as it could leave some other pages too..
	 * we should have a way to ask the cgroup to force_empty
	 * on last process exit but I did not find any for now
	 */
	//xcgroup_set_params(file_path,"memory.force_empty=1");

	return XCGROUP_SUCCESS;
}

int xcgroup_add_pids(char* cpath,pid_t* pids,int npids)
{
	int fstatus;
	char file_path[PATH_MAX];
	
	fstatus = XCGROUP_ERROR;

	if ( snprintf(file_path,PATH_MAX,"%s/tasks",
		      cpath) >= PATH_MAX ) {
		debug2("unable to add pids to '%s' : %m",cpath);
		return fstatus;
	}

	fstatus = _file_write_uint32s(file_path,(uint32_t*)pids,npids);
	if ( fstatus != XCGROUP_SUCCESS )
		debug2("unable to add pids to '%s'",cpath);
	return fstatus;
}

int
xcgroup_get_pids(char* cpath, pid_t **pids, int *npids)
{
	int fstatus;
	char file_path[PATH_MAX];

	fstatus = XCGROUP_ERROR;

	if ( pids == NULL || npids == NULL )
		return SLURM_ERROR;

	if ( snprintf(file_path,PATH_MAX,"%s/tasks",
		      cpath) >= PATH_MAX ) {
		debug2("unable to get pids of '%s' : %m",cpath);
		return fstatus;
	}
	
	fstatus = _file_read_uint32s(file_path,(uint32_t**)pids,npids);
	if ( fstatus != XCGROUP_SUCCESS )
		debug2("unable to get pids of '%s'",cpath);
	return fstatus;
}

int
xcgroup_find_by_pid(char* cpath, pid_t pid)
{
	int fstatus = SLURM_ERROR;
	char file_path[PATH_MAX];
	char* buf;
	size_t fsize;
	char* p;
	char* e;
	char* entry;

	/* build pid cgroup meta filepath */
	if ( snprintf(file_path,PATH_MAX,"/proc/%u/cgroup",
		      pid) >= PATH_MAX ) {
		debug2("unable to build cgroup meta filepath for pid=%u : %m",
		       pid);
		return XCGROUP_ERROR;
	}	

	/* read file content */
	fstatus = _file_read_content(file_path,&buf,&fsize);
	if ( fstatus == XCGROUP_SUCCESS ) {
		fstatus = XCGROUP_ERROR;
		p = buf;
		if ( index(p,'\n') != NULL ) {
			e = index(p,'\n');
			*e='\0';
			entry = rindex(p,':');
			if ( entry != NULL ) {
				entry++;
				snprintf(cpath,PATH_MAX,"%s%s",
					 CGROUP_BASEDIR,entry);
				fstatus = XCGROUP_SUCCESS;
			}
		}
		xfree(buf);
	}

	return fstatus;
}

int xcgroup_set_memlimit(char* cpath,uint32_t memlimit)
{
	int fstatus;
	char file_path[PATH_MAX];
	uint64_t ml;
	
	fstatus = XCGROUP_ERROR;
	
	if ( snprintf(file_path,PATH_MAX,"%s/memory.limit_in_bytes",
		      cpath) >= PATH_MAX ) {
		debug2("unable to set memory limit of '%s' : %m",cpath);
		return fstatus;
	}
	
	ml = (uint64_t) memlimit * 1024 * 1024;
	fstatus = _file_write_uint64s(file_path,&ml,1);
	if ( fstatus != XCGROUP_SUCCESS )
		debug2("unable to set memory limit of '%s' : %m",cpath);
	else
		debug3("memory limit set to %uMB for '%s'",memlimit,cpath);

	return fstatus;
}

int xcgroup_get_memlimit(char* cpath,uint32_t* memlimit)
{
	int fstatus;
	char file_path[PATH_MAX];
	uint64_t* ml;
	int i;
	
	fstatus = XCGROUP_ERROR;
	
	if ( snprintf(file_path,PATH_MAX,"%s/memory.limit_in_bytes",
		      cpath) >= PATH_MAX ) {
		debug2("unable to get memory limit of '%s' : %m",cpath);
		return fstatus;
	}
	
	fstatus = _file_read_uint64s(file_path,&ml,&i);
	if ( fstatus != XCGROUP_SUCCESS || 
	     i == 0 )
		debug2("unable to get memory limit of '%s' : %m",cpath);
	else {
		if ( *ml == 0 ) {
			*memlimit = 0;	
		}
		else {
			/* convert into MB */
			*ml /= 1024 * 1024;
			/* memlimit is stored into a uint32_t */
			/* so cap the memlimit value to the max value */
			/* of an uint32_t */
			*memlimit = -1 ;
			if ( *ml < *memlimit ) {
				*memlimit = *ml;
			}
		}
		debug3("memory limit of '%s' is %uMB",cpath,*memlimit);
		xfree(ml);
	}

	return fstatus;
}

int xcgroup_set_memswlimit(char* cpath,uint32_t memlimit)
{
	int fstatus;
	char file_path[PATH_MAX];
	uint64_t ml;
	
	fstatus = XCGROUP_ERROR;

	if ( snprintf(file_path,PATH_MAX,"%s/memory.memsw.limit_in_bytes",
		      cpath) >= PATH_MAX ) {
		debug2("unable to set memsw limit of '%s' : %m",cpath);
		return fstatus;
	}

	ml = (uint64_t) memlimit * 1024 * 1024;
	fstatus = _file_write_uint64s(file_path,&ml,1);
	if ( fstatus != XCGROUP_SUCCESS )
		debug2("unable to set memsw limit of '%s' : %m",cpath);
	else
		debug3("mem+swap limit set to %uMB for '%s'",memlimit,cpath);

	return fstatus;
}

int xcgroup_get_memswlimit(char* cpath,uint32_t* memlimit)
{
	int fstatus;
	char file_path[PATH_MAX];
	uint64_t *ml;
	int i;
	
	fstatus = XCGROUP_ERROR;

	if ( snprintf(file_path,PATH_MAX,"%s/memory.memsw.limit_in_bytes",
		      cpath) >= PATH_MAX ) {
		debug2("unable to get memsw limit of '%s' : %m",cpath);
		return fstatus;
	}

	fstatus = _file_read_uint64s(file_path,&ml,&i);
	if ( fstatus != XCGROUP_SUCCESS || 
	     i ==0 )
		debug2("unable to get memsw limit of '%s' : %m",cpath);
	else {
		if ( *ml == 0 ) {
			*memlimit = 0;	
		}
		else {
			/* convert into MB */
			*ml /= 1024 * 1024;
			/* memlimit is stored into a uint32_t */
			/* so cap the memlimit value to the max value */
			/* of an uint32_t */
			*memlimit = -1 ;
			if ( *ml < *memlimit ) {
				*memlimit = *ml;
			}
		}
		debug3("mem+swap limit of '%s' is %uMB",cpath,*memlimit);
		xfree(ml);
	}

	return fstatus;
}

int xcgroup_set_mem_use_hierarchy(char* cpath,int flag)
{
	if ( flag )
		return xcgroup_set_params(cpath,"memory.use_hierarchy=1");
	else
		return xcgroup_set_params(cpath,"memory.use_hierarchy=0");
}

int xcgroup_set_cpuset_cpus(char* cpath,char* range)
{
	int fstatus;
	char file_path[PATH_MAX];
	
	fstatus = XCGROUP_ERROR;
	
	if ( snprintf(file_path,PATH_MAX,"%s/cpuset.cpus",
		      cpath) >= PATH_MAX ) {
		debug2("unable to set cpuset.cpus to '%s' for '%s' : %m",
		       range,cpath);
		return fstatus;
	}
	
	fstatus = _file_write_content(file_path,range,strlen(range));
	if ( fstatus != XCGROUP_SUCCESS )
		debug2("unable to set cpuset.cpus to '%s' for '%s' : %m",
		       range,cpath);
	else
		debug3("cpuset.cpus set to '%s' for '%s'",range,cpath);
	
	return fstatus;
}

int xcgroup_set_params(char* cpath,char* parameters)
{
	int fstatus;
	char file_path[PATH_MAX];
	char* params;
	char* value;
	char* p;
	char* next;

	fstatus = XCGROUP_ERROR;

	params = (char*) xstrdup(parameters);

	p = params;
	while ( p != NULL && *p != '\0' ) {
		next = index(p,' ');
		if ( next ) {
			*next='\0';
			next++;
			while ( *next == ' ' )
				next++;
		}
		value = index(p,'=');
		if ( value != NULL ) {
			*value='\0';
			value++;
			if ( snprintf(file_path,PATH_MAX,"%s/%s",cpath,p) 
			     >= PATH_MAX ) {
				debug2("unable to build filepath for '%s' and"
				       " parameter '%s' : %m",cpath,p);
				goto next_loop;
			}
			fstatus = _file_write_content(file_path,value,
						      strlen(value));
			if ( fstatus != XCGROUP_SUCCESS )
				debug2("unable to set parameter '%s' to "
				       "'%s' for '%s'",p,value,cpath);
			else
				debug3("parameter '%s' set to '%s' for '%s'",
				       p,value,cpath);
		}
		else
			debug2("bad paramters format for entry '%s'",p);
	next_loop:
		p = next;
	}

	xfree(params);

	return fstatus;
}

int xcgroup_get_param(char* cpath,char* parameter,char **content,size_t *csize)
{
	int fstatus;
	char file_path[PATH_MAX];

	fstatus = XCGROUP_ERROR;

	if ( snprintf(file_path,PATH_MAX,"%s/%s",cpath,parameter) 
	     >= PATH_MAX ) {
		debug2("unable to build filepath for '%s' and"
		       " parameter '%s' : %m",cpath,parameter);
	}
	else {
		fstatus = _file_read_content(file_path,content,csize);
		if ( fstatus != XCGROUP_SUCCESS )
			debug2("unable to get parameter '%s'");
	}

	return fstatus;
}


size_t _file_getsize(int fd)
{
	int rc;
	size_t fsize;
	off_t offset;
	char c;

	/* store current position and rewind */
	offset = lseek(fd,0,SEEK_CUR);
	if ( offset < 0 )
		return -1;
	lseek(fd,0,SEEK_SET);

	/* get file size */
	fsize=0;
	do {
		rc = read(fd,(void*)&c,1);
		if ( rc > 0 )
			fsize++;
	}
	while ( (rc < 0 && errno == EINTR) || rc > 0 );

	/* restore position */
	lseek(fd,offset,SEEK_SET);

	if ( rc < 0 )
		return -1;
	else
		return fsize;
}

int
_file_write_uint64s(char* file_path,uint64_t* values,int nb)
{
	int fstatus;
	int rc;
	int fd;
	char tstr[256];
	uint64_t value;
	int i;
	
	/* open file for writing */
	fd = open(file_path, O_WRONLY, 0700);
	if (fd < 0) {
		debug2("unable to open '%s' for writing : %m",
		       file_path);
		return XCGROUP_ERROR;
	}

	/* add one value per line */
	fstatus = XCGROUP_SUCCESS;
	for ( i=0 ; i < nb ; i++ ) {
		
		value = values[i];
		
		rc = snprintf(tstr, sizeof(tstr), "%llu",
			      (long long unsigned int)value);
		if ( rc < 0 ) {
			debug2("unable to build %llu string value, skipping",
			       value);
			fstatus = XCGROUP_ERROR;
			continue;
		}

		do {
			rc = write(fd, tstr, strlen(tstr)+1);
		}
		while ( rc != 0 && errno == EINTR); 
		if (rc < 1) {
			debug2("unable to add value '%s' to file '%s' : %m",
			       tstr,file_path);
			fstatus = XCGROUP_ERROR;
		}

	}

	/* close file */
	close(fd);

	return fstatus;
}


int
_file_read_uint64s(char* file_path,uint64_t** pvalues,int* pnb)
{
	int rc;
	int fd;

	size_t fsize;
	char* buf;
	char* p;

	uint64_t* pa=NULL;
	int i;
	
	/* check input pointers */
	if ( pvalues == NULL || pnb == NULL )
		return XCGROUP_ERROR;
	
	/* open file for reading */
	fd = open(file_path, O_RDONLY, 0700);
	if (fd < 0) {
		debug2("unable to open '%s' for reading : %m",
		       file_path);
		return XCGROUP_ERROR;
	}

	/* get file size */
	fsize=_file_getsize(fd);
	if ( fsize == -1 ) {
		close(fd);
		return XCGROUP_ERROR;
	}

	/* read file contents */
	buf = (char*) xmalloc((fsize+1)*sizeof(char));
	do {
		rc = read(fd,buf,fsize);
	}
	while ( rc < 0 && errno == EINTR );
	close(fd);
	buf[fsize]='\0';

	/* count values (splitted by \n) */
	i=0;
	if ( rc > 0 ) {
		p = buf;
		while ( index(p,'\n') != NULL ) {
			i++;
			p = index(p,'\n') + 1;
		}
	}

	/* build uint32_t list */
	if ( i > 0 ) {
		pa = (uint64_t*) xmalloc(sizeof(uint64_t) * i);
		p = buf;
		i = 0;
		while ( index(p,'\n') != NULL ) {
			long long unsigned int ll_tmp;
			sscanf(p,"%llu",&ll_tmp);
			pa[i++] = ll_tmp;
			p = index(p,'\n') + 1;
		}
	}
	
	/* free buffer */
	xfree(buf);

	/* set output values */
	*pvalues = pa;
	*pnb = i;

	return XCGROUP_SUCCESS;
}

int
_file_write_uint32s(char* file_path,uint32_t* values,int nb)
{
	int fstatus;
	int rc;
	int fd;
	char tstr[256];
	uint32_t value;
	int i;
	
	/* open file for writing */
	fd = open(file_path, O_WRONLY, 0700);
	if (fd < 0) {
		debug2("unable to open '%s' for writing : %m",
		       file_path);
		return XCGROUP_ERROR;
	}

	/* add one value per line */
	fstatus = XCGROUP_SUCCESS;
	for ( i=0 ; i < nb ; i++ ) {
		
		value = values[i];
		
		rc = snprintf(tstr, sizeof(tstr), "%u",value);
		if ( rc < 0 ) {
			debug2("unable to build %u string value, skipping",
			       value);
			fstatus = XCGROUP_ERROR;
			continue;
		}

		do {
			rc = write(fd, tstr, strlen(tstr)+1);
		}
		while ( rc != 0 && errno == EINTR); 
		if (rc < 1) {
			debug2("unable to add value '%s' to file '%s' : %m",
			       tstr,file_path);
			fstatus = XCGROUP_ERROR;
		}

	}

	/* close file */
	close(fd);

	return fstatus;
}


int
_file_read_uint32s(char* file_path,uint32_t** pvalues,int* pnb)
{
	int rc;
	int fd;

	size_t fsize;
	char* buf;
	char* p;

	uint32_t* pa=NULL;
	int i;
	
	/* check input pointers */
	if ( pvalues == NULL || pnb == NULL )
		return XCGROUP_ERROR;
	
	/* open file for reading */
	fd = open(file_path, O_RDONLY, 0700);
	if (fd < 0) {
		debug2("unable to open '%s' for reading : %m",
		       file_path);
		return XCGROUP_ERROR;
	}

	/* get file size */
	fsize=_file_getsize(fd);
	if ( fsize == -1 ) {
		close(fd);
		return XCGROUP_ERROR;
	}

	/* read file contents */
	buf = (char*) xmalloc((fsize+1)*sizeof(char));
	do {
		rc = read(fd,buf,fsize);
	}
	while ( rc < 0 && errno == EINTR );
	close(fd);
	buf[fsize]='\0';

	/* count values (splitted by \n) */
	i=0;
	if ( rc > 0 ) {
		p = buf;
		while ( index(p,'\n') != NULL ) {
			i++;
			p = index(p,'\n') + 1;
		}
	}

	/* build uint32_t list */
	if ( i > 0 ) {
		pa = (uint32_t*) xmalloc(sizeof(uint32_t) * i);
		p = buf;
		i = 0;
		while ( index(p,'\n') != NULL ) {
			sscanf(p,"%u",pa+i);
			p = index(p,'\n') + 1;
			i++;
		}
	}
	
	/* free buffer */
	xfree(buf);

	/* set output values */
	*pvalues = pa;
	*pnb = i;

	return XCGROUP_SUCCESS;
}

int
_file_write_content(char* file_path,char* content,size_t csize)
{
	int fstatus;
	int rc;
	int fd;
	
	/* open file for writing */
	fd = open(file_path, O_WRONLY, 0700);
	if (fd < 0) {
		debug2("unable to open '%s' for writing : %m",
		       file_path);
		return XCGROUP_ERROR;
	}

	/* write content */
	do {
		rc = write(fd,content,csize);
	}
	while ( rc != 0 && errno == EINTR); 

	/* check read size */
	if (rc < csize) {
		debug2("unable to write %u bytes to file '%s' : %m",
		       csize,file_path);
		fstatus = XCGROUP_ERROR;
	}
	else
		fstatus = XCGROUP_SUCCESS;

	/* close file */
	close(fd);
	
	return fstatus;
}

int
_file_read_content(char* file_path,char** content,size_t *csize)
{
	int fstatus;
	int rc;
	int fd;

	size_t fsize;
	char* buf;

	fstatus = XCGROUP_ERROR;

	/* check input pointers */
	if ( content == NULL || csize == NULL )
		return fstatus;
	
	/* open file for reading */
	fd = open(file_path, O_RDONLY, 0700);
	if (fd < 0) {
		debug2("unable to open '%s' for reading : %m",
		       file_path);
		return fstatus;
	}

	/* get file size */
	fsize=_file_getsize(fd);
	if ( fsize == -1 ) {
		close(fd);
		return fstatus;
	}

	/* read file contents */
	buf = (char*) xmalloc((fsize+1)*sizeof(char));
	buf[fsize]='\0';
	do {
		rc = read(fd,buf,fsize);
	}
	while ( rc < 0 && errno == EINTR );

	/* set output values */
	if ( rc >= 0 ) {
		*content = buf;
		*csize = rc;
		fstatus = XCGROUP_SUCCESS;
	}

	/* close file */
	close(fd);

	return fstatus;
}


int _xcgroup_cpuset_init(char* file_path)
{
	int fstatus;
	char path[PATH_MAX];

	char* cpuset_metafiles[] = { 
		"cpuset.cpus",
		"cpuset.mems"
	};
	char* cpuset_meta;
	char* cpuset_conf;
	size_t csize;

	int i;

	fstatus = XCGROUP_ERROR;

	/* when cgroups are configured with cpuset, at least 
	 * cpuset.cpus and cpuset.mems must be set or the cgroup 
	 * will not be available at all.
	 * we duplicate the ancestor configuration in the init step */
	for ( i = 0 ; i < 2 ; i++ ) { 

		cpuset_meta = cpuset_metafiles[i];

		/* try to read ancestor configuration */
		if ( snprintf(path,PATH_MAX,"%s/../%s",
			      file_path,cpuset_meta) >= PATH_MAX ) {
			debug2("unable to get ancestor %s for cgroup '%s' : %m",
			       cpuset_meta,file_path);
			return fstatus;
		}
		if ( _file_read_content(path,&cpuset_conf,&csize) !=
		     XCGROUP_SUCCESS ) {
			debug3("assuming no cpuset support for '%s'",path);
			return XCGROUP_SUCCESS;
		}
		
		/* duplicate ancestor conf in current cgroup */
		if ( snprintf(path,PATH_MAX,"%s/%s",
			      file_path,cpuset_meta) >= PATH_MAX ) {
			debug2("unable to set %s for cgroup '%s' : %m",
			       cpuset_meta,file_path);
			return fstatus;
		}
		if ( _file_write_content(path,cpuset_conf,csize) !=
		     XCGROUP_SUCCESS ) {
			debug2("unable to write %s configuration (%s) of '%s'",
			       cpuset_meta,cpuset_conf,file_path);
			return fstatus;
		}

	}

	return XCGROUP_SUCCESS;
}
