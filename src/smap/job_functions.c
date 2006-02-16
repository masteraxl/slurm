/*****************************************************************************\
 *  job_functions.c - Functions related to job display mode of smap.
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
 *
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

#include "src/common/uid.h"
#include "src/common/node_select.h"
#include "src/smap/smap.h"

static void _print_header_job(void);
static int  _print_text_job(job_info_t * job_ptr);

extern void get_job()
{
	int error_code = -1, i, j, recs;
	static int printed_jobs = 0;
	static int count = 0;
	static job_info_msg_t *job_info_ptr = NULL, *new_job_ptr = NULL;
	job_info_t job;

	if (job_info_ptr) {
		error_code = slurm_load_jobs(job_info_ptr->last_update,
				&new_job_ptr, 0);
		if (error_code == SLURM_SUCCESS)
			slurm_free_job_info_msg(job_info_ptr);
		else if (slurm_get_errno() == SLURM_NO_CHANGE_IN_DATA) {
			error_code = SLURM_SUCCESS;
			new_job_ptr = job_info_ptr;
		}
	} else
		error_code = slurm_load_jobs((time_t) NULL, &new_job_ptr, 0);

	if (error_code) {
		if (quiet_flag != 1) {
			mvwprintw(ba_system_ptr->text_win,
				ba_system_ptr->ycord, 1,
				"slurm_load_job: %s", 
				slurm_strerror(slurm_get_errno()));
			ba_system_ptr->ycord++;
		}
	}

	if (!params.no_header)
		_print_header_job();

	if (new_job_ptr)
		recs = new_job_ptr->record_count;
	else
		recs = 0;
	
	if(!params.commandline)
		if((text_line_cnt+printed_jobs) > count) 
			text_line_cnt--;
	printed_jobs = 0;
	count = 0;
	for (i = 0; i < recs; i++) {
		job = new_job_ptr->job_array[i];
		
		if ((job.job_state == JOB_COMPLETE)
		    || (job.job_state == JOB_END)
		    || (job.job_state == JOB_FAILED))
			continue;	/* job has completed */

		if (job.node_inx[0] != -1) {
			job.num_nodes = 0;
			j = 0;
			while (job.node_inx[j] >= 0) {
				job.num_nodes +=
				    (job.node_inx[j + 1] + 1) -
				    job.node_inx[j];
				set_grid(job.node_inx[j],
					 job.node_inx[j + 1], count);
				j += 2;
			}
			
			if(!params.commandline) {
				if((count>=text_line_cnt)
				   && (printed_jobs 
				       < (ba_system_ptr->text_win->_maxy-3))) {
					job.num_procs = (int)letters[count%62];
					wattron(ba_system_ptr->text_win,
						COLOR_PAIR(colors[count%6]));
					_print_text_job(&job);
					wattroff(ba_system_ptr->text_win,
						 COLOR_PAIR(colors[count%6]));
					printed_jobs++;
				} 
			} else {
				job.num_procs = (int)letters[count%62];
				wattron(ba_system_ptr->text_win,
					COLOR_PAIR(colors[count%6]));
				_print_text_job(&job);
				wattroff(ba_system_ptr->text_win,
					 COLOR_PAIR(colors[count%6]));
			}
			count++;			
		}
		if(count==128)
			count=0;
	}
		
	for (i = 0; i < recs; i++) {
		job = new_job_ptr->job_array[i];
		
		if (job.job_state != JOB_PENDING)
			continue;	/* job has completed */

		if(!params.commandline) {
			if((count>=text_line_cnt)
			   && (printed_jobs 
			       < (ba_system_ptr->text_win->_maxy-3))) {
				job.nodes = "waiting...";
				job.num_procs = (int) letters[count%62];
				wattron(ba_system_ptr->text_win,
					COLOR_PAIR(colors[count%6]));
				_print_text_job(&job);
				wattroff(ba_system_ptr->text_win,
					 COLOR_PAIR(colors[count%6]));
				printed_jobs++;
			} 
		} else {
			job.nodes = "waiting...";
			job.num_procs = (int) letters[count%62];
			wattron(ba_system_ptr->text_win,
				COLOR_PAIR(colors[count%6]));
			_print_text_job(&job);
			wattroff(ba_system_ptr->text_win,
				 COLOR_PAIR(colors[count%6]));
			printed_jobs++;
		}
		count++;			
		
		if(count==128)
			count=0;
	}

	if (params.commandline && params.iterate)
		printf("\n");

	ba_system_ptr->ycord++;
	
	job_info_ptr = new_job_ptr;
	return;
}

static void _print_header_job(void)
{
	if(!params.commandline) {
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "ID");
		ba_system_ptr->xcord += 3;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "JOBID");
		ba_system_ptr->xcord += 6;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "PARTITION");
		ba_system_ptr->xcord += 10;
#ifdef HAVE_BG
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "BG_BLOCK");
		ba_system_ptr->xcord += 18;
#endif
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "USER");
		ba_system_ptr->xcord += 9;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "NAME");
		ba_system_ptr->xcord += 10;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "ST");
		ba_system_ptr->xcord += 8;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "TIME");
		ba_system_ptr->xcord += 5;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "NODES");
		ba_system_ptr->xcord += 6;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "NODELIST");
		ba_system_ptr->xcord = 1;
		ba_system_ptr->ycord++;
	} else {
		printf("JOBID ");
		printf("PARTITION ");
#ifdef HAVE_BG
		printf("        BG_BLOCK ");
#endif
		printf("    USER ");
		printf("  NAME ");
		printf("ST ");
		printf("      TIME ");
		printf("NODES ");
		printf("NODELIST\n");
	}
}

static int _print_text_job(job_info_t * job_ptr)
{
	time_t time;
	int printed = 0;
	int tempxcord;
	int prefixlen = 0;
	int i = 0;
	int width = 0;
	char time_buf[20];

#ifdef HAVE_BG
	uint16_t quarter = (uint16_t) NO_VAL;
	uint16_t segment = (uint16_t) NO_VAL;
	uint32_t node_cnt = 0;
	select_g_get_jobinfo(job_ptr->select_jobinfo, 
			     SELECT_DATA_QUARTER, 
			     &quarter);
	select_g_get_jobinfo(job_ptr->select_jobinfo, 
			     SELECT_DATA_SEGMENT, 
			     &segment);
	select_g_get_jobinfo(job_ptr->select_jobinfo, 
			     SELECT_DATA_NODE_CNT, 
			     &node_cnt);
	if(!strcasecmp(job_ptr->nodes,"waiting...")) 
		quarter = (uint16_t) NO_VAL;
#endif
	if(!params.commandline) {
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%c", job_ptr->num_procs);
		ba_system_ptr->xcord += 3;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%d", job_ptr->job_id);
		ba_system_ptr->xcord += 6;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%.10s", job_ptr->partition);
		ba_system_ptr->xcord += 10;
#ifdef HAVE_BG
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%.16s", 
			  select_g_sprint_jobinfo(job_ptr->select_jobinfo, 
						  time_buf, 
						  sizeof(time_buf), 
						  SELECT_PRINT_BG_ID));
		ba_system_ptr->xcord += 18;
#endif
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%.8s", 
			  uid_to_string((uid_t) job_ptr->user_id));
		ba_system_ptr->xcord += 9;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%.9s", job_ptr->name);
		ba_system_ptr->xcord += 10;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%.2s",
			  job_state_string_compact(job_ptr->job_state));
		ba_system_ptr->xcord += 2;
		if(!strcasecmp(job_ptr->nodes,"waiting...")) {
			sprintf(time_buf,"0:00:00");
		} else {
			time = ba_system_ptr->now_time - job_ptr->start_time;
			snprint_time(time_buf, sizeof(time_buf), time);
		}
		width = strlen(time_buf);
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord + (10 - width), "%s",
			  time_buf);
		ba_system_ptr->xcord += 11;

#ifdef HAVE_BG
		if(node_cnt >= 1024) {
			i = node_cnt % 1024;
			if(i > 0) {
				i *= 10;
				i /= 1024;
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "%2d.%dk", 
					  node_cnt/1024, i);
			} else {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "%4dk", 
					  node_cnt/1024);
			}	
		} else
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "%5d", 
				  node_cnt);
#else
		mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "%5d", 
				  job_ptr->num_nodes);
#endif
		ba_system_ptr->xcord += 6;

		tempxcord = ba_system_ptr->xcord;
		
		i=0;
		while (job_ptr->nodes[i] != '\0') {
			if ((printed = mvwaddch(ba_system_ptr->text_win,
						ba_system_ptr->ycord, 
						ba_system_ptr->xcord,
						job_ptr->nodes[i])) < 0)
				return printed;
			ba_system_ptr->xcord++;
			width = ba_system_ptr->text_win->_maxx 
				- ba_system_ptr->xcord;
			if (job_ptr->nodes[i] == '[')
				prefixlen = i + 1;
			else if (job_ptr->nodes[i] == ',' 
				 && (width - 9) <= 0) {
				ba_system_ptr->ycord++;
				ba_system_ptr->xcord = tempxcord + prefixlen;
			}
			i++;
		}
		if(quarter != (uint16_t) NO_VAL) {
			if(segment != (uint16_t) NO_VAL) {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, ".%d.%d", 
					  quarter,
					  segment);
				ba_system_ptr->xcord += 4;
			} else {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, ".%d", 
					  quarter);
				ba_system_ptr->xcord += 2;
			}
		}

		ba_system_ptr->xcord = 1;
		ba_system_ptr->ycord++;
	} else {
		printf("%5d ", job_ptr->job_id);
		printf("%9.9s ", job_ptr->partition);
#ifdef HAVE_BG
		printf("%16.16s ", 
		       select_g_sprint_jobinfo(job_ptr->select_jobinfo, 
					       time_buf, 
					       sizeof(time_buf), 
					       SELECT_PRINT_BG_ID));
#endif
		printf("%8.8s ", uid_to_string((uid_t) job_ptr->user_id));
		printf("%6.6s ", job_ptr->name);
		printf("%2.2s ",
		       job_state_string_compact(job_ptr->job_state));
		if(!strcasecmp(job_ptr->nodes,"waiting...")) {
			sprintf(time_buf,"0:00:00");
		} else {
			time = ba_system_ptr->now_time - job_ptr->start_time;
			snprint_time(time_buf, sizeof(time_buf), time);
		}
		
		printf("%10.10s ", time_buf);
#ifdef HAVE_BG		
		if(node_cnt >= 1024) {
			i = node_cnt % 1024;
			if(i > 0) {
				i *= 10;
				i /= 1024;
				printf("%2d.%dk", node_cnt/1024, i);
			} else {
				printf("%4dk ", node_cnt/1024);
			}
		} else
			printf("%5d ", node_cnt);
#else
		printf("%5d ", job_ptr->num_nodes);
#endif
		printf("%s", job_ptr->nodes);
		if(quarter != (uint16_t) NO_VAL) {
			if(segment != (uint16_t) NO_VAL)
				printf(".%d.%d", quarter, segment);
			else
				printf(".%d", quarter);
		}

		printf("\n");
		
	}
	return printed;
}
