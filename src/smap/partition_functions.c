/*****************************************************************************\
 *  partition_functions.c - Functions related to partition display 
 *  mode of smap.
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
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

#include "src/smap/smap.h"
#include "src/common/node_select.h"
#include "src/api/node_select_info.h"

#define _DEBUG 0
DEF_TIMERS;

typedef struct {
	char *bg_user_name;
	char *bg_block_name;
	char *slurm_part_name;
	char *nodes;
	enum connection_type bg_conn_type;
	enum node_use_type bg_node_use;
	rm_partition_state_t state;
	int letter_num;
	List nodelist;
	int size;
	int quarter;	
	int segment;	
	int node_cnt;	
	bool printed;

} db2_block_info_t;

#ifdef HAVE_BG
static List block_list = NULL;
#endif

static char* _convert_conn_type(enum connection_type conn_type);
static char* _convert_node_use(enum node_use_type node_use);
static int _marknodes(db2_block_info_t *block_ptr, int count);
static void _print_header_part(void);
static char *_part_state_str(rm_partition_state_t state);
static int  _print_text_part(partition_info_t *part_ptr, 
			     db2_block_info_t *db2_info_ptr);
#ifdef HAVE_BG
static void _block_list_del(void *object);
static void _nodelist_del(void *object);
static int _list_match_all(void *object, void *key);
static int _in_slurm_partition(List slurm_nodes, List bg_nodes);
static int _print_rest(db2_block_info_t *block_ptr);
static int _addto_node_list(List nodelist, int *start, int *end);
static int _make_nodelist(char *nodes, List nodelist);
#endif

extern void get_slurm_part()
{
	int error_code, i, j, recs, count = 0;
	static partition_info_msg_t *part_info_ptr = NULL;
	static partition_info_msg_t *new_part_ptr = NULL;
	partition_info_t part;

	if (part_info_ptr) {
		error_code = slurm_load_partitions(part_info_ptr->last_update, 
						   &new_part_ptr, SHOW_ALL);
		if (error_code == SLURM_SUCCESS)
			slurm_free_partition_info_msg(part_info_ptr);
		else if (slurm_get_errno() == SLURM_NO_CHANGE_IN_DATA) {
			error_code = SLURM_SUCCESS;
			new_part_ptr = part_info_ptr;
		}
	} else {
		error_code = slurm_load_partitions((time_t) NULL, 
						   &new_part_ptr, SHOW_ALL);
	}
	if (error_code) {
		if (quiet_flag != 1) {
			if(!params.commandline) {
				mvwprintw(ba_system_ptr->text_win,
					  ba_system_ptr->ycord, 1,
					  "slurm_load_partitions: %s",
					  slurm_strerror(slurm_get_errno()));
				ba_system_ptr->ycord++;
			} else {
				printf("slurm_load_partitions: %s",
				       slurm_strerror(slurm_get_errno()));
			}
		}
		return;
	}
	
	if (!params.no_header)
		_print_header_part();
	
	if (new_part_ptr)
		recs = new_part_ptr->record_count;
	else
		recs = 0;
	if (!params.commandline)
		if((recs - text_line_cnt) < (ba_system_ptr->text_win->_maxy-3))
			text_line_cnt--;

	for (i = 0; i < recs; i++) {
		j = 0;
		part = new_part_ptr->partition_array[i];
		
		if (!part.nodes || (part.nodes[0] == '\0'))
			continue;	/* empty partition */
		
				
		while (part.node_inx[j] >= 0) {
			set_grid(part.node_inx[j],
				 part.node_inx[j + 1], count);
			j += 2;
		}
		if(i>=text_line_cnt) {
			part.root_only =
				(int) letters[count%62];
			wattron(ba_system_ptr->text_win,
				COLOR_PAIR(colors[count%6]));
			_print_text_part(&part, NULL);
			wattroff(ba_system_ptr->text_win,
				 COLOR_PAIR(colors[count%6]));
		}
		count++;
			
	}
	if(count==128)
		count=0;
	if (params.commandline && params.iterate)
		printf("\n");
	
	part_info_ptr = new_part_ptr;
	return;
}

extern void get_bg_part()
{
#ifdef HAVE_BG
	int error_code, i, j, recs=0, count = 0, last_count = -1;
	static partition_info_msg_t *part_info_ptr = NULL;
	static partition_info_msg_t *new_part_ptr = NULL;
	static node_select_info_msg_t *bg_info_ptr = NULL;
	static node_select_info_msg_t *new_bg_ptr = NULL;

	partition_info_t part;
	int number, start[BA_SYSTEM_DIMENSIONS], end[BA_SYSTEM_DIMENSIONS];
	db2_block_info_t *block_ptr = NULL;
	ListIterator itr;
	List nodelist = NULL;

	if (part_info_ptr) {
		error_code = slurm_load_partitions(part_info_ptr->last_update, 
						   &new_part_ptr, SHOW_ALL);
		if (error_code == SLURM_SUCCESS)
			slurm_free_partition_info_msg(part_info_ptr);
		else if (slurm_get_errno() == SLURM_NO_CHANGE_IN_DATA) {
			error_code = SLURM_SUCCESS;
			new_part_ptr = part_info_ptr;
		}
	} else {
		error_code = slurm_load_partitions((time_t) NULL, 
						   &new_part_ptr, SHOW_ALL);
	}

	if (error_code) {
		if (quiet_flag != 1) {
			if(!params.commandline) {
				mvwprintw(ba_system_ptr->text_win,
					  ba_system_ptr->ycord, 1,
					  "slurm_load_partitions: %s",
					  slurm_strerror(slurm_get_errno()));
				ba_system_ptr->ycord++;
			} else {
				printf("slurm_load_partitions: %s",
				       slurm_strerror(slurm_get_errno()));
			}
		}
		return;
	}
	if (bg_info_ptr) {
		error_code = slurm_load_node_select(bg_info_ptr->last_update, 
						   &new_bg_ptr);
		if (error_code == SLURM_SUCCESS)
			select_g_free_node_info(&bg_info_ptr);
		else if (slurm_get_errno() == SLURM_NO_CHANGE_IN_DATA) {
			error_code = SLURM_SUCCESS;
			new_bg_ptr = bg_info_ptr;
		}
	} else {
		error_code = slurm_load_node_select((time_t) NULL, 
						    &new_bg_ptr);
	}
	if (error_code) {
		if (quiet_flag != 1) {
			if(!params.commandline) {
				mvwprintw(ba_system_ptr->text_win,
					  ba_system_ptr->ycord, 1,
					  "slurm_load_partitions: %s",
					  slurm_strerror(slurm_get_errno()));
				ba_system_ptr->ycord++;
			} else {
				printf("slurm_load_partitions: %s",
					  slurm_strerror(slurm_get_errno()));
			}
		}
		return;
	}
	if (block_list) {
		/* clear the old list */
		list_delete_all(block_list, _list_match_all, NULL);
	} else {
		block_list = list_create(_block_list_del);
		if (!block_list) {
			fprintf(stderr, "malloc error\n");
			return;
		}
	}
	if (!params.commandline)
		if((new_bg_ptr->record_count - text_line_cnt) 
		   < (ba_system_ptr->text_win->_maxy-3))
			text_line_cnt--;
	
	for (i=0; i<new_bg_ptr->record_count; i++) {
		block_ptr = xmalloc(sizeof(db2_block_info_t));
		list_append(block_list, block_ptr);
		
		block_ptr->bg_block_name 
			= xstrdup(new_bg_ptr->bg_info_array[i].bg_block_id);
		block_ptr->nodes 
			= xstrdup(new_bg_ptr->bg_info_array[i].nodes);
		block_ptr->nodelist = list_create(_nodelist_del);
		_make_nodelist(block_ptr->nodes,block_ptr->nodelist);
		
		block_ptr->bg_user_name 
			= xstrdup(new_bg_ptr->bg_info_array[i].owner_name);
		block_ptr->state 
			= new_bg_ptr->bg_info_array[i].state;
		block_ptr->bg_conn_type 
			= new_bg_ptr->bg_info_array[i].conn_type;
		block_ptr->bg_node_use 
			= new_bg_ptr->bg_info_array[i].node_use;
		block_ptr->quarter 
			= new_bg_ptr->bg_info_array[i].quarter;
		block_ptr->segment 
			= new_bg_ptr->bg_info_array[i].segment;
		block_ptr->node_cnt 
			= new_bg_ptr->bg_info_array[i].node_cnt;
		if(block_ptr->quarter < 1 && block_ptr->segment < 1) {
			last_count++;
			_marknodes(block_ptr, last_count);
		} else 
			block_ptr->letter_num = last_count;
		
		if(block_ptr->bg_conn_type == SELECT_SMALL)
			block_ptr->size = 0;
		
	}
	
	if (!params.no_header)
		_print_header_part();

	if (new_part_ptr)
		recs = new_part_ptr->record_count;
	else
		recs = 0;

	for (i = 0; i < recs; i++) {
		j = 0;
		part = new_part_ptr->partition_array[i];
		
		if (!part.nodes || (part.nodes[0] == '\0'))
			continue;	/* empty partition */
		nodelist = list_create(_nodelist_del);
		_make_nodelist(part.nodes,nodelist);	
		
		if (block_list) {
			itr = list_iterator_create(block_list);
			while ((block_ptr = (db2_block_info_t*) 
				list_next(itr)) != NULL) {
				if(_in_slurm_partition(nodelist,
						       block_ptr->nodelist)) {
					block_ptr->slurm_part_name 
						= xstrdup(part.name);
				}
			}
			list_iterator_destroy(itr);
		}
		list_destroy(nodelist);
	}

	/* Report the BG Blocks */
	if (block_list) {
		itr = list_iterator_create(block_list);
		while ((block_ptr = (db2_block_info_t*) 
			list_next(itr)) != NULL) {
			if (params.commandline)
				block_ptr->printed = 1;
			else {
				if(count>=text_line_cnt)
					block_ptr->printed = 1;
			}
			_print_rest(block_ptr);
			count++;			
		}
		
		list_iterator_destroy(itr);
	}

	
	if (params.commandline && params.iterate)
		printf("\n");

	part_info_ptr = new_part_ptr;
	bg_info_ptr = new_bg_ptr;
#endif /* HAVE_BG */
	return;
}

static int _marknodes(db2_block_info_t *block_ptr, int count)
{
	int j=0;
	int start[BA_SYSTEM_DIMENSIONS];
	int end[BA_SYSTEM_DIMENSIONS];
	int number = 0;
	
	block_ptr->letter_num = count;
	while (block_ptr->nodes[j] != '\0') {
		if ((block_ptr->nodes[j] == '['
		     || block_ptr->nodes[j] == ',')
		    && (block_ptr->nodes[j+8] == ']' 
			|| block_ptr->nodes[j+8] == ',')
		    && (block_ptr->nodes[j+4] == 'x'
			|| block_ptr->nodes[j+4] == '-')) {
			j++;
			number = atoi(block_ptr->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j += 4;
			number = atoi(block_ptr->nodes + j);
			end[X] = number / 100;
			end[Y] = (number % 100) / 10;
			end[Z] = (number % 10);
			j += 3;
			
			if(block_ptr->state != RM_PARTITION_FREE) 
				block_ptr->size += set_grid_bg(start,
							       end,
							       count,
							       1);
			else
				block_ptr->size += set_grid_bg(start, 
							       end, 
							       count, 
							       0);
			if(block_ptr->nodes[j] != ',')
				break;
			j--;
		} else if((block_ptr->nodes[j] < 58 
			   && block_ptr->nodes[j] > 47)) {
					
			number = atoi(block_ptr->nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j+=3;
			block_ptr->size += set_grid_bg(start, 
							start, 
							count, 
							0);
			if(block_ptr->nodes[j] != ',')
				break;
		}
		j++;
	}
	return SLURM_SUCCESS;
}

static void _print_header_part(void)
{
	if(!params.commandline) {
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "ID");
		ba_system_ptr->xcord += 4;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "PARTITION");
		ba_system_ptr->xcord += 10;
	
		if (params.display != BGPART) {
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "AVAIL");
			ba_system_ptr->xcord += 7;
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "TIMELIMIT");
			ba_system_ptr->xcord += 11;
		} else {
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "BG_BLOCK");
			ba_system_ptr->xcord += 18;
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "STATE");
			ba_system_ptr->xcord += 8;
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "USER");
			ba_system_ptr->xcord += 12;
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "CONN");
			ba_system_ptr->xcord += 7;
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "NODE_USE");
			ba_system_ptr->xcord += 10;
		}

		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "NODES");
		ba_system_ptr->xcord += 7;
		mvwprintw(ba_system_ptr->text_win, ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "NODELIST");
		ba_system_ptr->xcord = 1;
		ba_system_ptr->ycord++;
	} else {
		printf("PARTITION ");
		if (params.display != BGPART) {
			printf("AVAIL ");
			printf("TIMELIMIT ");
		} else {
			printf("        BG_BLOCK ");
			printf("STATE ");
			printf("    USER ");
			printf(" CONN ");
			printf(" NODE_USE ");
		}

		printf("NODES ");
		printf("NODELIST\n");	
	}	
}

static char *_part_state_str(rm_partition_state_t state)
{
	static char tmp[16];

#ifdef HAVE_BG
	switch (state) {
		case RM_PARTITION_BUSY: 
			return "BUSY";
		case RM_PARTITION_CONFIGURING:
			return "CONFIG";
		case RM_PARTITION_DEALLOCATING:
			return "DEALLOC";
		case RM_PARTITION_ERROR:
			return "ERROR";
		case RM_PARTITION_FREE:
			return "FREE";
		case RM_PARTITION_NAV:
			return "NAV";
		case RM_PARTITION_READY:
			return "READY";
	}
#endif

	snprintf(tmp, sizeof(tmp), "%d", state);
	return tmp;
}

static int _print_text_part(partition_info_t *part_ptr, 
			    db2_block_info_t *db2_info_ptr)
{
	int printed = 0;
	int tempxcord;
	int prefixlen;
	int i = 0;
	int width = 0;
	char *nodes = NULL, time_buf[20];

	if(!params.commandline) {
		
		mvwprintw(ba_system_ptr->text_win, 
			  ba_system_ptr->ycord,
			  ba_system_ptr->xcord, "%c", 
			  part_ptr->root_only);
		ba_system_ptr->xcord += 4;

		if (part_ptr->name) {
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "%.9s", 
				  part_ptr->name);
			ba_system_ptr->xcord += 10;
			if (params.display != BGPART) {
				if (part_ptr->state_up)
					mvwprintw(ba_system_ptr->text_win, 
						  ba_system_ptr->ycord,
						  ba_system_ptr->xcord, 
						  "up");
				else
					mvwprintw(ba_system_ptr->text_win, 
						  ba_system_ptr->ycord,
						  ba_system_ptr->xcord, 
						  "down");
				ba_system_ptr->xcord += 7;
			
				if (part_ptr->max_time == INFINITE)
					snprintf(time_buf, sizeof(time_buf), 
						 "infinite");
				else {
					snprint_time(time_buf, 
						     sizeof(time_buf), 
						     (part_ptr->max_time 
						      * 60));
				}
			
				width = strlen(time_buf);
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord + (9 - width), 
					  "%s", 
					  time_buf);
				ba_system_ptr->xcord += 11;
			} 
		} else
			ba_system_ptr->xcord += 10;

		if (params.display == BGPART) {
			if (db2_info_ptr) {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "%.16s", 
					  db2_info_ptr->bg_block_name);
				ba_system_ptr->xcord += 18;
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, 
					  _part_state_str(
						  db2_info_ptr->state));
				ba_system_ptr->xcord += 8;
				
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "%.11s", 
					  db2_info_ptr->bg_user_name);
				ba_system_ptr->xcord += 12;
			
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "%.5s", 
					  _convert_conn_type(
						  db2_info_ptr->
						  bg_conn_type));
				ba_system_ptr->xcord += 7;
				mvwprintw(ba_system_ptr->text_win,
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "%.9s",
					  _convert_node_use(
						  db2_info_ptr->bg_node_use));
				ba_system_ptr->xcord += 10;
			} else {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "?");
				ba_system_ptr->xcord += 12;
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "?");
				ba_system_ptr->xcord += 8;
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "?");
				ba_system_ptr->xcord += 12;
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "?");
				ba_system_ptr->xcord += 6;
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, "?");
				ba_system_ptr->xcord += 10;
			}
		}
		if(part_ptr->total_nodes >= 1024) {
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "%4dk", 
				  part_ptr->total_nodes/1024);
		} else	
			mvwprintw(ba_system_ptr->text_win, 
				  ba_system_ptr->ycord,
				  ba_system_ptr->xcord, "%5d", 
				  part_ptr->total_nodes);
		ba_system_ptr->xcord += 7;

		tempxcord = ba_system_ptr->xcord;
		
		if (params.display == BGPART)
			nodes = part_ptr->allow_groups;
		else
			nodes = part_ptr->nodes;
		prefixlen = i;
		while (nodes && nodes[i]) {
			width = ba_system_ptr->text_win->_maxx 
				- ba_system_ptr->xcord;

			if (!prefixlen && nodes[i] == '[' 
			    && nodes[i - 1] == ',')
				prefixlen = i + 1;

			if (nodes[i - 1] == ',' && (width - 12) <= 0) {
				ba_system_ptr->ycord++;
				ba_system_ptr->xcord = tempxcord + prefixlen;
			} else if (ba_system_ptr->xcord >
				   ba_system_ptr->text_win->_maxx) {
				ba_system_ptr->ycord++;
				ba_system_ptr->xcord = tempxcord + prefixlen;
			}


			if ((printed = mvwaddch(ba_system_ptr->text_win,
						ba_system_ptr->ycord, 
						ba_system_ptr->xcord,
						nodes[i])) < 0)
				return printed;
			ba_system_ptr->xcord++;

			i++;
		}
		if((params.display == BGPART) 
		   && (db2_info_ptr->quarter != -1)) {
			if(db2_info_ptr->segment != -1) {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, ".%d.%d", 
					  db2_info_ptr->quarter,
					  db2_info_ptr->segment);
			} else {
				mvwprintw(ba_system_ptr->text_win, 
					  ba_system_ptr->ycord,
					  ba_system_ptr->xcord, ".%d", 
					  db2_info_ptr->quarter);
			}
		}
			
		ba_system_ptr->xcord = 1;
		ba_system_ptr->ycord++;
	} else {
		if (part_ptr->name) {
			printf("%9.9s ", part_ptr->name);
			
			if (params.display != BGPART) {
				if (part_ptr->state_up)
					printf("   up ");
				else
					printf(" down ");
							
				if (part_ptr->max_time == INFINITE)
					snprintf(time_buf, sizeof(time_buf), 
						 "infinite");
				else {
					snprint_time(time_buf, 
						     sizeof(time_buf), 
						     (part_ptr->max_time 
						      * 60));
				}
			
				width = strlen(time_buf);
				printf("%9.9s ", time_buf);		
			} 
		}

		if (params.display == BGPART) {
			if (db2_info_ptr) {
				printf("%16.16s ",
				       db2_info_ptr->bg_block_name);
				printf("%5.5s ", 
				       _part_state_str(db2_info_ptr->state));
				
				printf("%8.8s ", db2_info_ptr->bg_user_name);
				
				printf("%5.5s ", _convert_conn_type(
					       db2_info_ptr->bg_conn_type));
				printf("%9.9s ",  _convert_node_use(
					       db2_info_ptr->bg_node_use));
			} 
		}
		
		if(part_ptr->total_nodes >= 1024)
			printf("%4dk ", part_ptr->total_nodes/1024);
		else	
			printf("%5d ", part_ptr->total_nodes);
		
		tempxcord = ba_system_ptr->xcord;
		
		if (params.display == BGPART)
			nodes = part_ptr->allow_groups;
		else
			nodes = part_ptr->nodes;
		
		if((params.display == BGPART) 
		   && (db2_info_ptr->quarter != -1)) {
			if(db2_info_ptr->segment != -1)
				printf("%s.%d.%d\n", nodes, 
				       db2_info_ptr->quarter,
				       db2_info_ptr->segment);
			else 
				printf("%s.%d\n", nodes, 
				       db2_info_ptr->quarter);
		} else
			printf("%s\n",nodes);
	}
	return printed;
}

#ifdef HAVE_BG
static void _block_list_del(void *object)
{
	db2_block_info_t *block_ptr = (db2_block_info_t *)object;

	if (block_ptr) {
		xfree(block_ptr->bg_user_name);
		xfree(block_ptr->bg_block_name);
		xfree(block_ptr->slurm_part_name);
		xfree(block_ptr->nodes);
		if(block_ptr->nodelist)
			list_destroy(block_ptr->nodelist);
		
		xfree(block_ptr);
		
	}
}

static void _nodelist_del(void *object)
{
	int *coord = (int *)object;
	xfree(coord);
	return;
}

static int _list_match_all(void *object, void *key)
{
	
	return 1;
}

static int _in_slurm_partition(List slurm_nodes, List bg_nodes)
{
	ListIterator slurm_itr;
	ListIterator bg_itr;
	int *coord = NULL;
	int *slurm_coord = NULL;
	int found = 0;
	
	bg_itr = list_iterator_create(bg_nodes);
	slurm_itr = list_iterator_create(slurm_nodes);
	while ((coord = list_next(bg_itr)) != NULL) {
		list_iterator_reset(slurm_itr);
		found = 0;
		while ((slurm_coord = list_next(slurm_itr)) != NULL) {
			if((coord[X] == slurm_coord[X])
			   && (coord[Y] == slurm_coord[Y])
			   && (coord[Z] == slurm_coord[Z])) {
				found=1;
				break;
			}			
		}
		if(!found) {
			break;
		}
	}
	list_iterator_destroy(slurm_itr);
	list_iterator_destroy(bg_itr);
			
	if(found)
		return 1;
	else
		return 0;
	
}

static int _print_rest(db2_block_info_t *block_ptr)
{
	partition_info_t part;
	db2_block_info_t *db2_info_ptr = NULL;
	ListIterator itr;
	int set = 0;
		
	if(block_ptr->node_cnt == 0)
		block_ptr->node_cnt = block_ptr->size;
	part.total_nodes = block_ptr->node_cnt;
	if(block_ptr->slurm_part_name)
		part.name = block_ptr->slurm_part_name;
	else
		part.name = "no part";

	if (!block_ptr->printed)
		return SLURM_SUCCESS;
	part.allow_groups = block_ptr->nodes;
	part.root_only = (int) letters[block_ptr->letter_num%62];
	wattron(ba_system_ptr->text_win, 
		COLOR_PAIR(colors[block_ptr->letter_num%6]));
	_print_text_part(&part, block_ptr);
	wattroff(ba_system_ptr->text_win,
		 COLOR_PAIR(colors[block_ptr->letter_num%6]));
	
	return SLURM_SUCCESS;
}

static int _addto_nodelist(List nodelist, int *start, int *end)
{
	int *coord = NULL;
	int x,y,z;
	
	assert(end[X] < DIM_SIZE[X]);
	assert(start[X] >= 0);
	assert(end[Y] < DIM_SIZE[Y]);
	assert(start[Y] >= 0);
	assert(end[Z] < DIM_SIZE[Z]);
	assert(start[Z] >= 0);
	
	for (x = start[X]; x <= end[X]; x++) {
		for (y = start[Y]; y <= end[Y]; y++) {
			for (z = start[Z]; z <= end[Z]; z++) {
				coord = xmalloc(sizeof(int)*3);
				coord[X] = x;
				coord[Y] = y;
				coord[Z] = z;
				list_append(nodelist, coord);
			}
		}
	}
	return 1;
}

static int _make_nodelist(char *nodes, List nodelist)
{
	int j = 0;
	int number;
	int start[BA_SYSTEM_DIMENSIONS];
	int end[BA_SYSTEM_DIMENSIONS];
	
	if(!nodelist)
		nodelist = list_create(_nodelist_del);
	while (nodes[j] != '\0') {
		if ((nodes[j] == '['
		     || nodes[j] == ',')
		    && (nodes[j+8] == ']' 
			|| nodes[j+8] == ',')
		    && (nodes[j+4] == 'x'
			|| nodes[j+4] == '-')) {
			j++;
			number = atoi(nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j += 4;
			number = atoi(nodes + j);
			end[X] = number / 100;
			end[Y] = (number % 100) / 10;
			end[Z] = (number % 10);
			j += 3;
			_addto_nodelist(nodelist, start, end);
			if(nodes[j] != ',')
				break;
			j--;
		} else if((nodes[j] < 58 
			   && nodes[j] > 47)) {
					
			number = atoi(nodes + j);
			start[X] = number / 100;
			start[Y] = (number % 100) / 10;
			start[Z] = (number % 10);
			j+=3;
			_addto_nodelist(nodelist, start, start);
			if(nodes[j] != ',')
				break;
		}
		j++;
	}
	return 1;
}

#endif

static char* _convert_conn_type(enum connection_type conn_type)
{
	switch (conn_type) {
	case (SELECT_MESH):
		return "MESH";
	case (SELECT_TORUS):
		return "TORUS";
	case (SELECT_SMALL):
		return "SMALL";
	case (SELECT_NAV):
		return "NAV";
	}
	return "?";
}

static char* _convert_node_use(enum node_use_type node_use)
{
	switch (node_use) {
	case (SELECT_COPROCESSOR_MODE):
		return "COPROCESSOR";
	case (SELECT_VIRTUAL_NODE_MODE):
		return "VIRTUAL";
	case (SELECT_NAV_MODE):
		return "NAV";
	}
	return "?";
}

