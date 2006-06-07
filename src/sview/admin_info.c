/*****************************************************************************\
 *  admin_info.c - Functions related to admin display 
 *  mode of sview.
 *****************************************************************************
 *  Copyright (C) 2004-2006 The Regents of the University of California.
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

#include "src/sview/sview.h"

#define _DEBUG 0
DEF_TIMERS;

enum { 
	SORTID_POS = POS_LOC,
	SORTID_PARTITION, 
	SORTID_AVAIL, 
	SORTID_TIMELIMIT, 
	SORTID_NODES, 
	SORTID_NODELIST, 
	SORTID_PARTITION_CNT
};

static display_data_t display_data_admin[] = {
	{G_TYPE_INT, SORTID_POS, NULL, FALSE, -1},
	{G_TYPE_STRING, SORTID_PARTITION, "PARTITION", TRUE, -1},
	{G_TYPE_STRING, SORTID_AVAIL, "AVAIL", TRUE, -1},
	{G_TYPE_STRING, SORTID_TIMELIMIT, "TIMELIMIT", TRUE, -1},
	{G_TYPE_STRING, SORTID_NODES, "NODES", TRUE, -1},
#ifdef HAVE_BG
	{G_TYPE_STRING, SORTID_NODELIST, "BP_LIST", TRUE, -1},
#else
	{G_TYPE_STRING, SORTID_NODELIST, "NODELIST", TRUE, -1},
#endif
	{G_TYPE_NONE, -1, NULL, FALSE, -1}};
static display_data_t *local_display_data = NULL;

static void _set_up_button(GtkTreeView *tree_view, GdkEventButton *event, 
			    gpointer user_data)
{
	local_display_data->user_data = user_data;
	button_pressed(tree_view, event, local_display_data);
}

extern void get_info_admin(GtkTable *table, display_data_t *display_data)
{
	local_display_data = display_data;	
}


extern void set_fields_admin(GtkMenu *menu)
{
	make_fields_menu(menu, display_data_admin);
}

extern void row_clicked_admin(GtkTreeView *tree_view,
			      GtkTreePath *path,
			      GtkTreeViewColumn *column,
			      gpointer user_data)
{
	/* job_info_msg_t *job_info_ptr = (job_info_msg_t *)user_data; */
/* 	job_info_t *job_ptr = NULL; */
	int line = get_row_number(tree_view, path);
	GtkWidget *popup = NULL;
	GtkWidget *label = NULL;
	char *info = NULL;
	if(line == -1) {
		g_error("problem getting line number");
		return;
	}
	
/* 	part_ptr = &new_part_ptr->partition_array[line]; */
	/* if(!(info = slurm_sprint_partition_info(part_ptr, 0))) { */
/* 		info = xmalloc(100); */
/* 		sprintf(info, "Problem getting partition info for %s",  */
/* 			part_ptr->name); */
/* 	}  */

	popup = gtk_dialog_new();

	label = gtk_label_new(info);
	gtk_box_pack_end(GTK_BOX(GTK_DIALOG(popup)->vbox), 
			   label, TRUE, TRUE, 0);
	xfree(info);
	gtk_widget_show(label);
	
	gtk_widget_show(popup);
	
}

