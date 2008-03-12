/****************************************************************************\
 *  grid.c - put display grid info here
 *****************************************************************************
 *  Copyright (C) 2002-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>, et. al.
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
#include "sview.h"

List grid_button_list = NULL;
List blinking_button_list = NULL;

char *sview_colors[] = {"#0000FF", "#00FF00", "#00FFFF", "#FFFF00", 
			"#FF0000", "#4D4DC6", "#F09A09", "#BDFA19",
			"#715627", "#6A8CA2", "#4C7127", "#25B9B9",
			"#A020F0", "#8293ED", "#FFA500", "#FFC0CB",
			"#8B6914", "#18A24E", "#F827FC", "#B8A40C"};
int sview_colors_cnt = 20;

GStaticMutex blinking_mutex = G_STATIC_MUTEX_INIT;


static void _open_node(GtkWidget *widget, GdkEventButton *event, 
		       grid_button_t *grid_button)
{
	GError *error = NULL;
	char title[100];
	ListIterator itr = NULL;
	popup_info_t *popup_win = NULL;

#ifdef HAVE_BG
	snprintf(title, 100, 
		 "Info about base partition %s", grid_button->node_name);
#else
	snprintf(title, 100, "Info about node %s", grid_button->node_name);
#endif
	itr = list_iterator_create(popup_list);
	while((popup_win = list_next(itr))) {
		if(popup_win->spec_info)
			if(!strcmp(popup_win->spec_info->title, title)) {
				break;
			} 
	}
	list_iterator_destroy(itr);

	if(!popup_win) {
		popup_win = create_popup_info(INFO_PAGE, NODE_PAGE, title);
		popup_win->spec_info->search_info->gchar_data =
			g_strdup(grid_button->node_name);
		if (!g_thread_create((gpointer)popup_thr, popup_win,
				     FALSE, &error))
		{
			g_printerr ("Failed to create grid popup thread: "
				    "%s\n", 
				    error->message);
			return;
		}	
	} else
		gtk_window_present(GTK_WINDOW(popup_win->popup));
	return;
}

#ifdef HAVE_BG
static void _open_block(GtkWidget *widget, GdkEventButton *event, 
			grid_button_t *grid_button)
{
	GError *error = NULL;
	char title[100];
	ListIterator itr = NULL;
	popup_info_t *popup_win = NULL;

	snprintf(title, 100, 
		 "Info about block containing %s", grid_button->node_name);

	itr = list_iterator_create(popup_list);
	while((popup_win = list_next(itr))) {
		if(popup_win->spec_info)
			if(!strcmp(popup_win->spec_info->title, title)) {
				break;
			} 
	}
	list_iterator_destroy(itr);

	if(!popup_win) {
		popup_win = create_popup_info(INFO_PAGE, BLOCK_PAGE, title);
		popup_win->spec_info->search_info->search_type =
			SEARCH_BLOCK_NODENAME;
		popup_win->spec_info->search_info->gchar_data =
			g_strdup(grid_button->node_name);
		if (!g_thread_create((gpointer)popup_thr, popup_win,
				     FALSE, &error))
		{
			g_printerr ("Failed to create block "
				    "grid popup thread: %s\n", 
				    error->message);
			return;
		}	
	} else
		gtk_window_present(GTK_WINDOW(popup_win->popup));
	return;
}
#endif

/* 
 * Comparator used for sorting buttons
 * 
 * returns: -1: button_a->inx > button_b->inx   
 *           0: rec_a == rec_b
 *           1: rec_a < rec_b
 * 
 */
static int _sort_button_inx(grid_button_t *button_a, grid_button_t *button_b)
{
	int inx_a = button_a->inx;
	int inx_b = button_b->inx;
	
	if (inx_a < inx_b)
		return -1;
	else if (inx_a > inx_b)
		return 1;
	return 0;
}

void _put_button_as_down(grid_button_t *grid_button, int state)
{
	GtkWidget *image = NULL;
	GdkColor color;

	if(GTK_IS_EVENT_BOX(grid_button->button)) {
		return;
	}
	gtk_widget_destroy(grid_button->button);		
		
	grid_button->button = gtk_event_box_new();
	gtk_tooltips_set_tip(grid_button->tip,
			     grid_button->button,
			     grid_button->node_name,
			     "click for node stats");
	gtk_widget_set_size_request(grid_button->button, 10, 10);
	gtk_event_box_set_above_child(GTK_EVENT_BOX(grid_button->button),
				      FALSE);
	g_signal_connect(G_OBJECT(grid_button->button),
			 "button-press-event",
			 G_CALLBACK(_open_node),
			 grid_button);
	if(grid_button->table) 
		gtk_table_attach(grid_button->table, grid_button->button,
				 grid_button->table_x,
				 (grid_button->table_x+1), 
				 grid_button->table_y,
				 (grid_button->table_y+1),
				 GTK_SHRINK, GTK_SHRINK,
				 1, 1);
	
	gdk_color_parse("black", &color);
	gtk_widget_modify_bg(grid_button->button, 
			     GTK_STATE_NORMAL, &color);
	gdk_color_parse("white", &color);
	gtk_widget_modify_bg(grid_button->button, 
			     GTK_STATE_PRELIGHT, &color);
	if(state == NODE_STATE_DRAIN)
		image = gtk_image_new_from_stock(GTK_STOCK_DIALOG_ERROR,
						 GTK_ICON_SIZE_SMALL_TOOLBAR);
	else 
		image = gtk_image_new_from_stock(GTK_STOCK_CANCEL,
						 GTK_ICON_SIZE_SMALL_TOOLBAR);
	gtk_container_add(GTK_CONTAINER(grid_button->button), image);
	gtk_widget_show_all(grid_button->button);
	return;	
}

void _put_button_as_up(grid_button_t *grid_button)
{
	if(GTK_IS_BUTTON(grid_button->button)) {
		return;
	}
	gtk_widget_destroy(grid_button->button);		
	grid_button->button = gtk_button_new();
	gtk_widget_set_size_request(grid_button->button, 10, 10);
	gtk_tooltips_set_tip(grid_button->tip,
			     grid_button->button,
			     grid_button->node_name,
			     "click for node stats");
	g_signal_connect(G_OBJECT(grid_button->button), 
			 "button-press-event",
			 G_CALLBACK(_open_node),
			 grid_button);
	if(grid_button->table) 
		gtk_table_attach(grid_button->table, grid_button->button,
				 grid_button->table_x,
				 (grid_button->table_x+1), 
				 grid_button->table_y,
				 (grid_button->table_y+1),
				 GTK_SHRINK, GTK_SHRINK,
				 1, 1);
	gtk_widget_show_all(grid_button->button);
	return;
}

#ifdef HAVE_BG
static int _block_in_node(int *bp_inx, int inx)
{
	int j=0;
	if(bp_inx[j] >= 0) {
		if ((bp_inx[j] == inx)
		    && (bp_inx[j+1] == inx))
			return 1;
	}
	return 0;
}
#endif

extern void destroy_grid_button(void *arg)
{
	grid_button_t *grid_button = (grid_button_t *)arg;
	if(grid_button) {
		if(grid_button->button) {
			gtk_widget_destroy(grid_button->button);
			grid_button->button = NULL;
		}
		xfree(grid_button->node_name);
		xfree(grid_button);
	}
}

/* we don't set the call back for the button here because sometimes we
 * need to get a different call back based on what we are doing with
 * the button, an example of this would be in
 * add_extra_bluegene_buttons were the small block buttons do
 * something different than they do regularly
 */

extern grid_button_t *create_grid_button_from_another(
	grid_button_t *grid_button, char *name, int color_inx)
{
	grid_button_t *send_grid_button = NULL;
	GdkColor color;
	uint16_t node_base_state;

	if(!grid_button || !name)
		return NULL;
	if(color_inx >= 0)
		color_inx %= sview_colors_cnt;
			
	send_grid_button = xmalloc(sizeof(grid_button_t));
	memcpy(send_grid_button, grid_button, sizeof(grid_button_t));
	node_base_state = send_grid_button->state & NODE_STATE_BASE;
	/* need to set the table to empty because we will want to fill
	   this into the new table later */
	send_grid_button->table = NULL;
	if((color_inx >= 0) && node_base_state == NODE_STATE_DOWN) {
		GtkWidget *image = gtk_image_new_from_stock(
			GTK_STOCK_CANCEL,
			GTK_ICON_SIZE_SMALL_TOOLBAR);
		send_grid_button->button = gtk_event_box_new();
		gtk_event_box_set_above_child(
			GTK_EVENT_BOX(send_grid_button->button),
			FALSE);
		gdk_color_parse("black", &color);
		gtk_widget_modify_bg(send_grid_button->button, 
				     GTK_STATE_NORMAL, &color);
		gdk_color_parse("white", &color);
		gtk_widget_modify_bg(send_grid_button->button, 
				     GTK_STATE_PRELIGHT, &color);
		gtk_container_add(
			GTK_CONTAINER(send_grid_button->button),
			image);
	} else if((color_inx >= 0)
		  && send_grid_button->state & NODE_STATE_DRAIN) {
		GtkWidget *image = gtk_image_new_from_stock(
			GTK_STOCK_DIALOG_ERROR,
			GTK_ICON_SIZE_SMALL_TOOLBAR);
		
		send_grid_button->button = gtk_event_box_new();
		gtk_event_box_set_above_child(
			GTK_EVENT_BOX(send_grid_button->button),
			FALSE);
		gdk_color_parse("black", &color);
		gtk_widget_modify_bg(send_grid_button->button, 
				     GTK_STATE_NORMAL, &color);
		gdk_color_parse("white", &color);
		gtk_widget_modify_bg(send_grid_button->button, 
				     GTK_STATE_PRELIGHT, &color);
		gtk_container_add(
			GTK_CONTAINER(send_grid_button->button),
			image);		
	} else {
		send_grid_button->button = gtk_button_new();
		if(color_inx >= 0)
			gdk_color_parse(sview_colors[color_inx], &color);
		else
			gdk_color_parse("white", &color);
		gtk_widget_modify_bg(send_grid_button->button, 
				     GTK_STATE_NORMAL, &color);
	}
	gtk_widget_set_size_request(send_grid_button->button, 10, 10);
	send_grid_button->tip = gtk_tooltips_new();
	
	send_grid_button->node_name = xstrdup(name);
	
	gtk_tooltips_set_tip(send_grid_button->tip,
			     send_grid_button->button,
			     send_grid_button->node_name,
			     "click for node stats");
	if(color_inx >= 0)
		send_grid_button->color = sview_colors[color_inx];
	else
		send_grid_button->color = "white";
	return send_grid_button;
}

extern char *change_grid_color(List button_list, int start, int end,
			       int color_inx)
{
	ListIterator itr = NULL;
	grid_button_t *grid_button = NULL;
	uint16_t node_base_state;
	GdkColor color;

	if(!button_list)
		return NULL;

	itr = list_iterator_create(button_list);
	color_inx %= sview_colors_cnt;
	gdk_color_parse(sview_colors[color_inx], &color);
	while((grid_button = list_next(itr))) {
		if ((grid_button->inx < start)
		    ||  (grid_button->inx > end)) 
			continue;
		node_base_state = grid_button->state & NODE_STATE_BASE;
		if (node_base_state == NODE_STATE_DOWN) {
			_put_button_as_down(grid_button, NODE_STATE_DOWN);
		} else if (grid_button->state & NODE_STATE_DRAIN) {
			_put_button_as_down(grid_button, NODE_STATE_DRAIN);
		} else {
			_put_button_as_up(grid_button);
			grid_button->color = sview_colors[color_inx];
			gtk_widget_modify_bg(grid_button->button, 
					     GTK_STATE_NORMAL, &color);
		}
	}
	list_iterator_destroy(itr);
	return sview_colors[color_inx];
}

extern void get_button_list_from_main(List *button_list, int start, int end,
				      int color_inx)
{
	ListIterator itr = NULL;
	ListIterator button_itr = NULL;
	grid_button_t *grid_button = NULL;
	grid_button_t *send_grid_button = NULL;
	
	if(!*button_list)
		*button_list = list_create(destroy_grid_button);
	
	color_inx %= sview_colors_cnt;
	itr = list_iterator_create(grid_button_list);
	while((grid_button = list_next(itr))) {
		if ((grid_button->inx < start)
		    ||  (grid_button->inx > end)) 
			continue;
		button_itr = list_iterator_create(*button_list);
		while((send_grid_button = list_next(button_itr))) {
			if(send_grid_button->inx == grid_button->inx)
				break;
		}
		list_iterator_destroy(button_itr);
		if(send_grid_button)
			continue;
		
		send_grid_button = create_grid_button_from_another(
			grid_button, grid_button->node_name, color_inx);
		if(send_grid_button) {
			g_signal_connect(G_OBJECT(send_grid_button->button),
			 "button-press-event",
			 G_CALLBACK(_open_node),
			 send_grid_button);
			list_append(*button_list, send_grid_button);
		}
	}
	list_iterator_destroy(itr);
	return;
}

extern List copy_main_button_list()
{
	ListIterator itr = NULL;
	grid_button_t *grid_button = NULL;
	grid_button_t *send_grid_button = NULL;
	List button_list = list_create(destroy_grid_button);
	
	itr = list_iterator_create(grid_button_list);
	while((grid_button = list_next(itr))) {
		send_grid_button = create_grid_button_from_another(
			grid_button, grid_button->node_name, -1);
		if(send_grid_button) {
			g_signal_connect(G_OBJECT(send_grid_button->button),
					 "button-press-event",
					 G_CALLBACK(_open_node),
					 send_grid_button);
			list_append(button_list, send_grid_button);
		}
	}
	list_iterator_destroy(itr);
	return button_list;
}

#ifdef HAVE_BG
extern void add_extra_bluegene_buttons(List *button_list, int inx, 
				       int *color_inx)
{
	node_select_info_msg_t *node_select_ptr = NULL;
	bg_info_record_t *bg_info_ptr = NULL;
	int error_code = SLURM_SUCCESS;
	ListIterator itr = NULL;
	grid_button_t *grid_button = NULL;
	grid_button_t *send_grid_button = NULL;
	int i=0;
	char *nodes = NULL;
	char tmp_nodes[256];
	int found = 0;
	int coord_y=0;
	
	error_code = get_new_info_node_select(&node_select_ptr, 0);
	
	if (error_code != SLURM_SUCCESS 
	    && error_code != SLURM_NO_CHANGE_IN_DATA) {
		return;
	}
	
	if(!*button_list) 
		*button_list = list_create(NULL);
	
	*color_inx %= sview_colors_cnt;
	
	itr = list_iterator_create(grid_button_list);
	while((grid_button = list_next(itr))) {
		if (grid_button->inx == inx)
			break;
	}
	list_iterator_destroy(itr);
	
	if(!grid_button)
		return;

	/* remove all (if any) buttons pointing to this node since we
	   will be creating all of them */

	itr = list_iterator_create(*button_list);
	while((send_grid_button = list_next(itr))) {
		if(send_grid_button->inx == grid_button->inx)
			list_remove(itr);
	}
	list_iterator_destroy(itr);
	
	for (i=0; i < node_select_ptr->record_count; i++) {
		bg_info_ptr = &node_select_ptr->bg_info_array[i];
		if(!_block_in_node(bg_info_ptr->bp_inx, inx))
			continue;
		found = 1;
		nodes = bg_info_ptr->nodes;			
		if(bg_info_ptr->ionodes) {
			sprintf(tmp_nodes, "%s[%s]", nodes,
				bg_info_ptr->ionodes);
			nodes = tmp_nodes;
		}
		send_grid_button = create_grid_button_from_another(
			grid_button, nodes, *color_inx);
		if(send_grid_button) {
			send_grid_button->table_x = 0;
			send_grid_button->table_y = coord_y++;
			g_signal_connect(
				G_OBJECT(send_grid_button->button),
				"button-press-event",
				G_CALLBACK(_open_block),
				send_grid_button);
			
			list_append(*button_list, send_grid_button);
			(*color_inx)++;
		}
	}
	if(!found) {
		send_grid_button = create_grid_button_from_another(
			grid_button, grid_button->node_name, *color_inx);
		if(send_grid_button) {
			send_grid_button->table_x = 0;
			send_grid_button->table_y = coord_y++;
			g_signal_connect(
				G_OBJECT(send_grid_button->button),
				"button-press-event",
				G_CALLBACK(_open_node),
				send_grid_button);
			
			list_append(*button_list, send_grid_button);
			(*color_inx)++;
		}
	}

}
#endif

extern void add_extra_cr_buttons(List *button_list, node_info_t *node_ptr)
{
	/* FIXME: this is here for consumable resources "multi-core"
	   and what not to add buttons for each.  This needs to be added
	   when HP is done with the multi-core code. */
	return;
}

extern void put_buttons_in_table(GtkTable *table, List button_list)
{
	int table_x=0, table_y=0;
#ifndef HAVE_BG
	int coord_x=0, coord_y=0;
#endif
	grid_button_t *grid_button = NULL;
	ListIterator itr = NULL;
	int node_count = list_count(button_list);
	
	list_sort(button_list, (ListCmpF) _sort_button_inx);
	
#ifdef HAVE_BG
	node_count = DIM_SIZE[X];
	table_x = DIM_SIZE[X] + DIM_SIZE[Z];
	table_y = (DIM_SIZE[Z] * DIM_SIZE[Y]) + DIM_SIZE[Y];
#else
	if(node_count < 50) {
		table_x = 1;
	} else if(node_count < 500) {
		table_x = 10;
	} else {
		table_x=20;
	}
	table_y = node_count/table_x;
	table_y++;
#endif
	//g_print("the table size is y=%d x=%d\n", table_y, table_x);
	gtk_table_resize(table, table_y, table_x);
	itr = list_iterator_create(button_list);
	while((grid_button = list_next(itr))) {
#ifdef HAVE_BG
		grid_button->table = table;
		gtk_table_attach(table, grid_button->button,
				 grid_button->table_x, 
				 (grid_button->table_x+1),
				 grid_button->table_y,
				 (grid_button->table_y+1),
				 GTK_SHRINK, GTK_SHRINK,
				 1, 1);
		if(!grid_button->table_x)
			gtk_table_set_row_spacing(table, 
						  grid_button->table_y, 5);
#else
		grid_button->table = table;
		grid_button->table_x = coord_x;
		grid_button->table_y = coord_y;
		gtk_table_attach(table, grid_button->button,
				 coord_x, (coord_x+1), coord_y, (coord_y+1),
				 GTK_SHRINK, GTK_SHRINK,
				 1, 1);
		
		coord_x++;
		
		if(coord_x == table_x) {
			coord_x = 0;
			coord_y++;
			if(!(coord_y%10)) {
				gtk_table_set_row_spacing(
					table, coord_y-1, 5);
			}
			
		}
		
		if(coord_y == table_y)
			break;
		
		if(coord_x && !(coord_x%10)) {
			gtk_table_set_col_spacing(table,
						  coord_x-1,
						  5);
		}
#endif
	}
	list_iterator_destroy(itr);

	gtk_widget_show_all(GTK_WIDGET(table));
}

extern int get_system_stats(GtkTable *table)
{
	int error_code = SLURM_SUCCESS;
	node_info_msg_t *node_info_ptr = NULL;
	List node_list = NULL;
	int changed = 1;
	static GtkWidget *label = NULL;
	char error_char[100];
	
	if(label)
		gtk_widget_destroy(label);

	if((error_code = get_new_info_node(&node_info_ptr, force_refresh))
	   == SLURM_NO_CHANGE_IN_DATA) { 
		changed = 0;
	} else if (error_code != SLURM_SUCCESS) {		
		snprintf(error_char, 100, "slurm_load_node: %s\n",
			 slurm_strerror(error_code));
		label = gtk_label_new(error_char);
		gtk_table_attach_defaults(table, label, 0, 1, 0, 1);
		gtk_widget_show(label);
		
		return SLURM_ERROR;
	}

	ba_init(node_info_ptr);	
	if(grid_button_list)
		return SLURM_SUCCESS;
	grid_button_list = list_create(destroy_grid_button);
	node_list = create_node_info_list(node_info_ptr, changed);
	setup_grid_table(main_grid_table, grid_button_list, node_list);
	gtk_widget_show_all(GTK_WIDGET(main_grid_table));
	return SLURM_SUCCESS;
}

extern int setup_grid_table(GtkTable *table, List button_list, List node_list)
{
	int error_code = SLURM_SUCCESS;
	int x=0, table_x=0, table_y=0;
	int coord_x=0, coord_y=0, i=0;
	grid_button_t *grid_button = NULL;
	int node_count = 0;
	ListIterator itr = NULL;
	sview_node_info_t *sview_node_info_ptr = NULL;
#ifdef HAVE_BG
	int y=0, z=0, x_offset=0, y_offset=0;
#endif

	if(!node_list) {
		g_print("setup_grid_table: no node_list given\n");
		return SLURM_ERROR;
	}
	
#ifdef HAVE_BG
	node_count = DIM_SIZE[X];
	table_x = DIM_SIZE[X] + DIM_SIZE[Z];
	table_y = (DIM_SIZE[Z] * DIM_SIZE[Y]) + DIM_SIZE[Y];
#else
	node_count = list_count(node_list);
	if(node_count < 50) {
		table_x = 1;
	} else if(node_count < 500) {
		table_x = 10;
	} else {
		table_x=20;
	}
	table_y = node_count/table_x;
	table_y++;
#endif

	gtk_table_resize(table, table_y, table_x);
	itr = list_iterator_create(node_list);
#ifdef HAVE_BG
	/* ok this is going to look much different than smap since we
	 * get the nodes from the controller going up from the Z dim
	 * instead of laying these out in a nice X fashion
	 */
	for (x=0; x<DIM_SIZE[X]; x++) {
		y_offset = (DIM_SIZE[Z] * DIM_SIZE[Y]);
		for (y=0; y<DIM_SIZE[Y]; y++) {
			coord_y = y_offset - y;
			x_offset = DIM_SIZE[Z] - 1;
			for (z=0; z<DIM_SIZE[Z]; z++){
				coord_x = x + x_offset;
				
				grid_button = xmalloc(sizeof(grid_button_t));
				grid_button->inx = i++;
				grid_button->table = table;
				grid_button->table_x = coord_x;
				grid_button->table_y = coord_y;
				grid_button->button = gtk_button_new();
				grid_button->tip = gtk_tooltips_new();
				if(!(sview_node_info_ptr = list_next(itr))) {
					g_print("no node for this "
						"inx %d!!!!\n",
						grid_button->inx);
					goto end_it;
				}
			
				grid_button->node_name = xstrdup(
					sview_node_info_ptr->node_ptr->name);
				
				gtk_tooltips_set_tip(grid_button->tip,
						     grid_button->button,
						     grid_button->node_name,
						     "click for node stats");
				gtk_widget_set_size_request(
					grid_button->button, 10, 10);
				g_signal_connect(G_OBJECT(grid_button->button),
						 "button-press-event",
						 G_CALLBACK(_open_node),
						 grid_button);
				list_append(button_list, grid_button);
								
				gtk_table_attach(table, grid_button->button,
						 coord_x, (coord_x+1),
						 coord_y, (coord_y+1),
						 GTK_SHRINK, GTK_SHRINK,
						 1, 1);
				
				coord_y++;
				x_offset--;
			}
			y_offset -= DIM_SIZE[Z];			
		}
		gtk_table_set_row_spacing(table, coord_y-1, 5);
	}
#else
	for (x=0; x<node_count; x++) {
		grid_button = xmalloc(sizeof(grid_button_t));
		grid_button->inx = i++;
		grid_button->table = table;
		grid_button->table_x = coord_x;
		grid_button->table_y = coord_y;
		
		grid_button->button = gtk_button_new();
		grid_button->tip = gtk_tooltips_new();
		if(!(sview_node_info_ptr = list_next(itr))) {
			g_print("no node for this "
				"inx %d!!!!\n",
				grid_button->inx);
			goto end_it;
		}
		grid_button->node_name = xstrdup(
			sview_node_info_ptr->node_ptr->name);
		gtk_tooltips_set_tip(grid_button->tip,
				     grid_button->button,
				     grid_button->node_name,
				     "click for node stats");
		gtk_widget_set_size_request(grid_button->button, 10, 10);
		g_signal_connect(G_OBJECT(grid_button->button), 
				 "button-press-event",
				 G_CALLBACK(_open_node),
				 grid_button);
		list_append(button_list, grid_button);
		
		gtk_table_attach(table, grid_button->button,
				 coord_x, (coord_x+1), coord_y, (coord_y+1),
				 GTK_SHRINK, GTK_SHRINK,
				 1, 1);
		
		coord_x++;
			
		if(coord_x == table_x) {
			coord_x = 0;
			coord_y++;
			if(!(coord_y%10)) {
				gtk_table_set_row_spacing(table,
							  coord_y-1, 5);
			}
			
		}
		
		if(coord_y == table_y)
			break;
	
		if(coord_x && !(coord_x%10)) {
			gtk_table_set_col_spacing(table,
						  coord_x-1, 5);
		}
	}
#endif
	
end_it:
	list_iterator_destroy(itr);
	list_sort(button_list, (ListCmpF) _sort_button_inx);
	
	return error_code;
}

extern void sview_init_grid()
{
	node_info_msg_t *node_info_ptr = NULL;
	int error_code = SLURM_SUCCESS;
	node_info_t *node_ptr = NULL;
	int i = 0;
	uint16_t node_base_state;
	ListIterator itr = NULL;
	grid_button_t *grid_button = NULL;
	GdkColor color;

	if((error_code = get_new_info_node(&node_info_ptr, force_refresh))
	   == SLURM_NO_CHANGE_IN_DATA) { 
		return;
	} else if (error_code != SLURM_SUCCESS) {
		return;
	}

	if(!grid_button_list) {
		g_print("you need to run get_system_stats() first\n");
		exit(0);
	}
	
	gdk_color_parse("white", &color);
	
	itr = list_iterator_create(grid_button_list);
	for(i=0; i<node_info_ptr->record_count; i++) {
		node_ptr = &node_info_ptr->node_array[i];
		list_iterator_reset(itr);
		while((grid_button = list_next(itr))) {
			if (grid_button->inx != i)
				continue;
			node_base_state = node_ptr->node_state 
				& NODE_STATE_BASE;
			if (node_base_state == NODE_STATE_DOWN) {
				_put_button_as_down(grid_button,
						    NODE_STATE_DOWN);
			} else if (node_ptr->node_state & NODE_STATE_DRAIN) {
				_put_button_as_down(grid_button,
						    NODE_STATE_DRAIN);
			} else {
				_put_button_as_up(grid_button);
				grid_button->color = "white";
				gtk_widget_modify_bg(grid_button->button, 
						     GTK_STATE_NORMAL, &color);
			}
			grid_button->state = node_ptr->node_state;
			break;
		}
	}
	list_iterator_destroy(itr);
}

extern void sview_reset_grid()
{
	grid_button_t *grid_button = NULL;
	uint16_t node_base_state;
	ListIterator itr = NULL;
	GdkColor color;
	
	if(!grid_button_list) {
		g_print("you need to run get_system_stats() first\n");
		exit(0);
	}
	gdk_color_parse("white", &color);
		
	itr = list_iterator_create(grid_button_list);
	while((grid_button = list_next(itr))) {
		node_base_state = grid_button->state & NODE_STATE_BASE;
		if ((node_base_state == NODE_STATE_DOWN)
		    || (grid_button->state & NODE_STATE_DRAIN)) {
			continue;
		}
		gtk_widget_modify_bg(grid_button->button, 
				     GTK_STATE_NORMAL, &color);
	}
	list_iterator_destroy(itr);
}
