/*****************************************************************************\
 *  common.c - common functions used by tabs in sview
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

typedef struct {
	GtkTreeModel *model;
	GtkTreeIter iter;
} treedata_t;


static int _sort_iter_compare_func_char(GtkTreeModel *model,
					GtkTreeIter  *a,
					GtkTreeIter  *b,
					gpointer      userdata)
{
	int sortcol = GPOINTER_TO_INT(userdata);
	int ret = 0;
	gchar *name1 = NULL, *name2 = NULL;

	gtk_tree_model_get(model, a, sortcol, &name1, -1);
	gtk_tree_model_get(model, b, sortcol, &name2, -1);
	
	if (name1 == NULL || name2 == NULL)
	{
		if (name1 == NULL && name2 == NULL)
			goto cleanup; /* both equal => ret = 0 */
		
		ret = (name1 == NULL) ? -1 : 1;
	}
	else
	{
		ret = g_utf8_collate(name1,name2);
	}
cleanup:
	g_free(name1);
	g_free(name2);
	
	return ret;
}

static int _sort_iter_compare_func_int(GtkTreeModel *model,
				       GtkTreeIter  *a,
				       GtkTreeIter  *b,
				       gpointer      userdata)
{
	int sortcol = GPOINTER_TO_INT(userdata);
	int ret = 0;
	gint int1, int2;

	gtk_tree_model_get(model, a, sortcol, &int1, -1);
	gtk_tree_model_get(model, b, sortcol, &int2, -1);
	
	if (int1 != int2)
		ret = (int1 > int2) ? 1 : -1;
	
	return ret;
}

static void _add_col_to_treeview(GtkTreeView *tree_view, 
				 display_data_t *display_data)
{
	GtkTreeViewColumn   *col;
	GtkCellRenderer     *renderer;
	renderer = gtk_cell_renderer_text_new();
	col = gtk_tree_view_column_new();	
	gtk_tree_view_column_pack_start (col, renderer, TRUE);
	gtk_tree_view_column_add_attribute (col, renderer, 
					    "text", display_data->id);
	gtk_tree_view_column_set_title (col, display_data->name);
	gtk_tree_view_column_set_reorderable(col, true);
	gtk_tree_view_column_set_resizable(col, true);
	gtk_tree_view_column_set_expand(col, true);
	gtk_tree_view_append_column(tree_view, col);
	gtk_tree_view_column_set_sort_column_id(col, display_data->id);

}

static void _toggle_state_changed(GtkCheckMenuItem *menuitem, 
				  display_data_t *display_data)
{
	if(display_data->show)
		display_data->show = FALSE;
	else
		display_data->show = TRUE;
	toggled = TRUE;
	refresh_main(NULL, NULL);
}

static void _popup_state_changed(GtkCheckMenuItem *menuitem, 
				 display_data_t *display_data)
{
	popup_info_t *popup_win = (popup_info_t *) display_data->user_data;
	popup_win->toggled = 1;
	if(display_data->show)
		display_data->show = FALSE;
	else
		display_data->show = TRUE;
	
	(display_data->refresh)(NULL, display_data->user_data);
}

static void _selected_page(GtkMenuItem *menuitem, 
			   display_data_t *display_data)
{
	treedata_t *treedata = (treedata_t *)display_data->user_data;

	switch(display_data->extra) {
	case PART_PAGE:
		popup_all_part(treedata->model, &treedata->iter, 
			       display_data->id);
		break;
	case JOB_PAGE:
		popup_all_job(treedata->model, &treedata->iter, 
			      display_data->id);
		break;
	case NODE_PAGE:
		popup_all_node(treedata->model, &treedata->iter, 
			       display_data->id);
		break;
	case BLOCK_PAGE: 
		popup_all_block(treedata->model, &treedata->iter, 
				display_data->id);
		break;
	default:
		g_print("common got %d %d\n", display_data->extra,
			display_data->id);
	}
	xfree(treedata);
}

static void _set_up_button(GtkTreeView *tree_view, GdkEventButton *event, 
			   const display_data_t *display_data)
{
	//local_display_data->user_data = user_data;
	row_clicked(tree_view, event, display_data);
}

extern void snprint_time(char *buf, size_t buf_size, time_t time)
{
	if (time == INFINITE) {
		snprintf(buf, buf_size, "UNLIMITED");
	} else {
		long days, hours, minutes, seconds;
		seconds = time % 60;
		minutes = (time / 60) % 60;
		hours = (time / 3600) % 24;
		days = time / 86400;

		if (days)
			snprintf(buf, buf_size,
				"%ld-%2.2ld:%2.2ld:%2.2ld",
				days, hours, minutes, seconds);
		else if (hours)
			snprintf(buf, buf_size,
				"%ld:%2.2ld:%2.2ld", 
				hours, minutes, seconds);
		else
			snprintf(buf, buf_size,
				"%ld:%2.2ld", minutes,seconds);
	}
}

extern int get_row_number(GtkTreeView *tree_view, GtkTreePath *path)
{
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreeIter iter;
	int line = 0;
	
	if(!model) {
		g_error("error getting the model from the tree_view");
		return -1;
	}
	
	if (!gtk_tree_model_get_iter(model, &iter, path)) {
		g_error("error getting iter from model");
		return -1;
	}	
	gtk_tree_model_get(model, &iter, POS_LOC, &line, -1);
	return line;
}

extern void load_header(GtkTreeView *tree_view, display_data_t *display_data)
{
	while(display_data++) {
		if(display_data->id == -1)
			break;
		else if(!display_data->show) 
			continue;
		_add_col_to_treeview(tree_view, display_data);
	}
}

extern void make_fields_menu(GtkMenu *menu, display_data_t *display_data)
{
	GtkWidget *menuitem = NULL;
	
	while(display_data++) {
		if(display_data->id == -1)
			break;
		if(!display_data->name)
			continue;
		menuitem = gtk_check_menu_item_new_with_label(
			display_data->name); 
	
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menuitem),
					       display_data->show);
		g_signal_connect(menuitem, "toggled",
				 G_CALLBACK(_toggle_state_changed), 
				 display_data);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	}
}

extern void make_options_menu(GtkTreeView *tree_view, GtkTreePath *path,
			      GtkMenu *menu, display_data_t *display_data)
{
	GtkWidget *menuitem = NULL;
	treedata_t *treedata = xmalloc(sizeof(treedata_t));
	treedata->model = gtk_tree_view_get_model(tree_view);
	if (!gtk_tree_model_get_iter(treedata->model, &treedata->iter, path)) {
		g_error("error getting iter from model\n");
		return;
	}	
	if(display_data->user_data)
		xfree(display_data->user_data);
		
	while(display_data++) {
		if(display_data->id == -1)
			break;
		if(!display_data->name)
			continue;
		
		display_data->user_data = treedata;
		menuitem = gtk_menu_item_new_with_label(display_data->name); 
		g_signal_connect(menuitem, "activate",
				 G_CALLBACK(_selected_page), 
				 display_data);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	}
}

extern void make_popup_fields_menu(popup_info_t *popup_win, GtkMenu *menu)
{
	GtkWidget *menuitem = NULL;
	display_data_t *display_data = popup_win->display_data;
	
	while(display_data++) {
		if(display_data->id == -1)
			break;
		if(!display_data->name)
			continue;
		display_data->user_data = popup_win;
		menuitem = 
			gtk_check_menu_item_new_with_label(display_data->name);
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menuitem),
					       display_data->show);
		g_signal_connect(menuitem, "toggled",
				 G_CALLBACK(_popup_state_changed), 
				 display_data);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	}
}


extern GtkScrolledWindow *create_scrolled_window()
{
	GtkScrolledWindow *scrolled_window = NULL;
	GtkWidget *table = NULL;
	table = gtk_table_new(1, 1, FALSE);

	gtk_container_set_border_width(GTK_CONTAINER(table), 10);	

	scrolled_window = GTK_SCROLLED_WINDOW(gtk_scrolled_window_new(
						      NULL, NULL));	
	gtk_container_set_border_width(GTK_CONTAINER(scrolled_window), 10);
    
	gtk_scrolled_window_set_policy(scrolled_window,
				       GTK_POLICY_AUTOMATIC,
				       GTK_POLICY_AUTOMATIC);
    
	gtk_scrolled_window_add_with_viewport(scrolled_window, table);

	return scrolled_window;
}
extern void create_page(GtkNotebook *notebook, display_data_t *display_data)
{
	GtkScrolledWindow *scrolled_window = create_scrolled_window();
	GtkWidget *event_box = NULL;
	GtkWidget *label = NULL;
	int err;
		
	event_box = gtk_event_box_new();
	gtk_event_box_set_above_child(GTK_EVENT_BOX(event_box), FALSE);
	g_signal_connect(G_OBJECT(event_box), "button-press-event",
			 G_CALLBACK(tab_pressed),
			 display_data);
	
	label = gtk_label_new(display_data->name);
	gtk_container_add(GTK_CONTAINER(event_box), label);
	gtk_widget_show(label);
	//(display_data->set_fields)(GTK_MENU(menu));
	if((err = gtk_notebook_append_page(GTK_NOTEBOOK(notebook), 
					   GTK_WIDGET(scrolled_window), 
					   event_box)) == -1) {
		g_error("Couldn't add page to notebook\n");
	}
	
	display_data->extra = err;

}

extern GtkTreeView *create_treeview(display_data_t *local, gpointer user_data)
{
	GtkTreeView *tree_view = GTK_TREE_VIEW(gtk_tree_view_new());
	/* g_signal_connect(G_OBJECT(tree_view), "row-activated", */
/* 			 G_CALLBACK(local->row_clicked), */
/* 			 user_data); */

	local->user_data = user_data;
	g_signal_connect(G_OBJECT(tree_view), "button_press_event",
			 G_CALLBACK(_set_up_button),
			 local);
	

	gtk_widget_show(GTK_WIDGET(tree_view));
	
	return tree_view;

}

extern GtkTreeStore *create_treestore(GtkTreeView *tree_view, 
				      display_data_t *display_data,
				      int count)
{
	GtkTreeStore *treestore = NULL;
	GType types[count];
	int i=0;
	
	/*set up the types defined in the display_data_t */
	for(i=0; i<count; i++)
		types[i] = display_data[i].type;

	treestore = gtk_tree_store_newv(count, types);
	if(!treestore)
		return NULL;

	for(i=1; i<count; i++) {
		if(display_data[i].show) {
			switch(display_data[i].type) {
			case G_TYPE_INT:
				gtk_tree_sortable_set_sort_func(
					GTK_TREE_SORTABLE(treestore), 
					i, 
					_sort_iter_compare_func_int,
					GINT_TO_POINTER(i), 
					NULL); 
				
				break;
			case G_TYPE_STRING:
				gtk_tree_sortable_set_sort_func(
				GTK_TREE_SORTABLE(treestore), 
				i, 
				_sort_iter_compare_func_char,
				GINT_TO_POINTER(i), 
				NULL); 
				break;
			default:
				g_print("unknown type %d",
					(int)display_data[i].type);
			}
		}
	}
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(treestore), 
					     1, 
					     GTK_SORT_ASCENDING);

	load_header(tree_view, display_data);
	gtk_tree_view_set_model(tree_view, GTK_TREE_MODEL(treestore));
	g_object_unref(GTK_TREE_MODEL(treestore));

	return treestore;
}

extern void right_button_pressed(GtkTreeView *tree_view, 
				 GtkTreePath *path,
				 GdkEventButton *event, 
				 const display_data_t *display_data,
				 int type)
{
	if(event->button == 3) {
		GtkMenu *menu = GTK_MENU(gtk_menu_new());
	
		(display_data->set_menu)(tree_view, path, menu, type);
				
		gtk_widget_show_all(GTK_WIDGET(menu));
		gtk_menu_popup(menu, NULL, NULL, NULL, NULL,
			       (event != NULL) ? event->button : 0,
			       gdk_event_get_time((GdkEvent*)event));
	}
}

extern void row_clicked(GtkTreeView *tree_view, GdkEventButton *event, 
			const display_data_t *display_data)
{
	GtkTreePath *path = NULL;
	GtkTreeSelection *selection = NULL;
	
        if(!gtk_tree_view_get_path_at_pos(tree_view,
					  (gint) event->x, 
					  (gint) event->y,
					  &path, NULL, NULL, NULL)) {
		return;
	}
	selection = gtk_tree_view_get_selection(tree_view);
	gtk_tree_selection_unselect_all(selection);
	gtk_tree_selection_select_path(selection, path);
			 	
	/* expand/collapse row 
	   or right mouse button
	   or double click?
	*/
	if(event->x <= 20) {	
		if(!gtk_tree_view_expand_row(tree_view, path, FALSE))
			gtk_tree_view_collapse_row(tree_view, path);
	} else if(event->button == 3) {
		right_button_pressed(tree_view, path, event, 
				     display_data, ROW_CLICKED);
	} else if(event->type==GDK_2BUTTON_PRESS ||
		  event->type==GDK_3BUTTON_PRESS) {
		(display_data->row_clicked)(tree_view, path, 
					    NULL, display_data->user_data);
	}
	gtk_tree_path_free(path);
}

extern popup_info_t *create_popup_info(int type, int dest_type, char *title)
{
	GtkScrolledWindow *window = NULL;
	GtkBin *bin = NULL;
	GtkViewport *view = NULL;
	GtkBin *bin2 = NULL;
	GtkWidget *table = NULL;
	GtkWidget *popup = NULL;
	popup_info_t *popup_win = xmalloc(sizeof(popup_info_t));

	list_push(popup_list, popup_win);
	
	popup_win->spec_info = xmalloc(sizeof(specific_info_t));
	popup_win->popup = gtk_dialog_new();
	popup_win->toggled = 0;
	popup_win->type = dest_type;

	gtk_window_set_default_size(GTK_WINDOW(popup_win->popup), 
				    600, 400);
	gtk_window_set_title(GTK_WINDOW(popup_win->popup), title);
	
	popup = popup_win->popup;

	table = gtk_table_new(1, 2, FALSE);
	gtk_container_set_border_width(GTK_CONTAINER(table), 10);

	popup_win->event_box = gtk_event_box_new();
	gtk_event_box_set_above_child(
		GTK_EVENT_BOX(popup_win->event_box), 
		FALSE);
	popup_win->button = gtk_button_new_with_label("Refresh");
	gtk_table_attach_defaults(GTK_TABLE(table), 
				  popup_win->event_box,
				  0, 1, 0, 1); 
	gtk_table_attach(GTK_TABLE(table), 
			 popup_win->button,
			 1, 2, 0, 1,
			 GTK_SHRINK, GTK_EXPAND | GTK_FILL,
			 0, 0); 
		
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(popup)->vbox), 
			   table, FALSE, FALSE, 0);

	window = create_scrolled_window();
	bin = GTK_BIN(&window->container);
	view = GTK_VIEWPORT(bin->child);
	bin2 = GTK_BIN(&view->bin);
	popup_win->table = GTK_TABLE(bin2->child);
	
	gtk_box_pack_end(GTK_BOX(GTK_DIALOG(popup)->vbox), 
			 GTK_WIDGET(window), TRUE, TRUE, 0);
	popup_win->spec_info->type = type;
	popup_win->spec_info->title = xstrdup(title);
	g_signal_connect(G_OBJECT(popup_win->popup), "delete_event",
			 G_CALLBACK(delete_popup), 
			 popup_win->spec_info->title);
	gtk_widget_show_all(popup_win->popup);
	return popup_win;
}

extern void setup_popup_info(popup_info_t *popup_win, 
			     display_data_t *display_data, 
			     int cnt)
{
	GtkWidget *label = NULL;
	int i = 0;
	specific_info_t *spec_info = popup_win->spec_info;
	
	popup_win->display_data = xmalloc(sizeof(display_data_t)*(cnt+2));
	for(i=0; i<cnt+1; i++) {
		memcpy(&popup_win->display_data[i], 
		       &display_data[i], 
		       sizeof(display_data_t));
	}
	
	g_signal_connect(G_OBJECT(popup_win->event_box), 
			 "button-press-event",
			 G_CALLBACK(redo_popup),
			 popup_win);
	
	g_signal_connect(G_OBJECT(popup_win->button), 
			 "pressed",
			 G_CALLBACK(popup_win->display_data->refresh),
			 popup_win);
	
	label = gtk_label_new(spec_info->title);
	gtk_container_add(GTK_CONTAINER(popup_win->event_box), label);
	gtk_widget_show(label);
}

extern void redo_popup(GtkWidget *widget, GdkEventButton *event, 
		       popup_info_t *popup_win)
{
	if(event->button == 3) {
		GtkMenu *menu = GTK_MENU(gtk_menu_new());
		
		(popup_win->display_data->set_menu)(popup_win, 
						    NULL, 
						    menu, POPUP_CLICKED);
		
		gtk_widget_show_all(GTK_WIDGET(menu));
		gtk_menu_popup(menu, NULL, NULL, NULL, NULL,
			       (event != NULL) ? event->button : 0,
			       gdk_event_get_time((GdkEvent*)event));
	}
}

extern void destroy_specific_info(void *arg)
{
	specific_info_t *spec_info = (specific_info_t *)arg;
	if(spec_info) {
		xfree(spec_info->title);
		if(spec_info->data) {
			g_free(spec_info->data);
			spec_info->data = NULL;
		}
		if(spec_info->display_widget) {
			gtk_widget_destroy(spec_info->display_widget);
			spec_info->display_widget = NULL;
		}
		xfree(spec_info);
	}
}

extern void destroy_popup_info(void *arg)
{
	popup_info_t *popup_win = (popup_info_t *)arg;
	if(popup_win) {
		*popup_win->running = 0;
		/* these are all childern of each other so must 
		   be freed in this order */
		if(popup_win->table) {
			gtk_widget_destroy(GTK_WIDGET(popup_win->table));
			popup_win->table = NULL;
		}
		if(popup_win->button) {
			gtk_widget_destroy(popup_win->button);
			popup_win->button = NULL;
		}
		if(popup_win->event_box) {
			gtk_widget_destroy(popup_win->event_box);
			popup_win->event_box = NULL;
		}
		if(popup_win->popup) {
			gtk_widget_destroy(popup_win->popup);
			popup_win->popup = NULL;
		}
		
		destroy_specific_info(popup_win->spec_info);
		xfree(popup_win->display_data);
		xfree(popup_win);
	}
}

extern gboolean delete_popup(GtkWidget *widget, GtkWidget *event, char *title)
{
	ListIterator itr = list_iterator_create(popup_list);
	popup_info_t *popup_win = NULL;
	
	while((popup_win = list_next(itr))) {
		if(popup_win->spec_info) {
			if(!strcmp(popup_win->spec_info->title, title)) {
				//g_print("removing %s\n", title);
				list_remove(itr);
				destroy_popup_info(popup_win);
				break;
			}
		}
	}
	list_iterator_destroy(itr);
	

	return FALSE;
}

extern void *popup_thr(popup_info_t *popup_win)
{
	void (*specifc_info) (popup_info_t *popup_win) = NULL;
	int running = 1;
	switch(popup_win->type) {
	case PART_PAGE:
		specifc_info = specific_info_part;
		break;
	case JOB_PAGE:
		specifc_info = specific_info_job;
		break;
	case NODE_PAGE:
		specifc_info = specific_info_node;
		break;
	case BLOCK_PAGE: 
		specifc_info = specific_info_block;
		break;
	case ADMIN_PAGE: 
	case SUBMIT_PAGE: 
	default:
		g_print("thread got unknown type %d\n", popup_win->type);
		return NULL;
	}
	/* this will switch to 0 when popup is closed. */
	popup_win->running = &running;
	/* when popup is killed toggled will be set to -1 */
	while(running) {
		gdk_threads_enter();
		(specifc_info)(popup_win);
		gdk_flush();
		gdk_threads_leave();
		sleep(5);
	}	
	return NULL;
}

extern void remove_old(GtkTreeModel *model, int updated)
{
	GtkTreePath *path = gtk_tree_path_new_first();
	GtkTreeIter iter;
	int i;
	
	/* remove all old partitions */
	if (gtk_tree_model_get_iter(model, &iter, path)) {
		while(1) {
			gtk_tree_model_get(model, &iter, updated, &i, -1);
			if(!i) {
				if(!gtk_tree_store_remove(
					   GTK_TREE_STORE(model), 
					   &iter))
					break;
				else
					continue;
			}
			if(!gtk_tree_model_iter_next(model, &iter)) {
				break;
			}
		}
	}
	gtk_tree_path_free(path);
}
