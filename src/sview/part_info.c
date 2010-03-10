/*****************************************************************************\
 *  part_info.c - Functions related to partition display
 *  mode of sview.
 *****************************************************************************
 *  Copyright (C) 2004-2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Danny Auble <da@llnl.gov>
 *
 *  CODE-OCEC-09-009. All rights reserved.
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
#include "src/common/parse_time.h"

#define _DEBUG 0

typedef struct {
	uint32_t cpu_alloc_cnt;
	uint32_t cpu_error_cnt;
	uint32_t cpu_idle_cnt;
	uint32_t disk_total;
	char *features;
	hostlist_t hl;
	uint32_t mem_total;
	uint32_t node_cnt;
	List node_ptr_list;
	uint16_t node_state;
	partition_info_t* part_ptr;
	char *reason;
} sview_part_sub_t;

/* Collection of data for printing reports. Like data is combined here */
typedef struct {
	int color_inx;
	/* part_info contains partition, avail, max_time, job_size,
	 * root, share, groups */
	partition_info_t* part_ptr;
	List sub_list;
} sview_part_info_t;

enum {
	EDIT_PART_STATE = 1,
	EDIT_EDIT
};

/* These need to be in alpha order (except POS and CNT) */
enum {
	SORTID_POS = POS_LOC,
#ifdef HAVE_BG
	SORTID_NODELIST,
	SORTID_NODES_ALLOWED,
#endif
	SORTID_COLOR,
	SORTID_CPUS,
	SORTID_DEFAULT,
	SORTID_FEATURES,
	SORTID_GROUPS,
	SORTID_HIDDEN,
	SORTID_JOB_SIZE,
	SORTID_MEM,
	SORTID_NAME,
#ifndef HAVE_BG
	SORTID_NODELIST,
	SORTID_NODES_ALLOWED,
#endif
	SORTID_NODE_INX,
	SORTID_NODE_STATE,
	SORTID_NODE_STATE_NUM,
	SORTID_NODES,
	SORTID_NODES_MAX,
	SORTID_NODES_MIN,
	SORTID_ONLY_LINE,
	SORTID_PART_STATE,
	SORTID_PRIORITY,
	SORTID_REASON,
	SORTID_ROOT,
	SORTID_SHARE,
	SORTID_TMP_DISK,
	SORTID_TIMELIMIT,
	SORTID_UPDATED,
	SORTID_CNT
};

static display_data_t display_data_part[] = {
	{G_TYPE_INT, SORTID_POS, NULL, FALSE, EDIT_NONE, refresh_part},
	{G_TYPE_STRING, SORTID_NAME, "Partition", TRUE,
	 EDIT_NONE, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_COLOR, NULL, TRUE, EDIT_NONE, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_DEFAULT, "Default", TRUE,
	 EDIT_MODEL, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_HIDDEN, "Hidden", FALSE,
	 EDIT_MODEL, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_PART_STATE, "Part State", TRUE,
	 EDIT_MODEL, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_TIMELIMIT, "Time Limit",
	 TRUE, EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_NODES, "Node Count",
	 TRUE, EDIT_NONE, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_CPUS, "CPU Count",
	 FALSE, EDIT_NONE, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_NODE_STATE, "Node State",
	 TRUE, EDIT_MODEL, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_JOB_SIZE, "Job Size", FALSE,
	 EDIT_NONE, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_PRIORITY, "Priority", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_NODES_MIN, "Nodes Min", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_NODES_MAX, "Nodes Max", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_ROOT, "Root", FALSE, EDIT_MODEL, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_SHARE, "Share", FALSE, EDIT_MODEL, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_GROUPS, "Groups Allowed", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
#ifdef HAVE_BG
	{G_TYPE_STRING, SORTID_NODES_ALLOWED, "BPs Allowed Allocating", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
#else
	{G_TYPE_STRING, SORTID_NODES_ALLOWED, "Nodes Allowed Allocating", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
#endif
	{G_TYPE_STRING, SORTID_TMP_DISK, "Temp Disk", FALSE,
	 EDIT_NONE, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_MEM, "Memory", FALSE, EDIT_NONE, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_FEATURES, "Features", FALSE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_STRING, SORTID_REASON, "Reason", FALSE,
	 EDIT_NONE, refresh_part, create_model_part, admin_edit_part},
#ifdef HAVE_BG
	{G_TYPE_STRING, SORTID_NODELIST, "BP List", TRUE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
#else
	{G_TYPE_STRING, SORTID_NODELIST, "NodeList", TRUE,
	 EDIT_TEXTBOX, refresh_part, create_model_part, admin_edit_part},
#endif
	{G_TYPE_INT, SORTID_NODE_STATE_NUM, NULL, FALSE, EDIT_NONE, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_INT, SORTID_ONLY_LINE, NULL, FALSE, EDIT_NONE, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_POINTER, SORTID_NODE_INX, NULL, FALSE, EDIT_NONE,
	 refresh_part, create_model_part, admin_edit_part},
	{G_TYPE_INT, SORTID_UPDATED, NULL, FALSE, EDIT_NONE, refresh_part,
	 create_model_part, admin_edit_part},
	{G_TYPE_NONE, -1, NULL, FALSE, EDIT_NONE}
};

static display_data_t options_data_part[] = {
	{G_TYPE_INT, SORTID_POS, NULL, FALSE, EDIT_NONE},
	{G_TYPE_STRING, INFO_PAGE, "Full Info", TRUE, PART_PAGE},
#ifdef HAVE_BG
	{G_TYPE_STRING, PART_PAGE, "Drain Base Partitions", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Resume Base Partitions", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Put Base Partitions Down",
	 TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Make Base Partitions Idle",
	 TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Update Base Partition Features",
	 TRUE, ADMIN_PAGE},
#else
	{G_TYPE_STRING, PART_PAGE, "Drain Nodes", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Resume Nodes", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Put Nodes Down", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Make Nodes Idle", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Update Node Features", TRUE, ADMIN_PAGE},
#endif
	{G_TYPE_STRING, PART_PAGE, "Change Part State Up/Down",
	 TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, PART_PAGE, "Edit Part", TRUE, ADMIN_PAGE},
	{G_TYPE_STRING, JOB_PAGE, "Jobs", TRUE, PART_PAGE},
#ifdef HAVE_BG
	{G_TYPE_STRING, BLOCK_PAGE, "Blocks", TRUE, PART_PAGE},
	{G_TYPE_STRING, NODE_PAGE, "Base Partitions", TRUE, PART_PAGE},
#else
	{G_TYPE_STRING, NODE_PAGE, "Nodes", TRUE, PART_PAGE},
#endif
	//{G_TYPE_STRING, SUBMIT_PAGE, "Job Submit", FALSE, PART_PAGE},
	{G_TYPE_STRING, RESV_PAGE, "Reservations", TRUE, PART_PAGE},
	{G_TYPE_NONE, -1, NULL, FALSE, EDIT_NONE}
};

static display_data_t *local_display_data = NULL;

static char *got_edit_signal = NULL;
static char *got_features_edit_signal = NULL;

static void _update_part_sub_record(sview_part_sub_t *sview_part_sub,
				    GtkTreeStore *treestore,
				    GtkTreeIter *iter);
static void _append_part_sub_record(sview_part_sub_t *sview_part_sub,
				    GtkTreeStore *treestore, GtkTreeIter *iter,
				    int line);

static int _build_min_max_32_string(char *buffer, int buf_size,
				    uint32_t min, uint32_t max, bool range)
{
	char tmp_min[8];
	char tmp_max[8];
	convert_num_unit((float)min, tmp_min, sizeof(tmp_min), UNIT_NONE);
	convert_num_unit((float)max, tmp_max, sizeof(tmp_max), UNIT_NONE);

	if (max == min)
		return snprintf(buffer, buf_size, "%s", tmp_max);
	else if (range) {
		if (max == (uint32_t) INFINITE)
			return snprintf(buffer, buf_size, "%s-infinite",
					tmp_min);
		else
			return snprintf(buffer, buf_size, "%s-%s",
					tmp_min, tmp_max);
	} else
		return snprintf(buffer, buf_size, "%s+", tmp_min);
}

static void _set_active_combo_part(GtkComboBox *combo,
				   GtkTreeModel *model, GtkTreeIter *iter,
				   int type)
{
	char *temp_char = NULL;
	int action = 0;
	int i = 0, unknown_found = 0;
	char *upper = NULL;

	gtk_tree_model_get(model, iter, type, &temp_char, -1);
	if(!temp_char)
		goto end_it;
	switch(type) {
	case SORTID_DEFAULT:
	case SORTID_HIDDEN:
	case SORTID_ROOT:
		if(!strcmp(temp_char, "yes"))
			action = 0;
		else if(!strcmp(temp_char, "no"))
			action = 1;
		else
			action = 0;

		break;
	case SORTID_SHARE:
		if(!strncmp(temp_char, "force", 5))
			action = 0;
		else if(!strcmp(temp_char, "no"))
			action = 1;
		else if(!strncmp(temp_char, "yes", 3))
			action = 2;
		else if(!strcmp(temp_char, "exclusive"))
			action = 3;
		else
			action = 0;
		break;
	case SORTID_PART_STATE:
		if(!strcmp(temp_char, "up"))
			action = 0;
		else if(!strcmp(temp_char, "down"))
			action = 1;
		else
			action = 0;
		break;
	case SORTID_NODE_STATE:
		if(!strcasecmp(temp_char, "drain"))
			action = 0;
		else if(!strcasecmp(temp_char, "resume"))
			action = 1;
		else
			for(i = 0; i < NODE_STATE_END; i++) {
				upper = node_state_string(i);
				if(!strcmp(upper, "UNKNOWN")) {
					unknown_found++;
					continue;
				}

				if(!strcasecmp(temp_char, upper)) {
					action = i + 2 - unknown_found;
					break;
				}
			}

		break;
	default:
		break;
	}
	g_free(temp_char);
end_it:
	gtk_combo_box_set_active(combo, action);

}

static uint16_t _set_part_share_popup()
{
	GtkWidget *table = gtk_table_new(1, 2, FALSE);
	GtkWidget *label = NULL;
	GtkObject *adjustment = gtk_adjustment_new(4,
						   1, 1000,
						   1, 60,
						   0);
	GtkWidget *spin_button =
		gtk_spin_button_new(GTK_ADJUSTMENT(adjustment), 1, 0);
	GtkWidget *popup = gtk_dialog_new_with_buttons(
		"Count",
		GTK_WINDOW (main_window),
		GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
		NULL);
	int response = 0;
	uint16_t count = 4;

	label = gtk_dialog_add_button(GTK_DIALOG(popup),
				      GTK_STOCK_OK, GTK_RESPONSE_OK);
	gtk_window_set_default(GTK_WINDOW(popup), label);

	label = gtk_label_new("Shared Job Count ");

	gtk_container_set_border_width(GTK_CONTAINER(table), 10);

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(popup)->vbox),
			   table, FALSE, FALSE, 0);

	gtk_table_attach_defaults(GTK_TABLE(table), label, 0, 1, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), spin_button, 1, 2, 0, 1);

	gtk_widget_show_all(popup);
	response = gtk_dialog_run (GTK_DIALOG(popup));

	if (response == GTK_RESPONSE_OK) {
		count = gtk_spin_button_get_value_as_int(
			GTK_SPIN_BUTTON(spin_button));
	}

	gtk_widget_destroy(popup);

	return count;
}

/* don't free this char */
static const char *_set_part_msg(update_part_msg_t *part_msg,
				 const char *new_text,
				 int column)
{
	char *type = "";
	int temp_int = 0;

	global_edit_error = 0;

	if(!part_msg)
		return NULL;

	switch(column) {
	case SORTID_DEFAULT:
		if (!strcasecmp(new_text, "yes"))
			part_msg->default_part = 1;
		else
			part_msg->default_part = 0;

		type = "default";
		break;
	case SORTID_HIDDEN:
		if (!strcasecmp(new_text, "yes"))
			part_msg->hidden = 1;
		else
			part_msg->hidden = 0;

		type = "hidden";
		break;
	case SORTID_TIMELIMIT:
		if ((strcasecmp(new_text,"infinite") == 0))
			temp_int = INFINITE;
		else
			temp_int = time_str2mins((char *)new_text);

		type = "timelimit";
		if((temp_int <= 0) && (temp_int != INFINITE))
			goto return_error;
		part_msg->max_time = (uint32_t)temp_int;
		break;
	case SORTID_PRIORITY:
		temp_int = strtol(new_text, (char **)NULL, 10);
		type = "priority";
		part_msg->priority = (uint16_t)temp_int;
		break;
	case SORTID_NODES_MIN:
		temp_int = strtol(new_text, (char **)NULL, 10);
		type = "min_nodes";

		if(temp_int <= 0)
			goto return_error;
		part_msg->min_nodes = (uint32_t)temp_int;
		break;
	case SORTID_NODES_MAX:
		if (!strcasecmp(new_text, "infinite")) {
			temp_int = INFINITE;
		} else {
			temp_int = strtol(new_text, (char **)NULL, 10);
		}

		type = "max_nodes";
		if(temp_int <= 0 && temp_int != INFINITE)
			goto return_error;
		part_msg->max_nodes = (uint32_t)temp_int;
		break;
	case SORTID_ROOT:
		if (!strcasecmp(new_text, "yes")) {
			part_msg->root_only = 1;
		} else {
			part_msg->root_only = 0;
		}

		type = "root";
		break;
	case SORTID_SHARE:
		if (!strcasecmp(new_text, "yes")) {
			part_msg->max_share = _set_part_share_popup();
		} else if (!strcasecmp(new_text, "exclusive")) {
			part_msg->max_share = 0;
		} else if (!strcasecmp(new_text, "force")) {
			part_msg->max_share =
				_set_part_share_popup() | SHARED_FORCE;
		} else if (!strcasecmp(new_text, "no))
			part_msg->max_share = 1;
		else
			goto return_error
		type = "share";
		break;
	case SORTID_GROUPS:
		type = "groups";
		part_msg->allow_groups = xstrdup(new_text);
		break;
	case SORTID_NODES_ALLOWED:
		type = "allowed alloc nodes";
		part_msg->allow_alloc_nodes = xstrdup(new_text);
		break;
	case SORTID_NODELIST:
		part_msg->nodes = xstrdup(new_text);
		type = "nodelist";
		break;
	case SORTID_PART_STATE:
		if (!strcasecmp(new_text, "up"))
			part_msg->state_up = PARTITION_UP;
		else if (!strcasecmp(new_text, "down"))
			part_msg->state_up = PARTITION_DOWN;
		else if (!strcasecmp(new_text, "inactive"))
			part_msg->state_up = PARTITION_INACTIVE;
		else if (!strcasecmp(new_text, "drain"))
			part_msg->state_up = PARTITION_DRAIN;
		else
			goto return_error;
		type = "availability";

		break;
	case SORTID_NODE_STATE:
		type = (char *)new_text;
		got_edit_signal = xstrdup(new_text);
		break;
	case SORTID_FEATURES:
		type = "Update Features";
		got_features_edit_signal = xstrdup(new_text);
		break;
	default:
		type = "unknown";
		break;
	}

	if(strcmp(type, "unknown"))
		global_send_update_msg = 1;

	return type;

return_error:
	global_edit_error = 1;
	return type;

}

static void _admin_edit_combo_box_part(GtkComboBox *combo,
				       update_part_msg_t *part_msg)
{
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	int column = 0;
	char *name = NULL;

	if(!part_msg)
		return;

	if(!gtk_combo_box_get_active_iter(combo, &iter)) {
		g_print("nothing selected\n");
		return;
	}
	model = gtk_combo_box_get_model(combo);
	if(!model) {
		g_print("nothing selected\n");
		return;
	}

	gtk_tree_model_get(model, &iter, 0, &name, -1);
	gtk_tree_model_get(model, &iter, 1, &column, -1);

	_set_part_msg(part_msg, name, column);

	g_free(name);
}

static gboolean _admin_focus_out_part(GtkEntry *entry,
				      GdkEventFocus *event,
				      update_part_msg_t *part_msg)
{
	if(global_entry_changed) {
		const char *col_name = NULL;
		int type = gtk_entry_get_max_length(entry);
		const char *name = gtk_entry_get_text(entry);
		type -= DEFAULT_ENTRY_LENGTH;
		col_name = _set_part_msg(part_msg, name, type);
		if(global_edit_error) {
			if(global_edit_error_msg)
				g_free(global_edit_error_msg);
			global_edit_error_msg = g_strdup_printf(
				"Partition %s %s can't be set to %s",
				part_msg->name,
				col_name,
				name);
		}

		global_entry_changed = 0;
	}
	return false;
}

static GtkWidget *_admin_full_edit_part(update_part_msg_t *part_msg,
					GtkTreeModel *model, GtkTreeIter *iter)
{
	GtkScrolledWindow *window = create_scrolled_window();
	GtkBin *bin = NULL;
	GtkViewport *view = NULL;
	GtkTable *table = NULL;
	int i = 0, row = 0;
	display_data_t *display_data = display_data_part;

	gtk_scrolled_window_set_policy(window,
				       GTK_POLICY_NEVER,
				       GTK_POLICY_AUTOMATIC);
	bin = GTK_BIN(&window->container);
	view = GTK_VIEWPORT(bin->child);
	bin = GTK_BIN(&view->bin);
	table = GTK_TABLE(bin->child);
	gtk_table_resize(table, SORTID_CNT, 2);

	gtk_table_set_homogeneous(table, FALSE);

	for(i = 0; i < SORTID_CNT; i++) {
		while(display_data++) {
			if(display_data->id == -1)
				break;
			if(!display_data->name)
				continue;
			if(display_data->id != i)
				continue;
			display_admin_edit(
				table, part_msg, &row, model, iter,
				display_data,
				G_CALLBACK(_admin_edit_combo_box_part),
				G_CALLBACK(_admin_focus_out_part),
				_set_active_combo_part);
			break;
		}
		display_data = display_data_part;
	}
	gtk_table_resize(table, row, 2);

	return GTK_WIDGET(window);
}

static void _subdivide_part(sview_part_info_t *sview_part_info,
			    GtkTreeModel *model,
			    GtkTreeIter *sub_iter,
			    GtkTreeIter *iter)
{
	GtkTreeIter first_sub_iter;
	ListIterator itr = NULL;
	uint16_t state;
	int i = 0, line = 0;
	sview_part_sub_t *sview_part_sub = NULL;
	int set = 0;

	memset(&first_sub_iter, 0, sizeof(GtkTreeIter));

	/* make sure all the steps are still here */
	if (sub_iter) {
		first_sub_iter = *sub_iter;
		while(1) {
			gtk_tree_store_set(GTK_TREE_STORE(model), sub_iter,
					   SORTID_UPDATED, 0, -1);
			if(!gtk_tree_model_iter_next(model, sub_iter)) {
				break;
			}
		}
		memcpy(sub_iter, &first_sub_iter, sizeof(GtkTreeIter));
		set = 1;
	}
	itr = list_iterator_create(sview_part_info->sub_list);
	if(list_count(sview_part_info->sub_list) == 1) {
		gtk_tree_store_set(GTK_TREE_STORE(model), iter,
				   SORTID_ONLY_LINE, 1, -1);
		sview_part_sub = list_next(itr);
		_update_part_sub_record(sview_part_sub,
					GTK_TREE_STORE(model),
					iter);
	} else {
		while((sview_part_sub = list_next(itr))) {
			if (!sub_iter) {
				i = NO_VAL;
				goto adding;
			} else {
				memcpy(sub_iter, &first_sub_iter,
				       sizeof(GtkTreeIter));
			}
			line = 0;
			while(1) {
				/* search for the state number and
				   check to see if it is in the list */
				gtk_tree_model_get(model, sub_iter,
						   SORTID_NODE_STATE_NUM,
						   &state, -1);
				if(state == sview_part_sub->node_state) {
					/* update with new info */
					_update_part_sub_record(
						sview_part_sub,
						GTK_TREE_STORE(model),
						sub_iter);
					goto found;
				}

				line++;
				if(!gtk_tree_model_iter_next(model,
							     sub_iter)) {
					break;
				}
			}
		adding:
			_append_part_sub_record(sview_part_sub,
						GTK_TREE_STORE(model),
						iter, line);
/* 			if(i == NO_VAL) */
/* 				line++; */
		found:
			;
		}
	}
	list_iterator_destroy(itr);

	if(set) {
		sub_iter = &first_sub_iter;
		/* clear all steps that aren't active */
		while(1) {
			gtk_tree_model_get(model, sub_iter,
					   SORTID_UPDATED, &i, -1);
			if(!i) {
				if(!gtk_tree_store_remove(
					   GTK_TREE_STORE(model),
					   sub_iter))
					break;
				else
					continue;
			}
			if(!gtk_tree_model_iter_next(model, sub_iter)) {
				break;
			}
		}
	}
	return;
}

static void _layout_part_record(GtkTreeView *treeview,
				sview_part_info_t *sview_part_info,
				int update)
{
	GtkTreeIter iter;
	ListIterator itr = NULL;
	char time_buf[20], tmp_buf[20];
	char tmp_cnt[8];
	char tmp_cnt1[8];
	char tmp_cnt2[8];
	partition_info_t *part_ptr = sview_part_info->part_ptr;
	sview_part_sub_t *sview_part_sub = NULL;
	sview_part_sub_t *temp_part_sub = NULL;
	sview_part_sub_t alloc_part_sub;
	sview_part_sub_t idle_part_sub;
	sview_part_sub_t other_part_sub;
	char ind_cnt[1024];
	char *temp_char = NULL;
	int global_set = 0, i;
	int yes_no = -1;
	int up_down = -1;
	uint32_t limit_set = NO_VAL;
	GtkTreeStore *treestore =
		GTK_TREE_STORE(gtk_tree_view_get_model(treeview));

	memset(&alloc_part_sub, 0, sizeof(sview_part_sub_t));
	memset(&idle_part_sub, 0, sizeof(sview_part_sub_t));
	memset(&other_part_sub, 0, sizeof(sview_part_sub_t));

	itr = list_iterator_create(sview_part_info->sub_list);
	while((sview_part_sub = list_next(itr))) {
		if(sview_part_sub->node_state == NODE_STATE_IDLE)
			temp_part_sub = &idle_part_sub;
		else if(sview_part_sub->node_state == NODE_STATE_ALLOCATED)
			temp_part_sub = &alloc_part_sub;
		else
			temp_part_sub = &other_part_sub;
		temp_part_sub->node_cnt += sview_part_sub->node_cnt;
		temp_part_sub->cpu_alloc_cnt += sview_part_sub->cpu_alloc_cnt;
		temp_part_sub->cpu_error_cnt += sview_part_sub->cpu_error_cnt;
		temp_part_sub->cpu_idle_cnt += sview_part_sub->cpu_idle_cnt;
		/* temp_part_sub->disk_total += sview_part_sub->disk_total; */
/* 		temp_part_sub->mem_total += sview_part_sub->mem_total; */

		if(!global_set) {
			global_set = 1;
			/* store features and reasons in the others
			   group */
			other_part_sub.features = sview_part_sub->features;
			other_part_sub.reason = sview_part_sub->reason;
			other_part_sub.disk_total = sview_part_sub->disk_total;
			other_part_sub.mem_total = sview_part_sub->mem_total;
		} else {
			other_part_sub.disk_total += sview_part_sub->disk_total;
			other_part_sub.mem_total += sview_part_sub->mem_total;
		}
	}
	list_iterator_destroy(itr);

	convert_num_unit((float)alloc_part_sub.node_cnt,
			 tmp_cnt, sizeof(tmp_cnt), UNIT_NONE);
	convert_num_unit((float)idle_part_sub.node_cnt,
			 tmp_cnt1, sizeof(tmp_cnt1), UNIT_NONE);
	convert_num_unit((float)other_part_sub.node_cnt,
			 tmp_cnt2, sizeof(tmp_cnt2), UNIT_NONE);
	snprintf(ind_cnt, sizeof(ind_cnt), "%s/%s/%s",
		 tmp_cnt, tmp_cnt1, tmp_cnt2);

	for(i = 0; i < SORTID_CNT; i++) {
		switch(i) {
		case SORTID_PART_STATE:
			up_down = part_ptr->state_up;
			break;
		case SORTID_CPUS:
			convert_num_unit((float)part_ptr->total_cpus,
					 tmp_cnt, sizeof(tmp_cnt),
					 UNIT_NONE);
			temp_char = tmp_cnt;
			break;
		case SORTID_DEFAULT:
			yes_no = part_ptr->default_part;
			break;
		case SORTID_FEATURES:
			sview_part_sub = list_peek(sview_part_info->sub_list);
			temp_char = sview_part_sub->features;
			break;
		case SORTID_GROUPS:
			if(part_ptr->allow_groups)
				temp_char = part_ptr->allow_groups;
			else
				temp_char = "all";
			break;
		case SORTID_HIDDEN:
			yes_no = part_ptr->hidden;
			break;
		case SORTID_JOB_SIZE:
			_build_min_max_32_string(time_buf, sizeof(time_buf),
			      part_ptr->min_nodes,
			      part_ptr->max_nodes, true);
			temp_char = time_buf;
			break;
		case SORTID_MEM:
			convert_num_unit((float)other_part_sub.mem_total,
					 tmp_cnt, sizeof(tmp_cnt),
					 UNIT_MEGA);
			temp_char = tmp_cnt;
			break;
		case SORTID_NODELIST:
			temp_char = part_ptr->nodes;
			break;
		case SORTID_NODES_ALLOWED:
			temp_char = part_ptr->allow_alloc_nodes;
			break;
		case SORTID_NODES:
#ifdef HAVE_BG
			convert_num_unit((float)part_ptr->total_nodes, tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
#else
			sprintf(tmp_cnt, "%u", part_ptr->total_nodes);
#endif
			temp_char = tmp_cnt;
			break;
		case SORTID_NODES_MAX:
			limit_set = part_ptr->max_nodes;
			break;
		case SORTID_NODES_MIN:
			limit_set = part_ptr->min_nodes;
			break;
		case SORTID_NODE_INX:
			break;
		case SORTID_ONLY_LINE:
			break;
		case SORTID_PRIORITY:
			convert_num_unit((float)part_ptr->priority,
					 time_buf, sizeof(time_buf), UNIT_NONE);
			temp_char = time_buf;
			break;
		case SORTID_REASON:
			sview_part_sub = list_peek(sview_part_info->sub_list);
			temp_char = sview_part_sub->features;
			break;
		case SORTID_ROOT:
			yes_no = part_ptr->root_only;
			break;
		case SORTID_SHARE:
			if(part_ptr->max_share & SHARED_FORCE) {
				snprintf(tmp_buf, sizeof(tmp_buf), "force:%u",
					 (part_ptr->max_share
					  & ~(SHARED_FORCE)));
				temp_char = tmp_buf;
			} else if(part_ptr->max_share == 0)
				temp_char = "exclusive";
			else if(part_ptr->max_share > 1) {
				snprintf(tmp_buf, sizeof(tmp_buf), "yes:%u",
					 part_ptr->max_share);
				temp_char = tmp_buf;
			} else
				temp_char = "no";
			break;
		case SORTID_TMP_DISK:
			convert_num_unit(
				(float)other_part_sub.disk_total,
				time_buf, sizeof(time_buf), UNIT_NONE);
			temp_char = time_buf;
			break;
		case SORTID_TIMELIMIT:
			limit_set = part_ptr->max_time;
			break;
		default:
			break;
		}

		if(up_down != -1) {
			if(up_down)
				temp_char = "up";
			else
				temp_char = "down";
			up_down = -1;
		} if(yes_no != -1) {
			if(yes_no)
				temp_char = "yes";
			else
				temp_char = "no";
			yes_no = -1;
		} else if(limit_set != NO_VAL) {
			if (limit_set == (uint32_t) INFINITE)
				temp_char = "infinite";
			else {
				convert_num_unit(
					(float)limit_set,
					time_buf, sizeof(time_buf), UNIT_NONE);
				temp_char = time_buf;
			}
			limit_set = NO_VAL;
		}

		if(temp_char) {
			add_display_treestore_line(
				update, treestore, &iter,
				find_col_name(display_data_part,
					      i),
				temp_char);
			if(i == SORTID_NODES) {
				add_display_treestore_line(
					update, treestore, &iter,
					"Nodes (Allocated/Idle/Other)",
					ind_cnt);
			}
			temp_char = NULL;
		}
	}
}

static void _update_part_record(sview_part_info_t *sview_part_info,
				GtkTreeStore *treestore,
				GtkTreeIter *iter)
{
	char time_buf[20], tmp_buf[20];
	char tmp_cnt[8];
	char *temp_char = NULL;
	partition_info_t *part_ptr = sview_part_info->part_ptr;
	GtkTreeIter sub_iter;
	int childern = 0;

	gtk_tree_store_set(treestore, iter, SORTID_COLOR,
			   sview_colors[sview_part_info->color_inx], -1);

	gtk_tree_store_set(treestore, iter, SORTID_NAME, part_ptr->name, -1);

	if(part_ptr->default_part)
		temp_char = "yes";
	else
		temp_char = "no";
	gtk_tree_store_set(treestore, iter, SORTID_DEFAULT, temp_char, -1);

	if(part_ptr->hidden)
		temp_char = "yes";
	else
		temp_char = "no";
	gtk_tree_store_set(treestore, iter, SORTID_HIDDEN, temp_char, -1);

	if (part_ptr->state_up == PARTITION_UP)
		temp_char = "up";
	else if (part_ptr->state_up == PARTITION_DOWN)
		temp_char = "down";
	else if (part_ptr->state_up == PARTITION_INACTIVE)
		temp_char = "inact";
	else if (part_ptr->state_up == PARTITION_DRAIN)
		temp_char = "drain";
	else
		temp_char = "unk";

	gtk_tree_store_set(treestore, iter, SORTID_PART_STATE, temp_char, -1);

	if (part_ptr->max_time == INFINITE)
		snprintf(time_buf, sizeof(time_buf), "infinite");
	else {
		secs2time_str((part_ptr->max_time * 60),
			      time_buf, sizeof(time_buf));
	}

	gtk_tree_store_set(treestore, iter, SORTID_TIMELIMIT, time_buf, -1);

	_build_min_max_32_string(time_buf, sizeof(time_buf),
			      part_ptr->min_nodes,
			      part_ptr->max_nodes, true);
	gtk_tree_store_set(treestore, iter, SORTID_JOB_SIZE, time_buf, -1);

	convert_num_unit((float)part_ptr->priority,
			 time_buf, sizeof(time_buf), UNIT_NONE);
	gtk_tree_store_set(treestore, iter, SORTID_PRIORITY,
			   time_buf, -1);

	if (part_ptr->min_nodes == (uint32_t) INFINITE)
		snprintf(time_buf, sizeof(time_buf), "infinite");
	else {
		convert_num_unit((float)part_ptr->min_nodes,
				 time_buf, sizeof(time_buf), UNIT_NONE);
	}
	gtk_tree_store_set(treestore, iter, SORTID_NODES_MIN,
			   time_buf, -1);
	if (part_ptr->max_nodes == (uint32_t) INFINITE)
		snprintf(time_buf, sizeof(time_buf), "infinite");
	else {
		convert_num_unit((float)part_ptr->max_nodes,
				 time_buf, sizeof(time_buf), UNIT_NONE);
	}
	gtk_tree_store_set(treestore, iter, SORTID_NODES_MAX,
			   time_buf, -1);

	if(part_ptr->root_only)
		temp_char = "yes";
	else
		temp_char = "no";
	gtk_tree_store_set(treestore, iter, SORTID_ROOT, temp_char, -1);

	if(part_ptr->max_share & SHARED_FORCE) {
		snprintf(tmp_buf, sizeof(tmp_buf), "force:%u",
			 (part_ptr->max_share & ~(SHARED_FORCE)));
		temp_char = tmp_buf;
	} else if(part_ptr->max_share == 0)
		temp_char = "exclusive";
	else if(part_ptr->max_share > 1) {
		snprintf(tmp_buf, sizeof(tmp_buf), "yes:%u",
			 part_ptr->max_share);
		temp_char = tmp_buf;
	} else
		temp_char = "no";
	gtk_tree_store_set(treestore, iter, SORTID_SHARE, temp_char, -1);

	if(part_ptr->allow_groups)
		temp_char = part_ptr->allow_groups;
	else
		temp_char = "all";
	gtk_tree_store_set(treestore, iter, SORTID_GROUPS, temp_char, -1);

#ifdef HAVE_BG
	convert_num_unit((float)part_ptr->total_nodes, tmp_cnt,
			 sizeof(tmp_cnt), UNIT_NONE);
#else
	sprintf(tmp_cnt, "%u", part_ptr->total_nodes);
#endif
	gtk_tree_store_set(treestore, iter, SORTID_NODES, tmp_cnt, -1);

#ifdef HAVE_BG
	convert_num_unit((float)part_ptr->total_cpus, tmp_cnt,
			 sizeof(tmp_cnt), UNIT_NONE);
#else
	sprintf(tmp_cnt, "%u", part_ptr->total_cpus);
#endif
	gtk_tree_store_set(treestore, iter, SORTID_CPUS, tmp_cnt, -1);

	gtk_tree_store_set(treestore, iter, SORTID_NODELIST,
			   part_ptr->nodes, -1);

	gtk_tree_store_set(treestore, iter,
			   SORTID_NODE_INX, part_ptr->node_inx, -1);

	gtk_tree_store_set(treestore, iter, SORTID_ONLY_LINE, 0, -1);
	/* clear out info for the main listing */
	gtk_tree_store_set(treestore, iter, SORTID_NODE_STATE, "", -1);
	gtk_tree_store_set(treestore, iter, SORTID_NODE_STATE_NUM, -1, -1);
	gtk_tree_store_set(treestore, iter, SORTID_TMP_DISK, "", -1);
	gtk_tree_store_set(treestore, iter, SORTID_MEM, "", -1);
	gtk_tree_store_set(treestore, iter, SORTID_UPDATED, 1, -1);
	gtk_tree_store_set(treestore, iter, SORTID_FEATURES, "", -1);
	gtk_tree_store_set(treestore, iter, SORTID_REASON, "", -1);

	childern = gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),
						&sub_iter, iter);
	if(gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),
					&sub_iter, iter))
		_subdivide_part(sview_part_info,
				GTK_TREE_MODEL(treestore), &sub_iter, iter);
	else
		_subdivide_part(sview_part_info,
				GTK_TREE_MODEL(treestore), NULL, iter);

	return;
}

static void _update_part_sub_record(sview_part_sub_t *sview_part_sub,
				    GtkTreeStore *treestore, GtkTreeIter *iter)
{
	char tmp_cnt[20];
	char *cpu_tmp = NULL;
	char *node_tmp = NULL;
	partition_info_t *part_ptr = sview_part_sub->part_ptr;
	char *upper = NULL, *lower = NULL;
	char tmp[MAXHOSTRANGELEN];

	gtk_tree_store_set(treestore, iter, SORTID_NAME, part_ptr->name, -1);

	upper = node_state_string(sview_part_sub->node_state);
	lower = str_tolower(upper);
	gtk_tree_store_set(treestore, iter, SORTID_NODE_STATE,
			   lower, -1);
	xfree(lower);

	gtk_tree_store_set(treestore, iter, SORTID_NODE_STATE_NUM,
			   sview_part_sub->node_state, -1);

	if((sview_part_sub->node_state & NODE_STATE_BASE) == NODE_STATE_MIXED) {
		if(sview_part_sub->cpu_alloc_cnt) {
			convert_num_unit((float)sview_part_sub->cpu_alloc_cnt,
					 tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
			xstrfmtcat(cpu_tmp, "Alloc:%s", tmp_cnt);
#ifdef HAVE_BG
			convert_num_unit((float)(sview_part_sub->cpu_alloc_cnt
						 / cpus_per_node),
					 tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
			xstrfmtcat(node_tmp, "Alloc:%s", tmp_cnt);
#endif
		}
		if(sview_part_sub->cpu_error_cnt) {
			convert_num_unit((float)sview_part_sub->cpu_error_cnt,
					 tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
			if(cpu_tmp)
				xstrcat(cpu_tmp, " ");
			xstrfmtcat(cpu_tmp, "Err:%s", tmp_cnt);
#ifdef HAVE_BG
			convert_num_unit((float)(sview_part_sub->cpu_error_cnt
						 / cpus_per_node),
					 tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
			if(node_tmp)
				xstrcat(node_tmp, " ");
			xstrfmtcat(node_tmp, "Err:%s", tmp_cnt);
#endif
		}
		if(sview_part_sub->cpu_idle_cnt) {
			convert_num_unit((float)sview_part_sub->cpu_idle_cnt,
					 tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
			if(cpu_tmp)
				xstrcat(cpu_tmp, " ");
			xstrfmtcat(cpu_tmp, "Idle:%s", tmp_cnt);
#ifdef HAVE_BG
			convert_num_unit((float)(sview_part_sub->cpu_idle_cnt
						 / cpus_per_node),
					 tmp_cnt,
					 sizeof(tmp_cnt), UNIT_NONE);
			if(node_tmp)
				xstrcat(node_tmp, " ");
			xstrfmtcat(node_tmp, "Idle:%s", tmp_cnt);
#endif
		}
	} else {
		cpu_tmp = xmalloc(20);
		convert_num_unit((float)sview_part_sub->cpu_idle_cnt,
				 cpu_tmp, 20, UNIT_NONE);
	}
	gtk_tree_store_set(treestore, iter, SORTID_CPUS, cpu_tmp, -1);
	xfree(cpu_tmp);

	convert_num_unit((float)sview_part_sub->disk_total, tmp_cnt,
			 sizeof(tmp_cnt), UNIT_NONE);
	gtk_tree_store_set(treestore, iter, SORTID_TMP_DISK, tmp_cnt, -1);

	convert_num_unit((float)sview_part_sub->mem_total, tmp_cnt,
			 sizeof(tmp_cnt), UNIT_MEGA);
	gtk_tree_store_set(treestore, iter, SORTID_MEM, tmp_cnt, -1);

	if(!node_tmp) {
		convert_num_unit((float)sview_part_sub->node_cnt, tmp_cnt,
				 sizeof(tmp_cnt), UNIT_NONE);
		node_tmp = xstrdup(tmp_cnt);
	}
	gtk_tree_store_set(treestore, iter, SORTID_NODES, node_tmp, -1);
	xfree(node_tmp);

	hostlist_ranged_string(sview_part_sub->hl, sizeof(tmp), tmp);
	gtk_tree_store_set(treestore, iter, SORTID_NODELIST,
			   tmp, -1);
	gtk_tree_store_set(treestore, iter, SORTID_UPDATED, 1, -1);

	gtk_tree_store_set(treestore, iter, SORTID_FEATURES,
			   sview_part_sub->features, -1);
	gtk_tree_store_set(treestore, iter, SORTID_REASON,
			   sview_part_sub->reason, -1);

	return;
}

static void _append_part_record(sview_part_info_t *sview_part_info,
				GtkTreeStore *treestore, GtkTreeIter *iter,
				int line)
{
	gtk_tree_store_append(treestore, iter, NULL);
	gtk_tree_store_set(treestore, iter, SORTID_POS, line, -1);
	_update_part_record(sview_part_info, treestore, iter);
}

static void _append_part_sub_record(sview_part_sub_t *sview_part_sub,
				    GtkTreeStore *treestore, GtkTreeIter *iter,
				    int line)
{
	GtkTreeIter sub_iter;

	gtk_tree_store_append(treestore, &sub_iter, iter);
	gtk_tree_store_set(treestore, &sub_iter, SORTID_POS, line, -1);
	_update_part_sub_record(sview_part_sub, treestore, &sub_iter);
}

static void _update_info_part(List info_list,
			      GtkTreeView *tree_view)
{
	GtkTreePath *path = gtk_tree_path_new_first();
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreeIter iter;
	partition_info_t *part_ptr = NULL;
	int line = 0;
	char *host = NULL, *part_name = NULL;
	ListIterator itr = NULL;
	sview_part_info_t *sview_part_info = NULL;

	/* get the iter, or find out the list is empty goto add */
	if (gtk_tree_model_get_iter(model, &iter, path)) {
		/* make sure all the partitions are still here */
		while(1) {
			gtk_tree_store_set(GTK_TREE_STORE(model), &iter,
					   SORTID_UPDATED, 0, -1);
			if(!gtk_tree_model_iter_next(model, &iter)) {
				break;
			}
		}
	}

	itr = list_iterator_create(info_list);
	while ((sview_part_info = (sview_part_info_t*) list_next(itr))) {
		part_ptr = sview_part_info->part_ptr;
		/* get the iter, or find out the list is empty goto add */
		if (!gtk_tree_model_get_iter(model, &iter, path)) {
			goto adding;
		}
		line = 0;
		while(1) {
			/* search for the jobid and check to see if
			   it is in the list */
			gtk_tree_model_get(model, &iter, SORTID_NAME,
					   &part_name, -1);
			if(!strcmp(part_name, part_ptr->name)) {
				/* update with new info */
				g_free(part_name);
				_update_part_record(sview_part_info,
						    GTK_TREE_STORE(model),
						    &iter);
				goto found;
			}
			g_free(part_name);

			line++;
			if(!gtk_tree_model_iter_next(model, &iter)) {
				break;
			}
		}
	adding:
		_append_part_record(sview_part_info, GTK_TREE_STORE(model),
				    &iter, line);
	found:
		;
	}
	list_iterator_destroy(itr);
	if(host)
		free(host);

	gtk_tree_path_free(path);
	/* remove all old partitions */
	remove_old(model, SORTID_UPDATED);
	return;
}

static void _part_info_list_del(void *object)
{
	sview_part_info_t *sview_part_info = (sview_part_info_t *)object;

	if (sview_part_info) {
		if(sview_part_info->sub_list)
			list_destroy(sview_part_info->sub_list);
		xfree(sview_part_info);
	}
}

static void _destroy_part_sub(void *object)
{
	sview_part_sub_t *sview_part_sub = (sview_part_sub_t *)object;

	if (sview_part_sub) {
		xfree(sview_part_sub->features);
		xfree(sview_part_sub->reason);
		if(sview_part_sub->hl)
			hostlist_destroy(sview_part_sub->hl);
		if(sview_part_sub->node_ptr_list)
			list_destroy(sview_part_sub->node_ptr_list);
		xfree(sview_part_sub);
	}
}

static void _update_sview_part_sub(sview_part_sub_t *sview_part_sub,
				   node_info_t *node_ptr,
				   int node_scaling)
{
	int cpus_per_node = 1;
	int idle_cpus = node_ptr->cpus;
	uint16_t err_cpus = 0, alloc_cpus = 0;

	if(node_scaling)
		cpus_per_node = node_ptr->cpus / node_scaling;

	xassert(sview_part_sub);
	xassert(sview_part_sub->node_ptr_list);
	xassert(sview_part_sub->hl);

	if (sview_part_sub->node_cnt == 0) {	/* first node added */
		sview_part_sub->node_state = node_ptr->node_state;
		sview_part_sub->features   = xstrdup(node_ptr->features);
		sview_part_sub->reason     = xstrdup(node_ptr->reason);
	} else if (hostlist_find(sview_part_sub->hl, node_ptr->name) != -1) {
		/* we already have this node in this record,
		 * just return, don't duplicate */
		g_print("already been here\n");
		return;
	}

	if((sview_part_sub->node_state & NODE_STATE_BASE) == NODE_STATE_MIXED) {
		slurm_get_select_nodeinfo(node_ptr->select_nodeinfo,
					  SELECT_NODEDATA_SUBCNT,
					  NODE_STATE_ALLOCATED,
					  &alloc_cpus);
#ifdef HAVE_BG
		if(!alloc_cpus
		   && (IS_NODE_ALLOCATED(node_ptr)
		       || IS_NODE_COMPLETING(node_ptr)))
			alloc_cpus = node_ptr->cpus;
		else
			alloc_cpus *= cpus_per_node;
#endif
		idle_cpus -= alloc_cpus;

		slurm_get_select_nodeinfo(node_ptr->select_nodeinfo,
					  SELECT_NODEDATA_SUBCNT,
					  NODE_STATE_ERROR,
					  &err_cpus);
#ifdef HAVE_BG
		err_cpus *= cpus_per_node;
#endif
		idle_cpus -= err_cpus;
	}

	sview_part_sub->cpu_alloc_cnt += alloc_cpus;
	sview_part_sub->cpu_error_cnt += err_cpus;
	sview_part_sub->cpu_idle_cnt += idle_cpus;
	sview_part_sub->disk_total += node_ptr->tmp_disk;
	sview_part_sub->mem_total  += node_ptr->real_memory;
	sview_part_sub->node_cnt   += node_scaling;
	list_append(sview_part_sub->node_ptr_list, node_ptr);
	hostlist_push(sview_part_sub->hl, node_ptr->name);
}

/*
 * _create_sview_part_sub - create an sview_part_sub record for
 *                          the given partition
 * sview_part_sub OUT     - ptr to an inited sview_part_sub_t
 */
static sview_part_sub_t *_create_sview_part_sub(partition_info_t *part_ptr,
						node_info_t *node_ptr,
						int node_scaling)
{
	sview_part_sub_t *sview_part_sub_ptr =
		xmalloc(sizeof(sview_part_sub_t));

	if (!part_ptr) {
		g_print("got no part_ptr!\n");
		xfree(sview_part_sub_ptr);
		return NULL;
	}
	if (!node_ptr) {
		g_print("got no node_ptr!\n");
		xfree(sview_part_sub_ptr);
		return NULL;
	}
	sview_part_sub_ptr->part_ptr = part_ptr;
	sview_part_sub_ptr->hl = hostlist_create("");
	sview_part_sub_ptr->node_ptr_list = list_create(NULL);

	_update_sview_part_sub(sview_part_sub_ptr, node_ptr, node_scaling);

	return sview_part_sub_ptr;
}

static int _insert_sview_part_sub(sview_part_info_t *sview_part_info,
				  partition_info_t *part_ptr,
				  node_info_t *node_ptr,
				  int node_scaling)
{
	sview_part_sub_t *sview_part_sub = NULL;
	ListIterator itr = list_iterator_create(sview_part_info->sub_list);

	while((sview_part_sub = list_next(itr))) {
		if(sview_part_sub->node_state
		   == node_ptr->node_state) {
			_update_sview_part_sub(sview_part_sub,
					       node_ptr,
					       node_scaling);
			break;
		}
	}
	list_iterator_destroy(itr);

	if(!sview_part_sub) {
		if((sview_part_sub = _create_sview_part_sub(
			   part_ptr, node_ptr, node_scaling)))
			list_push(sview_part_info->sub_list,
				  sview_part_sub);
	}
	return SLURM_SUCCESS;
}

/*
 * _create_sview_part_info - create an sview_part_info record for
 *                           the given partition
 * part_ptr IN             - pointer to partition record to add
 * sview_part_info OUT     - ptr to an inited sview_part_info_t
 */
static sview_part_info_t *_create_sview_part_info(partition_info_t* part_ptr)
{
	sview_part_info_t *sview_part_info =
		xmalloc(sizeof(sview_part_info_t));


	sview_part_info->part_ptr = part_ptr;
	sview_part_info->sub_list = list_create(_destroy_part_sub);
	return sview_part_info;
}

static int _sview_part_sort_aval_dec(sview_part_info_t* rec_a,
				     sview_part_info_t* rec_b)
{
	int size_a = rec_a->part_ptr->total_nodes;
	int size_b = rec_b->part_ptr->total_nodes;

	if (size_a < size_b)
		return -1;
	else if (size_a > size_b)
		return 1;

	if(rec_a->part_ptr->nodes && rec_b->part_ptr->nodes) {
		size_a = strcmp(rec_a->part_ptr->nodes, rec_b->part_ptr->nodes);
		if (size_a < 0)
			return -1;
		else if (size_a > 0)
			return 1;
	}
	return 0;
}

static int _sview_sub_part_sort(sview_part_sub_t* rec_a,
				sview_part_sub_t* rec_b)
{
	int size_a = rec_a->node_state & NODE_STATE_BASE;
	int size_b = rec_b->node_state & NODE_STATE_BASE;

	if (size_a < size_b)
		return -1;
	else if (size_a > size_b)
		return 1;
	return 0;
}

static List _create_part_info_list(partition_info_msg_t *part_info_ptr,
				   node_info_msg_t *node_info_ptr,
				   int changed)
{
	sview_part_info_t *sview_part_info = NULL;
	partition_info_t *part_ptr = NULL;
	node_info_t *node_ptr = NULL;
	static List info_list = NULL;
	int i, j2;

	if(!changed && info_list) {
		return info_list;
	}

	if(info_list) {
		list_destroy(info_list);
	}
	info_list = list_create(_part_info_list_del);
	if (!info_list) {
		g_print("malloc error\n");
		return NULL;
	}

	for (i=0; i<part_info_ptr->record_count; i++) {
		part_ptr = &(part_info_ptr->partition_array[i]);
		if (!part_ptr->nodes || (part_ptr->nodes[0] == '\0'))
			continue;	/* empty partition */

		sview_part_info = _create_sview_part_info(part_ptr);
		list_append(info_list, sview_part_info);
		sview_part_info->color_inx = i % sview_colors_cnt;

		j2 = 0;
		while(part_ptr->node_inx[j2] >= 0) {
			int i2 = 0;

			for(i2 = part_ptr->node_inx[j2];
			    i2 <= part_ptr->node_inx[j2+1];
			    i2++) {
				node_ptr = &(node_info_ptr->node_array[i2]);
				_insert_sview_part_sub(sview_part_info,
						       part_ptr,
						       node_ptr,
						       g_node_scaling);
			}
			j2 += 2;
		}
		list_sort(sview_part_info->sub_list,
			  (ListCmpF)_sview_sub_part_sort);
	}
	list_sort(info_list, (ListCmpF)_sview_part_sort_aval_dec);

	return info_list;
}

void _display_info_part(List info_list,	popup_info_t *popup_win)
{
	specific_info_t *spec_info = popup_win->spec_info;
	char *name = (char *)spec_info->search_info->gchar_data;
	int found = 0;
	partition_info_t *part_ptr = NULL;
	GtkTreeView *treeview = NULL;
	ListIterator itr = NULL;
	sview_part_info_t *sview_part_info = NULL;
	int update = 0;
	int j = 0;
	int first_time = 0;

	if(!spec_info->search_info->gchar_data) {
		//info = xstrdup("No pointer given!");
		goto finished;
	}
	if(!list_count(popup_win->grid_button_list))
		first_time = 1;

need_refresh:
	if(!spec_info->display_widget) {
		treeview = create_treeview_2cols_attach_to_table(
			popup_win->table);
		spec_info->display_widget =
			gtk_widget_ref(GTK_WIDGET(treeview));
	} else {
		treeview = GTK_TREE_VIEW(spec_info->display_widget);
		update = 1;
	}

	itr = list_iterator_create(info_list);
	while ((sview_part_info = (sview_part_info_t*) list_next(itr))) {
		part_ptr = sview_part_info->part_ptr;
		if(!strcmp(part_ptr->name, name)) {
			j=0;
			while(part_ptr->node_inx[j] >= 0) {
				change_grid_color(
					popup_win->grid_button_list,
					part_ptr->node_inx[j],
					part_ptr->node_inx[j+1],
					sview_part_info->color_inx,
					true, 0);
				j += 2;
			}
			_layout_part_record(treeview, sview_part_info, update);
			found = 1;
			break;
		}
	}
	list_iterator_destroy(itr);
	post_setup_popup_grid_list(popup_win);

	if(!found) {
		if(!popup_win->not_found) {
			char *temp = "PARTITION DOESN'T EXSIST\n";
			GtkTreeIter iter;
			GtkTreeModel *model = NULL;

			/* only time this will be run so no update */
			model = gtk_tree_view_get_model(treeview);
			add_display_treestore_line(0,
						   GTK_TREE_STORE(model),
						   &iter,
						   temp, "");
		}
		popup_win->not_found = true;
	} else {
		if(popup_win->not_found) {
			popup_win->not_found = false;
			gtk_widget_destroy(spec_info->display_widget);

			goto need_refresh;
		}
	}
	gtk_widget_show(spec_info->display_widget);

finished:

	return;
}

extern void refresh_part(GtkAction *action, gpointer user_data)
{
	popup_info_t *popup_win = (popup_info_t *)user_data;
	xassert(popup_win != NULL);
	xassert(popup_win->spec_info != NULL);
	xassert(popup_win->spec_info->title != NULL);
	popup_win->force_refresh = 1;
	specific_info_part(popup_win);
}

extern int get_new_info_part(partition_info_msg_t **part_ptr, int force)
{
	static partition_info_msg_t *part_info_ptr = NULL;
	static partition_info_msg_t *new_part_ptr = NULL;
	int error_code = SLURM_NO_CHANGE_IN_DATA;
	time_t now = time(NULL);
	static time_t last;
	static bool changed = 0;

	if(!force && ((now - last) < global_sleep_time)) {
		if(*part_ptr != part_info_ptr)
			error_code = SLURM_SUCCESS;
		*part_ptr = part_info_ptr;
		if(changed)
			return SLURM_SUCCESS;
		return error_code;
	}
	last = now;
	if (part_info_ptr) {
		error_code = slurm_load_partitions(part_info_ptr->last_update,
						   &new_part_ptr, SHOW_ALL);
		if (error_code == SLURM_SUCCESS) {
			slurm_free_partition_info_msg(part_info_ptr);
			changed = 1;
		} else if (slurm_get_errno() == SLURM_NO_CHANGE_IN_DATA) {
			error_code = SLURM_NO_CHANGE_IN_DATA;
			new_part_ptr = part_info_ptr;
				changed = 0;
		}
	} else {
		error_code = slurm_load_partitions((time_t) NULL,
						   &new_part_ptr, SHOW_ALL);
		changed = 1;
	}

	part_info_ptr = new_part_ptr;

	if(*part_ptr != part_info_ptr)
		error_code = SLURM_SUCCESS;

	*part_ptr = new_part_ptr;
	return error_code;
}

extern GtkListStore *create_model_part(int type)
{
	GtkListStore *model = NULL;
	GtkTreeIter iter;
	char *upper = NULL, *lower = NULL;
	int i=0;
	switch(type) {
	case SORTID_DEFAULT:
		model = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "yes",
				   1, SORTID_DEFAULT,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "no",
				   1, SORTID_DEFAULT,
				   -1);

		break;
	case SORTID_HIDDEN:
		model = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "yes",
				   1, SORTID_HIDDEN,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "no",
				   1, SORTID_HIDDEN,
				   -1);

		break;
	case SORTID_PRIORITY:
	case SORTID_TIMELIMIT:
	case SORTID_NODES_MIN:
	case SORTID_NODES_MAX:
		break;
	case SORTID_ROOT:
		model = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "yes",
				   1, SORTID_ROOT,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "no",
				   1, SORTID_ROOT,
				   -1);
		break;
	case SORTID_SHARE:
		model = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "force",
				   1, SORTID_SHARE,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "no",
				   1, SORTID_SHARE,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "yes",
				   1, SORTID_SHARE,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "exclusive",
				   1, SORTID_SHARE,
				   -1);
		break;
	case SORTID_GROUPS:
		break;
	case SORTID_NODELIST:
		break;
	case SORTID_PART_STATE:
		model = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "up",
				   1, SORTID_PART_STATE,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "down",
				   1, SORTID_PART_STATE,
				   -1);
		break;
	case SORTID_NODE_STATE:
		model = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "drain",
				   1, SORTID_NODE_STATE,
				   -1);
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter,
				   0, "resume",
				   1, SORTID_NODE_STATE,
				   -1);
		for(i = 0; i < NODE_STATE_END; i++) {
			upper = node_state_string(i);
			if(!strcmp(upper, "UNKNOWN"))
				continue;

			gtk_list_store_append(model, &iter);
			lower = str_tolower(upper);
			gtk_list_store_set(model, &iter,
					   0, lower,
					   1, SORTID_NODE_STATE,
					   -1);
			xfree(lower);
		}

		break;

	}
	return model;
}

extern void admin_edit_part(GtkCellRendererText *cell,
			    const char *path_string,
			    const char *new_text,
			    gpointer data)
{
	GtkTreeStore *treestore = GTK_TREE_STORE(data);
	GtkTreePath *path = gtk_tree_path_new_from_string(path_string);
	GtkTreeIter iter;
	update_part_msg_t *part_msg = xmalloc(sizeof(update_part_msg_t));

	char *temp = NULL;
	char *old_text = NULL;
	const char *type = NULL;
	int column = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(cell),
						       "column"));

	if(!new_text || !strcmp(new_text, ""))
		goto no_input;

	gtk_tree_model_get_iter(GTK_TREE_MODEL(treestore), &iter, path);

	if(column != SORTID_NODE_STATE) {
		slurm_init_part_desc_msg(part_msg);
		gtk_tree_model_get(GTK_TREE_MODEL(treestore), &iter,
				   SORTID_NAME, &temp,
				   column, &old_text,
				   -1);
		part_msg->name = xstrdup(temp);
		g_free(temp);
	}

	type = _set_part_msg(part_msg, new_text, column);
	if(global_edit_error)
		goto print_error;
	if(got_edit_signal) {
		temp = got_edit_signal;
		got_edit_signal = NULL;
		admin_part(GTK_TREE_MODEL(treestore), &iter, temp);
		xfree(temp);
		goto no_input;
	}

	if(got_features_edit_signal) {
		admin_part(GTK_TREE_MODEL(treestore), &iter, (char *)type);
		goto no_input;
	}

	if(column != SORTID_NODE_STATE && column != SORTID_FEATURES ) {
		if(old_text && !strcmp(old_text, new_text)) {
			temp = g_strdup_printf("No change in value.");
		} else if(slurm_update_partition(part_msg) == SLURM_SUCCESS) {
			gtk_tree_store_set(treestore, &iter, column,
					   new_text, -1);
			temp = g_strdup_printf("Partition %s %s changed to %s",
					       part_msg->name,
					       type,
					       new_text);
		} else {
		print_error:
			temp = g_strdup_printf("Partition %s %s can't be "
					       "set to %s",
					       part_msg->name,
					       type,
					       new_text);
		}
		display_edit_note(temp);
		g_free(temp);
	}
no_input:
	slurm_free_update_part_msg(part_msg);
	gtk_tree_path_free (path);
	g_free(old_text);
	g_static_mutex_unlock(&sview_mutex);
}

extern void get_info_part(GtkTable *table, display_data_t *display_data)
{
	int part_error_code = SLURM_SUCCESS;
	int node_error_code = SLURM_SUCCESS;
	static int view = -1;
	static partition_info_msg_t *part_info_ptr = NULL;
	static node_info_msg_t *node_info_ptr = NULL;
	char error_char[100];
	GtkWidget *label = NULL;
	GtkTreeView *tree_view = NULL;
	static GtkWidget *display_widget = NULL;
	List info_list = NULL;
	int changed = 1;
	int j=0;
	sview_part_info_t *sview_part_info = NULL;
	partition_info_t *part_ptr = NULL;
	ListIterator itr = NULL;

	if(display_data)
		local_display_data = display_data;
	if(!table) {
		display_data_part->set_menu = local_display_data->set_menu;
		return;
	}
	if(display_widget && toggled) {
		gtk_widget_destroy(display_widget);
		display_widget = NULL;
		goto display_it;
	}

	if((part_error_code = get_new_info_part(&part_info_ptr, force_refresh))
	   == SLURM_NO_CHANGE_IN_DATA) {
		// just goto the new info node
	} else 	if (part_error_code != SLURM_SUCCESS) {
		if(view == ERROR_VIEW)
			goto end_it;
		if(display_widget)
			gtk_widget_destroy(display_widget);
		view = ERROR_VIEW;
		snprintf(error_char, 100, "slurm_load_partitions: %s",
			slurm_strerror(slurm_get_errno()));
		label = gtk_label_new(error_char);
		display_widget = gtk_widget_ref(GTK_WIDGET(label));
		gtk_table_attach_defaults(table, label, 0, 1, 0, 1);
		gtk_widget_show(label);
		goto end_it;
	}

	if((node_error_code = get_new_info_node(&node_info_ptr, force_refresh))
	   == SLURM_NO_CHANGE_IN_DATA) {
		if((!display_widget || view == ERROR_VIEW)
		   || (part_error_code != SLURM_NO_CHANGE_IN_DATA))
			goto display_it;
		changed = 0;
	} else if (node_error_code != SLURM_SUCCESS) {
		if(view == ERROR_VIEW)
			goto end_it;
		if(display_widget)
			gtk_widget_destroy(display_widget);
		view = ERROR_VIEW;
		snprintf(error_char, 100, "slurm_load_node: %s",
			slurm_strerror(slurm_get_errno()));
		label = gtk_label_new(error_char);
		display_widget = gtk_widget_ref(GTK_WIDGET(label));
		gtk_table_attach_defaults(table, label, 0, 1, 0, 1);
		gtk_widget_show(label);
		goto end_it;
	}

display_it:

	info_list = _create_part_info_list(part_info_ptr,
					   node_info_ptr,
					   changed);
	if(!info_list)
		return;

	/* set up the grid */
	itr = list_iterator_create(info_list);
	while ((sview_part_info = list_next(itr))) {
		part_ptr = sview_part_info->part_ptr;
		j=0;
		while(part_ptr->node_inx[j] >= 0) {
				change_grid_color(grid_button_list,
						  part_ptr->node_inx[j],
						  part_ptr->node_inx[j+1],
						  sview_part_info->color_inx,
						  true, 0);
			j += 2;
		}
	}
	list_iterator_destroy(itr);
	change_grid_color(grid_button_list, -1, -1, MAKE_WHITE, true, 0);
	if(grid_speedup) {
		gtk_widget_set_sensitive(GTK_WIDGET(main_grid_table), 0);
		gtk_widget_set_sensitive(GTK_WIDGET(main_grid_table), 1);
	}

	if(view == ERROR_VIEW && display_widget) {
		gtk_widget_destroy(display_widget);
		display_widget = NULL;
	}
	if(!display_widget) {
		tree_view = create_treeview(local_display_data,
					    &grid_button_list);

		display_widget = gtk_widget_ref(GTK_WIDGET(tree_view));
		gtk_table_attach_defaults(table,
					  GTK_WIDGET(tree_view),
					  0, 1, 0, 1);
		/* since this function sets the model of the tree_view
		   to the treestore we don't really care about
		   the return value */
		create_treestore(tree_view, display_data_part,
				 SORTID_CNT, SORTID_NAME, SORTID_COLOR);
	}

	if(gtk_tree_selection_count_selected_rows(
		   gtk_tree_view_get_selection(
			   GTK_TREE_VIEW(display_widget)))) {
		GtkTreePath *path = NULL;
		GtkTreeViewColumn *focus_column = NULL;
		/* highlight the correct nodes from the last selection */
		gtk_tree_view_get_cursor(GTK_TREE_VIEW(display_widget),
					 &path, &focus_column);
		if(path)
			highlight_grid(GTK_TREE_VIEW(display_widget), path,
				       SORTID_NODE_INX, grid_button_list);
	}

	view = INFO_VIEW;
	_update_info_part(info_list, GTK_TREE_VIEW(display_widget));
end_it:
	toggled = FALSE;
	force_refresh = FALSE;

	return;
}

extern void specific_info_part(popup_info_t *popup_win)
{
	int part_error_code = SLURM_SUCCESS;
	int node_error_code = SLURM_SUCCESS;
	static partition_info_msg_t *part_info_ptr = NULL;
	static node_info_msg_t *node_info_ptr = NULL;
	specific_info_t *spec_info = popup_win->spec_info;
	char error_char[100];
	GtkWidget *label = NULL;
	GtkTreeView *tree_view = NULL;
	List info_list = NULL;
	List send_info_list = NULL;
	int changed = 1;
	int j=0, i=-1;
	sview_part_info_t *sview_part_info_ptr = NULL;
	partition_info_t *part_ptr = NULL;
	ListIterator itr = NULL;
	hostset_t hostset = NULL;

	if(!spec_info->display_widget)
		setup_popup_info(popup_win, display_data_part, SORTID_CNT);

	if(spec_info->display_widget && popup_win->toggled) {
		gtk_widget_destroy(spec_info->display_widget);
		spec_info->display_widget = NULL;
		goto display_it;
	}

	if((part_error_code = get_new_info_part(&part_info_ptr,
						popup_win->force_refresh))
	   == SLURM_NO_CHANGE_IN_DATA)  {

	} else if (part_error_code != SLURM_SUCCESS) {
		if(spec_info->view == ERROR_VIEW)
			goto end_it;
		if(spec_info->display_widget) {
			gtk_widget_destroy(spec_info->display_widget);
			spec_info->display_widget = NULL;
		}
		spec_info->view = ERROR_VIEW;
		snprintf(error_char, 100, "slurm_load_partitions: %s",
			 slurm_strerror(slurm_get_errno()));
		label = gtk_label_new(error_char);
		spec_info->display_widget = gtk_widget_ref(GTK_WIDGET(label));
		gtk_table_attach_defaults(popup_win->table, label, 0, 1, 0, 1);
		gtk_widget_show(label);
		goto end_it;
	}

	if((node_error_code = get_new_info_node(&node_info_ptr,
						popup_win->force_refresh))
	   == SLURM_NO_CHANGE_IN_DATA) {
		if((!spec_info->display_widget
		    || spec_info->view == ERROR_VIEW)
		   || (part_error_code != SLURM_NO_CHANGE_IN_DATA))
			goto display_it;
		changed = 0;
	} else if (node_error_code != SLURM_SUCCESS) {
		if(spec_info->view == ERROR_VIEW)
			goto end_it;
		if(spec_info->display_widget)
			gtk_widget_destroy(spec_info->display_widget);
		spec_info->view = ERROR_VIEW;
		snprintf(error_char, 100, "slurm_load_node: %s",
			slurm_strerror(slurm_get_errno()));
		label = gtk_label_new(error_char);
		spec_info->display_widget = gtk_widget_ref(GTK_WIDGET(label));
		gtk_table_attach_defaults(popup_win->table, label, 0, 1, 0, 1);
		gtk_widget_show(label);
		goto end_it;
	}

display_it:

	info_list = _create_part_info_list(part_info_ptr,
					   node_info_ptr,
					   changed);
	if(!info_list)
		return;

	if(spec_info->view == ERROR_VIEW && spec_info->display_widget) {
		gtk_widget_destroy(spec_info->display_widget);
		spec_info->display_widget = NULL;
	}

	if(spec_info->type != INFO_PAGE && !spec_info->display_widget) {
		tree_view = create_treeview(local_display_data,
					    &popup_win->grid_button_list);

		spec_info->display_widget =
			gtk_widget_ref(GTK_WIDGET(tree_view));
		gtk_table_attach_defaults(popup_win->table,
					  GTK_WIDGET(tree_view),
					  0, 1, 0, 1);

		/* since this function sets the model of the tree_view
		   to the treestore we don't really care about
		   the return value */
		create_treestore(tree_view, popup_win->display_data,
				 SORTID_CNT, SORTID_NAME, SORTID_COLOR);
	}

	setup_popup_grid_list(popup_win);

	spec_info->view = INFO_VIEW;
	if(spec_info->type == INFO_PAGE) {
		_display_info_part(info_list, popup_win);
		goto end_it;
	}

	/* just linking to another list, don't free the inside, just
	   the list */
	send_info_list = list_create(NULL);

	itr = list_iterator_create(info_list);
	i = -1;
	while ((sview_part_info_ptr = list_next(itr))) {
		i++;
		part_ptr = sview_part_info_ptr->part_ptr;
		switch(spec_info->type) {
		case RESV_PAGE:
		case NODE_PAGE:
			if(!part_ptr->nodes)
				continue;

			if(!(hostset = hostset_create(
				     spec_info->search_info->gchar_data)))
				continue;
			if(!hostset_intersects(hostset, part_ptr->nodes)) {
				hostset_destroy(hostset);
				continue;
			}
			hostset_destroy(hostset);
			break;
		case PART_PAGE:
			switch(spec_info->search_info->search_type) {
			case SEARCH_PARTITION_NAME:
				if(!spec_info->search_info->gchar_data)
					continue;

				if(strcmp(part_ptr->name,
					  spec_info->search_info->gchar_data))
					continue;
				break;
			case SEARCH_PARTITION_STATE:
				if(spec_info->search_info->int_data == NO_VAL)
					continue;
				if(part_ptr->state_up !=
				   spec_info->search_info->int_data)
					continue;
				break;
			default:
				continue;
				break;
			}
			break;
		case BLOCK_PAGE:
		case JOB_PAGE:
			if(!spec_info->search_info->gchar_data)
				continue;

			if(strcmp(part_ptr->name,
				  spec_info->search_info->gchar_data))
				continue;
			break;
		default:
			g_print("Unknown type %d\n", spec_info->type);
			list_iterator_destroy(itr);
			goto end_it;
		}
		list_push(send_info_list, sview_part_info_ptr);
		j=0;
		while(part_ptr->node_inx[j] >= 0) {
			change_grid_color(
				popup_win->grid_button_list,
				part_ptr->node_inx[j],
				part_ptr->node_inx[j+1],
				sview_part_info_ptr->color_inx, true, 0);
			j += 2;
		}
	}
	list_iterator_destroy(itr);
	post_setup_popup_grid_list(popup_win);

	_update_info_part(send_info_list,
			  GTK_TREE_VIEW(spec_info->display_widget));
	list_destroy(send_info_list);
end_it:
	popup_win->toggled = 0;
	popup_win->force_refresh = 0;

	return;
}

extern void set_menus_part(void *arg, void *arg2, GtkTreePath *path, int type)
{
	GtkTreeView *tree_view = (GtkTreeView *)arg;
	popup_info_t *popup_win = (popup_info_t *)arg;
	GtkMenu *menu = (GtkMenu *)arg2;
	List button_list = (List)arg2;

	switch(type) {
	case TAB_CLICKED:
		make_fields_menu(NULL, menu, display_data_part, SORTID_CNT);
		break;
	case ROW_CLICKED:
		make_options_menu(tree_view, path, menu, options_data_part);
		break;
	case ROW_LEFT_CLICKED:
		highlight_grid(tree_view, path, SORTID_NODE_INX, button_list);
		break;
	case FULL_CLICKED:
	{
		GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
		GtkTreeIter iter;
		if (!gtk_tree_model_get_iter(model, &iter, path)) {
			g_error("error getting iter from model\n");
			break;
		}

		popup_all_part(model, &iter, INFO_PAGE);

		break;
	}
	case POPUP_CLICKED:
		make_fields_menu(popup_win, menu,
				 popup_win->display_data, SORTID_CNT);
		break;
	default:
		g_error("UNKNOWN type %d given to set_fields\n", type);
	}
}

extern void popup_all_part(GtkTreeModel *model, GtkTreeIter *iter, int id)
{
	char *name = NULL;
	char *state = NULL;
	char title[100];
	int only_line = 0;
	ListIterator itr = NULL;
	popup_info_t *popup_win = NULL;
	GError *error = NULL;
	GtkTreeIter par_iter;

	gtk_tree_model_get(model, iter, SORTID_NAME, &name, -1);

	switch(id) {
	case JOB_PAGE:
		snprintf(title, 100, "Job(s) in partition %s", name);
		break;
	case RESV_PAGE:
		snprintf(title, 100, "Reservation(s) in partition %s", name);
		break;
	case NODE_PAGE:
		gtk_tree_model_get(model, iter, SORTID_ONLY_LINE,
				   &only_line, -1);
		if(!only_line)
			gtk_tree_model_get(model, iter,
					   SORTID_NODE_STATE, &state, -1);
#ifdef HAVE_BG
		if(!state || !strlen(state))
			snprintf(title, 100,
				 "Base partition(s) in partition %s",
				 name);
		else
			snprintf(title, 100,
				 "Base partition(s) in partition %s "
				 "that are in '%s' state",
				 name, state);
#else
		if(!state || !strlen(state))
			snprintf(title, 100, "Node(s) in partition %s ",
				 name);
		else
			snprintf(title, 100,
				 "Node(s) in partition %s that are in "
				 "'%s' state",
				 name, state);
#endif
		break;
	case BLOCK_PAGE:
		snprintf(title, 100, "Block(s) in partition %s", name);
		break;
	case SUBMIT_PAGE:
		snprintf(title, 100, "Submit job in partition %s", name);
		break;
	case INFO_PAGE:
		snprintf(title, 100, "Full info for partition %s", name);
		break;
	default:
		g_print("part got %d\n", id);
	}

	itr = list_iterator_create(popup_list);
	while((popup_win = list_next(itr))) {
		if(popup_win->spec_info)
			if(!strcmp(popup_win->spec_info->title, title)) {
				break;
			}
	}
	list_iterator_destroy(itr);

	if(!popup_win) {
		if(id == INFO_PAGE)
			popup_win = create_popup_info(id, PART_PAGE, title);
		else
			popup_win = create_popup_info(PART_PAGE, id, title);
	} else {
		g_free(name);
		g_free(state);
		gtk_window_present(GTK_WINDOW(popup_win->popup));
		return;
	}

	/* Pass the model and the structs from the iter so we can always get
	   the current node_inx.
	*/
	popup_win->model = model;
	popup_win->iter = *iter;
	popup_win->node_inx_id = SORTID_NODE_INX;

	switch(id) {
	case JOB_PAGE:
	case BLOCK_PAGE:
	case INFO_PAGE:
		popup_win->spec_info->search_info->gchar_data = name;
		//specific_info_job(popup_win);
		break;
	case RESV_PAGE:
	case NODE_PAGE:
		g_free(name);
		/* we want to include the parent's nodes here not just
		   the subset */
		if(gtk_tree_model_iter_parent(model, &par_iter, iter))
			gtk_tree_model_get(model, &par_iter,
					   SORTID_NODELIST, &name, -1);
		else
			gtk_tree_model_get(model, iter,
					   SORTID_NODELIST, &name, -1);
		popup_win->spec_info->search_info->gchar_data = name;
		if(state && strlen(state)) {
			popup_win->spec_info->search_info->search_type =
				SEARCH_NODE_STATE;
			gtk_tree_model_get(
				model, iter, SORTID_NODE_STATE_NUM,
				&popup_win->spec_info->search_info->int_data,
				-1);
		} else {
			popup_win->spec_info->search_info->search_type =
				SEARCH_NODE_NAME;
		}
		g_free(state);

		//specific_info_node(popup_win);
		break;
	case SUBMIT_PAGE:
		break;
	default:
		g_print("part got unknown type %d\n", id);
	}
	if (!g_thread_create((gpointer)popup_thr, popup_win, FALSE, &error))
	{
		g_printerr ("Failed to create part popup thread: %s\n",
			    error->message);
		return;
	}
}

extern void admin_part(GtkTreeModel *model, GtkTreeIter *iter, char *type)
{
	update_part_msg_t *part_msg = xmalloc(sizeof(update_part_msg_t));
	char *partid = NULL;
	char *nodelist = NULL;
	char tmp_char[100];
	char *temp = NULL;
	int edit_type = 0;
	int response = 0;
	GtkWidget *label = NULL;
	GtkWidget *entry = NULL;
	GtkWidget *popup = gtk_dialog_new_with_buttons(
		type,
		GTK_WINDOW(main_window),
		GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
		NULL);
	gtk_window_set_transient_for(GTK_WINDOW(popup), NULL);

	gtk_tree_model_get(model, iter, SORTID_NAME, &partid, -1);
	gtk_tree_model_get(model, iter, SORTID_NODELIST, &nodelist, -1);
	slurm_init_part_desc_msg(part_msg);

	part_msg->name = xstrdup(partid);

	if(!strcasecmp("Change Part State Up/Down", type)) {
		char *state = NULL;
		label = gtk_dialog_add_button(GTK_DIALOG(popup),
					      GTK_STOCK_YES, GTK_RESPONSE_OK);
		gtk_window_set_default(GTK_WINDOW(popup), label);
		gtk_dialog_add_button(GTK_DIALOG(popup),
				      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL);
		gtk_tree_model_get(model, iter, SORTID_PART_STATE, &state, -1);
		if(!strcasecmp("down", state)) {
			temp = "up";
			part_msg->state_up = PARTITION_UP;
		} else {
			temp = "down";
			part_msg->state_up = PARTITION_DOWN;
		}
		g_free(state);
		snprintf(tmp_char, sizeof(tmp_char),
			 "Are you sure you want to set partition %s %s?",
			 partid, temp);
		label = gtk_label_new(tmp_char);
		edit_type = EDIT_PART_STATE;
	} else if(!strcasecmp("Edit Part", type)) {
		label = gtk_dialog_add_button(GTK_DIALOG(popup),
					      GTK_STOCK_OK, GTK_RESPONSE_OK);
		gtk_window_set_default(GTK_WINDOW(popup), label);
		gtk_dialog_add_button(GTK_DIALOG(popup),
				      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL);

		gtk_window_set_default_size(GTK_WINDOW(popup), 200, 400);
		snprintf(tmp_char, sizeof(tmp_char),
			 "Editing partition %s think before you type",
			 partid);
		label = gtk_label_new(tmp_char);
		edit_type = EDIT_EDIT;
		entry = _admin_full_edit_part(part_msg, model, iter);
	} else if(!strncasecmp("Update", type, 6)) {
		char *old_features = NULL;
		if(got_features_edit_signal)
			old_features = got_features_edit_signal;
		else
			gtk_tree_model_get(model, iter, SORTID_FEATURES,
					   &old_features, -1);
		update_features_node(GTK_DIALOG(popup),
				     nodelist, old_features);
		if(got_features_edit_signal) {
			got_features_edit_signal = NULL;
			xfree(old_features);
		} else
			g_free(old_features);
		goto end_it;
	} else {
		/* something that has to deal with a node state change */
		update_state_node(GTK_DIALOG(popup), nodelist, type);
		goto end_it;
	}

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(popup)->vbox),
			   label, FALSE, FALSE, 0);
	if(entry)
		gtk_box_pack_start(GTK_BOX(GTK_DIALOG(popup)->vbox),
				   entry, TRUE, TRUE, 0);
	gtk_widget_show_all(popup);
	response = gtk_dialog_run (GTK_DIALOG(popup));

	if (response == GTK_RESPONSE_OK) {
		if(global_edit_error)
			temp = global_edit_error_msg;
		else if(!global_send_update_msg) {
			temp = g_strdup_printf("No change detected.");
		} else if(slurm_update_partition(part_msg) == SLURM_SUCCESS) {
			temp = g_strdup_printf(
				"Partition %s updated successfully",
				partid);
		} else {
			temp = g_strdup_printf(
				"Problem updating partition %s.",
				partid);
		}
		display_edit_note(temp);
		g_free(temp);
	}
end_it:

	g_free(partid);
	g_free(nodelist);
	global_entry_changed = 0;
	slurm_free_update_part_msg(part_msg);
	gtk_widget_destroy(popup);
	if(got_edit_signal) {
		type = got_edit_signal;
		got_edit_signal = NULL;
		admin_part(model, iter, type);
		xfree(type);
	}
	if(got_features_edit_signal) {
		type = "Update Features";
		admin_part(model, iter, type);
	}
	return;
}

