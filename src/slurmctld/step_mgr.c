/*****************************************************************************\
 *  step_mgr.c - manage the job step information of slurm
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by moe jette <jette1@llnl.gov>, Joseph Ekstrom <ekstrom1@llnl.gov>
 *  UCRL-CODE-2002-040.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <unistd.h>

#include <slurm/slurm_errno.h>

#include "src/common/bitstring.h"
#include "src/common/switch.h"
#include "src/common/xstring.h"
#include "src/slurmctld/agent.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/slurmctld.h"

#define MAX_RETRIES 10

static void _pack_ctld_job_step_info(struct step_record *step, Buf buffer);
static bitstr_t * _pick_step_nodes (struct job_record  *job_ptr, 
				    step_specs *step_spec );
/* 
 * create_step_record - create an empty step_record for the specified job.
 * IN job_ptr - pointer to job table entry to have step record added
 * RET a pointer to the record or NULL if error
 * NOTE: allocates memory that should be xfreed with delete_step_record
 */
struct step_record * 
create_step_record (struct job_record *job_ptr) 
{
	struct step_record *step_record_point;

	xassert(job_ptr);
	step_record_point = 
		(struct step_record *) xmalloc (sizeof (struct step_record));

	last_job_update = time(NULL);
	step_record_point->job_ptr = job_ptr; 
	step_record_point->step_id = (job_ptr->next_step_id)++;
	step_record_point->start_time = time ( NULL ) ;

	if (list_append (job_ptr->step_list, step_record_point) == NULL)
		fatal ("create_step_record: unable to allocate memory");

	return step_record_point;
}


/* 
 * delete_all_step_records - delete all step record for specified job_ptr
 * IN job_ptr - pointer to job table entry to have step record added
 */
void 
delete_all_step_records (struct job_record *job_ptr) 
{
	ListIterator step_iterator;
	struct step_record *step_ptr;

	xassert(job_ptr);
	step_iterator = list_iterator_create (job_ptr->step_list);

	last_job_update = time(NULL);
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		list_remove (step_iterator);
		g_switch_free_jobinfo(step_ptr->switch_job);
		xfree(step_ptr->host);
		xfree(step_ptr->step_node_list);
		FREE_NULL_BITMAP(step_ptr->step_node_bitmap);
		xfree(step_ptr);
	}		

	list_iterator_destroy (step_iterator);
}


/* 
 * delete_step_record - delete record for job step for specified job_ptr 
 *	and step_id
 * IN job_ptr - pointer to job table entry to have step record removed
 * IN step_id - id of the desired job step
 * RET 0 on success, errno otherwise
 */
int 
delete_step_record (struct job_record *job_ptr, uint32_t step_id) 
{
	ListIterator step_iterator;
	struct step_record *step_ptr;
	int error_code;

	xassert(job_ptr);
	error_code = ENOENT;
	step_iterator = list_iterator_create (job_ptr->step_list);		

	last_job_update = time(NULL);
	while ((step_ptr = (struct step_record *) list_next (step_iterator))) {
		if (step_ptr->step_id == step_id) {
			list_remove (step_iterator);
			g_switch_free_jobinfo (step_ptr->switch_job);
			xfree(step_ptr->host);
			xfree(step_ptr->step_node_list);
			FREE_NULL_BITMAP(step_ptr->step_node_bitmap);
			xfree(step_ptr);
			error_code = 0;
			break;
		}
	}		

	list_iterator_destroy (step_iterator);
	return error_code;
}


/*
 * dump_step_desc - dump the incoming step initiate request message
 * IN step_spec - job step request specification from RPC
 */
void
dump_step_desc(step_specs *step_spec)
{
	if (step_spec == NULL) 
		return;

	debug3("StepDesc: user_id=%u job_id=%u node_count=%u, cpu_count=%u", 
		step_spec->user_id, step_spec->job_id, 
		step_spec->node_count, step_spec->cpu_count);
	debug3("   num_tasks=%u relative=%u task_dist=%u node_list=%s", 
		step_spec->num_tasks, step_spec->relative, 
		step_spec->task_dist, step_spec->node_list);
	debug3("   host=%s port=%u", 
		step_spec->host, step_spec->port);
}


/* 
 * find_step_record - return a pointer to the step record with the given 
 *	job_id and step_id
 * IN job_ptr - pointer to job table entry to have step record added
 * IN step_id - id of the desired job step
 * RET pointer to the job step's record, NULL on error
 */
struct step_record *
find_step_record(struct job_record *job_ptr, uint16_t step_id) 
{
	ListIterator step_record_iterator;
	struct step_record *step_ptr;

	if (job_ptr == NULL)
		return NULL;

	step_record_iterator = list_iterator_create (job_ptr->step_list);		
	while ((step_ptr = (struct step_record *) 
			   list_next (step_record_iterator))) {
		if (step_ptr->step_id == step_id) {
			break;
		}
	}		

	list_iterator_destroy (step_record_iterator);
	return step_ptr;
}


/* 
 * job_step_signal - signal the specified job step
 * IN job_id - id of the job to be cancelled
 * IN step_id - id of the job step to be cancelled
 * IN signal - user id of user issuing the RPC
 * IN uid - user id of user issuing the RPC
 * RET 0 on success, otherwise ESLURM error code 
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 */
int job_step_signal(uint32_t job_id, uint32_t step_id, 
		    uint16_t signal, uid_t uid)
{
	struct job_record *job_ptr;
	struct step_record *step_ptr;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {

		error("job_step_cancel: invalid job id %u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}

	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;

	if ((job_ptr->user_id != uid) && (uid != 0) && (uid != getuid())) {
		error("Security violation, JOB_CANCEL RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}

	step_ptr = find_step_record(job_ptr, (uint16_t)step_id);
	if (step_ptr == NULL) {
		info("job_step_cancel step %u.%u not found",
		     job_id, step_id);
		return ESLURM_INVALID_JOB_ID;
	}

	signal_step_tasks(step_ptr, signal);
	return SLURM_SUCCESS;
}

/*
 * signal_step_tasks - send specific signal to specific job step
 * IN step_ptr - step record pointer
 * IN signal - signal to send
 */
void signal_step_tasks(struct step_record *step_ptr, uint16_t signal)
{
	int i;
	kill_tasks_msg_t *kill_tasks_msg;
	agent_arg_t *agent_args;
	pthread_attr_t attr_agent;
	pthread_t thread_agent;
	int buf_rec_size = 0;

	xassert(step_ptr);
	agent_args = xmalloc(sizeof(agent_arg_t));
	agent_args->msg_type = REQUEST_KILL_TASKS;
	agent_args->retry = 1;
	kill_tasks_msg = xmalloc(sizeof(kill_tasks_msg_t));
	kill_tasks_msg->job_id      = step_ptr->job_ptr->job_id;
	kill_tasks_msg->job_step_id = step_ptr->step_id;
	kill_tasks_msg->signal      = signal;

	for (i = 0; i < node_record_count; i++) {
		if (bit_test(step_ptr->step_node_bitmap, i) == 0)
			continue;
		if ((agent_args->node_count + 1) > buf_rec_size) {
			buf_rec_size += 32;
			xrealloc((agent_args->slurm_addr),
				 (sizeof(struct sockaddr_in) *
				  buf_rec_size));
			xrealloc((agent_args->node_names),
				 (MAX_NAME_LEN * buf_rec_size));
		}
		agent_args->slurm_addr[agent_args->node_count] =
		    node_record_table_ptr[i].slurm_addr;
		strncpy(&agent_args->
			node_names[MAX_NAME_LEN * agent_args->node_count],
			node_record_table_ptr[i].name, MAX_NAME_LEN);
		agent_args->node_count++;
	}

	if (agent_args->node_count == 0) {
		xfree(kill_tasks_msg);
		xfree(agent_args);
		return;
	}

	agent_args->msg_args = kill_tasks_msg;
	debug2("Spawning signal agent");
	if (pthread_attr_init(&attr_agent))
		fatal("pthread_attr_init error %m");
	if (pthread_attr_setdetachstate(&attr_agent, PTHREAD_CREATE_DETACHED))
		error("pthread_attr_setdetachstate error %m");
#ifdef PTHREAD_SCOPE_SYSTEM
	if (pthread_attr_setscope(&attr_agent, PTHREAD_SCOPE_SYSTEM))
		error("pthread_attr_setscope error %m");
#endif
	if (pthread_create(&thread_agent, &attr_agent, agent, 
			(void *) agent_args)) {
		error("pthread_create error %m");
		agent_queue_request(agent_args);
	}
	return;
}


/* 
 * job_step_complete - note normal completion the specified job step
 * IN job_id - id of the job to be completed
 * IN step_id - id of the job step to be completed
 * IN uid - user id of user issuing the RPC
 * IN requeue - job should be run again if possible
 * IN job_return_code - job's return code, if set then set state to JOB_FAILED
 * RET 0 on success, otherwise ESLURM error code 
 * global: job_list - pointer global job list
 *	last_job_update - time of last job table update
 */
int job_step_complete(uint32_t job_id, uint32_t step_id, uid_t uid,
		      bool requeue, uint32_t job_return_code)
{
	struct job_record *job_ptr;
	int error_code;

	job_ptr = find_job_record(job_id);
	if (job_ptr == NULL) {
		info("job_step_complete: invalid job id %u", job_id);
		return ESLURM_INVALID_JOB_ID;
	}

	if ((job_ptr->kill_on_step_done) &&
	    (list_count(job_ptr->step_list) <= 1))
		return job_complete(job_id, uid, requeue, job_return_code);

	if (IS_JOB_FINISHED(job_ptr))
		return ESLURM_ALREADY_DONE;

	if ((job_ptr->user_id != uid) && (uid != 0) && (uid != getuid())) {
		error("Security violation, JOB_COMPLETE RPC from uid %d",
		      uid);
		return ESLURM_USER_ID_MISSING;
	}

	last_job_update = time(NULL);
	error_code = delete_step_record(job_ptr, step_id);
	if (error_code == ENOENT) {
		info("job_step_complete step %u.%u not found", job_id,
		     step_id);
		return ESLURM_ALREADY_DONE;
	}
	return SLURM_SUCCESS;
}


/* 
 * _pick_step_nodes - select nodes for a job step that satify its requirements
 *	we satify the super-set of constraints.
 * IN job_ptr - pointer to job to have new step started
 * IN step_spec - job step specification
 * global: node_record_table_ptr - pointer to global node table
 * NOTE: returns all of a job's nodes if step_spec->node_count == INFINITE
 * NOTE: returned bitmap must be freed by the caller using bit_free()
 */
static bitstr_t *
_pick_step_nodes (struct job_record  *job_ptr, step_specs *step_spec ) {

	bitstr_t *nodes_avail = NULL, *nodes_picked = NULL, *node_tmp = NULL;
	int error_code, nodes_picked_cnt = 0, cpus_picked_cnt, i;

	if (job_ptr->node_bitmap == NULL)
		return NULL;
	
	nodes_avail = bit_copy (job_ptr->node_bitmap);
	bit_and (nodes_avail, avail_node_bitmap);

	if ( step_spec->node_count == INFINITE)	/* use all nodes */
		return nodes_avail;

	if (step_spec->node_list) {
		error_code = node_name2bitmap (step_spec->node_list, false, 
						&nodes_picked);
		if (error_code) {
			info ("_pick_step_nodes: invalid node list %s", 
				step_spec->node_list);
			goto cleanup;
		}
		if (bit_super_set (nodes_picked, job_ptr->node_bitmap) == 0) {
			info ("_pick_step_nodes: requested nodes %s not part "
				"of job %u", 
				step_spec->node_list, job_ptr->job_id);
			goto cleanup;
		}
		if (bit_super_set (nodes_picked, avail_node_bitmap) == 0) {
			info ("_pick_step_nodes: some requested node %s down",
				step_spec->node_list);
			goto cleanup;
		}
	}
	else if (step_spec->relative) {
		/* Remove first (step_spec->relative) nodes from  
		 * available list */
		bitstr_t *relative_nodes = NULL;
		relative_nodes = 
			bit_pick_cnt (nodes_avail, step_spec->relative);
		if (relative_nodes == NULL) {
			info ("_pick_step_nodes: Invalid relative value (%u) for job %u",
				step_spec->relative, job_ptr->job_id);
			goto cleanup;
		}
		bit_not (relative_nodes);
		bit_and (nodes_avail, relative_nodes);
		bit_free (relative_nodes);
	}
	else
		nodes_picked = bit_alloc (bit_size (nodes_avail) );

	/* if user specifies step needs a specific processor count and  */
	/* all nodes have the same processor count, just translate this to */
	/* a node count */
	if (step_spec->cpu_count && (job_ptr->num_cpu_groups == 1)) {
		i = (step_spec->cpu_count + (job_ptr->cpus_per_node[0] - 1) ) 
				/ job_ptr->cpus_per_node[0];
		step_spec->node_count = (i > step_spec->node_count) ? 
						i : step_spec->node_count ;
		step_spec->cpu_count = 0;
	}

	if (step_spec->node_count) {
		nodes_picked_cnt = bit_set_count(nodes_picked);
		if (step_spec->node_count > nodes_picked_cnt) {
			node_tmp = bit_pick_cnt(nodes_avail, 
						(step_spec->node_count - 
						 nodes_picked_cnt));
			if (node_tmp == NULL)
				goto cleanup;
			bit_or  (nodes_picked, node_tmp);
			bit_not (node_tmp);
			bit_and (nodes_avail, node_tmp);
			bit_free (node_tmp);
			node_tmp = NULL;
			nodes_picked_cnt = step_spec->node_count;
		}
	}

	if (step_spec->cpu_count) {
		cpus_picked_cnt = count_cpus(nodes_picked);
		if (step_spec->cpu_count > cpus_picked_cnt) {
			int first_bit, last_bit;
			first_bit = bit_ffs(nodes_avail);
			last_bit  = bit_fls(nodes_avail);
 			for (i = first_bit; i <= last_bit; i++) {
				if (bit_test (nodes_avail, i) != 1)
					continue;
				bit_set (nodes_picked, i);
				cpus_picked_cnt += 
					node_record_table_ptr[i].cpus;
				if (cpus_picked_cnt >= step_spec->cpu_count)
					break;
			}
			if (step_spec->cpu_count > cpus_picked_cnt)
				goto cleanup;
		}
	}

	FREE_NULL_BITMAP(nodes_avail);
	return nodes_picked;

cleanup:
	FREE_NULL_BITMAP(nodes_avail);
	FREE_NULL_BITMAP(nodes_picked);
	return NULL;
}


/*
 * step_create - creates a step_record in step_specs->job_id, sets up the
 *	according to the step_specs.
 * IN step_specs - job step specifications
 * OUT new_step_record - pointer to the new step_record (NULL on error)
 * IN kill_job_when_step_done - if set kill the job on step completion
 * RET - 0 or error code
 * NOTE: don't free the returned step_record because that is managed through
 * 	the job.
 */
int
step_create ( step_specs *step_specs, struct step_record** new_step_record,
	      bool kill_job_when_step_done )
{
	struct step_record *step_ptr;
	struct job_record  *job_ptr;
	bitstr_t *nodeset;
	int node_count;
	time_t now = time(NULL);

	*new_step_record = NULL;
	job_ptr = find_job_record (step_specs->job_id);
	if (job_ptr == NULL)
		return ESLURM_INVALID_JOB_ID ;

	if ((step_specs->user_id != job_ptr->user_id) &&
	    (step_specs->user_id != 0))
		return ESLURM_ACCESS_DENIED ;

	if (IS_JOB_PENDING(job_ptr))
		return ESLURM_INVALID_JOB_ID ;

	if (IS_JOB_FINISHED(job_ptr) || 
	    (job_ptr->end_time <= time(NULL)))
		return ESLURM_ALREADY_DONE;

	if ((step_specs->task_dist != SLURM_DIST_CYCLIC) &&
	    (step_specs->task_dist != SLURM_DIST_BLOCK))
		return ESLURM_BAD_DIST;

	if (job_ptr->kill_on_step_done)
		/* Don't start more steps, job already being cancelled */
		return ESLURM_ALREADY_DONE;
	job_ptr->kill_on_step_done = kill_job_when_step_done;

	job_ptr->time_last_active = now;
	nodeset = _pick_step_nodes (job_ptr, step_specs);
	if (nodeset == NULL)
		return ESLURM_REQUESTED_NODE_CONFIG_UNAVAILABLE ;
	node_count = bit_set_count(nodeset);

	if (step_specs->num_tasks == NO_VAL) {
		if (step_specs->cpu_count != NO_VAL)
			step_specs->num_tasks = step_specs->cpu_count;
		else
			step_specs->num_tasks = node_count;
	}
	if ((step_specs->num_tasks < 1) ||
	    (step_specs->num_tasks > (node_count*MAX_TASKS_PER_NODE)))
		return ESLURM_BAD_TASK_COUNT;

	step_ptr = create_step_record (job_ptr);
	if (step_ptr == NULL)
		fatal ("create_step_record failed with no memory");

	/* set the step_record values */
	step_ptr->step_node_list = bitmap2node_name(nodeset);
	step_ptr->step_node_bitmap = nodeset;
	step_ptr->cyclic_alloc = 
		(uint16_t) (step_specs->task_dist == SLURM_DIST_CYCLIC);
	step_ptr->num_tasks = step_specs->num_tasks;
	step_ptr->time_last_active = now;
	step_ptr->port = step_specs->port;
	step_ptr->host = xstrdup(step_specs->host);

	if (g_switch_build_jobinfo(&step_ptr->switch_job, 
				step_ptr->step_node_list,
				step_specs->num_tasks, 
				step_ptr->cyclic_alloc) < 0) {
		error("g_switch_build_jobinfo: %m");
		delete_step_record (job_ptr, step_ptr->step_id);
		return ESLURM_INTERCONNECT_FAILURE;
	}

	*new_step_record = step_ptr;
	return SLURM_SUCCESS;
}

/* Pack the data for a specific job step record
 * IN step - pointer to a job step record
 * IN/OUT buffer - location to store data, pointers automatically advanced
 */
static void _pack_ctld_job_step_info(struct step_record *step, Buf buffer)
{
	pack_job_step_info_members(step->job_ptr->job_id,
				   step->step_id,
				   step->job_ptr->user_id,
				   step->num_tasks,
				   step->start_time,
				   step->job_ptr->partition,
				   step->step_node_list, buffer);
}

/* 
 * pack_ctld_job_step_info_response_msg - packs job step info
 * IN - job_id and step_id - zero for all
 * OUT buffer - location to store data, pointers automatically advanced 
 * RET - 0 or error code
 * NOTE: MUST free_buf buffer
 */
int pack_ctld_job_step_info_response_msg(uint32_t job_id, 
				         uint32_t step_id, Buf buffer)
{
	ListIterator job_record_iterator;
	ListIterator step_record_iterator;
	int error_code = 0;
	uint32_t steps_packed = 0, tmp_offset;
	struct step_record *step_ptr;
	struct job_record *job_ptr;
	time_t now = time(NULL);

	pack_time(now, buffer);
	pack32(steps_packed, buffer);	/* steps_packed placeholder */

	if (job_id == 0) {
		/* Return all steps for all jobs */
		job_record_iterator = list_iterator_create(job_list);
		while ((job_ptr =
			(struct job_record *)
			list_next(job_record_iterator))) {
			step_record_iterator =
			    list_iterator_create(job_ptr->step_list);
			while ((step_ptr =
				(struct step_record *)
				list_next(step_record_iterator))) {
				_pack_ctld_job_step_info(step_ptr, buffer);
				steps_packed++;
			}
			list_iterator_destroy(step_record_iterator);
		}
		list_iterator_destroy(job_record_iterator);

	} else if (step_id == 0) {
		/* Return all steps for specific job_id */
		job_ptr = find_job_record(job_id);
		if (job_ptr) {
			step_record_iterator =
			    list_iterator_create(job_ptr->step_list);
			while ((step_ptr =
				(struct step_record *)
				list_next(step_record_iterator))) {
				_pack_ctld_job_step_info(step_ptr, buffer);
				steps_packed++;
			}
			list_iterator_destroy(step_record_iterator);
		} else
			error_code = ESLURM_INVALID_JOB_ID;
	} else {
		/* Return  step with give step_id/job_id */
		job_ptr = find_job_record(job_id);
		step_ptr = find_step_record(job_ptr, step_id);
		if (step_ptr == NULL)
			error_code = ESLURM_INVALID_JOB_ID;
		else {
			_pack_ctld_job_step_info(step_ptr, buffer);
			steps_packed++;
		}
	}

	/* put the real record count in the message body header */
	tmp_offset = get_buf_offset(buffer);
	set_buf_offset(buffer, 0);
	pack_time(now, buffer);
	pack32(steps_packed, buffer);
	set_buf_offset(buffer, tmp_offset);

	return error_code;
}

/* 
 * step_on_node - determine if the specified job has any job steps allocated to 
 * 	the specified node 
 * IN job_ptr - pointer to an active job record
 * IN node_ptr - pointer to a node record
 * RET true of job has step on the node, false otherwise 
 */
bool step_on_node(struct job_record  *job_ptr, struct node_record *node_ptr)
{
	ListIterator step_record_iterator;
	struct step_record *step_record_point;
	bool found = false;
	int bit_position;

	if ((job_ptr == NULL) || (node_ptr == NULL))
		return false;

	bit_position = node_ptr - node_record_table_ptr;
	step_record_iterator = list_iterator_create (job_ptr->step_list);		
	while ((step_record_point = 
		(struct step_record *) list_next (step_record_iterator))) {
		if (bit_test(step_record_point->step_node_bitmap, 
			     bit_position)) {
			found = true;
			break;
		}
	}		

	list_iterator_destroy (step_record_iterator);
	return found;
}
