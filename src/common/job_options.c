/*****************************************************************************\
 *  src/common/job_options.c  - Extra job options 
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark Grondona <grondona1@llnl.gov>.
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <src/common/xassert.h>
#include <src/common/xmalloc.h>
#include <src/common/xstring.h>
#include <src/common/list.h>
#include <src/common/pack.h>

#include "src/common/job_options.h"

#define JOB_OPTIONS_PACK_TAG "job_options"

struct job_options {
#ifndef NDEBUG
#define JOB_OPTIONS_MAGIC 0xa1a2a3a4
	int magic;
#endif  /* !NDEBUG */
	List options;
	ListIterator iterator;
};


static struct job_option_info * 
job_option_info_create (int type, const char *opt, const char *optarg)
{
	struct job_option_info *ji = xmalloc (sizeof (*ji));

	ji->type =   type;
	ji->option = xstrdup (opt);
	ji->optarg = optarg ? xstrdup (optarg) : NULL;

	return (ji);
}

static void job_option_info_destroy (struct job_option_info *ji)
{
	xfree (ji->option);
	xfree (ji->optarg);
	ji->type = -1;
	xfree (ji);
	return;
}

static void job_option_info_pack (struct job_option_info *ji, Buf buf)
{
	pack32  (ji->type, buf);
	packstr (ji->option, buf);
	packstr (ji->optarg, buf); /* packstr() handles NULL optarg */
	return;
}

static struct job_option_info * job_option_info_unpack (Buf buf)
{
	struct job_option_info *ji = xmalloc (sizeof (*ji));
	uint32_t type;
	uint16_t len;

	unpack32 (&type, buf);
	unpackstr_xmalloc (&ji->option, &len, buf);
	unpackstr_xmalloc (&ji->optarg, &len, buf);

	ji->type = (int) type;

	return (ji);
}


/*
 *  Create generic job options container.
 */
job_options_t job_options_create (void) 
{
	job_options_t j = xmalloc (sizeof (*j));

	xassert (j->magic = JOB_OPTIONS_MAGIC);

	j->options = list_create ((ListDelF) job_option_info_destroy);
	j->iterator = list_iterator_create (j->options);

	return (j);
}

/*
 *  Destroy container, freeing all data associated with options.
 */
void job_options_destroy (job_options_t opts) 
{
	xassert (opts != NULL);
	xassert (opts->magic == JOB_OPTIONS_MAGIC);

	if (opts->options)
		list_destroy (opts->options);

	xassert (opts->magic = ~JOB_OPTIONS_MAGIC);
	xfree (opts);
	return;
}

/*
 *  Append option of type `type' and its argument to job options
 */
int job_options_append (job_options_t opts, int type, const char *opt, 
		        const char *optarg)
{
	xassert (opts != NULL);
	xassert (opts->magic == JOB_OPTIONS_MAGIC);
	xassert (opts->options != NULL);

	list_append (opts->options, job_option_info_create (type, opt, optarg));

	return (0);
}

/*
 *  Pack all accumulated options into Buffer "buf"
 */
int job_options_pack (job_options_t opts, Buf buf)
{
	uint32_t count = 0;
	ListIterator i;
	struct job_option_info *opt;

	packstr (JOB_OPTIONS_PACK_TAG, buf);

	if (opts == NULL) {
		pack32  (0, buf);
		return (0);
	}

	xassert (opts->magic == JOB_OPTIONS_MAGIC);
	xassert (opts->options != NULL);
	xassert (opts->iterator != NULL);

	count = list_count (opts->options);
	pack32  (count, buf);

	i = list_iterator_create (opts->options);

	while ((opt = list_next (i))) 
		job_option_info_pack (opt, buf);
	list_iterator_destroy (i);

	return (count);
}

/*
 *  Unpack options from buffer "buf" into options container opts.
 */
int job_options_unpack (job_options_t opts, Buf buf)
{
	uint32_t count;
	uint16_t len;
	char *   tag;
	int      i;

	unpackstr_xmalloc (&tag, &len, buf);

	if (strncmp (tag, JOB_OPTIONS_PACK_TAG, len) != 0)
		return (-1);

	unpack32 (&count, buf);

	for (i = 0; i < count; i++) {
		list_append (opts->options, job_option_info_unpack (buf));
	}

	return (0);
}

/*
 *  Iterate over all job options
 */
const struct job_option_info * job_options_next (job_options_t opts)
{
	xassert (opts != NULL);
	xassert (opts->magic == JOB_OPTIONS_MAGIC);
	xassert (opts->options != NULL);
	xassert (opts->iterator != NULL);
	
	return (list_next (opts->iterator));
}

void job_options_iterator_reset (job_options_t opts)
{
	xassert (opts != NULL);
	xassert (opts->magic == JOB_OPTIONS_MAGIC);
	xassert (opts->options != NULL);
	xassert (opts->iterator != NULL);
	
	list_iterator_reset (opts->iterator);
}
