
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "src/srun/job.h"
#include "src/srun/fname.h"
#include "src/srun/opt.h"

#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/xassert.h"

/* 
 * Max zero-padding width allowed
 */
#define MAX_WIDTH 10


/*
 * Fill in as much of filename as possible from srun, update
 * filename type to one of the io types ALL, NONE, PER_TASK, ONE
 */
io_filename_t *
fname_create(job_t *job, char *format)
{
	unsigned long int wid     = 0;
	unsigned long int taskid  = 0;
	io_filename_t *fname = NULL;
	char *p, *q, *name;

	fname = xmalloc(sizeof(*fname));
	fname->type = IO_ALL;
	fname->name = NULL;

	/* Handle special  cases
	 */

	if ((format == NULL)
	    || (strncasecmp(format, "all", (size_t) 3) == 0)
	    || (strncmp(format, "-", (size_t) 1) == 0)       ) {
		/* fname->type = IO_ALL; */
		/* fname->name = NULL; */
		return fname;
	}

	if (strncasecmp(format, "none", (size_t) 4) == 0) {
		/* fname->type = IO_ALL; */
		fname->name = "/dev/null";
		return fname;
	}

	taskid = strtoul(format, &p, 10);
	if ((*p == '\0') && ((int) taskid < opt.nprocs)) {
		fname->type   = IO_ONE;
		fname->taskid = (uint32_t) taskid;
		fname->name   = NULL;
		return fname;
	}

	name = NULL;
	q = p = format;
	while (*p != '\0') {
		if (*p == '%') {
			if (isdigit(*(++p))) {
				xmemcat(name, q, p - 1);
				if ((wid = strtoul(p, &p, 10)) > MAX_WIDTH)
					wid = MAX_WIDTH;
				q = p - 1;
				if (*p == '\0')
					break;
			}

			switch (*p) {
			 case 't':  /* '%t' => taskid         */
			 case 'n':  /* '%n' => nodeid         */
			 case 'N':  /* '%N' => node name      */

				 fname->type = IO_PER_TASK;
				 if (wid)
					 xstrcatchar(name, '%');
				 p++;
				 break;

			 case 'J':  /* '%J' => "jobid.stepid" */
			 case 'j':  /* '%j' => jobid          */

				 xmemcat(name, q, p - 1);
				 xstrfmtcat(name, "%0*d", wid, job->jobid);

				 if ((*p == 'J') && (job->stepid != NO_VAL)) 
					 xstrfmtcat(name, ".%d", job->stepid);
				 q = ++p;
				 break;

			 case 's':  /* '%s' => stepid         */
				 xmemcat(name, q, p - 1);
				 xstrfmtcat(name, "%0*d", wid, job->stepid);
				 q = ++p;
				 break;

			 default:
				 break;
			}
			wid = 0;
		} else
			p++;
	}

	if (q != p) 
		xmemcat(name, q, p);

	fname->name = name;
	return fname;
}

void 
fname_destroy(io_filename_t *f)
{
	if (f->name)
		xfree(f->name);
	xfree(f);
}
