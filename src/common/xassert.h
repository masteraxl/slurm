/*****************************************************************************\
 *  xassert.h: assert type macro with configurable handling
 *             If NDEBUG is defined, do nothing.
 *             If not, and expression is zero, log an error message and abort.
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2004 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Morris Jette <jette@llnl.gov>.
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

#ifndef _XASSERT_H
#define _XASSERT_H	1

#if HAVE_CONFIG_H
#  include "config.h"
#endif

#include "macros.h"

#ifdef NDEBUG

#  define xassert(expr)	((void) (0))

#else /* !NDEBUG */

#  define xassert(__ex)  _STMT_START { \
     (__ex) ? ((void)0) : \
     __xassert_failed(__STRING(__ex), __FILE__,  __LINE__, __CURRENT_FUNC__);\
     } _STMT_END 

/*  This prints the assertion failed message to the slurm log facility
 *  (see log.h) and aborts the calling program
 *  (messages go to stderr if log is not initialized)
 */
extern void __xassert_failed(char *, const char *, int, const char *) 
	    __NORETURN_ATTR;

#endif /* NDEBUG. */

#endif /* !__XASSERT_H */
