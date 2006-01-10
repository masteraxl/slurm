/*****************************************************************************\
 *  slurm_protocol_socket_common.h - slurm communications interface
 *	 definitions based upon sockets
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Kevin Tew <tew1@llnl.gov>, et. al.
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

#ifndef _SLURM_PROTOCOL_SOCKET_COMMON_H
#define _SLURM_PROTOCOL_SOCKET_COMMON_H

#if HAVE_CONFIG_H
#  include "config.h"
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  else
#    if HAVE_STDINT_H
#      include <stdint.h>
#    endif
#  endif  /* HAVE_INTTYPES_H */
#else   /* !HAVE_CONFIG_H */
#  include <inttypes.h>
#endif  /*  HAVE_CONFIG_H */

#include <netinet/in.h>

#define AF_SLURM AF_INET
#define SLURM_INADDR_ANY 0x00000000

/* LINUX SPECIFIC */
/* this is the slurm equivalent of the operating system file descriptor, 
 * which in linux is just an int */
typedef int32_t slurm_fd ;

/* this is the slurm equivalent of the BSD sockets sockaddr
 * also defined in slurm/slurm.h for users */
#ifndef __slurm_addr_defined
#  define  __slurm_addr_defined
   typedef struct sockaddr_in slurm_addr ;
#endif

/* this is the slurm equivalent of the BSD sockets fd_set */
typedef fd_set slurm_fd_set ;
typedef fd_set _slurm_fd_set ;
/*{
	int16_t family ;
	uint16_t port ;
	uint32_t address ;
	char pad[16 - sizeof ( int16_t ) - sizeof (uint16_t) - 
	         sizeof (uint32_t) ] ;
} ;
*/

#endif
