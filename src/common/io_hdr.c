/*****************************************************************************\
 * src/common/io_hdr.c - IO connection header functions
 * $Id$
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Mark A. Grondona <mgrondona@llnl.gov>.
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

#include "src/common/io_hdr.h"
#include "src/common/slurm_protocol_defs.h"

#define IO_PROTOCOL_VERSION 0xb001

/*
static void
_print_data(char *data, int datalen)
{
	char buf[1024];
	size_t len = 0;
	int i;

	for (i = 0; i < datalen; i += sizeof(char))
		len += sprintf(buf+len, "%02x", data[i]);

	info("data: %s", buf);
}
*/

void
io_hdr_pack(io_hdr_t *hdr, Buf buffer)
{
	pack16(hdr->type, buffer);
	pack16(hdr->gtaskid, buffer);
	pack16(hdr->ltaskid, buffer);
	pack32(hdr->length,  buffer);
}

int
io_hdr_unpack(io_hdr_t *hdr, Buf buffer)
{
	uint16_t val;

	safe_unpack16(&hdr->type, buffer);
	safe_unpack16(&hdr->gtaskid, buffer);
	safe_unpack16(&hdr->ltaskid, buffer);
	safe_unpack32(&hdr->length, buffer);

	return SLURM_SUCCESS;

    unpack_error:
	return SLURM_ERROR;
}

int 
io_hdr_packed_size()
{
	return sizeof(uint32_t) + 3*sizeof(uint16_t);
}

int 
io_init_msg_validate(struct slurm_io_init_msg *msg, const char *sig)
{
	debug2("Entering io_init_msg_validate");

	debug3("  msg->version = %x", msg->version);
	debug3("  msg->nodeid = %u", msg->nodeid);

	if (msg->version != IO_PROTOCOL_VERSION) {
		error("Invalid IO init header version");
		return SLURM_ERROR;
	}

	if (memcmp((void *)sig, (void *)msg->cred_signature,
		   SLURM_CRED_SIGLEN)) {
		error("Invalid IO init header signature");
		return SLURM_ERROR;
	}

	debug2("Leaving  io_init_msg_validate");
	return SLURM_SUCCESS;
}


static int
io_init_msg_packed_size(void)
{
	int len;

	len = sizeof(uint16_t)        /* version */
		+ sizeof(uint32_t)    /* nodeid */
		+ (SLURM_CRED_SIGLEN + sizeof(uint16_t));  /* signature */

	return len;
}

static void
io_init_msg_pack(struct slurm_io_init_msg *hdr, Buf buffer)
{
	pack16(hdr->version, buffer);
       	pack32(hdr->nodeid, buffer);
	packmem((char *) hdr->cred_signature,
		(uint16_t) SLURM_CRED_SIGLEN, buffer);
}


static int
io_init_msg_unpack(struct slurm_io_init_msg *hdr, Buf buffer)
{
	uint16_t val;

	safe_unpack16(&hdr->version, buffer);
	safe_unpack32(&hdr->nodeid, buffer);
	safe_unpackmem((char *) hdr->cred_signature, &val, buffer);
	if (val != SLURM_CRED_SIGLEN)
		goto unpack_error;

	return SLURM_SUCCESS;

    unpack_error:
	error("unpack error in io_init_msg_unpack");
	return SLURM_ERROR;
}


int
io_init_msg_write_to_fd(int fd, struct slurm_io_init_msg *msg)
{
	Buf buf;
	void *ptr;
	int n;
	
	xassert(msg);

	debug2("Entering io_init_msg_write_to_fd");
	msg->version = IO_PROTOCOL_VERSION;
	buf = init_buf(io_init_msg_packed_size());
	debug2("  msg->nodeid = %d", msg->nodeid);
	io_init_msg_pack(msg, buf);
	
	ptr = get_buf_data(buf);
again:
	if ((n = write(fd, ptr, io_init_msg_packed_size())) < 0) {
		if (errno == EINTR)
			goto again;
		free_buf(buf);
		return SLURM_ERROR;
	}
	if (n != io_init_msg_packed_size()) {
		error("io init msg write too small");
		free_buf(buf);
		return SLURM_ERROR;
	}

	free_buf(buf);
	debug2("Leaving  io_init_msg_write_to_fd");
	return SLURM_SUCCESS;
}


int
io_init_msg_read_from_fd(int fd, struct slurm_io_init_msg *msg)
{
	Buf buf;
	void *ptr;
	int n;
	
	xassert(msg);

	debug2("Entering io_init_msg_read_from_fd");
	buf = init_buf(io_init_msg_packed_size());
	ptr = get_buf_data(buf);
again:
	if ((n = read(fd, ptr, io_init_msg_packed_size())) < 0) {
		if (errno == EINTR)
			goto again;
		free_buf(buf);
		return SLURM_ERROR;
	}
	if (n != io_init_msg_packed_size()) {
		error("io init msg read too small");
		free_buf(buf);
		return SLURM_ERROR;
	}
	debug3("  read %d bytes", n);
	io_init_msg_unpack(msg, buf);

	free_buf(buf);

	debug2("Leaving  io_init_msg_read_from_fd");
	return SLURM_SUCCESS;
}
