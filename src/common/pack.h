/****************************************************************************\
 *  pack.h - definitions for lowest level un/pack functions. all functions 
 *	utilize a Buf structure. Call init_buf, un/pack, and free_buf
 *****************************************************************************
 *  Copyright (C) 2002 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Kevin Tew <tew1@llnl.gov>, Moe Jette <jette1@llnl.gov>, et. al.
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
\****************************************************************************/

#ifndef _PACK_INCLUDED
#define _PACK_INCLUDED

#if HAVE_CONFIG_H
#  include "config.h"
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#else  /* !HAVE_INTTYPES_H */
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#endif  /* HAVE_INTTYPES_H */
#else   /* !HAVE_CONFIG_H */
#include <stdint.h>
#endif  /* HAVE_CONFIG_H */

#include <assert.h>
#include <time.h>
#include <string.h>

#define BUF_MAGIC 0x42554545
#define BUF_SIZE 4096

struct slurm_buf {
	uint32_t magic;
	char *head;
	uint32_t size;
	uint32_t processed;
};

typedef struct slurm_buf * Buf;

#define get_buf_data(__buf)		(__buf->head)
#define get_buf_offset(__buf)		(__buf->processed)
#define set_buf_offset(__buf,__val)	(__buf->processed = __val)
#define remaining_buf(__buf)		(__buf->size - __buf->processed)
#define size_buf(__buf)			(__buf->size)

Buf	create_buf (char *data, int size);
void	free_buf(Buf my_buf);
Buf	init_buf(int size);
void	*xfer_buf_data(Buf my_buf);

void	pack_time(time_t val, Buf buffer);
int	unpack_time(time_t *valp, Buf buffer);

void 	pack32(uint32_t val, Buf buffer);
int	unpack32(uint32_t *valp, Buf buffer);

void	pack16(uint16_t val, Buf buffer);
int	unpack16(uint16_t *valp, Buf buffer);

void	pack8(uint8_t val, Buf buffer);
int	unpack8(uint8_t *valp, Buf buffer);

void	pack32_array(uint32_t *valp, uint32_t size_val, Buf buffer);
int	unpack32_array( uint32_t **valp, uint32_t* size_val, Buf buffer);

void	packmem(char *valp, uint16_t size_val, Buf buffer);
int	unpackmem(char *valp, uint16_t *size_valp, Buf buffer);
int	unpackmem_ptr(char **valp, uint16_t *size_valp, Buf buffer);
int	unpackmem_xmalloc(char **valp, uint16_t *size_valp, Buf buffer);
int	unpackmem_malloc(char **valp, uint16_t *size_valp, Buf buffer);

void	packstr_array(char **valp, uint16_t size_val, Buf buffer);
int	unpackstr_array(char ***valp, uint16_t* size_val, Buf buffer);

void	packmem_array(char *valp, uint32_t size_val, Buf buffer);
int	unpackmem_array(char *valp, uint32_t size_valp, Buf buffer);

#define safe_pack_time(val,buf) do {			\
	assert(sizeof(val) == sizeof(time_t)); 		\
	assert(buf->magic == BUF_MAGIC);		\
	pack_time(val,buf);				\
} while (0)

#define safe_unpack_time(valp,buf) do {			\
	assert((valp) != NULL); 			\
	assert(sizeof(*valp) == sizeof(time_t));	\
	assert(buf->magic == BUF_MAGIC);		\
        if (unpack_time(valp,buf))			\
		goto unpack_error;			\
} while (0)

#define safe_pack32(val,buf) do {			\
	assert(sizeof(val) == sizeof(uint32_t)); 	\
	assert(buf->magic == BUF_MAGIC);		\
	pack32(val,buf);				\
} while (0)

#define safe_unpack32(valp,buf) do {			\
	assert((valp) != NULL); 			\
	assert(sizeof(*valp) == sizeof(uint32_t));      \
	assert(buf->magic == BUF_MAGIC);		\
        if (unpack32(valp,buf))				\
		goto unpack_error;			\
} while (0)

#define safe_pack16(val,buf) do {			\
	assert(sizeof(val) == sizeof(uint16_t)); 	\
	assert(buf->magic == BUF_MAGIC);		\
	pack16(val,buf);				\
} while (0)

#define safe_unpack16(valp,buf) do {			\
	assert((valp) != NULL); 			\
	assert(sizeof(*valp) == sizeof(uint16_t)); 	\
	assert(buf->magic == BUF_MAGIC);		\
        if (unpack16(valp,buf))				\
		goto unpack_error;			\
} while (0)

#define safe_pack8(val,buf) do {			\
	assert(sizeof(val) == sizeof(uint8_t)); 	\
	assert(buf->magic == BUF_MAGIC);		\
	pack8(val,buf);					\
} while (0)

#define safe_unpack8(valp,buf) do {			\
	assert((valp) != NULL); 			\
	assert(sizeof(*valp) == sizeof(uint8_t)); 	\
	assert(buf->magic == BUF_MAGIC);		\
        if (unpack8(valp,buf))				\
		goto unpack_error;			\
} while (0)

#define safe_pack32_array(array,size_val,buf) do {	\
	assert(size_val == 0 || array != NULL);		\
	assert(buf->magic == BUF_MAGIC);		\
	pack32_array(array,size_val,buf);		\
} while (0)				

#define safe_unpack32_array(valp,size_valp,buf) do {	\
	assert(valp != NULL);				\
	assert(sizeof(*size_valp) == sizeof(uint32_t)); \
	assert(buf->magic == BUF_MAGIC);		\
	if (unpack32_array(valp,size_valp,buf))		\
		goto unpack_error;			\
} while (0)

#define safe_packmem(valp,size_val,buf) do {		\
	assert(sizeof(size_val) == sizeof(uint16_t)); 	\
	assert(size_val == 0 || valp != NULL);		\
	assert(buf->magic == BUF_MAGIC);		\
	packmem(valp,size_val,buf);			\
} while (0)

#define safe_unpackmem(valp,size_valp,buf) do {		\
	assert(valp != NULL);		                \
	assert(sizeof(*size_valp) == sizeof(uint16_t)); \
	assert(buf->magic == BUF_MAGIC);		\
	if (unpackmem(valp,size_valp,buf))		\
		goto unpack_error;			\
} while (0)

#define safe_unpackmem_ptr(valp,size_valp,buf) do {	\
	assert(valp != NULL);				\
	assert(sizeof(*size_valp) == sizeof(uint16_t)); \
	assert(buf->magic == BUF_MAGIC);		\
	if (unpackmem_ptr(valp,size_valp,buf))		\
		goto unpack_error;			\
} while (0)

#define safe_unpackmem_xmalloc(valp,size_valp,buf) do {	\
	assert(valp != NULL);				\
	assert(sizeof(*size_valp) == sizeof(uint16_t)); \
	assert(buf->magic == BUF_MAGIC);		\
	if (unpackmem_xmalloc(valp,size_valp,buf))	\
		goto unpack_error;			\
} while (0)

#define safe_unpackmem_malloc(valp,size_valp,buf) do {	\
	assert(valp != NULL);				\
	assert(sizeof(*size_valp) == sizeof(uint16_t)); \
	assert(buf->magic == BUF_MAGIC);		\
	if (unpackmem_malloc(valp,size_valp,buf))	\
		goto unpack_error;			\
} while (0)

#define safe_pack_bit_fmt(bitmap,max_len,buf) do {	\
	assert(buf->magic == BUF_MAGIC);		\
	assert(max_len < 0xffff);			\
	if (bitmap) {					\
		char _tmp_str[max_len];			\
		uint32_t _size;				\
		bit_fmt(_tmp_str,max_len,bitmap);	\
		_size = strlen(_tmp_str)+1;		\
		packmem(_tmp_str,(uint16_t)_size,buf);	\
	} else						\
		packmem(NULL,(uint16_t)0,buf);		\
} while (0)				

#define safe_packstr(str,max_len,buf) do {		\
	uint32_t _size;					\
	assert(buf->magic == BUF_MAGIC);		\
	assert(max_len <= 0xffff);			\
	_size = (str ? strlen(str)+1 : 0);		\
	assert(_size == 0 || str != NULL);		\
	if (_size <= max_len)				\
		packmem(str,(uint16_t)_size,buf);	\
	else {						\
		char tmp_str[max_len];			\
		strncpy(tmp_str, str, max_len-1);	\
		tmp_str[max_len - 1] = (char) NULL;	\
		packmem(tmp_str,(uint16_t)max_len,buf);	\
	}						\
} while (0)				

#define packstr(str,buf) do {				\
	uint32_t _size;					\
	_size = (uint32_t)(str ? strlen(str)+1 : 0);	\
        assert(_size == 0 || str != NULL);		\
	assert(_size <= 0xffff);			\
	assert(buf->magic == BUF_MAGIC);		\
	packmem(str,(uint16_t)_size,buf);		\
} while (0)				

#define pack_bit_fmt(bitmap,buf) do {	\
	assert(buf->magic == BUF_MAGIC);		\
	if (bitmap) {					\
		char _tmp_str[0xfffe];			\
		uint32_t _size;				\
		bit_fmt(_tmp_str,0xfffe,bitmap);	\
		_size = strlen(_tmp_str)+1;		\
		packmem(_tmp_str,(uint16_t)_size,buf);	\
	} else						\
		packmem(NULL,(uint16_t)0,buf);		\
} while (0)				

#define unpackstr_ptr		                        \
        unpackmem_ptr

#define unpackstr_malloc	                        \
        unpackmem_malloc

#define unpackstr_xmalloc	                        \
        unpackmem_xmalloc

#define safe_unpackstr_malloc	                        \
        safe_unpackmem_malloc

#define safe_unpackstr_xmalloc	                        \
        safe_unpackmem_xmalloc

#define safe_packstr_array(array,size_val,buf) do {	\
	assert(size_val == 0 || array != NULL);		\
	assert(buf->magic == BUF_MAGIC);		\
	packstr_array(array,size_val,buf);		\
} while (0)				

#define safe_unpackstr_array(valp,size_valp,buf) do {	\
	assert(valp != NULL);				\
	assert(size_valp != NULL);			\
	assert(sizeof(*size_valp) == sizeof(uint16_t)); \
	assert(buf->magic == BUF_MAGIC);		\
	if (unpackstr_array(valp,size_valp,buf))	\
		goto unpack_error;			\
} while (0)

#define safe_packmem_array(valp,size,buf) do {		\
	assert(size == 0 || valp != NULL);		\
	assert(sizeof(size) == sizeof(uint32_t));	\
	assert(buf->magic == BUF_MAGIC);		\
	packmem_array(valp,size,buf);			\
} while (0)

#define safe_unpackmem_array(valp,size,buf) do {	\
	assert(valp != NULL);				\
	assert(sizeof(size) == sizeof(uint32_t)); 	\
	assert(buf->magic == BUF_MAGIC);		\
	if (unpackmem_array(valp,size,buf))		\
		goto unpack_error;			\
} while (0)

#endif /* _PACK_INCLUDED */
