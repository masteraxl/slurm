/*****************************************************************************\
 *  parse_config.c - parse any slurm.conf-like configuration file
 *
 *  NOTE: when you see the prefix "s_p_", think "slurm parser".
 *
 *  $Id$
 *****************************************************************************
 *  Copyright (C) 2006 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  Written by Christopher J. Morrone <morrone2@llnl.gov>.
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

#include <string.h>
#include <sys/types.h>
#include <regex.h>
#include <stdint.h>

/* #include "src/common/slurm_protocol_defs.h" */
#include "src/common/log.h"
#include "src/common/macros.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/xassert.h"
/* #include "src/common/slurm_rlimits_info.h" */
#include "src/common/parse_config.h"

#include <slurm/slurm.h>

#define BUFFER_SIZE 4096

#define CONF_HASH_LEN 26

static regex_t keyvalue_re;
static char *keyvalue_pattern =
	"(^|[[:space:]])([[:alpha:]]+)"
	"[[:space:]]*=[[:space:]]*"
	"([[:graph:]]+)([[:space:]]|$)";
static bool keyvalue_initialized = false;

struct s_p_values {
	char *key;
	int type;
	int data_count;
	void *data;
	int (*handler)(void **, slurm_parser_enum_t,
		       const char *, const char *, const char *);
	void (*destroy)(void *);
	s_p_values_t *next;
};

/*
 * NOTE - "key" is case insensitive.
 */
static int _conf_hashtbl_index(const char *key)
{
	int i;
	int idx = 0;

	xassert(key);
	for (i = 0; i < 10; i++) {
		if (key[i] == '\0')
			break;
		idx += tolower(key[i]);
	}
	return idx % CONF_HASH_LEN;
}

static void _conf_hashtbl_insert(s_p_hashtbl_t *hashtbl,
				 s_p_values_t *value)
{
	int idx;

	xassert(value);
	idx = _conf_hashtbl_index(value->key);
	value->next = hashtbl[idx];
	hashtbl[idx] = value;
}

/*
 * NOTE - "key" is case insensitive.
 */
static s_p_values_t *_conf_hashtbl_lookup(
	const s_p_hashtbl_t *hashtbl, const char *key)
{
	int idx;
	s_p_values_t *p;

	xassert(key);
	if (hashtbl == NULL)
		return NULL;

	idx = _conf_hashtbl_index(key);
	for (p = hashtbl[idx]; p != NULL; p = p->next) {
		if (strcasecmp(p->key, key) == 0)
			return p;
	}
	return NULL;
}

s_p_hashtbl_t *s_p_hashtbl_create(
	s_p_options_t options[])
{
	s_p_options_t *op = NULL;
	s_p_values_t *value;
	s_p_hashtbl_t *hashtbl;
	int len;

	len = CONF_HASH_LEN * sizeof(s_p_values_t *);
	hashtbl = (s_p_hashtbl_t *)xmalloc(len);
	memset(hashtbl, 0, len);
					      
	for (op = options; op->key != NULL; op++) {
		value = xmalloc(sizeof(s_p_values_t));
		value->key = xstrdup(op->key);
		value->type = op->type;
		value->data_count = 0;
		value->data = NULL;
		value->next = NULL;
		value->handler = op->handler;
		value->destroy = op->destroy;
		_conf_hashtbl_insert(hashtbl, value);
	}

	return hashtbl;
}

static void _conf_file_values_free(s_p_values_t *p)
{
	int i;

	if (p->data_count > 0) {
		switch(p->type) {
		case S_P_ARRAY:
			for (i = 0; i < p->data_count; i++) {
				void **ptr_array = (void **)p->data;
				if (p->destroy != NULL) {
					p->destroy(ptr_array[i]);
				} else {
					xfree(ptr_array[i]);
				}
			}
			xfree(p->data);
			break;
		default:
			if (p->destroy != NULL) {
				p->destroy(p->data);
			} else {
				xfree(p->data);
			}
			break;
		}
	}
	xfree(p->key);
	xfree(p);
}

void s_p_hashtbl_destroy(s_p_hashtbl_t *hashtbl) {
	int i;
	s_p_values_t *p, *next;

	for (i = 0; i < CONF_HASH_LEN; i++) {
		for (p = hashtbl[i]; p != NULL; p = next) {
			next = p->next;
			_conf_file_values_free(p);
		}
	}
	xfree(hashtbl);
}

static void _keyvalue_regex_init(void)
{
	if (!keyvalue_initialized) {
		if (regcomp(&keyvalue_re, keyvalue_pattern,
			    REG_EXTENDED) != 0) {
			/* FIXME - should be fatal */
			error("keyvalue regex compilation failed\n");
		}
		keyvalue_initialized = true;
	}
}

/*
 * IN line - string to be search for a key=value pair
 * OUT key - pointer to the key string (caller must free with free())
 * OUT value - pointer to the value string (caller must free with free())
 * OUT remaining - pointer into the "line" string denoting the start
 *                 of the unsearched portion of the string
 * Return 0 when a key-value pair is found, and -1 otherwise.
 */
static int _keyvalue_regex(const char *line,
			   char **key, char **value, char **remaining)
{
        size_t nmatch = 5;
        regmatch_t pmatch[5];
	char *start;
	size_t len;
	char *match;

	memset(pmatch, 0, sizeof(regmatch_t)*nmatch);
	if (regexec(&keyvalue_re, line, nmatch, pmatch, 0)
	    == REG_NOMATCH) {
		return -1;
	}
	
	*key = (char *)(xstrndup(line + pmatch[2].rm_so,
				 pmatch[2].rm_eo - pmatch[2].rm_so));
	*value = (char *)(xstrndup(line + pmatch[3].rm_so,
				   pmatch[3].rm_eo - pmatch[3].rm_so));
	*remaining = (char *)(line + pmatch[3].rm_eo);
	return 0;
}

static int _strip_continuation(char *buf, int len)
{
	char *ptr;
	int bs = 0;

	for (ptr = buf+len-1; ptr >= buf; ptr--) {
		if (*ptr == '\\')
			bs++;
		else if (isspace(*ptr) && bs == 0)
			continue;
		else
			break;
	}
	/* Check for an odd number of contiguous backslashes at
	   the end of the line */
	if (bs % 2 == 1) {
		ptr = ptr + bs;
		*ptr = '\0';
		return (ptr - buf);
	} else {
		return len; /* no continuation */
	}
}

/*
 * Strip out trailing carriage returns and newlines
 */
static void _strip_cr_nl(char *line)
{
	int len = strlen(line);
	char *ptr;

	for (ptr = line+len-1; ptr >= line; ptr--) {
		if (*ptr=='\r' || *ptr=='\n') {
			*ptr = '\0';
		} else {
			return;
		}
	}
}

/* Strip comments from a line by terminating the string
 * where the comment begins.
 * Everything after a non-escaped "#" is a comment.
 */
static void _strip_comments(char *line)
{
	int i;
	int len = strlen(line);
	int bs_count = 0;

	for (i = 0; i < len; i++) {
		/* if # character is preceded by an even number of
		 * escape characters '\' */
		if (line[i] == '#' && (bs_count%2) == 0) {
			line[i] = '\0';
 			break;
		} else if (line[i] == '\\') {
			bs_count++;
		} else {
			bs_count = 0;
		}
	}
}

/*
 * Strips any escape characters, "\".  If you WANT a back-slash,
 * it must be escaped, "\\".
 */
static void _strip_escapes(char *line)
{
	int i, j;
	int len = strlen(line);

	for (i = 0, j = 0; i < len+1; i++, j++) {
		if (line[i] == '\\')
			i++;
		line[j] = line[i];
	}
}

/*
 * Reads the next line from the "file" into buffer "buf".
 *
 * Concatonates together lines that are continued on
 * the next line by a trailing "\".  Strips out comments,
 * replaces escaped "\#" with "#", and replaces "\\" with "\".
 */
static int _get_next_line(char *buf, int buf_size, FILE *file)
{
	char *ptr = buf;
	int leftover = buf_size;
	int read_size, new_size;
	int eof = 1;

	while (fgets(ptr, leftover, file)) {
		eof = 0;
		_strip_comments(ptr);
		read_size = strlen(ptr);
		new_size = _strip_continuation(ptr, read_size);
		if (new_size < read_size) {
			ptr += new_size;
			leftover -= new_size;
		} else { /* no continuation */
			break;
		}
	}
	/*_strip_cr_nl(buf);*/ /* not necessary */
	_strip_escapes(buf);
	
	return !eof;
}

static int _handle_string(s_p_values_t *v,
			  const char *value, const char *line)
{
	if (v->data_count != 0) {
		error("%s specified more than once", v->key);
		return -1;
	}

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&v->data, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		v->data = xstrdup(value);
	}

	v->data_count = 1;
	return 1;
}

static int _handle_long(s_p_values_t *v,
		       const char *value, const char *line)
{
	if (v->data_count != 0) {
		error("%s specified more than once", v->key);
		return -1;
	}

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&v->data, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		char *endptr;
		long num;
		errno = 0;
		num = strtol(value, &endptr, 0);
		if ((num == 0 && errno == EINVAL)
		    || (*endptr != '\0')) {
			if (strcasecmp(value, "INFINITE") == 0) {
				num = (long)-1;
			} else {
				error("\"%s\" is not a valid number", value);
				return -1;
			}
		} else if (errno == ERANGE) {
			error("\"%s\" is out of range", value);
			return -1;
		}
		v->data = xmalloc(sizeof(long));
		*(long *)v->data = num;
	}

	v->data_count = 1;
	return 1;
}

static int _handle_uint16(s_p_values_t *v,
			  const char *value, const char *line)
{
	if (v->data_count != 0) {
		error("%s specified more than once", v->key);
		return -1;
	}

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&v->data, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		char *endptr;
		long num;

		errno = 0;
		num = strtol(value, &endptr, 0);
		if ((num == 0 && errno == EINVAL)
		    || (*endptr != '\0')) {
			if (strcasecmp(value, "INFINITE") == 0) {
				num = (uint16_t)-1;
			} else {
				error("\"%s\" is not a valid number", value);
				return -1;
			}
		} else if (errno == ERANGE) {
			error("\"%s\" is out of range", value);
			return -1;
		} else if (num < 0) {
			error("\"%s\" is less than zero", value);
			return -1;
		} else if (num > 0xffff) {
			error("\"%s\" is greater than 65535", value);
			return -1;
		}
		v->data = xmalloc(sizeof(uint16_t));
		*(uint16_t *)v->data = (uint16_t)num;
	}

	v->data_count = 1;
	return 1;
}

static int _handle_uint32(s_p_values_t *v,
			  const char *value, const char *line)
{
	if (v->data_count != 0) {
		error("%s specified more than once", v->key);
		return -1;
	}

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&v->data, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		char *endptr;
		long long num;

		errno = 0;
		num = strtoll(value, &endptr, 0);
		if ((num == 0 && errno == EINVAL)
		    || (*endptr != '\0')) {
			if (strcasecmp(value, "INFINITE") == 0) {
				num = (uint32_t)-1;
			} else {
				error("\"%s\" is not a valid number", value);
				return -1;
			}
		} else if (errno == ERANGE) {
			error("\"%s\" is out of range", value);
			return -1;
		} else if (num < 0) {
			error("\"%s\" is less than zero", value);
			return -1;
		} else if (num > 0xffffffff) {
			error("\"%s\" is greater than 4294967295", value);
			return -1;
		}
		v->data = xmalloc(sizeof(uint32_t));
		*(uint32_t *)v->data = (uint32_t)num;
	}

	v->data_count = 1;
	return 1;
}

static int _handle_pointer(s_p_values_t *v,
			   const char *value, const char *line)
{
	if (v->data_count != 0) {
		error("%s specified more than once", v->key);
		return -1;
	}

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&v->data, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		v->data = xstrdup(value);
	}

	v->data_count = 1;
	return 1;
}

static int _handle_array(s_p_values_t *v,
			 const char *value, const char *line)
{
	void *new_ptr;
	void **data;

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&new_ptr, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		new_ptr = xstrdup(value);
	}
	v->data_count += 1;
	v->data = xrealloc(v->data, (v->data_count)*sizeof(void *));
	data = &((void**)v->data)[v->data_count-1];
	*data = new_ptr;

	return 1;
}

static int _handle_boolean(s_p_values_t *v,
			   const char *value, const char *line)
{
	if (v->data_count != 0) {
		error("%s specified more than once", v->key);
		return -1;
	}

	if (v->handler != NULL) {
		/* call the handler function */
		int rc;
		rc = v->handler(&v->data, v->type, v->key, value, line);
		if (rc != 1)
			return rc == 0 ? 0 : -1;
	} else {
		bool flag;

		if (!strcasecmp(value, "yes")
		    || !strcasecmp(value, "up")
		    || !strcasecmp(value, "1")) {
			flag = true;
		} else if (!strcasecmp(value, "no")
			   || !strcasecmp(value, "down")
			   || !strcasecmp(value, "0")) {
			flag = false;
		} else {
			error("\"%s\" is not a valid option for \"%s\"",
			      value, v->key);
			return -1;
		}

		v->data = xmalloc(sizeof(bool));
		*(bool *)v->data = flag;
	}

	v->data_count = 1;
	return 1;
}

static void _handle_keyvalue_match(s_p_values_t *v,
				   const char *value, const char *line)
{
	/* debug3("key = %s, value = %s, line = \"%s\"", */
	/*        v->key, value, line); */
	switch (v->type) {
	case S_P_STRING:
		_handle_string(v, value, line);
		break;
	case S_P_LONG:
		_handle_long(v, value, line);
		break;
	case S_P_UINT16:
		_handle_uint16(v, value, line);
		break;
	case S_P_UINT32:
		_handle_uint32(v, value, line);
		break;
	case S_P_POINTER:
		_handle_pointer(v, value, line);
		break;
	case S_P_ARRAY:
		_handle_array(v, value, line);
		break;
	case S_P_BOOLEAN:
		_handle_boolean(v, value, line);
		break;
	}
}

void s_p_parse_file(s_p_hashtbl_t *hashtbl, char *filename)
{
	FILE *f;
	char line[BUFFER_SIZE];
	char *key, *value, *leftover;

	_keyvalue_regex_init();

	f = fopen(filename, "r");

	while(_get_next_line(line, BUFFER_SIZE, f)) {
		/* skip empty lines */
		if (line[0] == '\0')
			continue;
		/* debug3("line = \"%s\"", line); */

		if (_keyvalue_regex(line, &key, &value, &leftover) == 0) {
			s_p_values_t *p;

			if (p = _conf_hashtbl_lookup(hashtbl, key)) {
				_handle_keyvalue_match(p, value, line);
			} else {
				fatal("UNRECOGNIZED KEY %s!", key);
			}
			xfree(key);
			xfree(value);
		}
	}

	fclose(f);
}

void s_p_parse_line(s_p_hashtbl_t *hashtbl, const char *line)
{
	char *key, *value, *leftover;
	const char *ptr = line;
	s_p_values_t *p;

	_keyvalue_regex_init();

	while (_keyvalue_regex(ptr, &key, &value, &leftover) == 0) {
		if (p = _conf_hashtbl_lookup(hashtbl, key)) {
			_handle_keyvalue_match(p, value, leftover);
			ptr = leftover;
		} else {
			fatal("UNRECOGNIZED KEY %s!", key);
		}
		xfree(key);
		xfree(value);
	}
}

/*
 * s_p_get_string - Search for a key in a s_p_hashtbl_t with value of type
 *                  string.  If the key is found and has a set value, the
 *                  value is retuned in "str".
 *
 * IN hashtbl - hash table created by s_p_hashtbl_create()
 * IN key - hash table key.
 * OUT str - pointer to a copy of the string value
 *           (caller is resonsible for freeing str with xfree())
 *
 * Returns 1 when a value was set for "key" during parsing and "str"
 *   was successfully set, otherwise returns 0;
 */
int s_p_get_string(const s_p_hashtbl_t *hashtbl, const char *key, char **str)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_STRING) {
		error("Key \"%s\" is not a string\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*str = xstrdup((char *)p->data);

	return 1;
}

/*
 * s_p_get_long - Search for a key in a s_p_hashtbl_t with value of type
 *                  long.  If the key is found and has a set value, the
 *                  value is retuned in "num".
 *
 * IN hashtbl - hash table created by s_p_hashtbl_create()
 * IN key - hash table key
 * OUT num - pointer to a long where the value is returned
 *
 * Returns 1 when a value was set for "key" during parsing and "num"
 *   was successfully set, otherwise returns 0;
 */
int s_p_get_long(const s_p_hashtbl_t *hashtbl, const char *key, long *num)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_LONG) {
		error("Key \"%s\" is not a long\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*num = *(long *)p->data;

	return 1;
}

/*
 * s_p_get_uint16 - Search for a key in a s_p_hashtbl_t with value of type
 *                  uint16.  If the key is found and has a set value, the
 *                  value is retuned in "num".
 *
 * IN hashtbl - hash table created by s_p_hashtbl_create()
 * IN key - hash table key
 * OUT num - pointer to a uint16_t where the value is returned
 *
 * Returns 1 when a value was set for "key" during parsing and "num"
 *   was successfully set, otherwise returns 0;
 */
int s_p_get_uint16(const s_p_hashtbl_t *hashtbl, const char *key,
		   uint16_t *num)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_UINT16) {
		error("Key \"%s\" is not a uint16_t\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*num = *(uint16_t *)p->data;

	return 1;
}

/*
 * s_p_get_uint32 - Search for a key in a s_p_hashtbl_t with value of type
 *                  uint32.  If the key is found and has a set value, the
 *                  value is retuned in "num".
 *
 * IN hashtbl - hash table created by s_p_hashtbl_create()
 * IN key - hash table key
 * OUT num - pointer to a uint32_t where the value is returned
 *
 * Returns 1 when a value was set for "key" during parsing and "num"
 *   was successfully set, otherwise returns 0;
 */
int s_p_get_uint32(const s_p_hashtbl_t *hashtbl, const char *key,
		   uint32_t *num)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_UINT32) {
		error("Key \"%s\" is not a uint32_t\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*num = *(uint32_t *)p->data;

	return 1;
}

/*
 * s_p_get_pointer - Search for a key in a s_p_hashtbl_t with value of type
 *                   pointer.  If the key is found and has a set value, the
 *                   value is retuned in "ptr".
 *
 * IN hashtbl - hash table created by s_p_hashtbl_create()
 * IN key - hash table key
 * OUT num - pointer to a void pointer where the value is returned
 *
 * Returns 1 when a value was set for "key" during parsing and "ptr"
 *   was successfully set, otherwise returns 0;
 */
int s_p_get_pointer(const s_p_hashtbl_t *hashtbl, const char *key, void **ptr)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_POINTER) {
		error("Key \"%s\" is not a pointer\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*ptr = p->data;

	return 1;
}

int s_p_get_array(const s_p_hashtbl_t *hashtbl, const char *key,
		  void **ptr_array[], int *count)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_ARRAY) {
		error("Key \"%s\" is not an array\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*ptr_array = (void **)p->data;
	*count = p->data_count;

	return 1;
}

/*
 * s_p_get_boolean - Search for a key in a s_p_hashtbl_t with value of type
 *                   boolean.  If the key is found and has a set value, the
 *                   value is retuned in "flag".
 *
 * IN hashtbl - hash table created by s_p_hashtbl_create()
 * IN key - hash table key
 * OUT flag - pointer to a bool where the value is returned
 *
 * Returns 1 when a value was set for "key" during parsing and "num"
 *   was successfully set, otherwise returns 0;
 */
int s_p_get_boolean(const s_p_hashtbl_t *hashtbl, const char *key, bool *flag)
{
	s_p_values_t *p;

	p = _conf_hashtbl_lookup(hashtbl, key);
	if (p == NULL) {
		error("Invalid key \"%s\"", key);
		return 0;
	}
	if (p->type != S_P_BOOLEAN) {
		error("Key \"%s\" is not a boolean\n", key);
		return 0;
	}
	if (p->data_count == 0) {
		return 0;
	}

	*flag = *(bool *)p->data;

	return 1;
}


/*
 * Given an "options" array, print the current values of all
 * options in supplied hash table "hashtbl".
 *
 * Primarily for debugging purposes.
 */
void s_p_dump_values(const s_p_hashtbl_t *hashtbl,
		     const s_p_options_t options[])
{
	const s_p_options_t *op = NULL;
	long num;
	uint16_t num16;
	uint32_t num32;
	char *str;
	void *ptr;
	void **ptr_array;
	int count;
	int i;

	for (op = options; op->key != NULL; op++) {
		switch(op->type) {
		case S_P_STRING:
			if (s_p_get_string(hashtbl, op->key, &str)) {
			        debug("%s = %s", op->key, str);
				xfree(str);
			} else {
				debug("%s", op->key);
			}
			break;
		case S_P_LONG:
			if (s_p_get_long(hashtbl, op->key, &num))
				debug("%s = %ld", op->key, num);
			else
				debug("%s", op->key);
			break;
		case S_P_UINT16:
			if (s_p_get_uint16(hashtbl, op->key, &num16))
				debug("%s = %hu", op->key, num16);
			else
				debug("%s", op->key);
			break;
		case S_P_UINT32:
			if (s_p_get_uint32(hashtbl, op->key, &num32))
				debug("%s = %u", op->key, num32);
			else
				debug("%s", op->key);
			break;
		case S_P_POINTER:
			if (s_p_get_pointer(hashtbl, op->key, &ptr))
				debug("%s = %x", op->key, ptr);
			else
				debug("%s", op->key);
			break;
		case S_P_ARRAY:
			if (s_p_get_array(hashtbl, op->key,
					  &ptr_array, &count)) {
				debug("%s, count = %d", op->key, count);
			} else {
				debug("%s", op->key);
			}
			break;
		}
	}
}
