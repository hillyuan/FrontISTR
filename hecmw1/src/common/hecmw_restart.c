/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_config.h"
#include "hecmw_restart.h"

struct restart_list {
	void *data;
	size_t size;
	struct restart_list *next;
};

struct restart_ent {
	int n;
	struct restart_list *rst;
};


static struct restart_list *restart_list;
static FILE *restart_fp;

struct fortran_remainder {
	void *ptr;
	struct fortran_remainder *next;
};

static struct fortran_remainder *remainder;	/* for Fortran */


static void
clear(void) 
{
	struct restart_list *p,*q;

	for(p=restart_list; p; p=q) {
		q = p->next;
		HECMW_free(p);
	}
	restart_list = NULL;
}




int
HECMW_restart_open_by_name(char *name_ID)
{
	char filename[HECMW_FILENAME_LEN+1];

	if(name_ID) {
		if(HECMW_ctrl_get_restart_file(name_ID, filename, sizeof(filename)) == NULL) {
			return -1;
		}
	} else {
		/* io is bitmap */
		int io = HECMW_CTRL_FILE_IO_IN | HECMW_CTRL_FILE_IO_INOUT;
		if(HECMW_ctrl_get_restart_file_by_io(io, filename, sizeof(filename)) == NULL) {
			return -1;
		}
	}

	if((restart_fp = fopen(filename, "rb")) == NULL) {
		HECMW_set_error(HECMW_UTIL_E0101, "File: %s, %s", filename, HECMW_strmsg(errno));
		return -1;
	}
	return 0;
}


int
HECMW_restart_open(void)
{
	return HECMW_restart_open_by_name(NULL);
}


int
HECMW_restart_close(void)
{
	if(restart_fp == NULL) return 0;
	if(fclose(restart_fp)) {
		HECMW_set_error(HECMW_UTIL_E0102, HECMW_strmsg(errno));
		return -1;
	}
	restart_fp = NULL;
	return 0;
}


static int
restart_read(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	int rc,msgno;

	rc = fread(ptr, size, nmemb, stream);
	if(rc != nmemb) {
		if(feof(stream)) {
			msgno = HECMW_UTIL_E0104;
		} else if(ferror(stream)) {
			msgno = HECMW_UTIL_E0105;
		}
		return msgno;
	}
	return 0;
}


void *
HECMW_restart_read(void *addr)
{
	int rc;
	size_t size;
	void *data;
	
	if(restart_fp == NULL) {
		HECMW_set_error(HECMW_UTIL_E0103, "");
		return NULL;
	}

	/* read size */
	rc = restart_read(&size, sizeof(size), 1, restart_fp);
	if(rc) {
		HECMW_set_error(rc, "");
		return NULL;
	}
	HECMW_assert(size > 0);

	/* read data */
	if(addr == NULL) {
		data = HECMW_malloc(size);
		if(data == NULL) {
			HECMW_set_error(errno, "");
			return NULL;
		}
	} else {
		data = addr;
	}
	rc = restart_read(data, size, 1, restart_fp);
	if(rc) {
		HECMW_set_error(rc, "");
		return NULL;
	}

	return data;
}


int
HECMW_restart_add(void *data, size_t size)
{
	struct restart_list *p,*q,*rst;

	rst = HECMW_malloc(sizeof(*rst));
	if(rst == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	
	rst->data = data; 
	rst->size = size;
	rst->next = NULL;

	q = NULL;	
	for(p=restart_list; p; p=(q=p)->next) ;

	if(q == NULL) {
		restart_list = rst;
	} else {
		q->next = rst;
	}

	return 0;
}


int
HECMW_restart_add_int(int *data, int n_data)
{
	return HECMW_restart_add(data, sizeof(int) * n_data);
}


int
HECMW_restart_add_double(double *data, int n_data)
{
	return HECMW_restart_add(data, sizeof(double) * n_data);
}


static int
restart_write(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	int rc,msgno;

	rc = fwrite(ptr, size, nmemb, stream);
	if(rc != nmemb) {
		if(feof(stream)) {
			msgno = HECMW_UTIL_E0104;
		} else if(ferror(stream)) {
			msgno = HECMW_UTIL_E0105;
		}
		return msgno;
	}
	return 0;
}


int
HECMW_restart_write_by_name(char *name_ID)
{
	int rc,n;
	FILE *fp;
	struct restart_list *p;
	char filename[HECMW_FILENAME_LEN+1];

	n = 0;
	for(p=restart_list; p; p=p->next) {
		n++;
	}
	if(n == 0) return 0;

	if(name_ID) {
		if(HECMW_ctrl_get_restart_file(name_ID, filename, sizeof(filename)) == NULL) {
			return -1;
		}
	} else {
		/* io is bitmap */
		int io = HECMW_CTRL_FILE_IO_OUT | HECMW_CTRL_FILE_IO_INOUT;
		if(HECMW_ctrl_get_restart_file_by_io(io, filename, sizeof(filename)) == NULL) {
			return -1;
		}
	}

	if((fp = fopen(filename , "w")) == NULL) {
		HECMW_set_error(HECMW_UTIL_E0101, "File: %s, %s", filename, HECMW_strmsg(errno));
		return -1;
	}

	for(p=restart_list; p; p=p->next) {
		/* write size */
		rc = restart_write(&p->size, sizeof(p->size), 1, fp);
		if(rc) {
			HECMW_set_error(rc, "");
			return -1;
		}
		/* write data */
		rc = restart_write(p->data, p->size, 1, fp);
		if(rc) {
			HECMW_set_error(rc, "");
			return -1;
		}
	}

	if(fclose(fp)) {
		HECMW_set_error(errno, "");
		return -1;
	}

	clear();

	return 0;
}


int
HECMW_restart_write(void)
{
	return HECMW_restart_write_by_name(NULL);
}



void
hecmw_restart_open_by_name_if(char *name_ID, int *err, int len)
{
	char *name = NULL;
	char cname[HECMW_NAME_LEN+1];

	*err = 1;

	if(name_ID) {
		if(HECMW_strcpy_f2c_r(name_ID, len, cname, sizeof(cname)) == NULL) {
			return;
		}
		name = cname;
	}

	if(HECMW_restart_open_by_name(name)) {
		return;
	}

	*err = 0;
}



void
hecmw_restart_open_by_name_if_(char *name_ID, int *err, int len)
{
	hecmw_restart_open_by_name_if(name_ID, err, len);
}



void
hecmw_restart_open_by_name_if__(char *name_ID, int *err, int len)
{
	hecmw_restart_open_by_name_if(name_ID, err, len);
}



void
HECMW_RESTART_OPEN_BY_NAME_IF(char *name_ID, int *err, int len)
{
	hecmw_restart_open_by_name_if(name_ID, err, len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_open_if(int *err) 
{
	hecmw_restart_open_by_name_if(NULL, err, 0);
}



void
hecmw_restart_open_if_(int *err) 
{
	hecmw_restart_open_if(err);
}



void
hecmw_restart_open_if__(int *err) 
{
	hecmw_restart_open_if(err);
}



void
HECMW_RESTART_OPEN_IF(int *err) 
{
	hecmw_restart_open_if(err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_close_if(int *err)
{
	*err = HECMW_restart_close() ? 1 : 0;
}



void
hecmw_restart_close_if_(int *err)
{
	hecmw_restart_close_if(err);
}



void
hecmw_restart_close_if__(int *err)
{
	hecmw_restart_close_if(err);
}



void
HECMW_RESTART_CLOSE_IF(int *err)
{
	hecmw_restart_close_if(err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_read_int_if(int *dst, int *err)
{
	*err = 1;

	if(dst == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "");
		return;
	}
	if(HECMW_restart_read(dst) == NULL) {
		return;
	}

	*err = 0;
}



void
hecmw_restart_read_int_if_(int *dst, int *err)
{
	hecmw_restart_read_int_if(dst, err);
}



void
hecmw_restart_read_int_if__(int *dst, int *err)
{
	hecmw_restart_read_int_if(dst, err);
}



void
HECMW_RESTART_READ_INT_IF(int *dst, int *err)
{
	hecmw_restart_read_int_if(dst, err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_read_real_if(double *dst, int *err)
{
	*err = 1;

	if(dst == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "");
		return;
	}
	if(HECMW_restart_read(dst) == NULL) {
		return;
	}

	*err = 0;
}



void
hecmw_restart_read_real_if_(double *dst, int *err)
{
	hecmw_restart_read_real_if(dst, err);
}



void
hecmw_restart_read_real_if__(double *dst, int *err)
{
	hecmw_restart_read_real_if(dst, err);
}



void
HECMW_RESTART_READ_REAL_IF(double *dst, int *err)
{
	hecmw_restart_read_real_if(dst, err);
}


/*----------------------------------------------------------------------------*/

static void *
restart_add_alloc(void *data, int byte, int n_data)
{
	int size;
	void *cdata;
	struct fortran_remainder *remain;

	size = byte * n_data;
	cdata = HECMW_malloc(size);
	if(cdata == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	memcpy(cdata, data, size);

	remain = HECMW_malloc(sizeof(*remain));
	if(remain == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	remain->ptr = cdata;
	remain->next = remainder;
	remainder = remain;

	return cdata;
}

/*----------------------------------------------------------------------------*/


void
hecmw_restart_add_int_if(int *data, int *n_data, int *err)
{
	int *cdata;

	cdata = restart_add_alloc(data, sizeof(int), *n_data);
	if(cdata == NULL) {
		*err = 1;
		return;
	}

	*err = HECMW_restart_add_int(cdata, *n_data);
}



void
hecmw_restart_add_int_if_(int *data, int *n_data, int *err)
{
	hecmw_restart_add_int_if(data, n_data, err);
}



void
hecmw_restart_add_int_if__(int *data, int *n_data, int *err)
{
	hecmw_restart_add_int_if(data, n_data, err);
}



void
HECMW_RESTART_ADD_INT_IF(int *data, int *n_data, int *err)
{
	hecmw_restart_add_int_if(data, n_data, err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_add_real_if(double *data, int *n_data, int *err)
{
	void *cdata;

	cdata = restart_add_alloc(data, sizeof(double), *n_data);
	if(cdata == NULL) {
		*err = 1;
		return;
	}
	*err = HECMW_restart_add_double(cdata, *n_data);
}



void
hecmw_restart_add_real_if_(double *data, int *n_data, int *err)
{
	hecmw_restart_add_real_if(data, n_data, err);
}



void
hecmw_restart_add_real_if__(double *data, int *n_data, int *err)
{
	hecmw_restart_add_real_if(data, n_data, err);
}



void
HECMW_RESTART_ADD_REAL_IF(double *data, int *n_data, int *err)
{
	hecmw_restart_add_real_if(data, n_data, err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_write_by_name_if(char *name_ID, int *err, int len)
{
	char *name = NULL;
	char cname[HECMW_NAME_LEN+1];
	struct fortran_remainder *p,*q;

	*err = 1;

	if(name_ID) {
		if(HECMW_strcpy_f2c_r(name_ID, len, cname, sizeof(cname)) == NULL) {
			return;
		}
		name = cname;
	}

	if(HECMW_restart_write_by_name(name)) {
		return;
	}

	for(p=remainder; p; p=q) {
		q = p->next;
		HECMW_free(p->ptr);
		HECMW_free(p);
	}
	remainder = NULL;

	*err = 0;
}



void
hecmw_restart_write_by_name_if_(char *name_ID, int *err, int len)
{
	hecmw_restart_write_by_name_if(name_ID, err, len);
}



void
hecmw_restart_write_by_name_if__(char *name_ID, int *err, int len)
{
	hecmw_restart_write_by_name_if(name_ID, err, len);
}



void
HECMW_RESTART_WRITE_BY_NAME_IF(char *name_ID, int *err, int len)
{
	hecmw_restart_write_by_name_if(name_ID, err, len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_restart_write_if(int *err)
{
	hecmw_restart_write_by_name_if(NULL, err, 0);
}



void
hecmw_restart_write_if_(int *err)
{
	hecmw_restart_write_if(err);
}



void
hecmw_restart_write_if__(int *err)
{
	hecmw_restart_write_if(err);
}



void
HECMW_RESTART_WRITE_IF(int *err)
{
	hecmw_restart_write_if(err);
}
