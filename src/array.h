#ifndef ARRAY_H
#define ARRAY_H

#include <stdlib.h>
#include "types.h"

struct _Array
{
	void          **pdata;
	size_t          len;
	size_t          alloc;
	DestroyNotify   element_free_func;
};

typedef struct _Array Array;

Array  * array_new  (DestroyNotify element_free_func);
void  ** array_free (Array *array, int free_segment);
void     array_add  (Array *array, void *ptr);
void     array_sort (Array *array, CompareFunc compare_func);
void     array_uniq (Array *array, CompareFunc compare_func);

int     array_find                (Array *array, const void *needle, int *index_);
int     array_find_with_equal_fun (Array *array, const void *needle,
		EqualFun equal_fun, int *index_);

int    array_remove       (Array *array, void *data);
void * array_remove_index (Array *array, int index_);

#define array_len(array)        ((array)->len)
#define array_get(array, index) ((array)->pdata[index])

#endif /* array.h */
