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

#define array_get(array, index) ((array)->pdata)[index]

Array  * array_new  (DestroyNotify element_free_func);
void  ** array_free (Array *array, int free_segment);
void     array_add  (Array *array, void *ptr);
void     array_sort (Array *array, CompareFunc compare_func);
void     array_uniq (Array *array, CompareFunc compare_func);

#endif /* array.h */
