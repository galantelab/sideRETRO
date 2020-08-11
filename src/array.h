/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2020 Thiago L. A. Miller <tmiller@mochsl.org.br
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
#define array_data(array)       ((array)->pdata)
#define array_get(array, index) ((array)->pdata[index])

#endif /* array.h */
