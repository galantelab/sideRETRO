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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "wrapper.h"
#include "array.h"

Array *
array_new (DestroyNotify element_free_func)
{
	Array *array = xcalloc (1, sizeof (Array));
	array->element_free_func = element_free_func;
	return array;
}

void **
array_free (Array *array, int free_segment)
{
	if (array == NULL)
		return NULL;

	void **segment = NULL;

	if (free_segment)
		{
			void **data = array->pdata;
			array->pdata = NULL;

			if (array->element_free_func != NULL)
				{
					for (int i = 0; i < array->len; i++)
						array->element_free_func (data[i]);
				}

			xfree (data);
		}
	else
		segment = array->pdata;

	xfree (array);
	return segment;
}

static inline size_t
nearest_pow (size_t num)
{
	size_t n = 1;

	while (n < num && n > 0)
		n <<= 1;

	return n ? n : num;
}

static void
array_maybe_expand (Array *array, int len)
{
	if ((array->len + len) > array->alloc)
		{
			size_t old_alloc = array->alloc;
			array->alloc = nearest_pow (array->len + len);

			array->pdata = xrealloc (array->pdata,
					sizeof (void *) * array->alloc);

			for ( ; old_alloc < array->alloc; old_alloc++)
				array->pdata [old_alloc] = NULL;
		}
}

void
array_add (Array *array, void *ptr)
{
	assert (array != NULL);
	assert (array->len == 0 || (array->len != 0 && array->pdata != NULL));

	array_maybe_expand (array, 1);
	array->pdata [array->len++] = ptr;
}

int
array_find_with_equal_fun (Array *array, const void *needle,
		EqualFun equal_fun, int *index_)
{
	assert (array != NULL && equal_fun != NULL && index_ != NULL);

	int i = 0;
	for (; i < array->len; i++)
		{
			if (equal_fun (needle, array->pdata[i]))
				{
					*index_ = i;
					return 1;
				}
		}

	return 0;
}

static int
pointer_equal (const void *a, const void *b)
{
	return a == b;
}

int
array_find (Array *array, const void *needle, int *index_)
{
	return array_find_with_equal_fun (array, needle,
			pointer_equal, index_);
}

void *
array_remove_index (Array *array, int index_)
{
	assert (array != NULL && index_ >= 0);

	if (index_ >= array->len)
		return NULL;

	int i = 0;
	void *data = array->pdata[index_];

	if (array->element_free_func != NULL)
		array->element_free_func (data);

	for (i = index_ + 1; i < array->len; i++)
		array->pdata[i - 1] = array->pdata[i];

	array->pdata[i - 1] = NULL;
	array->len--;

	return data;
}

int
array_remove (Array *array, void *data)
{
	assert (array != NULL);

	int index_ = 0;

	if (array_find (array, data, &index_))
		{
			array_remove_index (array, index_);
			return 1;
		}

	return 0;
}

void
array_sort (Array *array, CompareFunc compare_func)
{
	assert (array != NULL);
	assert (compare_func != NULL);

	qsort (array->pdata, array->len, sizeof (void *), compare_func);
}

void
array_uniq (Array *array, CompareFunc compare_func)
{
	assert (array != NULL);
	assert (compare_func != NULL);

	void **data;
	size_t beacon, len;

	data = array->pdata;
	len = array->len;
	beacon = 0;

	array_sort (array, compare_func);

	for (size_t i = 1; i < array->len; i++)
		{
			if (!compare_func (&data[beacon], &data[i]))
				{
					if (array->element_free_func != NULL)
						{
							array->element_free_func (data[i]);
							data[i] = NULL;
						}
					len--;
					continue;
				}

			data[++beacon] = data[i];
		}

	array->len = len;
}
