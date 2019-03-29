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
