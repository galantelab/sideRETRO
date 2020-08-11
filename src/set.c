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

#include <assert.h>
#include "wrapper.h"
#include "set.h"

struct _Set
{
	List          *list;
	Hash          *index;
	HashFunc       hash_fun;
	EqualFun       match_fun;
	DestroyNotify  destroy_fun;
};

Set *
set_new_full (HashFunc hash_fun, EqualFun match_fun, DestroyNotify destroy_fun)
{
	assert (hash_fun != NULL && match_fun != NULL);

	Set *set = xcalloc (1, sizeof (Set));

	set->list = list_new (destroy_fun);
	set->index = hash_new_full (hash_fun, match_fun, NULL, NULL);
	set->hash_fun = hash_fun;
	set->match_fun = match_fun;
	set->destroy_fun = destroy_fun;

	return set;
}

Set *
set_new (DestroyNotify destroy_fun)
{
	return set_new_full (direct_hash, direct_equal, destroy_fun);
}

void
set_free (Set *set)
{
	if (set == NULL)
		return;

	hash_free (set->index);
	list_free (set->list);
	xfree (set);
}

List *
set_list (const Set *set)
{
	return set->list;
}

static inline int
__set_is_member (const Set *set, const void *data)
{
	return hash_contains (set->index, data);
}

int
set_is_member (const Set *set, const void *data)
{
	return __set_is_member (set, data);
}

int
set_insert (Set *set, const void *data)
{
	assert (set != NULL && data != NULL);

	if (__set_is_member (set, data))
		return 0;

	list_append (set->list, data);
	hash_insert (set->index, data, data);

	return 1;
}

int
set_remove (Set *set, void **data)
{
	assert (set != NULL && data != NULL && *data != NULL);

	if (hash_remove (set->index, *data))
		{
			ListElmt *cur = list_head (set->list);

			for (; cur != NULL; cur = list_next (cur))
				if (set->match_fun (*data, list_data (cur)))
					break;

			list_remove (set->list, cur, data);
			return 1;
		}

	return 0;
}

Set *
set_union (const Set *set1, const Set *set2)
{
	assert (set1 != NULL && set2 != NULL);

	Set *setu = NULL;
	ListElmt *cur = NULL;

	setu = set_new_full (set1->hash_fun, set1->match_fun, NULL);

	for (cur = list_head (set1->list); cur != NULL; cur = list_next (cur))
		set_insert (setu, list_data (cur));

	for (cur = list_head (set2->list); cur != NULL; cur = list_next (cur))
		set_insert (setu, list_data (cur));

	return setu;
}

Set *
set_intersection (const Set *set1, const Set *set2)
{
	assert (set1 != NULL && set2 != NULL);

	Set *seti = NULL;
	ListElmt *cur = NULL;

	seti = set_new_full (set1->hash_fun, set1->match_fun, NULL);

	for (cur = list_head (set1->list); cur != NULL; cur = list_next (cur))
		if (__set_is_member (set2, list_data (cur)))
			set_insert (seti, list_data (cur));

	return seti;
}

Set *
set_difference (const Set *set1, const Set *set2)
{
	assert (set1 != NULL && set2 != NULL);

	Set *setd = NULL;
	ListElmt *cur = NULL;

	setd = set_new_full (set1->hash_fun, set1->match_fun, NULL);

	for (cur = list_head (set1->list); cur != NULL; cur = list_next (cur))
		if (!__set_is_member (set2, list_data (cur)))
			set_insert (setd, list_data (cur));

	return setd;
}

int
set_is_subset (const Set *set1, const Set *set2)
{
	assert (set1 != NULL && set2 != NULL);

	ListElmt *cur = NULL;

	if (list_size (set1->list) > list_size (set2->list))
		return 0;

	for (cur = list_head (set1->list); cur != NULL; cur = list_next (cur))
		if (!__set_is_member (set2, list_data (cur)))
			return 0;

	return 1;
}

int
set_is_equal (const Set *set1, const Set *set2)
{
	assert (set1 != NULL && set2 != NULL);

	if (list_size (set1->list) != list_size (set2->list))
		return 0;

	return set_is_subset (set1, set2);
}

size_t
set_size (const Set *set)
{
	return list_size (set->list);
}
