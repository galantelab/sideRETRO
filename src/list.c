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
#include "list.h"

List *
list_new (DestroyNotify destroy_fun)
{
	List *list = xcalloc (1, sizeof (List));
	list->destroy_fun = destroy_fun;
	return list;
}

static inline ListElmt *
list_element_new (const void *data)
{
	assert (data != NULL);
	ListElmt *element = xcalloc (1, sizeof (ListElmt));
	element->data = (void *) data;
	return element;
}

void
list_ins_next (List *list, ListElmt *element, const void *data)
{
	ListElmt *new_element = list_element_new (data);
	list_ins_next_link (list, element, new_element);
}

void
list_ins_next_link (List *list, ListElmt *element, ListElmt *new_element)
{
	// Do not allow a NULL element unless the list is empty
	assert ((element == NULL && list->size == 0)
			|| (element != NULL && list->size > 0));

	// Do nor allow a NULL new_element
	assert (new_element != NULL);

	if (list->size == 0)
		{
			// Handle insertion when the list is empty
			list->head = new_element;
			list->tail = new_element;
		}
	else
		{
			// Handle insertion when the list is not empty
			new_element->next = element->next;
			new_element->prev = element;

			if (element->next == NULL)
				list->tail = new_element;
			else
				element->next->prev = new_element;

			element->next = new_element;
		}

	list->size++;
}

void
list_ins_prev (List *list, ListElmt *element, const void *data)
{
	ListElmt *new_element = list_element_new (data);
	list_ins_prev_link (list, element, new_element);
}

void
list_ins_prev_link (List *list, ListElmt *element, ListElmt *new_element)
{
	// Do not allow a NULL element unless the list is empty
	assert ((element == NULL && list->size == 0)
			|| (element != NULL && list->size > 0));

	// Do nor allow a NULL new_element
	assert (new_element != NULL);

	if (list->size == 0)
		{
			// Handle insertion when the list is empty
			list->head = new_element;
			list->tail = new_element;
		}
	else
		{
			// Handle insertion when the list is not empty
			new_element->next = element;
			new_element->prev = element->prev;

			if (element->prev == NULL)
				list->head = new_element;
			else
				element->prev->next = new_element;

			element->prev = new_element;
		}

	list->size++;
}

void
list_free (List *list)
{
	if (list == NULL)
		return;

	while (list->size > 0)
		list_remove (list, list->tail, NULL);

	xfree (list);
}

void
list_remove (List *list, ListElmt *element, void **data)
{
	assert (list != NULL && element != NULL);
	list_remove_link (list, element);

	if (data != NULL)
		*data = element->data;
	else
		{
			if (list->destroy_fun != NULL)
				list->destroy_fun (element->data);
		}

	xfree (element);
}

void
list_remove_link (List *list, ListElmt *element)
{
	assert (list != NULL && list->size > 0 && element != NULL);

	if (element == list->head)
		{
			// Handle removal from the head of the list
			list->head = element->next;
			if (list->head == NULL)
				list->tail = NULL;
			else
				element->next->prev = NULL;
		}
	else
		{
			// Handle removal from other than the head of the list
			element->prev->next = element->next;
			if (element->next == NULL)
				list->tail = element->prev;
			else
				element->next->prev = element->prev;
		}

	element->prev = NULL;
	element->next = NULL;
	list->size--;
}

void
list_foreach (List *list, Func func, void *user_data)
{
	assert (list != NULL && func != NULL);
	ListElmt *cur = list->head;

	for (; cur != NULL; cur = cur->next)
		func (cur->data, user_data);
}
