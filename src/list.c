#include <config.h>

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

	free (list);
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

	free (element);
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
