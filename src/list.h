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

#ifndef LIST_H
#define LIST_H

#include <stdlib.h>
#include "types.h"

struct _ListElmt
{
	void             *data;
	struct _ListElmt *next;
	struct _ListElmt *prev;
};

typedef struct _ListElmt ListElmt;

struct _List
{
	size_t         size;
	ListElmt      *head;
	ListElmt      *tail;
	DestroyNotify  destroy_fun;
};

typedef struct _List List;

List *list_new           (DestroyNotify destroy_fun);
void  list_free          (List *list);
void  list_ins_prev      (List *list, ListElmt *element, const void *data);
void  list_ins_prev_link (List *list, ListElmt *element, ListElmt *new_element);
void  list_ins_next      (List *list, ListElmt *element, const void *data);
void  list_ins_next_link (List *list, ListElmt *element, ListElmt *new_element);
void  list_remove        (List *list, ListElmt *element, void **data);
void  list_remove_link   (List *list, ListElmt *element);
void  list_foreach       (List *list, Func func, void *user_data);

#define list_append(list,data) \
	(list_ins_next ((list), (list)->tail, (data)))
#define list_append_link(list,element) \
	(list_ins_next_link ((list), (list)->tail, (element)))
#define list_prepend(list,data) \
	(list_ins_prev ((list), (list)->head, (data)))
#define list_prepend_link(list,element) \
	(list_ins_prev_link ((list), (list)->head, (element)))

#define list_size(list)       ((list)->size)
#define list_head(list)       ((list)->head)
#define list_tail(list)       ((list)->tail)
#define list_is_head(element) ((element)->prev == NULL ? 1 : 0)
#define list_is_tail(element) ((element)->next == NULL ? 1 : 0)
#define list_data(element)    ((element)->data)
#define list_prev(element)    ((element)->prev)
#define list_next(element)    ((element)->next)

#endif /* list.h */
