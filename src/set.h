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

#ifndef SET_H
#define SET_H

#include <stdlib.h>
#include "hash.h"
#include "list.h"
#include "types.h"

typedef struct _Set Set;

Set    * set_new_full     (HashFunc hash_fun, EqualFun match_fun, DestroyNotify destroy_fun);
Set    * set_new          (DestroyNotify destroy_fun);
void     set_free         (Set *set);
List   * set_list         (const Set *set);
int      set_is_member    (const Set *set, const void *data);
int      set_insert       (Set *set, const void *data);
int      set_remove       (Set *set, void **data);
Set    * set_union        (const Set *set1, const Set *set2);
Set    * set_intersection (const Set *set1, const Set *set2);
Set    * set_difference   (const Set *set1, const Set *set2);
int      set_is_subset    (const Set *set1, const Set *set2);
int      set_is_equal     (const Set *set1, const Set *set2);
size_t   set_size         (const Set *set);

#endif /* set.h */
