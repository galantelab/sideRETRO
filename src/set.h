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
int      set_insert       (Set *set, const void *data);
//void     set_insert_all   (Set *setu, const Set *set);
int      set_remove       (Set *set, void **data);
//void     set_remove_all   (Set *setd, const Set *set);
Set    * set_union        (const Set *set1, const Set *set2);
Set    * set_intersection (const Set *set1, const Set *set2);
Set    * set_difference   (const Set *set1, const Set *set2);
int      set_is_subset    (const Set *set1, const Set *set2);
int      set_is_equal     (const Set *set1, const Set *set2);
size_t   set_size         (const Set *set);

#endif /* set.h */
