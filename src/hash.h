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

#ifndef HASH_H
#define HASH_H

#include <stdlib.h>
#include <stdint.h>
#include "list.h"
#include "array.h"
#include "types.h"

typedef uint32_t (*HashFunc) (const void *key);
typedef void     (*HFunc)    (void *key, void *value, void *user_data);

typedef struct _Hash Hash;

struct _HashIter
{
	/* PRIVATE  */
	size_t  dummy1;
	size_t  dummy2;
	size_t  dummy3;
	void   *dummy4;
};

typedef struct _HashIter HashIter;

uint32_t str_hash     (const void *key);
uint32_t int_hash     (const void *key);
uint32_t direct_hash  (const void *key);
int      str_equal    (const void *key1, const void *key2);
int      int_equal    (const void *key1, const void *key2);
int      direct_equal (const void *key1, const void *key2);

Hash  * hash_new_full (HashFunc hash_fun, EqualFun match_fun,
		DestroyNotify destroy_key_fun, DestroyNotify destroy_value_fun);

Hash  * hash_new (DestroyNotify destroy_key_fun, DestroyNotify destroy_value_fun);

void    hash_free (Hash *hash);
size_t  hash_size (Hash *hash);

int     hash_insert (Hash *hash, const void *key, const void *value);
int     hash_remove (Hash *hash, const void *key);

void  * hash_lookup   (Hash *hash, const void *key);
int     hash_contains (Hash *hash, const void *key);

void    hash_foreach  (Hash *hash, HFunc func, void *user_data);

void    hash_iter_init (HashIter *iter, Hash *hash);
int     hash_iter_next (HashIter *iter, void **key, void **value);

List  * hash_get_keys_as_list    (Hash *hash);
List  * hash_get_values_as_list  (Hash *hash);
Array * hash_get_keys_as_array   (Hash *hash);
Array * hash_get_values_as_array (Hash *hash);

#endif /* hash.h */
