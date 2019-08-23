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

uint32_t str_hash  (const void *key);
uint32_t int_hash  (const void *key);
int      str_equal (const void *key1, const void *key2);
int      int_equal (const void *key1, const void *key2);

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
