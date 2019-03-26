#ifndef HASH_H
#define HASH_H

#include <stdlib.h>
#include <stdint.h>
#include "list.h"
#include "types.h"

typedef uint32_t (*HashFunc) (const void *key);

struct _Hash
{
	size_t          buckets;

	HashFunc        hash_fun;
	EqualFun        match_fun;
	DestroyNotify   destroy_key_fun;
	DestroyNotify   destroy_value_fun;

	size_t          size;
	List          **table;
};

typedef struct _Hash Hash;

#define HASH_SMALL_SIZE  1021
#define HASH_MEDIUM_SIZE 65521
#define HASH_LARGE_SIZE  1048573

uint32_t str_hash  (const void *key);
int      str_equal (const void *key1, const void *key2);

Hash * hash_new_full (size_t buckets, HashFunc hash_fun, EqualFun match_fun,
		DestroyNotify destroy_key_fun, DestroyNotify destroy_value_fun);
Hash * hash_new (size_t buckets, DestroyNotify destroy_key_fun,
		DestroyNotify destroy_value_fun);

void   hash_free     (Hash *hash);
int    hash_remove   (Hash *hash, const void *key);
int    hash_insert   (Hash *hash, const void *key, const void *value);
void * hash_lookup   (Hash *hash, const void *key);
int    hash_contains (Hash *hash, const void *key);

#define hash_size(hash) ((hash)->size)

#endif /* hash.h */
