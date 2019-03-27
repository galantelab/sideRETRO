#ifndef HASH_H
#define HASH_H

#include <stdlib.h>
#include <stdint.h>
#include "list.h"
#include "types.h"

typedef uint32_t (*HashFunc) (const void *key);
typedef void     (*HFunc)    (void *key, void *value, void *user_data);

struct _Hash
{
	size_t          version;
	size_t          buckets;

	HashFunc        hash_fun;
	EqualFun        match_fun;
	DestroyNotify   destroy_key_fun;
	DestroyNotify   destroy_value_fun;

	size_t          size;
	List          **table;
};

typedef struct _Hash Hash;

struct _HashIter
{
	size_t    version;
	size_t    bucket;
	Hash     *hash;
	ListElmt *cur;
};

typedef struct _HashIter HashIter;

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
void   hash_foreach  (Hash *hash, HFunc func, void *user_data);

void   hash_iter_init (HashIter *iter, Hash *hash);
int    hash_iter_next (HashIter *iter, void **key, void **value);

#define hash_size(hash) ((hash)->size)

#endif /* hash.h */
