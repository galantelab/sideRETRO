#include <config.h>

#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "wrapper.h"
#include "hash.h"

struct _HashElmt
{
	void *key;
	void *value;
};

typedef struct _HashElmt HashElmt;

static inline HashElmt *
hash_elmt_new (const void *key, const void *value)
{
	HashElmt *hash_elmt = xcalloc (1, sizeof (HashElmt));
	hash_elmt->key = (void *) key;
	hash_elmt->value = (void *) value;
	return hash_elmt;
}

static inline void
hash_elmt_clean (Hash *hash, HashElmt *hash_elmt)
{
	if (hash->destroy_key_fun != NULL)
		hash->destroy_key_fun (hash_elmt->key);
	if (hash->destroy_value_fun != NULL)
		hash->destroy_value_fun (hash_elmt->value);
}

static inline void
hash_elmt_free (Hash *hash, HashElmt *hash_elmt)
{
	hash_elmt_clean (hash, hash_elmt);
	free (hash_elmt);
}

size_t
str_hash (const void *key)
{
	/**
	*
	* copyright (c) 2014 joseph werle <joseph.werle@gmail.com>
	*/
	uint32_t len = strlen ((char *) key);
	uint32_t seed = 0;
	uint32_t c1 = 0xcc9e2d51;
	uint32_t c2 = 0x1b873593;
	uint32_t r1 = 15;
	uint32_t r2 = 13;
	uint32_t m = 5;
	uint32_t n = 0xe6546b64;
	uint32_t h = 0;
	uint32_t k = 0;
	uint8_t *d = (uint8_t *) key; // 32 bit extract from `key'
	const uint32_t *chunks = NULL;
	const uint8_t *tail = NULL; // tail - last 8 bytes
	int i = 0;
	int l = len / 4; // chunk length

	h = seed;

	chunks = (const uint32_t *) (d + l * 4); // body
	tail = (const uint8_t *) (d + l * 4); // last 8 byte chunk of `key'

	// for each 4 byte chunk of `key'
	for (i = -l; i != 0; ++i)
		{
			// next 4 byte chunk of `key'
			k = chunks[i];

			// encode next 4 byte chunk of `key'
			k *= c1;
			k = (k << r1) | (k >> (32 - r1));
			k *= c2;

			// append to hash
			h ^= k;
			h = (h << r2) | (h >> (32 - r2));
			h = h * m + n;
		}

	k = 0;

	// remainder
	switch (len & 3)
		{ // `len % 4'
		case 3: k ^= (tail[2] << 16);
		case 2: k ^= (tail[1] << 8);

		case 1:
			k ^= tail[0];
			k *= c1;
			k = (k << r1) | (k >> (32 - r1));
			k *= c2;
			h ^= k;
		}

	h ^= len;

	h ^= (h >> 16);
	h *= 0x85ebca6b;
	h ^= (h >> 13);
	h *= 0xc2b2ae35;
	h ^= (h >> 16);

	return (size_t) h;
}

int
str_equal (const void *key1, const void *key2)
{
	return !strcmp ((const char *) key1, (const char *) key2);
}

static inline size_t
hash_calculate_bucket (Hash *hash, const void *key)
{
	return (hash->hash_fun (key) * 11) % hash->buckets;
}

Hash *
hash_new_full (size_t buckets, HashFunc hash_fun, EqualFun match_fun,
		DestroyNotify destroy_key_fun, DestroyNotify destroy_value_fun)
{
	assert (buckets > 0 && hash_fun != NULL
			&& match_fun != NULL);

	Hash *hash = NULL;
	size_t i = 0;

	hash = xcalloc (1, sizeof (Hash));
	hash->table = xcalloc (buckets, sizeof (List *));

	for (; i < buckets; i++)
		hash->table[i] = list_new (NULL);

	hash->buckets = buckets;
	hash->hash_fun = hash_fun;
	hash->match_fun = match_fun;
	hash->destroy_key_fun = destroy_key_fun;
	hash->destroy_value_fun = destroy_value_fun;

	return hash;
}

Hash *
hash_new (size_t buckets, DestroyNotify destroy_key_fun,
		DestroyNotify destroy_value_fun)
{
	return hash_new_full (buckets, str_hash, str_equal,
			destroy_key_fun, destroy_value_fun);
}

void
hash_free (Hash *hash)
{
	if (hash == NULL)
		return;

	size_t i = 0;

	for (; i < hash->buckets; i++)
		{
			List *list = hash->table[i];

			while (list_size (list) > 0)
				{
					ListElmt *list_elmt = list_tail (list);
					HashElmt *hash_elmt = list_data (list_elmt);
					hash_elmt_free (hash, hash_elmt);
					list_remove (list, list_elmt, NULL);
				}

			list_free (list);
		}

	free (hash->table);
	free (hash);
}

static ListElmt *
hash_fetch (Hash *hash, const void *key)
{
	assert (hash != NULL && key != NULL);

	size_t bucket = 0;
	size_t i = 0;

	bucket = hash_calculate_bucket (hash, key);
	ListElmt *cur = list_head (hash->table[bucket]);

	while (cur != NULL)
		{
			HashElmt *element = list_data (cur);
			if (hash->match_fun (element->key, key))
				break;
			cur = list_next (cur);
		}

	return cur;
}

int
hash_remove (Hash *hash, const void *key)
{
	ListElmt *list_elmt = hash_fetch (hash, key);

	if (list_elmt != NULL)
		{
			HashElmt *hash_elmt = list_data (list_elmt);

			size_t bucket = hash_calculate_bucket (hash, key);
			List *list = hash->table[bucket];

			hash_elmt_free (hash, hash_elmt);
			list_remove (list, list_elmt, NULL);

			hash->size--;
			hash->version++;
			return 1;
		}

	return 0;
}

int
hash_insert (Hash *hash, const void *key, const void *value)
{
	ListElmt *list_elmt = hash_fetch (hash, key);

	if (list_elmt != NULL)
		{
			HashElmt *hash_elmt = list_data (list_elmt);
			hash_elmt_clean (hash, hash_elmt);

			hash_elmt->key = (void *) key;
			hash_elmt->value = (void *) value;

			return 0;
		}

	size_t bucket = hash_calculate_bucket (hash, key);
	list_append (hash->table[bucket], hash_elmt_new (key, value));

	hash->size++;
	hash->version++;
	return 1;
}

void *
hash_lookup (Hash *hash, const void *key)
{
	ListElmt *list_elmt = hash_fetch (hash, key);

	if (list_elmt != NULL)
		{
			HashElmt *hash_elmt = list_data (list_elmt);
			return hash_elmt->value;
		}

	return NULL;
}

int
hash_contains (Hash *hash, const void *key)
{
	return hash_fetch (hash, key) != NULL;
}

static inline int
hash_find_bucket_in_use (Hash *hash, size_t *from)
{
	assert (from != NULL);
	size_t i = *from;

	for (; i < hash->buckets; i++)
		{
			List *list = hash->table[i];
			if (list_size (list) > 0)
				break;
		}

	*from = i;
	return i < hash->buckets;
}

void
hash_foreach (Hash *hash, HFunc func, void *user_data)
{
	assert (hash != NULL && func != NULL);
	size_t i = 0;

	while (hash_find_bucket_in_use (hash, &i))
		{
			List *list = hash->table[i];
			ListElmt *cur = list_head (list);

			while (cur != NULL)
				{
					HashElmt *hash_elmt = list_data (cur);
					func (hash_elmt->key, hash_elmt->value, user_data);
					cur = list_next (cur);
				}

			i++;
		}
}

void
hash_iter_init (HashIter *iter, Hash *hash)
{
	assert (iter != NULL && hash != NULL);
	memset (iter, 0, sizeof (HashIter));

	size_t i = 0;

	if (hash_find_bucket_in_use (hash, &i))
		{
			List *list = hash->table[i];
			iter->cur = list_head (list);
		}

	iter->bucket = i;
	iter->version = hash->version;
	iter->hash = hash;
}

int
hash_iter_next (HashIter *iter, void **key, void **value)
{
	assert (iter != NULL);
	assert (iter->version == iter->hash->version);

	Hash *hash = iter->hash;
	HashElmt *hash_elmt = NULL;
	ListElmt *cur = iter->cur;

	if (cur == NULL)
		{
			size_t i = iter->bucket + 1;

			if (hash_find_bucket_in_use (hash, &i))
				{
					List *list = hash->table[i];
					iter->bucket = i;
					cur = list_head (list);
				}
			else
				return 0;
		}

	hash_elmt = list_data (cur);

	if (key != NULL)
		*key = hash_elmt->key;

	if (value != NULL)
		*value = hash_elmt->value;

	iter->cur = list_next (cur);
	return 1;
}

List *
hash_get_keys_as_list (Hash *hash)
{
	assert (hash != NULL);

	HashIter iter;
	List *keys = NULL;
	void *key = NULL;

	keys = list_new (NULL);
	hash_iter_init (&iter, hash);

	while (hash_iter_next (&iter, &key, NULL))
		list_append (keys, key);

	return keys;
}

List *
hash_get_values_as_list (Hash *hash)
{
	assert (hash != NULL);

	HashIter iter;
	List *values = NULL;
	void *value = NULL;

	values = list_new (NULL);
	hash_iter_init (&iter, hash);

	while (hash_iter_next (&iter, NULL, &value))
		list_append (values, value);

	return values;
}

Array *
hash_get_keys_as_array (Hash *hash)
{
	assert (hash != NULL);

	HashIter iter;
	Array *keys = NULL;
	void *key = NULL;

	keys = array_new (NULL);
	hash_iter_init (&iter, hash);

	while (hash_iter_next (&iter, &key, NULL))
		array_add (keys, key);

	return keys;
}

Array *
hash_get_values_as_array (Hash *hash)
{
	assert (hash != NULL);

	HashIter iter;
	Array *values = NULL;
	void *value = NULL;

	values = array_new (NULL);
	hash_iter_init (&iter, hash);

	while (hash_iter_next (&iter, NULL, &value))
		array_add (values, value);

	return values;
}
