#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "hash.h"

#define GLOBAL_DEPTH 3
#define LOCAL_DEPTH  3
#define BUCKET_SIZE  20

enum _HashColor
{
	HASH_BUCKET_WHITE,
	HASH_BUCKET_BLACK
};

typedef enum _HashColor HashColor;

struct _Record
{
	void     *key;
	void     *value;
	uint32_t  hash;
};

typedef struct _Record Record;

struct _Bucket
{
	size_t     depth;
	size_t     size;
	HashColor  color;
	Record    *records;
};

typedef struct _Bucket Bucket;

struct _Hash
{
	size_t          version;
	size_t          size;

	HashFunc        hash_fun;
	EqualFun        match_fun;
	DestroyNotify   destroy_key_fun;
	DestroyNotify   destroy_value_fun;

	size_t          global_depth;
	size_t          buckets;
	Bucket        **directory;
};

struct _RealIter
{
	size_t  version;
	size_t  bucket;
	size_t  record;
	Hash   *hash;
};

typedef struct _RealIter RealIter;

uint32_t
str_hash (const void *key)
{
	/*
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

	return h;
}

int
str_equal (const void *key1, const void *key2)
{
	return !strcmp ((const char *) key1, (const char *) key2);
}

uint32_t
int_hash (const void *key)
{
	return * (const int *) key;
}

int
int_equal (const void *key1, const void *key2)
{
	return * ((const int *) key1) == * ((const int *) key2);
}

static inline void
hash_record_destroy (Hash *hash, Record *record)
{
	if (hash->destroy_key_fun != NULL)
		hash->destroy_key_fun (record->key);
	if (hash->destroy_value_fun != NULL)
		hash->destroy_value_fun (record->value);
}

static Bucket *
hash_bucket_new (void)
{
	Bucket *bucket = NULL;
	size_t depth = LOCAL_DEPTH;
	size_t size = BUCKET_SIZE + 1;

	bucket = xcalloc (1, sizeof (Bucket));
	bucket->records = xcalloc (size, sizeof (Record));
	bucket->depth = depth;
	bucket->color = HASH_BUCKET_WHITE;

	return bucket;
}

static void
hash_bucket_free (Hash *hash, Bucket *bucket)
{
	Record *records = bucket->records;
	size_t i = 0;

	for (; i < bucket->size; i++)
		hash_record_destroy (hash, &records[i]);

	xfree (records);
	xfree (bucket);
}

static inline int
hash_bucket_is_full (Bucket *bucket)
{
	return bucket->size > BUCKET_SIZE;
}

Hash *
hash_new_full (HashFunc hash_fun, EqualFun match_fun,
		DestroyNotify destroy_key_fun, DestroyNotify destroy_value_fun)
{
	assert (hash_fun != NULL && match_fun != NULL);

	Hash *hash = NULL;
	size_t buckets = 1 << GLOBAL_DEPTH;
	size_t i = 0;

	hash = xcalloc (1, sizeof (Hash));
	hash->directory = xcalloc (buckets, sizeof (Bucket *));

	for (; i < buckets; i++)
		hash->directory[i] = hash_bucket_new ();

	hash->buckets = buckets;
	hash->global_depth = GLOBAL_DEPTH;

	hash->hash_fun = hash_fun;
	hash->match_fun = match_fun;
	hash->destroy_key_fun = destroy_key_fun;
	hash->destroy_value_fun = destroy_value_fun;

	return hash;
}

Hash *
hash_new (DestroyNotify destroy_key_fun, DestroyNotify destroy_value_fun)
{
	return hash_new_full (str_hash, str_equal,
			destroy_key_fun, destroy_value_fun);
}

size_t
hash_size (Hash *hash)
{
	return hash->size;
}

static inline void
hash_clean_buckets (Hash *hash)
{
	size_t i = 0;
	for (; i < hash->buckets; i++)
		hash->directory[i]->color = HASH_BUCKET_WHITE;
}

void
hash_free (Hash *hash)
{
	if (hash == NULL)
		return;

	Bucket **bucket = NULL;
	size_t i = 0;

	hash_clean_buckets (hash);

	for (i = 0; i < hash->buckets; i++)
		{
			bucket = &hash->directory[i];

			switch ((*bucket)->color)
				{
				case HASH_BUCKET_WHITE:
					{
						(*bucket)->color = HASH_BUCKET_BLACK;
						break;
					}
				case HASH_BUCKET_BLACK:
					{
						*bucket = NULL;
						break;
					}
				}
		}

	for (i = 0; i < hash->buckets; i++)
		{
			bucket = &hash->directory[i];

			if (*bucket != NULL)
				hash_bucket_free (hash, *bucket);
		}

	xfree (hash->directory);
	xfree (hash);
}

static inline uint32_t
hash_probing_func (Hash *hash, const void *key)
{
	return hash->hash_fun (key) * 11;
}

static inline Bucket *
hash_get_bucket (Hash *hash, const void *key, uint32_t *hash_key)
{
	*hash_key = hash_probing_func (hash, key);
	size_t index = *hash_key & ((1 << hash->global_depth) - 1);
	return hash->directory[index];
}

static inline Record *
hash_fetch_bucket_record (Hash *hash, Bucket *bucket,
		const void *key, uint32_t hash_key, size_t *i_)
{
	Record *records = bucket->records;
	size_t i = 0;

	for (; i < bucket->size; i++)
		{
			if (records[i].hash == hash_key
					&& hash->match_fun (records[i].key, key))
				{
					*i_ = i;
					return &records[i];
				}
		}

	return NULL;
}

static inline int
hash_insert_bucket_record (Hash *hash, Bucket *bucket,
		const void *key, const void *value, uint32_t hash_key)
{
	Record *record = NULL;
	size_t i = 0;
	int rc = 0;

	record = hash_fetch_bucket_record (hash, bucket,
			key, hash_key, &i);

	if (record == NULL)
		{
			record = &bucket->records[bucket->size];
			bucket->size++;
			rc = 1;
		}
	else
		hash_record_destroy (hash, record);

	record->key   = (void *) key;
	record->value = (void *) value;
	record->hash  = hash_key;

	return rc;
}

static inline int
hash_remove_bucket_record (Hash *hash, Bucket *bucket,
		const void *key, uint32_t hash_key)
{
	Record *record = NULL;
	size_t i = 0;

	record = hash_fetch_bucket_record (hash, bucket,
			key, hash_key, &i);

	if (record != NULL)
		{
			hash_record_destroy (hash, record);

			for (i++; i < bucket->size; i++)
				bucket->records[i - 1] = bucket->records[i];

			bucket->size--;
			return 1;
		}

	return 0;
}

static inline void
hash_duplicate_buckets (Hash *hash)
{
	size_t old_buckets = hash->buckets;
	size_t i = 0;

	hash->buckets *= 2;
	hash->directory = xrealloc (hash->directory,
			sizeof (Bucket *) * hash->buckets);

	for (i = old_buckets; i < hash->buckets; i++)
		hash->directory[i] = hash->directory[i - old_buckets];
}

static void
hash_maybe_expand (Hash *hash, Bucket *bucket,
		uint32_t hash_key)
{
	if (!hash_bucket_is_full (bucket))
		return;

	if (bucket->depth == hash->global_depth)
		{
			hash_duplicate_buckets (hash);
			hash->global_depth++;
		}

	Bucket *bucket0 = hash_bucket_new ();
	Bucket *bucket1 = hash_bucket_new ();

	bucket0->depth = bucket1->depth = bucket->depth + 1;
	size_t bit = 1 << bucket->depth;

	Bucket **new_bucket = NULL;
	Record *record = NULL;
	size_t i = 0;

	for (; i < bucket->size; i++)
		{
			record = &bucket->records[i];

			new_bucket = record->hash & bit
				? &bucket1
				: &bucket0;

			(*new_bucket)->records[(*new_bucket)->size] = *record;
			(*new_bucket)->size++;
		}

	for (i = hash_key & (bit - 1); i < hash->buckets; i += bit)
		hash->directory[i] = i & bit ? bucket1 : bucket0;

	bucket->size = 0;
	hash_bucket_free (hash, bucket);

	hash_maybe_expand (hash, bucket0, hash_key);
	hash_maybe_expand (hash, bucket1, hash_key);
}

int
hash_insert (Hash *hash, const void *key, const void *value)
{
	assert (hash != NULL && key != NULL);

	Bucket *bucket = NULL;
	uint32_t hash_key = 0;

	bucket = hash_get_bucket (hash, key, &hash_key);

	if (!hash_insert_bucket_record (hash, bucket,
				key, value, hash_key))
		return 0;

	hash_maybe_expand (hash, bucket, hash_key);

	hash->size++;
	hash->version++;
	return 1;
}

void *
hash_lookup (Hash *hash, const void *key)
{
	assert (hash != NULL && key != NULL);

	Bucket *bucket = NULL;
	Record *record = NULL;

	uint32_t hash_key = 0;
	size_t i = 0;

	bucket = hash_get_bucket (hash, key, &hash_key);
	record = hash_fetch_bucket_record (hash, bucket,
			key, hash_key, &i);

	return record != NULL ? record->value : NULL;
}

int
hash_contains (Hash *hash, const void *key)
{
	return hash_lookup (hash, key) != NULL;
}

int
hash_remove (Hash *hash, const void *key)
{
	assert (hash != NULL && key != NULL);

	Bucket *bucket = NULL;
	uint32_t hash_key = 0;

	bucket = hash_get_bucket (hash, key, &hash_key);

	if (!hash_remove_bucket_record (hash, bucket,
				key, hash_key))
		return 0;

	hash->size--;
	hash->version++;
	return 1;
}

static inline Bucket *
hash_find_bucket (Hash *hash, size_t *from)
{
	Bucket *bucket = NULL;
	size_t i = *from;

	for (; i < hash->buckets; i++)
		{
			bucket = hash->directory[i];
			if (bucket->color == HASH_BUCKET_WHITE
					&& bucket->size > 0)
				{
					*from = i;
					return bucket;
				}
		}

	return NULL;
}

void
hash_foreach (Hash *hash, HFunc func, void *user_data)
{
	assert (hash != NULL && func != NULL);

	Bucket *bucket = NULL;
	Record *record = NULL;
	size_t i = 0;
	size_t j = 0;

	hash_clean_buckets (hash);

	while (hash_find_bucket (hash, &i))
		{
			bucket = hash->directory[i];

			for (j = 0; j < bucket->size; j++)
				{
					record = &bucket->records[j];
					func (record->key, record->value, user_data);
				}

			bucket->color = HASH_BUCKET_BLACK;
		}
}

void
hash_iter_init (HashIter *iter, Hash *hash)
{
	assert (iter != NULL && hash != NULL);

	RealIter *ri = (RealIter *) iter;

	hash_clean_buckets (hash);

	ri->bucket = 0;
	ri->record = 0;
	ri->version = hash->version;
	ri->hash = hash;
}

int
hash_iter_next (HashIter *iter, void **key, void **value)
{
	assert (iter != NULL);

	RealIter *ri = (RealIter *) iter;
	assert (ri->version == ri->hash->version);

	Hash *hash = NULL;
	Bucket *bucket = NULL;
	Record *record = NULL;
	size_t i = 0;
	size_t j = 0;

	hash = ri->hash;

	i = ri->bucket;
	j = ri->record;

	bucket = hash->directory[i];

	if (j == bucket->size)
		{
			bucket->color = HASH_BUCKET_BLACK;

			if (hash_find_bucket (hash, &i) == NULL)
				return 0;

			bucket = hash->directory[i];
			j = 0;
		}

	record = &bucket->records[j];

	if (key != NULL)
		*key = record->key;

	if (value != NULL)
		*value = record->value;

	ri->record = j + 1;
	ri->bucket = i;
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
