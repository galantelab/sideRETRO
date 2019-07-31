#include "config.h"

#include <stdlib.h>
#include <check.h>
#include "check_sider.h"
#include "../src/wrapper.h"
#include "../src/hash.h"

START_TEST (test_hash_new)
{
	Hash *hash = hash_new_full (HASH_MEDIUM_SIZE, str_hash,
			str_equal, free, free);

	ck_assert_ptr_ne (hash, NULL);
	ck_assert_ptr_ne (hash->table, NULL);
	ck_assert_int_eq (hash->size, 0);
	ck_assert_int_eq (hash->buckets, HASH_MEDIUM_SIZE);

	hash_free (hash);
}
END_TEST

START_TEST (test_hash_insert)
{
	Hash *hash = hash_new (10, NULL, free);

	hash_insert (hash, "eu", xstrdup ("Thiago"));
	hash_insert (hash, "tu", xstrdup ("Fernanda"));
	hash_insert (hash, "ele", xstrdup ("Pedro"));

	ck_assert_int_eq (hash_size (hash), 3);

	hash_free (hash);
}
END_TEST

START_TEST (test_hash_lookup)
{
	Hash *hash = hash_new (10, NULL, free);

	hash_insert (hash, "eu", xstrdup ("Thiago"));
	hash_insert (hash, "tu", xstrdup ("Fernanda"));
	hash_insert (hash, "ele", xstrdup ("Pedro"));

	ck_assert_str_eq (hash_lookup (hash, "eu"), "Thiago");
	ck_assert_str_eq (hash_lookup (hash, "tu"), "Fernanda");
	ck_assert_str_eq (hash_lookup (hash, "ele"), "Pedro");

	hash_free (hash);
}
END_TEST

START_TEST (test_hash_replace)
{
	Hash *hash = hash_new (10, NULL, free);
	int *num_alloc = xcalloc (1, sizeof (int));
	int num;
	int rc;

	rc = hash_insert (hash, "eu", xstrdup ("Thiago"));
	ck_assert_int_eq (rc, 1);

	*num_alloc = 120;

	rc = hash_insert (hash, "eu", num_alloc);
	ck_assert_int_eq (rc, 0);

	num = * (int *) hash_lookup (hash, "eu");
	ck_assert_int_eq (num, *num_alloc);
	ck_assert_int_eq (num, 120);

	hash_free (hash);
}
END_TEST

START_TEST (test_hash_remove)
{
	Hash *hash = hash_new (10, free, NULL);
	int i = 0;
	int rc;
	char **str;

	str = xcalloc (20, sizeof (char *));

	for (i = 0; i < 20; i++)
		{
			int len = xasprintf (&str[i], "loop%d", i);
			hash_insert (hash, str[i], str[i]);
		}

	ck_assert_int_eq (hash_size (hash), 20);

	for (i = 0; i < 20; i++)
		{
			rc = hash_remove (hash, str[i]);
			ck_assert_int_eq (rc, 1);
			ck_assert_int_eq (hash_size (hash), 20 - i - 1);
		}

	rc = hash_remove (hash, "loop1");
	ck_assert_int_eq (rc, 0);
	ck_assert_int_eq (hash_size (hash), 0);

	xfree (str);
	hash_free (hash);
}
END_TEST

START_TEST (test_hash_contains)
{
	Hash *hash = hash_new (50, free, NULL);
	char **keys = xcalloc (100, sizeof (char *));
	int *values = xcalloc (100, sizeof (int));
	int i = 0;

	for (i = 0; i < 100; i++)
		{
			values[i] = i * 10;
			int len = xasprintf (&keys[i], "loop%d", i);
			hash_insert (hash, keys[i], &values[i]);
		}

	for (i = 0; i < 100; i++)
		{
			ck_assert_int_eq (hash_contains (hash, keys[i]), 1);
			int *ptr = hash_lookup (hash, keys[i]);
			ck_assert_int_eq (*ptr, values[i]);
		}

	hash_free (hash);
	xfree (keys);
	xfree (values);
}
END_TEST

static void
sum (void *key, void *value, void *user_data)
{
	char *keys[] = {"loop1", "loop10", "loop50"};
	int i = 0;

	for  (; i < sizeof (keys)/sizeof (char *); i++)
		{
			if (!strcmp (keys[i], (char *) key))
				{
					* (int *) user_data += * (int *) value;
					break;
				}
		}
}

START_TEST (test_hash_foreach)
{
	Hash *hash = hash_new (50, free, NULL);
	char **keys = xcalloc (100, sizeof (char *));
	int *values = xcalloc (100, sizeof (int));
	int i = 0;
	int total = 0;

	for (i = 0; i < 100; i++)
		{
			values[i] = i;
			int len = xasprintf (&keys[i], "loop%d", i);
			hash_insert (hash, keys[i], &values[i]);
		}

	hash_foreach (hash, sum, &total);
	ck_assert_int_eq (total, 61);

	hash_free (hash);
	xfree (keys);
	xfree (values);
}
END_TEST

START_TEST (test_hash_iter)
{
	Hash *hash = hash_new (50, free, NULL);
	char **keys = xcalloc (100, sizeof (char *));
	int *values = xcalloc (100, sizeof (int));
	int i = 0;

	for (i = 0; i < 100; i++)
		{
			values[i] = i;
			int len = xasprintf (&keys[i], "loop%d", i);
			hash_insert (hash, keys[i], &values[i]);
		}

	HashIter iter;
	int *value_ptr;
	int total = 0;

	hash_iter_init (&iter, hash);

	while (hash_iter_next (&iter, NULL, (void **) &value_ptr))
		total += *value_ptr;

	ck_assert_int_eq (total, 4950);

	hash_free (hash);
	xfree (keys);
	xfree (values);
}
END_TEST

START_TEST (test_hash_get_keys_as_list)
{
	Hash *hash = hash_new (10, NULL, NULL);
	char *pronouns[] = { "I", "you", "he", "she", "it", "we", "they" };
	int num_elem = sizeof (pronouns)/sizeof (char *);
	int i = 0;

	for (; i < num_elem; i++)
		hash_insert (hash, pronouns[i], NULL);

	List *list = hash_get_keys_as_list (hash);
	ck_assert_ptr_ne (list, NULL);
	ck_assert_int_eq (list_size (list), num_elem);

	ListElmt *cur = list_head (list);
	int match = 0;

	for (; cur != NULL; cur = list_next (cur))
		for (i = 0; i < num_elem; i++)
			if (list_data (cur) == pronouns[i])
				match ++;

	ck_assert_int_eq (match, num_elem);

	hash_free (hash);
	list_free (list);
}
END_TEST

START_TEST (test_hash_get_values_as_list)
{
	Hash *hash = hash_new (10, free, NULL);
	char *pronouns[] = { "I", "you", "he", "she", "it", "we", "they" };
	int num_elem = sizeof (pronouns)/sizeof (char *);
	int i = 0;

	for (; i < num_elem; i++)
		{
			char *key;
			int len = xasprintf (&key, "key%d", i);
			hash_insert (hash, key, pronouns[i]);
		}

	List *list = hash_get_values_as_list (hash);
	ck_assert_ptr_ne (list, NULL);
	ck_assert_int_eq (list_size (list), num_elem);

	ListElmt *cur = list_head (list);
	int match = 0;

	for (; cur != NULL; cur = list_next (cur))
		for (i = 0; i < num_elem; i++)
			if (list_data (cur) == pronouns[i])
				match ++;

	ck_assert_int_eq (match, num_elem);

	hash_free (hash);
	list_free (list);
}
END_TEST

START_TEST (test_hash_get_keys_as_array)
{
	Hash *hash = hash_new (10, NULL, NULL);
	char *pronouns[] = { "I", "you", "he", "she", "it", "we", "they" };
	int num_elem = sizeof (pronouns)/sizeof (char *);
	int i = 0;

	for (; i < num_elem; i++)
		hash_insert (hash, pronouns[i], NULL);

	Array *array = hash_get_keys_as_array (hash);
	ck_assert_ptr_ne (array, NULL);
	ck_assert_int_eq (array_len (array), num_elem);

	int match = 0;
	int j = 0;

	for (i = 0; i < num_elem; i++)
		for (j = 0; j < num_elem; j++)
			if (array_get (array, i) == pronouns[j])
				match ++;

	ck_assert_int_eq (match, num_elem);

	hash_free (hash);
	array_free (array, 1);
}
END_TEST

START_TEST (test_hash_get_values_as_array)
{
	Hash *hash = hash_new (10, free, NULL);
	char *pronouns[] = { "I", "you", "he", "she", "it", "we", "they" };
	int num_elem = sizeof (pronouns)/sizeof (char *);
	int i = 0;

	for (; i < num_elem; i++)
		{
			char *key;
			int len = xasprintf (&key, "key%d", i);
			hash_insert (hash, key, pronouns[i]);
		}

	Array *array = hash_get_values_as_array (hash);
	ck_assert_ptr_ne (array, NULL);
	ck_assert_int_eq (array_len (array), num_elem);

	int match = 0;
	int j = 0;

	for (i = 0; i < num_elem; i++)
		for (j = 0; j < num_elem; j++)
			if (array_get (array, i) == pronouns[j])
				match ++;

	ck_assert_int_eq (match, num_elem);

	hash_free (hash);
	array_free (array, 1);
}
END_TEST

START_TEST (test_hash_int)
{
	Hash *h = NULL;
	int i = 0;
	int *n = NULL;

	h = hash_new_full (HASH_SMALL_SIZE,
			int_hash, int_equal, NULL, xfree);

	for (i = 0; i < 100; i++)
		{
			n = xcalloc (1, sizeof (int));
			*n = i;
			hash_insert (h, n, n);
		}

	ck_assert_int_eq (hash_size (h), 100);

	for (i = 0; i < 100; i++)
		{
			n = hash_lookup (h, &i);
			ck_assert_int_eq (*n, i);
		}

	hash_free (h);
}
END_TEST

Suite *
make_hash_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_misc;

	s = suite_create ("Hash");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_hash_new);
	tcase_add_test (tc_core, test_hash_insert);
	tcase_add_test (tc_core, test_hash_lookup);
	tcase_add_test (tc_core, test_hash_replace);
	tcase_add_test (tc_core, test_hash_remove);
	tcase_add_test (tc_core, test_hash_contains);
	tcase_add_test (tc_core, test_hash_int);
	suite_add_tcase (s, tc_core);

	/* Misc test case */
	tc_misc = tcase_create ("Misc");
	tcase_add_test (tc_misc, test_hash_foreach);
	tcase_add_test (tc_misc, test_hash_iter);
	tcase_add_test (tc_misc, test_hash_get_keys_as_list);
	tcase_add_test (tc_misc, test_hash_get_values_as_list);
	tcase_add_test (tc_misc, test_hash_get_keys_as_array);
	tcase_add_test (tc_misc, test_hash_get_values_as_array);
	suite_add_tcase (s, tc_misc);

	return s;
}
