#include "config.h"

#include <stdlib.h>
#include <check.h>
#include "check_sider.h"
#include "../src/wrapper.h"
#include "../src/hash.h"

static Hash *hash = NULL;

static void
setup (void)
{
	hash = hash_new (xfree, NULL);
}

static void
teardown (void)
{
	hash_free (hash);
}

START_TEST (test_hash_free)
{
	hash_free (NULL);
}
END_TEST

START_TEST (test_hash_insert)
{
	hash_insert (hash, xstrdup ("Thiago"), "eu");
	hash_insert (hash, xstrdup ("Fernanda"), "tu");
	hash_insert (hash, xstrdup ("Pedro"), "ele");

	ck_assert_int_eq (hash_size (hash), 3);
}
END_TEST

START_TEST (test_hash_lookup)
{
	hash_insert (hash, xstrdup ("Thiago"), "eu");
	hash_insert (hash, xstrdup ("Fernanda"), "tu");
	hash_insert (hash, xstrdup ("Pedro"), "ele");

	ck_assert_str_eq (hash_lookup (hash, "Thiago"), "eu");
	ck_assert_str_eq (hash_lookup (hash, "Fernanda"), "tu");
	ck_assert_str_eq (hash_lookup (hash, "Pedro"), "ele");
}
END_TEST

START_TEST (test_hash_replace)
{
	int *num_alloc = xcalloc (1, sizeof (int));
	int num;
	int rc;

	rc = hash_insert (hash, xstrdup ("Thiago"), "eu");
	ck_assert_int_eq (rc, 1);

	*num_alloc = 120;

	rc = hash_insert (hash, xstrdup ("Thiago"), num_alloc);
	ck_assert_int_eq (rc, 0);

	num = * (int *) hash_lookup (hash, "Thiago");
	ck_assert_int_eq (num, *num_alloc);
	ck_assert_int_eq (num, 120);

	xfree (num_alloc);
}
END_TEST

START_TEST (test_hash_remove)
{
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
}
END_TEST

START_TEST (test_hash_contains)
{
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

	xfree (keys);
	xfree (values);
}
END_TEST

START_TEST (test_hash_iter)
{
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

	xfree (keys);
	xfree (values);
}
END_TEST

START_TEST (test_hash_get_keys_as_list)
{
	char *pronouns[] = { "I", "you", "he", "she", "it", "we", "they" };
	int num_elem = sizeof (pronouns)/sizeof (char *);
	int i = 0;

	for (; i < num_elem; i++)
		hash_insert (hash, xstrdup (pronouns[i]), NULL);

	List *list = hash_get_keys_as_list (hash);
	ck_assert_ptr_ne (list, NULL);
	ck_assert_int_eq (list_size (list), num_elem);

	ListElmt *cur = list_head (list);
	int match = 0;

	for (; cur != NULL; cur = list_next (cur))
		for (i = 0; i < num_elem; i++)
			if (!strcmp (list_data (cur), pronouns[i]))
				match ++;

	ck_assert_int_eq (match, num_elem);

	list_free (list);
}
END_TEST

START_TEST (test_hash_get_values_as_list)
{
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
			if (!strcmp (list_data (cur), pronouns[i]))
				match ++;

	ck_assert_int_eq (match, num_elem);

	list_free (list);
}
END_TEST

START_TEST (test_hash_get_keys_as_array)
{
	char *pronouns[] = { "I", "you", "he", "she", "it", "we", "they" };
	int num_elem = sizeof (pronouns)/sizeof (char *);
	int i = 0;

	for (; i < num_elem; i++)
		hash_insert (hash, xstrdup (pronouns[i]), NULL);

	Array *array = hash_get_keys_as_array (hash);
	ck_assert_ptr_ne (array, NULL);
	ck_assert_int_eq (array_len (array), num_elem);

	int match = 0;
	int j = 0;

	for (i = 0; i < num_elem; i++)
		for (j = 0; j < num_elem; j++)
			if (!strcmp (array_get (array, i), pronouns[j]))
				match ++;

	ck_assert_int_eq (match, num_elem);

	array_free (array, 1);
}
END_TEST

START_TEST (test_hash_get_values_as_array)
{
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
			if (!strcmp (array_get (array, i), pronouns[j]))
				match ++;

	ck_assert_int_eq (match, num_elem);

	array_free (array, 1);
}
END_TEST

START_TEST (test_hash_int)
{
	Hash *h = NULL;
	int i = 0;
	int *n = NULL;

	h = hash_new_full (int_hash, int_equal, NULL, xfree);

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

START_TEST (test_hash_direct)
{
	Hash *h = NULL;
	int k[100] = {};
	int i = 0;
	int *n = NULL;

	h = hash_new_full (direct_hash, direct_equal, NULL, NULL);

	for (i = 0; i < 100; i++)
		{
			k[i] = i;
			hash_insert (h, &k[i], &k[i]);
		}

	ck_assert_int_eq (hash_size (h), 100);

	for (i = 0; i < 100; i++)
		{
			n = hash_lookup (h, &k[i]);
			ck_assert_int_eq (*n, i);
		}

	hash_free (h);
}
END_TEST

START_TEST (test_hash_many_insertions)
{
	char *key = NULL;
	int size = 1e4;
	int i = 0;

	for (; i < size; i++)
		{
			xasprintf (&key, "ponga%d", i);
			hash_insert (hash, key, key);
		}

	ck_assert_int_eq (hash_size (hash), size);
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
	tcase_add_checked_fixture (tc_core, setup, teardown);

	tcase_add_test (tc_core, test_hash_free);
	tcase_add_test (tc_core, test_hash_insert);
	tcase_add_test (tc_core, test_hash_lookup);
	tcase_add_test (tc_core, test_hash_replace);
	tcase_add_test (tc_core, test_hash_remove);
	tcase_add_test (tc_core, test_hash_contains);
	tcase_add_test (tc_core, test_hash_int);
	tcase_add_test (tc_core, test_hash_direct);
	tcase_add_test (tc_core, test_hash_many_insertions);
	suite_add_tcase (s, tc_core);

	/* Misc test case */
	tc_misc = tcase_create ("Misc");
	tcase_add_checked_fixture (tc_misc, setup, teardown);

	tcase_add_test (tc_misc, test_hash_foreach);
	tcase_add_test (tc_misc, test_hash_iter);
	tcase_add_test (tc_misc, test_hash_get_keys_as_list);
	tcase_add_test (tc_misc, test_hash_get_values_as_list);
	tcase_add_test (tc_misc, test_hash_get_keys_as_array);
	tcase_add_test (tc_misc, test_hash_get_values_as_array);
	suite_add_tcase (s, tc_misc);

	return s;
}
