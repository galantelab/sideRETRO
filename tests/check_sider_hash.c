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
			xasprintf (&str[i], "loop%d", i);
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

	free (str);
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
			xasprintf (&keys[i], "loop%d", i);
			hash_insert (hash, keys[i], &values[i]);
		}

	for (i = 0; i < 100; i++)
		{
			ck_assert_int_eq (hash_contains (hash, keys[i]), 1);
			int *ptr = hash_lookup (hash, keys[i]);
			ck_assert_int_eq (*ptr, values[i]);
		}

	hash_free (hash);
	free (keys);
	free (values);
}
END_TEST

Suite *
make_hash_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Hash");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_hash_new);
	tcase_add_test (tc_core, test_hash_insert);
	tcase_add_test (tc_core, test_hash_lookup);
	tcase_add_test (tc_core, test_hash_replace);
	tcase_add_test (tc_core, test_hash_remove);
	tcase_add_test (tc_core, test_hash_contains);

	return s;
}
