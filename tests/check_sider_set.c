#include "config.h"

#include <check.h>
#include "check_sider.h"
#include "../src/wrapper.h"
#include "../src/set.h"

Set *set = NULL;

void
setup (void)
{
	set = set_new (xfree);
}

void
teardown (void)
{
	set_free (set);
}

START_TEST (test_set_free)
{
	Set *set2 = set_new_full (str_hash, str_equal, NULL);
	set_free (NULL);
	set_free (set2);
}
END_TEST

START_TEST (test_set_insert)
{
	int rc = 0;
	char *k = NULL;

	rc = set_insert (set, xstrdup ("ponga1"));
	ck_assert_int_eq (rc, 1);

	rc = set_insert (set, xstrdup ("ponga2"));
	ck_assert_int_eq (rc, 1);

	rc = set_insert (set, xstrdup ("ponga3"));
	ck_assert_int_eq (rc, 1);

	k = xstrdup ("ponga3");

	rc = set_insert (set, k);
	ck_assert_int_eq (rc, 1);

	rc = set_insert (set, k);
	ck_assert_int_eq (rc, 0);
}
END_TEST

START_TEST (test_set_remove)
{
	int i = 0;
	char *k[5];
	int size = sizeof (k) / sizeof (char *);

	for (i = 0; i < size; i++)
		{
			xasprintf (&k[i], "ponga%d", i);
			set_insert (set, k[i]);
		}

	ck_assert_int_eq (set_size (set), size);

	for (i = 0; i < size; i++)
		ck_assert_int_eq (set_remove (set, (void **) &k[i]), 1);

	ck_assert_int_eq (set_size (set), 0);

	for (i = 0; i < size; i++)
		xfree (k[i]);
}
END_TEST

START_TEST (test_set_union)
{
	Set *set1, *set2, *setu;

	set1 = set_new (NULL);
	set_insert (set1, "ponga0");
	set_insert (set1, "ponga1");

	set2 = set_new (NULL);
	set_insert (set2, "ponga1");
	set_insert (set2, "ponga2");

	setu = set_union (set1, set2);
	ck_assert_int_eq (set_size (setu), 3);

	List *l = set_list (setu);
	ListElmt *cur = list_head (l);

	ck_assert_str_eq (list_data (cur), "ponga0");
	cur = list_next (cur);
	ck_assert_str_eq (list_data (cur), "ponga1");
	cur = list_next (cur);
	ck_assert_str_eq (list_data (cur), "ponga2");

	set_free (set1);
	set_free (set2);
	set_free (setu);
}
END_TEST

START_TEST (test_set_intersection)
{
	Set *set1, *set2, *seti;

	set1 = set_new (NULL);
	set_insert (set1, "ponga0");
	set_insert (set1, "ponga1");

	set2 = set_new (NULL);
	set_insert (set2, "ponga1");
	set_insert (set2, "ponga2");

	seti = set_intersection (set1, set2);
	ck_assert_int_eq (set_size (seti), 1);

	List *l = set_list (seti);
	ListElmt *cur = list_head (l);

	ck_assert_str_eq (list_data (cur), "ponga1");

	set_free (set1);
	set_free (set2);
	set_free (seti);
}
END_TEST

START_TEST (test_set_difference)
{
	Set *set1, *set2, *setd;

	set1 = set_new (NULL);
	set_insert (set1, "ponga0");
	set_insert (set1, "ponga1");

	set2 = set_new (NULL);
	set_insert (set2, "ponga1");
	set_insert (set2, "ponga2");

	setd = set_difference (set1, set2);
	ck_assert_int_eq (set_size (setd), 1);

	List *l = set_list (setd);
	ListElmt *cur = list_head (l);

	ck_assert_str_eq (list_data (cur), "ponga0");

	set_free (set1);
	set_free (set2);
	set_free (setd);
}
END_TEST

START_TEST (test_set_is_subset)
{
	Set *set1, *set2;
	char *k[5];
	int size = sizeof (k) / sizeof (char *);
	int i = 0;

	for (i = 0; i < size; i++)
		{
			xasprintf (&k[i], "ponga%d", i);
			set_insert (set, k[i]);
		}

	set1 = set_new (NULL);
	for (i = size - 2; i < size; i++)
		set_insert (set1, k[i]);

	set2 = set_new (NULL);
	set_insert (set2, "ponga1000");
	set_insert (set2, "ponga1001");

	ck_assert_int_eq (set_is_subset (set, set1), 0);
	ck_assert_int_eq (set_is_subset (set1, set), 1);
	ck_assert_int_eq (set_is_subset (set2, set), 0);

	set_free (set1);
	set_free (set2);
}
END_TEST

START_TEST (test_set_is_equal)
{
	Set *set1, *set2;
	char *k[5];
	int size = sizeof (k) / sizeof (char *);
	int i = 0;

	for (i = 0; i < size; i++)
		{
			xasprintf (&k[i], "ponga%d", i);
			set_insert (set, k[i]);
		}

	set1 = set_new (NULL);
	for (i = 0; i < size; i++)
		set_insert (set1, k[i]);

	set2 = set_new (NULL);
	set_insert (set2, "ponga1000");
	set_insert (set2, "ponga1001");

	ck_assert_int_eq (set_is_equal (set, set2), 0);
	ck_assert_int_eq (set_is_equal (set, set1), 1);

	set_free (set1);
	set_free (set2);
}
END_TEST

Suite *
make_set_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Set");

	/* Core test case */
	tc_core = tcase_create ("Core");
	tcase_add_checked_fixture (tc_core, setup, teardown);

	tcase_add_test (tc_core, test_set_free);
	tcase_add_test (tc_core, test_set_insert);
	tcase_add_test (tc_core, test_set_remove);
	tcase_add_test (tc_core, test_set_union);
	tcase_add_test (tc_core, test_set_intersection);
	tcase_add_test (tc_core, test_set_difference);
	tcase_add_test (tc_core, test_set_is_subset);
	tcase_add_test (tc_core, test_set_is_equal);
	suite_add_tcase (s, tc_core);

	return s;
}
