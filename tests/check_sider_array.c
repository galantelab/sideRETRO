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

#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_sider.h"

#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/array.h"

static Array *arr = NULL;

static void
setup1 (void)
{
	arr = array_new (xfree);
}

static void
setup2 (void)
{
	arr = array_new (NULL);
}

static void
teardown (void)
{
	array_free (arr, 1);
}

START_TEST (test_array)
{
	int ponga_size = 5;

	char *ponga[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5"
	};

	for (int i = 0; i < ponga_size; i++)
		array_add (arr, xstrdup (ponga[i]));

	ck_assert_int_eq (arr->len, ponga_size);

	for (int i = 0; i < arr->len; i++)
		ck_assert_str_eq ((char *) array_get (arr, i), ponga[i]);
}
END_TEST

START_TEST (test_array_no_free_segment)
{
	Array *arr = array_new (xfree);
	int ponga_size = 5;

	char *ponga[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5"
	};

	for (int i = 0; i < ponga_size; i++)
		array_add (arr, xstrdup (ponga[i]));

	char **ponga2 = (char **)array_free (arr, 0);

	for (int i = 0; i < ponga_size; i++)
		ck_assert_str_eq (ponga2[i], ponga[i]);

	for (int i = 0; i < ponga_size; i++)
		xfree (ponga2[i]);

	xfree (ponga2);
}
END_TEST

START_TEST (test_array_uniq)
{
	int ponga_size = 15;

	char *ponga[] = {
		"Ponga1",
		"Ponga1",
		"Ponga5",
		"Ponga5",
		"Ponga5",
		"Ponga1",
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga4",
		"Ponga4",
		"Ponga4",
		"Ponga5",
		"Ponga5"
	};

	int uniq_size = 5;
	char *uniq[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5",
	};

	for (int i = 0; i < ponga_size; i++)
		array_add (arr, xstrdup (ponga[i]));

	array_uniq (arr, cmpstringp);
	ck_assert_int_eq (arr->len, uniq_size);

	for (int i = 0; i < arr->len; i++)
		ck_assert_str_eq ((char *) array_get (arr, i), uniq[i]);
}
END_TEST

START_TEST (test_array_find)
{
	int ponga_size = 5;
	int rc = 0;
	int i = 0;
	int index_ = 0;

	char *ponga[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5"
	};

	for (i = 0; i < ponga_size; i++)
		array_add (arr, ponga[i]);

	for (i = 0; i < ponga_size; i++)
		{
			rc = array_find (arr, ponga[i], &index_);
			ck_assert_int_eq (rc, 1);
			ck_assert_int_eq (index_, i);
		}

	rc = array_find (arr, "ponga66", &index_);
	ck_assert_int_eq (rc, 0);
}
END_TEST

START_TEST (test_array_find_with_equal_fun)
{
	int ponga_size = 5;
	int rc = 0;
	int i = 0;
	int index_ = 0;

	char *ponga[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5"
	};

	for (i = 0; i < ponga_size; i++)
		array_add (arr, xstrdup (ponga[i]));

	for (i = 0; i < ponga_size; i++)
		{
			rc = array_find_with_equal_fun (arr, ponga[i],
					equalstring, &index_);
			ck_assert_int_eq (rc, 1);
			ck_assert_int_eq (index_, i);
		}

	rc = array_find_with_equal_fun (arr, "ponga66",
			equalstring, &index_);
	ck_assert_int_eq (rc, 0);
}
END_TEST

START_TEST (test_array_remove)
{
	int ponga_size = 5;
	int rc = 0;
	int i = 0;

	char *ponga[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5"
	};

	for (i = 0; i < ponga_size; i++)
		array_add (arr, ponga[i]);

	rc = array_remove (arr, "Ponga1");
	ck_assert_int_eq (rc, 1);
	ck_assert_int_eq (array_len (arr), ponga_size - 1);

	rc = array_remove (arr, "Ponga3");
	ck_assert_int_eq (rc, 1);
	ck_assert_int_eq (array_len (arr), ponga_size - 2);

	rc = array_remove (arr, "Ponga5");
	ck_assert_int_eq (rc, 1);
	ck_assert_int_eq (array_len (arr), ponga_size - 3);

	ck_assert_str_eq (array_get (arr, 0), "Ponga2");
	ck_assert_str_eq (array_get (arr, 1), "Ponga4");
}
END_TEST

START_TEST (test_array_remove_index)
{
	int ponga_size = 5;
	int i = 0;

	char *ponga[] = {
		"Ponga1",
		"Ponga2",
		"Ponga3",
		"Ponga4",
		"Ponga5"
	};

	for (i = 0; i < ponga_size; i++)
		array_add (arr, ponga[i]);

	ck_assert (array_remove_index (arr, 10) == NULL);

	ck_assert_str_eq (array_remove_index (arr, 0), ponga[0]);
	ck_assert_int_eq (array_len (arr), ponga_size - 1);

	ck_assert_str_eq (array_remove_index (arr, 3), ponga[4]);
	ck_assert_int_eq (array_len (arr), ponga_size - 2);

	ck_assert_str_eq (array_remove_index (arr, 1), ponga[2]);
	ck_assert_int_eq (array_len (arr), ponga_size - 3);

	ck_assert_str_eq (array_remove_index (arr, 0), ponga[1]);
	ck_assert_int_eq (array_len (arr), ponga_size - 4);

	ck_assert_str_eq (array_remove_index (arr, 0), ponga[3]);
	ck_assert_int_eq (array_len (arr), ponga_size - 5);
}
END_TEST

Suite *
make_array_suite (void)
{
	Suite *s;
	TCase *tc_core_free;
	TCase *tc_core_null;
	TCase *tc_free;

	s = suite_create ("Array");

	/* Core free test case */
	tc_core_free = tcase_create ("CoreFree");
	tcase_add_checked_fixture (tc_core_free, setup1, teardown);

	tcase_add_test (tc_core_free, test_array);
	tcase_add_test (tc_core_free, test_array_uniq);
	tcase_add_test (tc_core_free, test_array_find_with_equal_fun);
	suite_add_tcase (s, tc_core_free);

	/* Core null test case */
	tc_core_null = tcase_create ("CoreNull");
	tcase_add_checked_fixture (tc_core_null, setup2, teardown);

	tcase_add_test (tc_core_null, test_array_find);
	tcase_add_test (tc_core_null, test_array_remove);
	tcase_add_test (tc_core_null, test_array_remove_index);
	suite_add_tcase (s, tc_core_null);

	/* Free test case */
	tc_free = tcase_create ("Free");

	tcase_add_test (tc_free, test_array_no_free_segment);
	suite_add_tcase (s, tc_free);

	return s;
}
