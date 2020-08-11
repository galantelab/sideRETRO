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

#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/list.h"
#include "../src/ibitree.h"

static IBiTree *tree = NULL;
static List *list = NULL;

static void
setup (void)
{
	tree = ibitree_new (xfree);
	list = list_new (NULL);
}

static void
teardown (void)
{
	ibitree_free (tree);
	list_free (list);
}

START_TEST (test_ibitree_insert)
{
	int i = 0;

	for (; i < 10; i++)
		ibitree_insert (tree, i, i + 1, xstrdup ("ponga"));

	ck_assert_int_eq (ibitree_size (tree), i);
}
END_TEST

static void
catcher (IBiTreeLookupData *ldata, void *user_data)
{
	List *list = user_data;
	list_append (list, ldata->data);
}

START_TEST (test_ibitree_lookup)
{
	ListElmt *cur = NULL;
	int acm = 0;
	int i = 0;

	ibitree_insert (tree, 1, 10, xstrdup ("1-10"));
	ibitree_insert (tree, 9, 12, xstrdup ("9-12"));
	ibitree_insert (tree, 5, 10, xstrdup ("5-10"));
	ibitree_insert (tree, 1, 5, xstrdup ("1-5"));
	ibitree_insert (tree, 1, 1, xstrdup ("1-1"));

	acm = ibitree_lookup (tree, 6, 20, -1, -1, 0, catcher, list);
	ck_assert_int_eq (acm, 3);
	ck_assert_int_eq (list_size (list), 3);

	cur = list_head (list);
	const char *overlap[] = {"1-10", "5-10", "9-12"};

	for (; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), overlap[i]);
}
END_TEST

START_TEST (test_ibitree_lookup_frac)
{
	ListElmt *cur = NULL;
	int acm = 0;
	int i = 0;

	ibitree_insert (tree, 1, 10, xstrdup ("1-10"));
	ibitree_insert (tree, 5, 15, xstrdup ("5-15"));
	ibitree_insert (tree, 10, 20, xstrdup ("10-20"));
	ibitree_insert (tree, 15, 25, xstrdup ("15-25"));

	acm = ibitree_lookup (tree, 1, 20, 1.0, 0.5, 0, catcher, list);
	ck_assert_int_eq (acm, 3);
	ck_assert_int_eq (list_size (list), 3);

	cur = list_head (list);
	const char *overlap[] = {"1-10", "5-15", "10-20"};

	for (; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), overlap[i]);
}
END_TEST

START_TEST (test_ibitree_lookup_either)
{
	ListElmt *cur = NULL;
	int acm = 0;
	int i = 0;

	ibitree_insert (tree, 1, 10, xstrdup ("1-10"));
	ibitree_insert (tree, 5, 15, xstrdup ("5-15"));
	ibitree_insert (tree, 10, 20, xstrdup ("10-20"));
	ibitree_insert (tree, 15, 25, xstrdup ("15-25"));

	acm = ibitree_lookup (tree, 1, 20, 0.5, 0.5, 1, catcher, list);
	ck_assert_int_eq (acm, 4);
	ck_assert_int_eq (list_size (list), 4);

	cur = list_head (list);
	const char *overlap[] = {"1-10", "5-15", "10-20", "15-25"};

	for (; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), overlap[i]);
}
END_TEST

Suite *
make_ibitree_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("IBiTree");

	/* Core test case */
	tc_core = tcase_create ("Core");
	tcase_add_checked_fixture (tc_core, setup, teardown);

	tcase_add_test (tc_core, test_ibitree_insert);
	tcase_add_test (tc_core, test_ibitree_lookup);
	tcase_add_test (tc_core, test_ibitree_lookup_frac);
	tcase_add_test (tc_core, test_ibitree_lookup_either);
	suite_add_tcase (s, tc_core);

	return s;
}
