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
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/list.h"

static List *l = NULL;

static void
setup (void)
{
	l = list_new (xfree);
}

static void
teardown (void)
{
	list_free (l);
}

START_TEST (test_list_new)
{
	ck_assert_ptr_ne (l, NULL);
	ck_assert_int_eq (l->size, 0);
	ck_assert_ptr_eq (l->head, NULL);
	ck_assert_ptr_eq (l->tail, NULL);
}
END_TEST

START_TEST (test_list_ins_prev)
{
	list_ins_prev (l, NULL, xstrdup ("eu"));
	list_ins_prev (l, list_head (l), xstrdup ("tu"));
	list_ins_prev (l, list_head (l), xstrdup ("ele"));

	ListElmt *cur = list_head (l);
	ck_assert_str_eq (list_data (cur), "ele"); cur = list_next (cur);
	ck_assert_str_eq (list_data (cur), "tu"); cur = list_next (cur);
	ck_assert_str_eq (list_data (cur), "eu");
}
END_TEST

START_TEST (test_list_ins_next)
{
	list_ins_next (l, NULL, xstrdup ("eu"));
	list_ins_next (l, list_tail (l), xstrdup ("tu"));
	list_ins_next (l, list_tail (l), xstrdup ("ele"));

	ListElmt *cur = list_head (l);
	ck_assert_str_eq (list_data (cur), "eu"); cur = list_next (cur);
	ck_assert_str_eq (list_data (cur), "tu"); cur = list_next (cur);
	ck_assert_str_eq (list_data (cur), "ele");
}
END_TEST

START_TEST (test_list_ins_prev_link)
{
	list_ins_prev (l, NULL, xstrdup ("eu"));
	list_ins_prev (l, list_head (l), xstrdup ("tu"));
	list_ins_prev (l, list_head (l), xstrdup ("ele"));

	ListElmt *tail = list_tail (l);
	list_remove_link (l, tail);
	list_ins_prev_link (l, list_head (l), tail);

	ck_assert_str_eq (list_data (list_head (l)), "eu");
}
END_TEST

START_TEST (test_list_ins_next_link)
{
	list_ins_next (l, NULL, xstrdup ("eu"));
	list_ins_next (l, list_tail (l), xstrdup ("tu"));
	list_ins_next (l, list_tail (l), xstrdup ("ele"));

	ListElmt *head = list_head (l);
	list_remove_link (l, head);
	list_ins_next_link (l, list_tail (l), head);

	ck_assert_str_eq (list_data (list_tail (l)), "eu");
}
END_TEST

START_TEST (test_list_remove)
{
	list_ins_next (l, NULL, xstrdup ("eu"));
	list_ins_next (l, list_tail (l), xstrdup ("tu"));
	list_ins_next (l, list_tail (l), xstrdup ("ele"));

	ListElmt *head = list_head (l);
	void *data;

	list_remove (l, head, &data);
	ck_assert_int_eq (list_size (l), 2);
	ck_assert_str_eq ((char *) data, "eu");

	list_remove (l, list_head (l), NULL);
	ck_assert_int_eq (list_size (l), 1);
	ck_assert_str_eq (list_data (list_head (l)), "ele");

	xfree (data);
}
END_TEST

START_TEST (test_list_remove_link)
{
	list_ins_next (l, NULL, xstrdup ("eu"));
	list_ins_next (l, list_tail (l), xstrdup ("tu"));
	list_ins_next (l, list_tail (l), xstrdup ("ele"));

	ListElmt *head = list_head (l);
	list_remove_link (l, head);

	ck_assert_int_eq (list_size (l), 2);
	ck_assert_str_eq (list_data (head), "eu");
	ck_assert_str_eq (list_data (list_head (l)), "tu");

	xfree (list_data (head));
	xfree (head);
}
END_TEST

static void
sum (void *a, void *b)
{
	int *ai = a;
	int *bi = b;
	*bi += *ai;
}

START_TEST (test_list_foreach)
{
	int prime[] = {1, 3, 5, 7, 11, 13, 17};
	int total = 0;
	int i = 0;

	for (; i < sizeof (prime)/sizeof (int); i++)
		{
			int *p = xcalloc (1, sizeof (int));
			*p = prime[i];
			list_append (l, p);
		}

	list_foreach (l, sum, &total);
	ck_assert_int_eq (total, 57);
}
END_TEST

Suite *
make_list_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("List");

	/* Core test case */
	tc_core = tcase_create ("Core");
	tcase_add_checked_fixture (tc_core, setup, teardown);

	tcase_add_test (tc_core, test_list_new);
	tcase_add_test (tc_core, test_list_ins_prev);
	tcase_add_test (tc_core, test_list_ins_next);
	tcase_add_test (tc_core, test_list_ins_prev_link);
	tcase_add_test (tc_core, test_list_ins_next_link);
	tcase_add_test (tc_core, test_list_remove);
	tcase_add_test (tc_core, test_list_remove_link);
	tcase_add_test (tc_core, test_list_foreach);
	suite_add_tcase (s, tc_core);

	return s;
}
