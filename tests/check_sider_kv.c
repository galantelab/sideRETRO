/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2023 Thiago L. A. Miller <tmiller@mochsl.org.br
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
#include "../src/kv.h"

static KV *kv = NULL;

static void
setup (void)
{
	char path[] = "/tmp/kv.db.XXXXXX";
	int fd;

	fd = xmkstemp (path);
	close (fd);

	kv = kv_new (path);
}

static void
teardown (void)
	{
	char *path = NULL;

	path = xstrdup (kv_path (kv));

	kv_free (kv);
	xunlink (path);

	xfree (path);
}

START_TEST (test_kv_path)
{
	ck_assert_ptr_nonnull (kv_path (kv));
}
END_TEST

START_TEST (test_kv_insert)
{
	int i = 0;
	int n = 1000;
	char *key = NULL;

	for (i = 0; i < n; i++)
		{
			xasprintf (&key, "key %i", i);

			kv_insert (kv, key, NULL);

			xfree (key);
			key = NULL;
		}
}
END_TEST

START_TEST (test_kv_count)
{
	int i = 0;
	int n = 10;
	char *key = NULL;

	for (i = 0; i < n; i++)
		{
			xasprintf (&key, "key %i", i);

			kv_insert (kv, key, key);
			ck_assert_int_eq (i + 1, kv_count (kv));

			xfree (key);
			key = NULL;
		}
}
END_TEST

START_TEST (test_kv_del_key)
{
	kv_insert (kv, "one key", "one value");
	ck_assert_int_eq (kv_contains (kv, "one key"), 1);

	kv_del_key (kv, "one key");
	ck_assert_int_eq (kv_contains (kv, "one key"), 0);
}
END_TEST

START_TEST (test_kv_get_value)
{
	int n = 7;
	int i = 0;

	char *all[] = {
		"I",  "you",  "he",  "she", "it",  "we",  "they",
		"my", "your", "his", "her", "its", "our", "their"
	};

	for (i = 0; i < n; i++)
		kv_insert (kv, all[i], all[i+n]);

	for (i = 0; i < n; i++)
		ck_assert_str_eq (kv_get_value (kv, all[i]), all[i+n]);
}
END_TEST

static void
foreach (const char *key, const char *value, void *user_data)
{
	char **data = user_data;
	int i = 0;

	for (i = 0; data[i] != NULL; i++)
		if (!strcmp (data[i], key))
			break;

	ck_assert_str_eq (data[i], key);
	ck_assert_str_eq (key, value);
}

START_TEST (test_kv_foreach)
{
	char *data[] = {
		"eu", "tu", "ele", "ela", "nos", "vos", "eles", "elas", NULL
	};

	int i = 0;

	for (i = 0; data[i] != NULL; i++)
		kv_insert (kv, data[i], NULL);

	kv_foreach (kv, foreach, data);
}
END_TEST

Suite *
make_kv_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("KV");

	/* Core test case */
	tc_core = tcase_create ("Core");
	tcase_add_checked_fixture (tc_core, setup, teardown);

	tcase_add_test (tc_core, test_kv_path);
	tcase_add_test (tc_core, test_kv_insert);
	tcase_add_test (tc_core, test_kv_count);
	tcase_add_test (tc_core, test_kv_del_key);
	tcase_add_test (tc_core, test_kv_get_value);
	tcase_add_test (tc_core, test_kv_foreach);
	suite_add_tcase (s, tc_core);

	return s;
}
