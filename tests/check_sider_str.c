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

#include <string.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/str.h"

START_TEST (test_string_new)
{
	const char *str = "PONGA";
	String *s = string_new (str);

	ck_assert_int_eq (s->len, strlen (str));
	ck_assert_int_eq (s->alloc, strlen (str) + 1);
	ck_assert_str_eq (s->str, str);

	char *s_str = string_free (s, 1);
	ck_assert (s_str == NULL);;
}
END_TEST

START_TEST (test_string_sized_new)
{
	size_t size = 1024;
	String *s = string_sized_new (size);

	ck_assert_int_eq (s->alloc, size);
	ck_assert_int_eq (s->len, 0);

	char *s_str = string_free (s, 1);
	ck_assert (s_str == NULL);;
}
END_TEST

START_TEST (test_string_set)
{
	const char *str = "PONGA";
	String *s = string_new (NULL);
	ck_assert_int_eq (s->alloc, 1);
	ck_assert_int_eq (s->len, 0);

	s = string_set (s, str);
	ck_assert_int_eq (s->len, strlen (str));
	ck_assert_int_eq (s->alloc, strlen (str) + 1);
	ck_assert_str_eq (s->str, str);

	char *s_str = string_free (s, 1);
	ck_assert (s_str == NULL);;
}
END_TEST

START_TEST (test_string_concat)
{
	const char *str = "PONGA1_PONGA2";
	const char *str1 = "PONGA1";
	const char *sep = "_";
	const char *str2 = "PONGA2";

	String *s = string_new (str1);
	s = string_concat (s, sep);
	s = string_concat (s, str2);

	ck_assert_int_eq (s->len, strlen (str));
	ck_assert_int_eq (s->alloc, strlen (str) + 1);
	ck_assert_str_eq (s->str, str);

	char *s_str = string_free (s, 1);
	ck_assert (s_str == NULL);;
}
END_TEST

START_TEST (test_string_printf)
{
	const char *str = "PONGA1_PONGA2";
	const char *str1 = "PONGA1";
	const char *str2 = "PONGA2";

	String *s = string_new (NULL);
	s = string_printf (s, "%s_%s", str1, str2);

	ck_assert_int_eq (s->len, strlen (str));
	ck_assert_int_eq (s->alloc, strlen (str) + 1);
	ck_assert_str_eq (s->str, str);

	char *s_str = string_free (s, 1);
	ck_assert (s_str == NULL);;
}
END_TEST

START_TEST (test_string_concat_printf)
{
	const char *str = "PONGA1_PONGA2";
	const char *str1 = "PONGA1";
	const char *str2 = "PONGA2";

	String *s = string_new (str1);
	s = string_concat_printf (s, "_%s", str2);

	ck_assert_int_eq (s->len, strlen (str));
	ck_assert_int_eq (s->alloc, strlen (str) + 1);
	ck_assert_str_eq (s->str, str);

	char *s_str = string_free (s, 1);
	ck_assert (s_str == NULL);;
}
END_TEST

START_TEST (test_string_free)
{
	const char *str = "PONGA";
	String *s = string_new (str);

	char *s_str = string_free (s, 0);
	ck_assert_str_eq (str, s_str);

	xfree (s_str);
}
END_TEST

START_TEST (test_string_clear)
{
	const char *str = "PONGA";
	String *s = string_new (str);

	s = string_clear (s);

	ck_assert_int_eq (s->len, 0);
	ck_assert_int_eq (s->str[0], 0);

	string_free (s, 1);
}
END_TEST

Suite *
make_str_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("String");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_string_new);
	tcase_add_test (tc_core, test_string_sized_new);
	tcase_add_test (tc_core, test_string_set);
	tcase_add_test (tc_core, test_string_concat);
	tcase_add_test (tc_core, test_string_printf);
	tcase_add_test (tc_core, test_string_concat_printf);
	tcase_add_test (tc_core, test_string_free);
	tcase_add_test (tc_core, test_string_clear);
	suite_add_tcase (s, tc_core);

	return s;
}
