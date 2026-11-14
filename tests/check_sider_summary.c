/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br
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

#include <time.h>
#include <string.h>
#include <check.h>
#include "check_sider.h"

#include "../src/summary.h"

START_TEST (test_summary_count)
{
	char buf[SUMMARY_NUM_LEN] = {};

	ck_assert_str_eq (summary_count (0, buf), "0");
	ck_assert_str_eq (summary_count (999, buf), "999");
	ck_assert_str_eq (summary_count (1000, buf), "1,000");
	ck_assert_str_eq (summary_count (1234567, buf), "1,234,567");
	ck_assert_str_eq (summary_count (-1234567, buf), "-1,234,567");
}
END_TEST

START_TEST (test_summary_percent)
{
	char buf[SUMMARY_NUM_LEN] = {};

	ck_assert_str_eq (summary_percent (1, 4, buf), "25.0%");
	ck_assert_str_eq (summary_percent (0, 10, buf), "0.0%");

	// No total: No percentage
	ck_assert_str_eq (summary_percent (1, 0, buf), "-");
}
END_TEST

START_TEST (test_summary_elapsed)
{
	char buf[SUMMARY_NUM_LEN] = {};
	time_t now = time (NULL);

	ck_assert_str_eq (summary_elapsed (now - 59, buf), "59s");
	ck_assert_str_eq (summary_elapsed (now - 61, buf), "1m 01s");
	ck_assert_str_eq (summary_elapsed (now - 3723, buf), "1h 02m 03s");

	// A start in the future is not negative
	ck_assert_str_eq (summary_elapsed (now + 100, buf), "0s");
}
END_TEST

START_TEST (test_summary_block)
{
	Summary *s = summary_new ("ponga summary", time (NULL));

	summary_section (s, "Section");
	summary_line (s, "Key", "%d", 42);
	summary_text (s, "%s", "free text");

	// The title is centered at the frame
	ck_assert (strstr (s->msg->str, "== ponga summary ==") != NULL);
	ck_assert (strstr (s->msg->str, "\n Section\n") != NULL);
	ck_assert (strstr (s->msg->str, "   Key") != NULL);
	ck_assert (strstr (s->msg->str, " 42\n") != NULL);
	ck_assert (strstr (s->msg->str, "   free text\n") != NULL);

	summary_log (s);
	ck_assert (strstr (s->msg->str, "Elapsed time") != NULL);

	summary_free (s);
	summary_free (NULL);
}
END_TEST

Suite *
make_summary_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Summary");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_summary_count);
	tcase_add_test (tc_core, test_summary_percent);
	tcase_add_test (tc_core, test_summary_elapsed);
	tcase_add_test (tc_core, test_summary_block);
	suite_add_tcase (s, tc_core);

	return s;
}
