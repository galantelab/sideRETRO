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
#include "../src/process_sample.h"

START_TEST (test_sider)
{
	ck_assert_int_eq (1, 1);
}
END_TEST


Suite *
make_process_sample_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("ProcessSample");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_sider);
	suite_add_tcase (s, tc_core);

	return s;
}
