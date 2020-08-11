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

#include "../src/correlation.h"

START_TEST (test_pearson)
{
	double x[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}; 
	double y[10] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};

	double r = pearson (x, y, 10);
	ck_assert (r >= 0.95);
}
END_TEST

START_TEST (test_spearman)
{
	double x[10] = {86, 97, 99, 100, 101, 103, 106, 110, 112, 113};
	double y[10] = {0, 20, 28, 27, 50, 29, 7, 17, 6, 12};
	double work1[20];

	double rho = spearman (x, y, 10, work1);
	ck_assert (rho < 0);
}
END_TEST

START_TEST (test_spearman_eq)
{
	double x[10] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
	double y[10] = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
	double work1[20];

	double rho = spearman (x, y, 10, work1);
	ck_assert (rho < 0);
}
END_TEST

START_TEST (test_spearman_permutation_test)
{
	double x[10] = {86, 97, 99, 100, 101, 103, 106, 110, 112, 113};
	double y[10] = {0, 20, 28, 27, 50, 29, 7, 17, 6, 12};
	double work1[20];
	double work2[20];
	unsigned int seed = 1;

	double rho = spearman (x, y, 10, work1);
	double p_value = spearman_permutation_test (x, y, 10,
			work1, work2, &seed, rho);

	ck_assert (p_value > 0.5);
}
END_TEST

Suite *
make_correlation_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Correlation");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_pearson);
	tcase_add_test (tc_core, test_spearman);
	tcase_add_test (tc_core, test_spearman_eq);
	tcase_add_test (tc_core, test_spearman_permutation_test);
	suite_add_tcase (s, tc_core);

	return s;
}
