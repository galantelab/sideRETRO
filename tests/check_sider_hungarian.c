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
#include "../src/hungarian.c"

static double **
array_to_matrix (double *m, int rows, int cols)
{
  int i,j;
  double **r;

  r = xcalloc(rows, sizeof (double *));

  for(i = 0; i < rows; i++)
		{
			r[i] = calloc (cols, sizeof (double));
			for(j = 0; j < cols; j++)
				r[i][j] = m[i * cols + j];
		}

  return r;
}

START_TEST (test_hungarian_solve)
{
	int n = 4;
	int m = 4;
	int i, j;

	double r[4*4] = {
		82, 83, 69, 92,
		77, 37, 49, 92,
		11, 69, 5, 86,
		8, 9, 98, 23
	};

	int assign[4][4] = {
		{0, 0, 1, 0},
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1}
	};

	int **assign_res;

	double **rr = array_to_matrix (r, n, m);

	Hungarian *p = hungarian_new (rr, n, m,
			HUNGARIAN_MODE_MINIMIZE_COST);

	hungarian_solve (p);
	assign_res = hungarian_assignment (p);

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			ck_assert_int_eq (assign[i][j], assign_res[i][j]);

  for (i = 0; i < n; i++)
		xfree (rr[i]);

  xfree (rr);
	hungarian_free (p);
}
END_TEST

Suite *
make_hungarian_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Hungarian");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_hungarian_solve);
	suite_add_tcase (s, tc_core);

	return s;
}
