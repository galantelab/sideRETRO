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
#include "../src/dbscan.h"

void
get_detail (Point *p, void *user_data)
{
	int (*point_detail)[3] = user_data;
	int i = * (int *) p->data;
	point_detail[i][0] = p->label;
	point_detail[i][1] = p->id;
	point_detail[i][2] = p->neighbors;
}

START_TEST (test_dbscan)
{
	// To test
	int pos[][2] = {
		{1000, 1100},
		{1050, 1150},
		{1300, 1400},
		{2000, 2100},
		{2500, 2600},
		{2560, 2660}
	};

	// Get label, id and neighbors
	int point_detail[6][3];

	int n_pos = 6;
	int i = 0;
	int j = 0;
	int *k = NULL;
	int acm = 0;

	DBSCAN *db = dbscan_new (xfree);

	for (i = 0; i < n_pos; i++)
		{
			k = xcalloc (1, sizeof (int));
			*k = i;
			dbscan_insert_point (db, pos[i][0],
					pos[i][1], k);
		}

	// Test eps = 300 and min_pts = 5
	// No grouping
	acm = dbscan_cluster (db, 300, 5, get_detail,
			point_detail);
	ck_assert_int_eq (acm, 0);

	// Test eps = 300 and min_pts = 3
	// label, id and neighbors
	int n_t300_3 = 3;
	int t300_3[][3] = {
		{CORE, 1, 3},
		{CORE, 1, 3},
		{CORE, 1, 3},
	};

	acm = dbscan_cluster (db, 300, 3, get_detail,
			point_detail);
	ck_assert_int_eq (acm, 1);

	for (i = 0; i < n_t300_3; i++)
		for (j = 0; j < 3; j++)
			ck_assert_int_eq (point_detail[i][j],
					t300_3[i][j]);

	// Test eps = 500 and min_pts = 3
	// label, id and neighbors
	int n_t500_3 = 6;
	int t500_3[][3] = {
		{CORE,      1, 3},
		{CORE,      1, 3},
		{CORE,      1, 3},
		{REACHABLE, 2, 2},
		{CORE,      2, 3},
		{REACHABLE, 2, 2},
	};

	acm = dbscan_cluster (db, 500, 3, get_detail,
			point_detail);
	ck_assert_int_eq (acm, 2);

	for (i = 0; i < n_t500_3; i++)
		for (j = 0; j < 3; j++)
			ck_assert_int_eq (point_detail[i][j],
					t500_3[i][j]);

	dbscan_free (db);
}
END_TEST

Suite *
make_dbscan_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("DBSCAN");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_dbscan);
	suite_add_tcase (s, tc_core);

	return s;
}
