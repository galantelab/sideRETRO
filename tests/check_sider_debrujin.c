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
#include "../src/log.h"
#include "../src/wrapper.h"
#include "../src/debrujin.c"

START_TEST (test_debrujin_new)
{
	DeBrujin *debrujin = debrujin_new (3);
	debrujin_free (debrujin);
	debrujin_free (NULL);
}
END_TEST

START_TEST (test_debrujin_insert)
{
	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;
	ListElmt *cur_set = NULL;
	DeBrujinVetex *vertex = NULL;

	DeBrujin *debrujin = debrujin_new (3);

	debrujin_insert (debrujin, "AAGACTC");
	debrujin_insert (debrujin, "ACTCCGACTG");
	debrujin_insert (debrujin, "ACTGGGAC");
	debrujin_insert (debrujin, "GGACTTT");

	cur = list_head (graph_adjlists (debrujin->graph));
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			vertex = adjlist->vertex;
			log_debug ("k_mer: %s", vertex->k_mer);

			cur_set = list_head (set_list (adjlist->adjacent));
			for (; cur_set != NULL; cur_set = list_next (cur_set))
				{
					vertex = list_data (cur_set);
					log_debug ("  adjacent: %s", vertex->k_mer);
				}
		}

	debrujin_free (debrujin);
}
END_TEST

Suite *
make_debrujin_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("DeBrujin");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_debrujin_new);
	tcase_add_test (tc_core, test_debrujin_insert);
	suite_add_tcase (s, tc_core);

	return s;
}
