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
#include "../src/graph_unipath.c"

static int num_vertices = 6;
static char *vertices[6] = {
	"TAC", "ACA", "CAC", "ACT", "CTC", "GAC"
};

static int num_edges = 5;
static int edges[5][2] = {
	{0, 1}, {1, 2}, {0, 3}, {3, 4}, {5, 1}
};

static int num_uni_vertices = 3;
static char *uni_vertices[3] = {
	"GAC", "ACA", "TAC"
};

static int num_uni_edges = 2;
static int uni_edges[2][2] = {
	{0, 1}, {2, 1}
};

START_TEST (test_graph_unipath_new)
{
	int i = 0;
	Graph *u = NULL;
	Graph *g = NULL;

	g = graph_new_full (str_hash, str_equal, NULL);

	for (i = 0; i < num_vertices; i++)
		graph_ins_vertex (g, vertices[i]);

	for (i = 0; i < num_edges; i++)
		graph_ins_edge (g, vertices[edges[i][0]],
				vertices[edges[i][1]]);

	u = graph_unipath_new (g, str_hash, str_equal);
	ck_assert (u != NULL);

	GraphIter ii = {};
	AdjList *aa;
	ListElmt *cur;
	graph_iter_init (&ii, u);
	while (graph_iter_next (&ii, &aa))
		{
			log_debug (":: %s", (char *) aa->vertex);
			cur = list_head (aa->adjacent);
			for (; cur != NULL; cur = list_next (cur))
				{
					log_debug ("  -> %s", (char *) list_data (cur));
				}
		}

	ck_assert_int_eq (graph_vcount (u), num_uni_vertices);
	ck_assert_int_eq (graph_ecount (u), num_uni_edges);

	for (i = 0; i < num_uni_edges; i++)
		ck_assert_int_eq (graph_is_adjacent (u,
					uni_vertices[uni_edges[i][0]],
					uni_vertices[uni_edges[i][1]]),
				1);

	graph_free (g);
	graph_free (u);
}
END_TEST

Suite *
make_graph_unipath_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Unipath");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_graph_unipath_new);
	suite_add_tcase (s, tc_core);

	return s;
}
