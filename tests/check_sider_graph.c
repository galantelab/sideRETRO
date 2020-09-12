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
#include "../src/graph.h"

START_TEST (test_graph_new)
{
	Graph *g = graph_new_full (str_hash, str_equal, NULL);
	graph_free (g);
	graph_free (NULL);
}
END_TEST

START_TEST (test_graph_ins_vertex)
{
	int i = 0;
	int *i_alloc = NULL;

	Graph *g = graph_new_full (int_hash, int_equal, xfree);

	for (i = 0; i < 10; i++)
		{
			i_alloc = xcalloc (1, sizeof (int));
			*i_alloc = i;
			ck_assert_int_eq (graph_ins_vertex (g, i_alloc), 1);
			ck_assert_int_eq (graph_ins_vertex (g, i_alloc), 0);
		}

	ck_assert_int_eq (graph_vcount (g), 10);
	ck_assert_int_eq (graph_ecount (g), 0);

	graph_free (g);
}
END_TEST

START_TEST (test_graph_ins_edge)
{
	int i = 0;
	int j = 0;
	int edges[4][10] = {};
	int degrees[4] = {};
	int num_vertex = 4;
	int num_edge = 0;
	int *i_alloc = NULL;
	AdjList *a = NULL;

	Graph *g = graph_new_full (int_hash, int_equal, xfree);

	for (i = 0; i < num_vertex; i++)
		{
			i_alloc = xcalloc (1, sizeof (int));
			*i_alloc = i;
			graph_ins_vertex (g, i_alloc);
		}

	/*
	* 0--------x1
	* |     /   x
	* |   /     |
	* x x       |
	* 2--------x3
	*/

	num_edge = 5;
	edges[0][0] = 1;
	edges[0][1] = 2;
	edges[0][2] = -1;
	edges[1][0] = 2;
	edges[1][1] = -1;
	edges[2][0] = 3;
	edges[2][1] = -1;
	edges[3][0] = 1;
	edges[3][1] = -1;

	degrees[0] = 2;
	degrees[1] = -1;
	degrees[2] = -1;
	degrees[3] = 0;

	for (i = 0; i < num_vertex; i++)
		{
			for (j = 0; edges[i][j] != -1; j++)
				{
					log_debug ("[%d] --> [%d]", i, edges[i][j]);
					ck_assert_int_eq (graph_ins_edge (g, &i, &edges[i][j]), 1);
					ck_assert_int_eq (graph_ins_edge (g, &i, &edges[i][j]), 0);
				}
		}

	ck_assert_int_eq (graph_vcount (g), num_vertex);
	ck_assert_int_eq (graph_ecount (g), num_edge);

	for (i = 0; i < num_vertex; i++)
		{
			a = graph_adjlist (g, &i);

			ck_assert_int_eq (* (int *) (a->vertex), i);
			ck_assert_int_eq (degrees[i], graph_diff_degree (a));

			for (j = 0; edges[i][j] != -1; j++)
				ck_assert_int_eq (graph_is_adjacent (g, &i, &edges[i][j]), 1);
		}

	// Unknown vetex for adjacency
	i = 0; j = 10;
	ck_assert_int_eq (graph_ins_edge (g, &i, &j), 0);
	ck_assert_int_eq (graph_is_adjacent (g, &i, &j), 0);

	// Unknown vetex
	i = 10; j = 0;
	ck_assert_int_eq (graph_ins_edge (g, &i, &j), 0);
	ck_assert_int_eq (graph_is_adjacent (g, &i, &j), 0);

	graph_free (g);
}
END_TEST

START_TEST (test_graph_ins_multi_edge)
{
	int i = 0;
	int n = 3;
	char *k_mers[] = {"ATCG", "TCGG", "CGGA"};
	Graph *graph = graph_new_full (str_hash, str_equal, NULL);

	for (i = 0; i < n; i++)
		graph_ins_vertex (graph, k_mers[i]);

	graph_ins_multi_edge (graph, k_mers[1], k_mers[2]);
	graph_ins_multi_edge (graph, k_mers[0], k_mers[1]);
	graph_ins_multi_edge (graph, k_mers[0], k_mers[1]);
	ck_assert_int_eq (graph_ins_multi_edge (graph, k_mers[0], k_mers[1]), 1);

	ck_assert_int_eq (graph_ecount (graph), 4);

	graph_free (graph);
}
END_TEST

START_TEST (test_graph_rem_vertex)
{
	int i = 0;
	int *i_alloc = NULL;
	int *j_alloc = NULL;
	void *ii_alloc = NULL;

	Graph *g = graph_new_full (int_hash, int_equal, xfree);

	for (i = 0; i < 10; i++)
		{
			i_alloc = xcalloc (1, sizeof (int));
			*i_alloc = i;
			graph_ins_vertex (g, i_alloc);
		}

	for (i = 0; i < 10; i++)
		{
			ii_alloc = &i;

			ck_assert_int_eq (graph_rem_vertex (g, &ii_alloc), 1);
			ck_assert_int_eq (* (int *) ii_alloc, i);
			ck_assert_int_eq (graph_vcount (g), 9 - i);

			xfree (ii_alloc);
		}

	// Coverage issues
	i_alloc = xcalloc (1, sizeof (int));
	j_alloc = xcalloc (1, sizeof (int));
	*i_alloc = 0;
	*j_alloc = 1;

	ck_assert_int_eq (graph_rem_vertex (g, (void **) &i_alloc), 0);

	graph_ins_vertex (g, i_alloc);
	graph_ins_vertex (g, j_alloc);

	graph_ins_edge (g, i_alloc, j_alloc);
	ck_assert_int_eq (graph_rem_vertex (g, (void **) &i_alloc), 0);
	ck_assert_int_eq (graph_rem_vertex (g, (void **) &j_alloc), 0);

	graph_free (g);
}
END_TEST

START_TEST (test_graph_rem_edge)
{
	int i = 0;
	void *edge = NULL;
	int num_vertex = 3;

	char *vertice[] =
		{
			"vertex1",
			"vertex2",
			"vertex3"
		};

	char *edges[] =
		{
			"vertex2",
			"vertex3",
			"vertex1"
		};

	Graph *g = graph_new_full (str_hash, str_equal, NULL);

	edge = vertice[0];
	ck_assert_int_eq (graph_rem_edge (g, "ponga", &edge), 0);

	for (i = 0; i < num_vertex; i++)
		ck_assert_int_eq (graph_ins_vertex (g, vertice[i]), 1);

	ck_assert_int_eq (graph_rem_edge (g, vertice[0], &edge), 0);

	for (i = 0; i < num_vertex; i++)
		ck_assert_int_eq (graph_ins_edge (g, vertice[i], edges[i]), 1);

	for (i = 0; i < num_vertex; i++)
		{
			edge = edges[i];
			ck_assert_int_eq (graph_rem_edge (g, vertice[i], &edge), 1);
			ck_assert (edge == edges[i]);
			ck_assert_int_eq (graph_ecount (g), num_vertex - 1 - i);
		}

	graph_free (g);
}
END_TEST

Suite *
make_graph_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Graph");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_graph_new);
	tcase_add_test (tc_core, test_graph_ins_vertex);
	tcase_add_test (tc_core, test_graph_ins_edge);
	tcase_add_test (tc_core, test_graph_ins_multi_edge);
	tcase_add_test (tc_core, test_graph_rem_vertex);
	tcase_add_test (tc_core, test_graph_rem_edge);
	suite_add_tcase (s, tc_core);

	return s;
}
