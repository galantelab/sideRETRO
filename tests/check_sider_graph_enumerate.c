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
#include "../src/hash.h"
#include "../src/graph.h"
#include "../src/graph_enumerate.h"

static const int k_mers_num = 42;

static const char *k_mers[42] = {
	"ATC", "TCG",
	"TTC", "TCG",
	"TCG", "CGG",
	"ACG", "CGG",
	"TCG", "CGT",
	"CGT", "GTA",
	"GTA", "TAG",
	"GTA", "TAG",
	"CGG", "GGA",
	"AGG", "GGA",
	"CGG", "GGG",
	"GGG", "GGG",
	"GGG", "GGT",
	"GGA", "GAT",
	"GGA", "GAC",
	"GAT", "ATT",
	"GAT", "ATG",
	"TAG", "AGA",
	"AGA", "GAA",
	"GAA", "AAA",
	"AAA", "AAA"
};

static const int edges_num = 21;

static const int edges[21][2] = {
	{0,  1},
	{2,  3},
	{4,  5},
	{6,  7},
	{8,  9},
	{10, 11},
	{12, 13},
	{14, 15},
	{16, 17},
	{18, 19},
	{20, 21},
	{22, 23},
	{24, 25},
	{26, 27},
	{28, 29},
	{30, 31},
	{32, 33},
	{34, 35},
	{36, 37},
	{38, 39},
	{40, 41}
};

static int enums_num = 10;

static struct TestGraphEnum {
	GraphLabel    label;
	int           sense;
	const char   *parent;
	const char   *child1;
	const char   *child2;
} enums[10] = {
	{GRAPH_LABEL_FORK,   1, "CGG", "GGA", "GGG"},
	{GRAPH_LABEL_TIP,    0, "CGG", "ACG", "TCG"},
	{GRAPH_LABEL_TIP,    1, "GGA", "GAC", "GAT"},
	{GRAPH_LABEL_TIP,    0, "GGA", "AGG", "CGG"},
	{GRAPH_LABEL_CYCLE,  1, "AAA", "AAA", NULL},
	{GRAPH_LABEL_FORK,   1, "TCG", "CGG", "GTA"},
	{GRAPH_LABEL_FORK,   0, "TCG", "ATC", "TTC"},
	{GRAPH_LABEL_BUBBLE, 1, "GTA", "TAG", NULL},
	{GRAPH_LABEL_FORK,   1, "GAT", "ATT", "ATG"},
	{GRAPH_LABEL_CYCLE,  1, "GGG", "GGG", NULL}
};

static int
graph_enum_exist (const GraphEnum *e)
{
	AdjList *p, *c1, *c2;
	int i = 0;

	for (; i < enums_num; i++)
		{
			p = e->parent;
			c1 = e->child1;
			c2 = e->child2;

			if (enums[i].label == e->label && enums[i].sense == e->sense
					&& str_equal (enums[i].parent, p->vertex)
					&& str_equal (enums[i].child1, c1->vertex)
					&& ((enums[i].child2 == NULL && c2 == NULL)
						|| (enums[i].child2 != NULL && c2 != NULL
							&& str_equal (enums[i].child2, c2->vertex))))
				{
					return 1;
				}
		}

	return 0;
}

static Graph *
build_graph (void)
{
	Graph *g = graph_new_full (str_hash, str_equal, NULL);
	int i = 0;

	for (i = 0; i < k_mers_num; i++)
		graph_ins_vertex (g, k_mers[i]);

	for (i = 0; i < edges_num; i++)
		graph_ins_multi_edge (g, k_mers[edges[i][0]],
				k_mers[edges[i][1]]);

	return g;
}

START_TEST (test_graph_enumerate)
{
	Graph *g = NULL;
	Unipath *u = NULL;
	List *e = NULL;
	ListElmt *cur = NULL;

	g = build_graph ();
	u = graph_unipath_new (g, str_hash, str_equal);

	e = graph_enumerate (u, str_equal);

	cur = list_head (e);
	for (; cur != NULL; cur = list_next (cur))
		{
			GraphEnum *ge = list_data (cur);

			log_debug ("Label: %d", ge->label);
			log_debug ("Sense: %c", ge->sense ? '+' : '-');

			AdjList *a = ge->parent;
			log_debug ("Parent: %s", (char *) a->vertex);

			a = ge->child1;
			log_debug ("Child1: %s", (char *) a->vertex);

			a = ge->child2;
			if (a != NULL)
				log_debug ("Child2: %s", (char *) a->vertex);

			ck_assert (graph_enum_exist (ge));
		}

	list_free (e);
	graph_free (u);
	graph_free (g);
}
END_TEST

Suite *
make_graph_enumerate_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Enumerate");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_graph_enumerate);
	suite_add_tcase (s, tc_core);

	return s;
}
