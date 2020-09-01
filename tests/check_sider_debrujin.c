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
	DeBrujinVetex *vertex = NULL;

	DeBrujin *debrujin = debrujin_new (4);

	debrujin_insert (debrujin, "AAGACTC");
	debrujin_insert (debrujin, "ACTCCGACTG");
	debrujin_insert (debrujin, "ACTGGGAC");
	debrujin_insert (debrujin, "GGACTTT");

	GraphIter giter;
	graph_iter_init (&giter, debrujin->graph);

	while (graph_iter_next (&giter, &adjlist))
		{
			vertex = adjlist->vertex;
			log_debug ("k_mer (%d) [%d,%d]: %s",
					vertex->depth, vertex->in_degree,
					vertex->out_degree, vertex->k_mer_affix);

			cur = list_head (adjlist->adjacent);
			for (; cur != NULL; cur = list_next (cur))
				{
					vertex = list_data (cur);
					log_debug ("  adjacent: %s",
							vertex->k_mer_affix);
				}
		}

	List *seqs = debrujin_assembly (debrujin);

	if (seqs != NULL)
		{
			ListElmt *cur = list_head (seqs);
			for (; cur != NULL; cur = list_next (cur))
				log_debug (":: PONGA: %s", (char *) list_data (cur));
		}

	list_free (seqs);
	debrujin_free (debrujin);
}
END_TEST

START_TEST (test_debrujin_assembly)
{
	List *seqs = NULL;
	ListElmt *cur = NULL;
	const char *seq = "AAGACTCCGACTGGGACTTT";

	DeBrujin *debrujin = debrujin_new (5);

	debrujin_insert (debrujin, "AAGACTC");
	debrujin_insert (debrujin, "ACTCCGACTG");
	debrujin_insert (debrujin, "ACTGGGAC");
	debrujin_insert (debrujin, "GGACTTT");

	seqs = debrujin_assembly (debrujin);

	log_debug (":: SEQ: %s", seq);

	if (seqs != NULL)
		{
			cur = list_head (seqs);
			for (; cur != NULL; cur = list_next (cur))
				log_debug (":: CONTIG: %s", (char *) list_data (cur));
		}

	list_free (seqs);
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
	tcase_add_test (tc_core, test_debrujin_assembly);
	suite_add_tcase (s, tc_core);

	return s;
}
