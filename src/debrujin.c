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

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "debrujin.h"

struct _DeBrujin
{
	size_t       k;
	Graph       *graph;
};

struct _DeBrujinVetex
{
	const char  *k_mer;
	VertexColor  color;
};

typedef struct _DeBrujinVetex DeBrujinVetex;

static void
debrujin_vertex_free (DeBrujinVetex *vertex)
{
	xfree ((char *) vertex->k_mer);
	xfree (vertex);
}

static uint32_t
debrujin_hash (const DeBrujinVetex *vertex)
{
	return str_hash (vertex->k_mer);
}

static int
debrujin_equal (const DeBrujinVetex *vertex1,
		const DeBrujinVetex *vertex2)
{
	return str_equal (vertex1->k_mer, vertex2->k_mer);
}

DeBrujin *
debrujin_new (size_t k)
{
	assert (k > 0);

	DeBrujin *debrujin = xcalloc (1, sizeof (DeBrujin));
	debrujin->graph = graph_new_full ((HashFunc) debrujin_hash,
			(EqualFun) debrujin_equal, (DestroyNotify) debrujin_vertex_free);
	debrujin->k = k;

	return debrujin;
}

void
debrujin_free (DeBrujin *debrujin)
{
	if (debrujin == NULL)
		return;

	graph_free (debrujin->graph);
	xfree (debrujin);
}

static inline const DeBrujinVetex *
debrujin_get_vertex (DeBrujin *debrujin, const char *k_mer)
{
	DeBrujinVetex temp = { .k_mer = k_mer };
	AdjList *adjlist = graph_adjlist (debrujin->graph, &temp);
	return adjlist != NULL ? adjlist->vertex : NULL;
}

static const DeBrujinVetex *
debrujin_insert_k_mer (DeBrujin *debrujin, const DeBrujinVetex *vertex,
		const char *k_mer)
{
	const DeBrujinVetex *cur_vertex = NULL;
	DeBrujinVetex *temp = NULL;

	cur_vertex = debrujin_get_vertex (debrujin, k_mer);

	if (cur_vertex == NULL)
		{
			temp = xcalloc (1, sizeof (DeBrujinVetex));
			temp->k_mer = xstrdup (k_mer);
			graph_ins_vertex (debrujin->graph, temp);
			cur_vertex = temp;
		}

	if (vertex != NULL)
		graph_ins_edge (debrujin->graph, vertex, cur_vertex);

	return cur_vertex;
}

void
debrujin_insert (DeBrujin *debrujin, const char *seq)
{
	assert (debrujin != NULL && seq != NULL);

	const DeBrujinVetex *vertex = NULL;
	size_t len = 0;
	size_t i = 0;
	char *k_mer = NULL;

	len = strlen (seq);
	if (len <= debrujin->k)
		return;

	// [k] == '\0'
	k_mer = xcalloc (debrujin->k + 1, sizeof (char));

	// There are len - k + 1 k_mers into a sequence
	for (i = 0; i < (len - debrujin->k + 1); i++)
		{
			strncpy (k_mer, seq + i, debrujin->k);
			vertex = debrujin_insert_k_mer (debrujin,
					vertex, k_mer);
		}

	xfree (k_mer);
}
