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
	int          k;
	char        *k_mer_buff;
	Graph       *graph;
};

struct _DeBrujinVetex
{
	int          in_degree;
	int          out_degree;

	const char  *k_mer;

	ListElmt    *cur_adj;
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
debrujin_new (int k)
{
	assert (k > 0);

	DeBrujin *debrujin = xcalloc (1, sizeof (DeBrujin));

	debrujin->graph = graph_new_full ((HashFunc) debrujin_hash,
			(EqualFun) debrujin_equal,
			(DestroyNotify) debrujin_vertex_free);

	// [k - 1] == '\0'
	debrujin->k_mer_buff = xcalloc (k, sizeof (char));
	debrujin->k = k;

	return debrujin;
}

void
debrujin_free (DeBrujin *debrujin)
{
	if (debrujin == NULL)
		return;

	graph_free (debrujin->graph);
	xfree (debrujin->k_mer_buff);

	xfree (debrujin);
}

static DeBrujinVetex *
debrujin_insert_k_mer (DeBrujin *debrujin)
{
	DeBrujinVetex *vertex = NULL;
	AdjList *adjlist = NULL;

	DeBrujinVetex temp = { .k_mer = debrujin->k_mer_buff };
	adjlist = graph_adjlist (debrujin->graph, &temp);

	if (adjlist == NULL)
		{
			vertex = xcalloc (1, sizeof (DeBrujinVetex));
			vertex->k_mer = xstrdup (debrujin->k_mer_buff);
			graph_ins_vertex (debrujin->graph, vertex);
		}
	else
		vertex = adjlist->vertex;

	return vertex;
}

void
debrujin_insert (DeBrujin *debrujin, const char *seq)
{
	assert (debrujin != NULL && seq != NULL);

	DeBrujinVetex *pref = NULL;
	DeBrujinVetex *suff = NULL;
	int len = 0;
	int i = 0;

	len = strlen (seq);
	if (len < debrujin->k)
		return;

	// There are len - k + 1 k_mers into a sequence
	for (i = 0; i < (len - debrujin->k + 1); i++)
		{
			// Set prefix k_mer
			strncpy (debrujin->k_mer_buff, seq + i,
					debrujin->k - 1);

			pref = debrujin_insert_k_mer (debrujin);

			// Set suffix k_mer
			strncpy (debrujin->k_mer_buff, seq + i + 1,
					debrujin->k - 1);

			suff = debrujin_insert_k_mer (debrujin);

			// Connect path for this k_mer
			if (graph_ins_edge (debrujin->graph, pref, suff))
				{
					pref->out_degree++;
					suff->in_degree++;
				}
		}
}

static void
debrujin_reset_path (DeBrujin *debrujin)
{
	AdjList *adjlist = NULL;
	DeBrujinVetex *vertex = NULL;
	GraphIter iter;

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			vertex = adjlist->vertex;
			vertex->cur_adj = list_head (adjlist->adjacent);
		}
}

static int
debrujin_has_eulerian_path (const DeBrujin *debrujin)
{
	AdjList *adjlist = NULL;
	DeBrujinVetex *vertex = NULL;
	GraphIter iter;
	int start_nodes = 0;
	int end_nodes = 0;

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			vertex = adjlist->vertex;
			if ((vertex->in_degree - vertex->out_degree) > 1
					|| (vertex->out_degree - vertex->in_degree) > 1)
				return 0;
			else if ((vertex->in_degree - vertex->out_degree) == 1)
				end_nodes++;
			else if ((vertex->out_degree - vertex->in_degree) == 1)
				start_nodes++;
		}

	return (end_nodes == 0 && start_nodes == 0)
		|| (end_nodes == 1 && start_nodes == 1);
}

AdjList *
debrujin_find_start_node (const DeBrujin *debrujin)
{
	AdjList *adjlist = NULL;
	AdjList *adj_start = NULL;
	DeBrujinVetex *vertex = NULL;
	GraphIter iter;

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			vertex = adjlist->vertex;
			if ((vertex->out_degree - vertex->in_degree) == 1)
				return adjlist;

			if (vertex->out_degree > 0)
				adj_start = adjlist;
		}

	return adj_start;
}

static void
debrujin_dfs (DeBrujin *debrujin, AdjList *adjlist, List *k_mers)
{
	AdjList *clr_adjlist = NULL;
	DeBrujinVetex *clr_vertex = NULL;
	DeBrujinVetex *vertex = NULL;

	vertex = adjlist->vertex;

	while (vertex->cur_adj != NULL)
		{
			// Determine the color of the next adjacent vertex
			clr_vertex = list_data (vertex->cur_adj);

			// Get adjacency list of the clr_vertex
			clr_adjlist = graph_adjlist (debrujin->graph, clr_vertex);

			// Move one vertex deeper
			vertex->cur_adj = list_next (vertex->cur_adj);
			debrujin_dfs (debrujin, clr_adjlist, k_mers);
		}

	// Sequence grows backwards
	list_prepend (k_mers, vertex->k_mer);
}

static const char *
debrujin_sequence (const DeBrujin *debrujin, const List *k_mers)
{
	ListElmt *cur = NULL;
	char *seq = NULL;
	int k = 0;
	int i = 0;

	cur = list_head (k_mers);
	if (cur == NULL)
		return NULL;

	k = debrujin->k - 1;

	// seq has length num_k-mer + k - 1
	seq = xcalloc (list_size (k_mers) + k, sizeof (char));

	strncpy (seq, list_data (cur), k);
	for (cur = list_next (cur); cur != NULL; cur = list_next (cur))
		seq[k + i++] = * ((char *) list_data (cur) + k - 1);

	return seq;
}

List *
debrujin_assembly (DeBrujin *debrujin)
{
	assert (debrujin != NULL);

	List *seqs = NULL;
	List *k_mers = NULL;
	const char *seq = NULL;

	if (!debrujin_has_eulerian_path (debrujin))
		return NULL;

	// Reset the vetice color
	debrujin_reset_path (debrujin);

	seqs = list_new (xfree);
	k_mers = list_new (NULL);

	debrujin_dfs (debrujin,
			debrujin_find_start_node (debrujin),
			k_mers);

	seq = debrujin_sequence (debrujin, k_mers);

	if (seq != NULL)
		list_append (seqs, seq);

	list_free (k_mers);
	return seqs;
}
