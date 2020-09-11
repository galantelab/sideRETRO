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
#include "hash.h"
#include "graph.h"
#include "log.h"
#include "debrujin.h"

struct _DeBrujin
{
	int                    k;
	int                    min_depth;

	char                  *buf;

	Hash                  *k_mers;
	Graph                 *graph;
};

struct _DeBrujinVetex
{
	// Degree
	int                    in_degree;
	int                    out_degree;

	// Unigraph
	int                    unipath;

	// Edge depth
	int                    depth;

	// The affix
	const char            *k_mer_affix;

	// For BFS/DFS
	VertexColor            color;
	int                    dist;
	struct _DeBrujinVetex *parent;

	// For DFS/Hierholzer
	ListElmt              *cur_adj;
};

typedef struct _DeBrujinVetex DeBrujinVetex;

static void
debrujin_vertex_free (DeBrujinVetex *vertex)
{
	xfree ((char *) vertex->k_mer_affix);
	xfree (vertex);
}

static uint32_t
debrujin_hash (const DeBrujinVetex *vertex)
{
	return str_hash (vertex->k_mer_affix);
}

static int
debrujin_equal (const DeBrujinVetex *vertex1,
		const DeBrujinVetex *vertex2)
{
	return str_equal (vertex1->k_mer_affix,
			vertex2->k_mer_affix);
}

DeBrujin *
debrujin_new (int k)
{
	assert (k > 0);

	DeBrujin *debrujin = xcalloc (1, sizeof (DeBrujin));

	debrujin->k_mers = hash_new_full (str_hash, str_equal,
			xfree, xfree);

	debrujin->graph = graph_new_full ((HashFunc) debrujin_hash,
			(EqualFun) debrujin_equal,
			(DestroyNotify) debrujin_vertex_free);

	// [k] == '\0'
	debrujin->buf = xcalloc (k + 1, sizeof (char));
	debrujin->k = k;

	return debrujin;
}

void
debrujin_free (DeBrujin *debrujin)
{
	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;
	GraphIter iter = {};

	if (debrujin == NULL)
		return;

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			cur = list_head (adjlist->adjacent);
			for (; cur != NULL; cur = list_next (cur))
				xfree (list_data (cur));
		}

	xfree (debrujin->buf);
	hash_free (debrujin->k_mers);
	graph_free (debrujin->graph);

	xfree (debrujin);
}

static inline AdjList *
debrujin_adjlist (const DeBrujin *debrujin,
		const char *k_mer_affix)
{
	DeBrujinVetex temp = { .k_mer_affix = k_mer_affix };
	return graph_adjlist (debrujin->graph, &temp);
}

static inline DeBrujinVetex *
debrujin_insert_k_mer_affix (DeBrujin *debrujin,
		const char *k_mer_affix)
{
	DeBrujinVetex *vertex = NULL;
	AdjList *adjlist = NULL;

	adjlist = debrujin_adjlist (debrujin, k_mer_affix);

	if (adjlist == NULL)
		{
			vertex = xcalloc (1, sizeof (DeBrujinVetex));
			vertex->k_mer_affix = xstrdup (k_mer_affix);
			graph_ins_vertex (debrujin->graph, vertex);
		}
	else
		vertex = adjlist->vertex;

	return vertex;
}

static inline void
debrujin_insert_edge (DeBrujin *debrujin,
		DeBrujinVetex *a, DeBrujinVetex *b)
{
	DeBrujinVetex *b_adj = NULL;
	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;

	adjlist = graph_adjlist (debrujin->graph, a);
	cur = list_head (adjlist->adjacent);

	for (; cur != NULL; cur = list_next (cur))
		{
			b_adj = list_data (cur);
			if (debrujin_equal (b, b_adj))
				break;
		}

	if (cur == NULL)
		{
			b_adj = xcalloc (1, sizeof (DeBrujinVetex));
			b_adj->k_mer_affix = b->k_mer_affix;
		}

	if (graph_ins_edge (debrujin->graph, a, b_adj))
		{
			a->out_degree++;
			b->in_degree++;
		}

	b_adj->depth++;
}

static inline void
debrujin_insert_multi_edge (DeBrujin *debrujin,
		DeBrujinVetex *a, DeBrujinVetex *b)
{
	DeBrujinVetex *b_adj = NULL;

	b_adj = xcalloc (1, sizeof (DeBrujinVetex));
	b_adj->k_mer_affix = b->k_mer_affix;

	if (graph_ins_multi_edge (debrujin->graph, a, b_adj))
		{
			a->out_degree++;
			b->in_degree++;
		}
}

void
debrujin_insert (DeBrujin *debrujin, const char *seq)
{
	assert (debrujin != NULL && seq != NULL);

	DeBrujinVetex *pref = NULL;
	DeBrujinVetex *suff = NULL;
	char non_overlap = 0;
	int *cov = NULL;
	int len = 0;
	int k = 0;
	int i = 0;

	k = debrujin->k;

	len = strlen (seq);
	if (len < k)
		return;

	// There are len - k + 1 k_mers into a sequence
	for (i = 0; i < (len - k + 1); i++)
		{
			// Set k_mer
			strncpy (debrujin->buf, seq + i, k);

			cov = hash_lookup (debrujin->k_mers, debrujin->buf);
			if (cov == NULL)
				{
					cov = xcalloc (1, sizeof (int));
					hash_insert (debrujin->k_mers,
							xstrdup (debrujin->buf), cov);
				}

			(*cov)++;

			// Backup non overlap base
			non_overlap = debrujin->buf[k - 1];

			// Set prefix k_mer
			debrujin->buf[k - 1] = '\0';
			pref = debrujin_insert_k_mer_affix (debrujin,
					debrujin->buf);

			// Get non overlap base back
			debrujin->buf[k - 1] = non_overlap;

			// Set suffix k_mer
			suff = debrujin_insert_k_mer_affix (debrujin,
					debrujin->buf + 1);

			// Connect path for this k_mer
			debrujin_insert_edge (debrujin, pref, suff);
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

static AdjList *
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

static void
debrujin_reset_bfs (DeBrujin *debrujin)
{
	AdjList *adjlist = NULL;
	DeBrujinVetex *vertex = NULL;
	GraphIter iter;

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			vertex = adjlist->vertex;
			vertex->color = VERTEX_WHITE;
			vertex->dist = -1;
			vertex->parent = NULL;
		}
}

static List *
debrujin_resolve_shortest_path (const DeBrujinVetex *start,
		const DeBrujinVetex *end)
{
	List *path = NULL;
	DeBrujinVetex *cur = NULL;

	path = list_new (NULL);
	cur = (DeBrujinVetex *) end;

	do
		{
			list_prepend (path, cur);

			if (debrujin_equal (start, cur))
				break;

			cur = cur->parent;
		}
	while (cur != NULL);

	return path;
}

static List *
debrujin_shortest_path (DeBrujin *debrujin,
		DeBrujinVetex *start, DeBrujinVetex *end)
{
	List *queue = NULL;
	AdjList *adjlist = NULL;
	AdjList *clr_adjlist = NULL;
	DeBrujinVetex *vertex = NULL;
	DeBrujinVetex *clr_vertex = NULL;
	DeBrujinVetex *adj_vertex = NULL;
	ListElmt *cur = NULL;

	// Init all vertices
	debrujin_reset_bfs (debrujin);

	// Init start node
	start->color = VERTEX_GRAY;
	start->dist = 0;

	// Init start adjlist
	queue = list_new (NULL);
	clr_adjlist = graph_adjlist (debrujin->graph, start);

	// Enqueue start adjlist
	list_append (queue, clr_adjlist);

	// Perform breadth-first search
	while (list_size (queue) > 0)
		{
			// Queue peek
			adjlist = list_data (list_head (queue));
			vertex = adjlist->vertex;

			// Found!
			if (debrujin_equal (vertex, end))
				break;

			// Traverse each vertex in the current adjacency list
			cur = list_head (adjlist->adjacent);
			for (; cur != NULL; cur = list_next (cur))
				{
					adj_vertex = list_data (cur);

					clr_adjlist = graph_adjlist (debrujin->graph, adj_vertex);
					clr_vertex = clr_adjlist->vertex;

					if (clr_vertex->color == VERTEX_WHITE)
						{
							// Dicover vertex and set parent and
							// distance from start node
							clr_vertex->color = VERTEX_GRAY;
							clr_vertex->dist = vertex->dist + 1;
							clr_vertex->parent = vertex;

							// Enqueue adjlist
							list_append (queue, clr_adjlist);
						}
				}

			// Dequeue the current adjacency list and color its vertex black
			list_remove (queue, list_head (queue), (void **) &adjlist);
			vertex->color = VERTEX_BLACK;
		}

	list_free (queue);

	return debrujin_equal (vertex, end)
		? debrujin_resolve_shortest_path (start, end)
		: NULL;
}

static void
debrujin_route_inspection (DeBrujin *debrujin)
{
	// PONGA
}

static void
debrujin_dfs (DeBrujin *debrujin, AdjList *adjlist, List *tour)
{
	AdjList *clr_adjlist = NULL;
	DeBrujinVetex *adj_vertex = NULL;
	DeBrujinVetex *vertex = NULL;

	vertex = adjlist->vertex;

	while (vertex->cur_adj != NULL)
		{
			// Determine the color of the next adjacent vertex
			adj_vertex = list_data (vertex->cur_adj);

			// Get adjacency list of the clr_vertex
			clr_adjlist = graph_adjlist (debrujin->graph, adj_vertex);

			// Move one vertex deeper
			vertex->cur_adj = list_next (vertex->cur_adj);
			debrujin_dfs (debrujin, clr_adjlist, tour);
		}

	// Sequence grows backwards
	list_prepend (tour, vertex->k_mer_affix);
}

static const char *
debrujin_sequence (const DeBrujin *debrujin, const List *tour)
{
	ListElmt *cur = NULL;
	char *seq = NULL;
	int k = 0;
	int i = 0;

	cur = list_head (tour);
	if (cur == NULL)
		return NULL;

	k = debrujin->k - 1;

	// seq has length num_k-mer + k - 1
	seq = xcalloc (list_size (tour) + k, sizeof (char));

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
	List *tour = NULL;
	const char *seq = NULL;

	/*while (!debrujin_has_eulerian_path (debrujin))*/
		/*debrujin_route_inspection (debrujin);*/

	debrujin_reset_path (debrujin);

	seqs = list_new (xfree);

	//
	AdjList *adjlist = NULL;
	AdjList *adj_start = NULL;
	DeBrujinVetex *vertex = NULL;
	GraphIter iter;

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			vertex = adjlist->vertex;
			if ((vertex->out_degree - vertex->in_degree) == 1)
				{
					tour = list_new (NULL);

					debrujin_dfs (debrujin,
							adjlist,
							tour);

					seq = debrujin_sequence (debrujin, tour);

					if (seq != NULL)
						list_append (seqs, seq);

					list_free (tour);
				}

			if (vertex->out_degree > 0)
				adj_start = adjlist;
		}

	if (list_size (seqs) == 0 && adj_start != NULL)
		{
			tour = list_new (NULL);

			debrujin_dfs (debrujin,
					adjlist,
					tour);

			seq = debrujin_sequence (debrujin, tour);

			if (seq != NULL)
				list_append (seqs, seq);

			list_free (tour);
		}
	//

	return seqs;
}
