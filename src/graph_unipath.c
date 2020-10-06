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

#include <assert.h>
#include "hash.h"
#include "wrapper.h"
#include "graph_unipath.h"

static inline int
graph_unipath_is_non_mergeable (const AdjList *adjlist)
{
	return graph_in_degree (adjlist) > 1
		|| graph_out_degree (adjlist) > 1
		|| graph_in_degree (adjlist) == 0
		|| graph_out_degree (adjlist) == 0;
}

static inline void *
graph_unipath_child (const AdjList *adjlist)
{
	return list_size (adjlist->adjacent) > 0
		? list_data (list_head (adjlist->adjacent))
		: NULL;
}

static void
graph_unipath_build (Graph *unigraph, const Graph *graph, Hash *connect)
{
	AdjList *adjlist = NULL;
	AdjList *adjlist_cpy = NULL;
	ListElmt *cur = NULL;
	void *vertex1 = NULL;
	void *vertex2 = NULL;
	GraphIter iter = {};

	// Add vertices
	graph_iter_init (&iter, graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			if (graph_unipath_is_non_mergeable (adjlist))
				graph_ins_vertex (unigraph, adjlist->vertex);
			else if (graph_out_degree (adjlist) > 0)
				hash_insert (connect, adjlist->vertex,
						graph_unipath_child (adjlist));
		}

	graph_iter_init (&iter, unigraph);
	while (graph_iter_next (&iter, &adjlist_cpy))
		{
			// Get the original vertex
			adjlist = graph_adjlist (graph,
					adjlist_cpy->vertex);

			// Search for all edges
			cur = list_head (adjlist->adjacent);
			for (; cur != NULL; cur = list_next (cur))
				{
					vertex1 = list_data (cur);
					vertex2 = hash_lookup (connect, vertex1);

					// Maybe there is a path from vertex1 to any
					// unigraph vertex. Search for the other
					// end of the path
					for (; vertex2 != NULL; vertex2 = hash_lookup (connect, vertex1))
						vertex1 = vertex2;

					// Maybe this child is a unipath vertex
					if (graph_adjlist (unigraph, vertex1) != NULL)
						graph_ins_multi_edge (unigraph, adjlist->vertex, vertex1);
				}
		}
}

Graph *
graph_unipath_new (const Graph *graph, HashFunc hash_fun,
		EqualFun equal_fun)
{
	assert (graph != NULL);

	Graph *unigraph = NULL;
	Hash *connect = NULL;

	unigraph = graph_new_full (hash_fun, equal_fun, NULL);
	connect = hash_new_full (hash_fun, equal_fun, NULL, NULL);

	graph_unipath_build (unigraph, graph, connect);

	hash_free (connect);
	return unigraph;
}
