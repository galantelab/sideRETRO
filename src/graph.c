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
#include "wrapper.h"
#include "graph.h"

Graph *
graph_new_full (HashFunc hash_fun, EqualFun match_fun,
		DestroyNotify destroy_fun)
{
	assert (hash_fun != NULL && match_fun != NULL);

	Graph *graph = xcalloc (1, sizeof (Graph));

	graph->adjlists = hash_new_full (hash_fun,
			match_fun, NULL, NULL);

	graph->destroy_fun = destroy_fun;
	graph->match_fun = match_fun;

	return graph;
}

void
graph_free (Graph *graph)
{
	HashIter iter;
	AdjList *adjlist;

	if (graph == NULL)
		return;

	hash_iter_init (&iter, graph->adjlists);

	while (hash_iter_next (&iter, NULL, (void **) &adjlist))
		{
			if (graph->destroy_fun != NULL)
				graph->destroy_fun (adjlist->vertex);

			list_free (adjlist->adjacent);
			xfree (adjlist);
		}

	hash_free (graph->adjlists);
	xfree (graph);
}

int
graph_ins_vertex (Graph *graph, const void *data)
{
	assert (graph != NULL && data != NULL);

	AdjList *adjlist = NULL;

	adjlist = hash_lookup (graph->adjlists, data);

	if (adjlist != NULL)
		return 0;

	adjlist = xcalloc (1, sizeof (AdjList));

	adjlist->vertex = (void *) data;
	adjlist->adjacent = list_new (NULL);

	hash_insert (graph->adjlists, data, adjlist);
	graph->vcount++;

	return 1;
}

static inline ListElmt *
__graph_adjlist_elmt (const Graph *graph,
		const AdjList *adjlist, const void *data)
{
	ListElmt *cur = NULL;

	cur = list_head (adjlist->adjacent);
	for (; cur != NULL; cur = list_next (cur))
		{
			if (graph->match_fun (list_data (cur), data))
				break;
		}

	return cur;
}

static inline int
__graph_ins_edge (Graph *graph, const void *data1,
		const void *data2, int is_multi)
{
	AdjList *adjlist1 = NULL;
	AdjList *adjlist2 = NULL;

	adjlist1 = hash_lookup (graph->adjlists, data1);
	if (adjlist1 == NULL)
		return 0;

	adjlist2 = hash_lookup (graph->adjlists, data2);
	if (adjlist2 == NULL)
		return 0;

	// In non multi-edge graph, do not allow the insertion of
	// repetitive vetice into the adjacency list
	if (!is_multi && __graph_adjlist_elmt (graph, adjlist1,
				adjlist2->vertex) != NULL)
		return 0;

	// Insert the second vertex into the adjacency
	// list of the first vertex.
	list_append (adjlist1->adjacent, adjlist2->vertex);
	graph->ecount++;

	return 1;
}

int
graph_ins_multi_edge (Graph *graph, const void *data1, const void *data2)
{
	assert (graph != NULL && data1 != NULL && data2 != NULL);
	return __graph_ins_edge (graph, data1, data2, 1);
}

int
graph_ins_edge (Graph *graph, const void *data1, const void *data2)
{
	assert (graph != NULL && data1 != NULL && data2 != NULL);
	return __graph_ins_edge (graph, data1, data2, 0);
}

int
graph_is_adjacent (const Graph *graph, const void *data1,
		const void *data2)
{
	assert (graph != NULL && data1 != NULL && data2 != NULL);

	AdjList *adjlist = NULL;

	// Locate the adjacency list for the first vertex
	adjlist = hash_lookup (graph->adjlists, data1);
	if (adjlist == NULL)
		return 0;

	// Return whether the second vertex is in the
	// adjacency list of the first
	return __graph_adjlist_elmt (graph, adjlist, data2) != NULL;
}

int
graph_rem_vertex (Graph *graph, void **data)
{
	assert (graph != NULL && data != NULL && *data != NULL);

	AdjList *adjlist = NULL;
	AdjList *iter_adjlist = NULL;
	HashIter iter;

	adjlist = hash_lookup (graph->adjlists, *data);
	if (adjlist == NULL)
		return 0;

	// Do not allow removal of the vertex
	// if its adjacency list is not empty
	if (list_size (adjlist->adjacent) > 0)
		return 0;

	// Do not allow removal of the vertex
	// if it is in an adjacency list
	hash_iter_init (&iter, graph->adjlists);
	while (hash_iter_next (&iter, NULL, (void **) &iter_adjlist))
		if (__graph_adjlist_elmt (graph,
					iter_adjlist, *data) != NULL)
			return 0;

	hash_remove (graph->adjlists, *data);
	*data = adjlist->vertex;

	list_free (adjlist->adjacent);
	xfree (adjlist);
	graph->vcount--;

	return 1;
}

int
graph_rem_edge (Graph *graph, const void *data1, void **data2)
{
	assert (graph != NULL && data1 != NULL && data2 != NULL
			&& *data2 != NULL);

	AdjList *adjlist = NULL;
	ListElmt *found = NULL;

	// Locate the adjacency list for the first vertex
	adjlist = hash_lookup (graph->adjlists, data1);
	if (adjlist == NULL)
		return 0;

	// Locate the second vertex into the adjacency
	// list of the first vertex
	found = __graph_adjlist_elmt (graph, adjlist, *data2);
	if (found == NULL)
		return 0;

	// Remove the second vertex from the adjacency
	// list of the first vertex
	list_remove (adjlist->adjacent, found, data2);
	graph->ecount--;

	return 1;
}
