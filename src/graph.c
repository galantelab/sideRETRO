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
	assert (match_fun != NULL);

	Graph *graph = xcalloc (1, sizeof (Graph));

	graph->adjlists = list_new (NULL);

	graph->hash_fun = hash_fun;
	graph->match_fun = match_fun;
	graph->destroy_fun = destroy_fun;

	return graph;
}

void
graph_free (Graph *graph)
{
	if (graph == NULL)
		return;

	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;

	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			set_free (adjlist->adjacent);

			if (graph->destroy_fun != NULL)
				graph->destroy_fun (adjlist->vertex);

			xfree (adjlist);
		}

	list_free (graph->adjlists);
	xfree (graph);
}

int
graph_ins_vertex (Graph *graph, const void *data)
{
	assert (graph != NULL && data != NULL);

	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;

	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			if (graph->match_fun (data, adjlist->vertex))
				return 0;
		}

	adjlist = xcalloc (1, sizeof (AdjList));
	adjlist->vertex = (void *) data;
	adjlist->adjacent = set_new_full (graph->hash_fun,
			graph->match_fun, NULL);

	list_append (graph->adjlists, adjlist);
	graph->vcount++;

	return 1;
}

int
graph_ins_edge (Graph *graph, const void *data1, const void *data2)
{
	assert (graph != NULL && data1 != NULL && data2 != NULL);

	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;

	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			if (graph->match_fun (data2, adjlist->vertex))
				break;
		}

	if (cur == NULL)
		return 0;

	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			if (graph->match_fun (data1, adjlist->vertex))
				break;
		}

	if (cur == NULL)
		return 0;

	// Insert the second vertex into the adjacency
	// list of the first vertex.
	if (!set_insert (adjlist->adjacent, data2))
		return 0;

	graph->ecount++;
	return 1;
}

int
graph_rem_vertex (Graph *graph, void **data)
{
	assert (graph != NULL && data != NULL && *data != NULL);

	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;
	ListElmt *found = NULL;

	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);

			// Do not allow removal of the vertex
			// if it is in an adjacency list
			if (set_is_member (adjlist->adjacent, *data))
				return 0;

			if (graph->match_fun (adjlist->vertex, *data))
				found = cur;
		}

	if (found == NULL)
		return 0;

	// Do not allow removal of the vertex
	// if its adjacency list is not empty
	adjlist = list_data (found);
	if (set_size (adjlist->adjacent) > 0)
		return 0;

	list_remove (graph->adjlists, found, NULL);

	*data = adjlist->vertex;
	set_free (adjlist->adjacent);
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
	ListElmt *cur = NULL;

	// Locate the adjacency list for the first vertex
	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			if (graph->match_fun (adjlist->vertex, data1))
				break;
		}

	if (cur == NULL)
		return 0;

	// Remove the second vertex from the adjacency
	// list of the first vertex
	if (!set_remove (adjlist->adjacent, data2))
		return 0;

	graph->ecount--;
	return 1;
}

AdjList *
graph_adjlist (const Graph *graph, const void *data)
{
	assert (graph != NULL && data != NULL);

	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;

	// Locate the adjacency list for the vertex
	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			if (graph->match_fun (adjlist->vertex, data))
				break;
		}

	return cur != NULL ? list_data (cur) : NULL;
}

int
graph_is_adjacent (const Graph *graph, const void *data1,
		const void *data2)
{
	assert (graph != NULL && data1 != NULL && data2 != NULL);

	AdjList *adjlist = NULL;
	ListElmt *cur = NULL;

	// Locate the adjacency list for the first vertex
	cur = list_head (graph->adjlists);
	for (; cur != NULL; cur = list_next (cur))
		{
			adjlist = list_data (cur);
			if (graph->match_fun (adjlist->vertex, data1))
				break;
		}

	if (cur == NULL)
		return 0;

	// Return whether the second vertex is in the
	// adjacency list of the first
	return set_is_member (adjlist->adjacent, data2);
}
