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

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "wrapper.h"
#include "floyd_warshall.h"

struct _FloydWarshal
{
	Graph      *graph;
	Hash       *nodes;

	HashFunc    hash_fun;
	EqualFun    match_fun;
	WFunc       weight_fun;

	size_t      size;
	double    **dist;
	void     ***next;
};

#define MTX_NEW(size, _type) ({ \
	_type **mtx = NULL; \
	do { \
		size_t i = 0; \
		mtx = xcalloc (size, sizeof (_type *)); \
		for (i = 0; i < size; i++) \
			mtx[i] = xcalloc (size, sizeof (_type)); \
	} while (0); \
	mtx; \
})

#define MTX_FREE(mtx,size) ({ \
	do { \
		size_t i = 0; \
		for (i = 0; i < size; i++) \
			xfree (mtx[i]); \
		xfree (mtx); \
	} while (0); \
})

static void
floyd_warshall_get_nodes (Hash *nodes, const Graph *graph)
{
	AdjList *adjlist = NULL;
	GraphIter iter = {};
	int i = 0;
	int *i_alloc = NULL;

	graph_iter_init (&iter, graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			i_alloc = xcalloc (1, sizeof (int));
			*i_alloc = i++;
			hash_insert (nodes, adjlist->vertex, i_alloc);
		}
}

FloydWarshal *
floyd_warshall_new  (const Graph *graph,
		HashFunc hash_fun, EqualFun match_fun,
		WFunc weight_fun)
{
	assert (graph != NULL && hash_fun != NULL
			&& match_fun != NULL
			&& weight_fun != NULL);

	FloydWarshal *fw = xcalloc (1, sizeof (FloydWarshal));
	fw->graph = (Graph *) graph;

	fw->hash_fun = hash_fun;
	fw->match_fun = match_fun;
	fw->weight_fun = weight_fun;

	fw->nodes = hash_new_full (hash_fun, match_fun,
			NULL, xfree);

	floyd_warshall_get_nodes (fw->nodes, graph);
	fw->size = hash_size (fw->nodes);

	fw->dist = MTX_NEW (fw->size, double);
	fw->next = MTX_NEW (fw->size, void *);

	return fw;
}

void
floyd_warshall_free (FloydWarshal *fw)
{
	if (fw == NULL)
		return;

	hash_free (fw->nodes);

	MTX_FREE (fw->dist, fw->size);
	MTX_FREE (fw->next, fw->size);

	xfree (fw);
}

static void
floyd_warshall_run_init (FloydWarshal *fw)
{
	Array *keys = NULL;
	int *i_alloc = NULL;
	int *j_alloc = NULL;
	void *a = NULL;
	void *b = NULL;
	double **dist = NULL;
	void ***next = NULL;
	size_t size, i, j;

	size = fw->size;
	dist = fw->dist;
	next = fw->next;

	keys = hash_get_keys_as_array (fw->nodes);

	for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
				{
					// Get nodes
					a = array_get (keys, i);
					b = array_get (keys, j);

					if (i == j)
						{
							dist[i][i] = 0;
							next[i][i] = a;
						}
					else
						{
							// The hash contains the index associated
							// with each node. So, correct the index
							i_alloc = hash_lookup (fw->nodes, a);
							j_alloc = hash_lookup (fw->nodes, b);

							if (graph_is_adjacent (fw->graph, a, b))
								{
									dist[*i_alloc][*j_alloc] = fw->weight_fun (a, b);
									next[*i_alloc][*j_alloc] = b;
								}
							else
								dist[*i_alloc][*j_alloc] = INFINITY;
						}
				}
		}

	array_free (keys, 1);
}

void
floyd_warshall_run (FloydWarshal *fw)
{
	assert (fw != NULL);

	double **dist = NULL;
	void ***next = NULL;
	size_t size, i, j, k;

	size = fw->size;
	dist = fw->dist;
	next = fw->next;

	// Init matrixes
	floyd_warshall_run_init (fw);

	// Standard Floyd-Warshall implementation
	for (k = 0; k < size; k++)
		{
			for (i = 0; i < size; i++)
				{
					for (j = 0; j < size; j++)
						{
							if (dist[i][j] > (dist[i][k] + dist[k][j]))
								{
									dist[i][j] = dist[i][k] + dist[k][j];
									next[i][j] = next[i][k];
								}
						}
				}
		}
}

double
floyd_warshall_dist (const FloydWarshal *fw,
		const void *from, const void *to)
{
	assert (fw != NULL && from != NULL && to != NULL);

	int *i_alloc = NULL;
	int *j_alloc = NULL;

	i_alloc = hash_lookup (fw->nodes, from);
	j_alloc = hash_lookup (fw->nodes, to);

	if (i_alloc == NULL || j_alloc == NULL)
		return NAN;

	return fw->dist[*i_alloc][*j_alloc];
}

List *
floyd_warshall_path (const FloydWarshal *fw,
		const void *from, const void *to)
{
	assert (fw != NULL && from != NULL && to != NULL);

	List *path = NULL;
	void *a = NULL;
	void *b = NULL;
	int *i_alloc = NULL;
	int *j_alloc = NULL;

	i_alloc = hash_lookup (fw->nodes, from);
	j_alloc = hash_lookup (fw->nodes, to);

	if (i_alloc == NULL || j_alloc == NULL)
		return NULL;

	a = (void *) from;
	b = (void *) to;

	path = list_new (NULL);
	list_append (path, a);

	while (!fw->match_fun (a, b))
		{
			a = fw->next[*i_alloc][*j_alloc];
			i_alloc = hash_lookup (fw->nodes, a);
			list_append (path, a);
		}

	return path;
}
