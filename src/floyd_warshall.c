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

#include <math.h>
#include <assert.h>
#include "wrapper.h"
#include "floyd_warshall.h"

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

FloydWarshal *
floyd_warshall_new (const Graph *graph, WFunc weight_fun)
{
	assert (graph != NULL && weight_fun != NULL);

	FloydWarshal *fw = xcalloc (1, sizeof (FloydWarshal));

	fw->graph = (Graph *) graph;
	fw->weight_fun = weight_fun;

	fw->nodes = graph_adjlists_as_array (graph);
	fw->size = array_len (fw->nodes);

	fw->dist = MTX_NEW (fw->size, double);
	fw->next = MTX_NEW (fw->size, int);

	return fw;
}

void
floyd_warshall_free (FloydWarshal *fw)
{
	if (fw == NULL)
		return;

	array_free (fw->nodes, 1);

	MTX_FREE (fw->dist, fw->size);
	MTX_FREE (fw->next, fw->size);

	xfree (fw);
}

void
floyd_warshall_run (FloydWarshal *fw)
{
	assert (fw != NULL);

	AdjList *a = NULL;
	AdjList *b = NULL;
	double **dist = NULL;
	int **next = NULL;
	size_t size, i, j, k;

	size = fw->size;
	dist = fw->dist;
	next = fw->next;

	// Init matrixes
	for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
				{
						if (i == j)
							{
								dist[i][i] = 0;
								next[i][i] = i;
							}
						else
							{
								a = array_get (fw->nodes, i);
								b = array_get (fw->nodes, j);

								if (graph_is_adjacent (fw->graph,
											a->vertex, b->vertex))
									{
										dist[i][j] = fw->weight_fun (a->vertex, b->vertex);
										next[i][j] = j;
									}
								else
									dist[i][j] = INFINITY;
							}
				}
		}

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
