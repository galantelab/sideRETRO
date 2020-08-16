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

#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include "list.h"
#include "set.h"
#include "types.h"

struct _AdjList
{
	void          *vertex;
	Set           *adjacent;
};

typedef struct _AdjList AdjList;

struct _Graph
{
	size_t         vcount;
	size_t         ecount;

	HashFunc       hash_fun;
	EqualFun       match_fun;
	DestroyNotify  destroy_fun;

	List          *adjlists;
};

typedef struct _Graph Graph;

enum _VertexColor
{
	VERTEX_WHITE,
	VERTEX_GRAY,
	VERTEX_BLACK
};

typedef enum _VertexColor VertexColor;

Graph   * graph_new_full    (HashFunc hash_fun, EqualFun equal_fun, DestroyNotify destroy_fun);
void      graph_free        (Graph *graph);
int       graph_ins_vertex  (Graph *graph, const void *data);
int       graph_ins_edge    (Graph *graph, const void *data1, const void *data2);
int       graph_rem_vertex  (Graph *graph, void **data);
int       graph_rem_edge    (Graph *graph, const void *data1, void **data2);
AdjList * graph_adjlist     (const Graph *graph, const void *data);
int       graph_is_adjacent (const Graph *graph, const void *data1, const void *data2);

#define graph_adjlists(graph) ((graph)->adjlists)
#define graph_vcount(graph)   ((graph)->vcount)
#define graph_ecount(graph)   ((graph)->ecount)

#endif /* graph.h */
