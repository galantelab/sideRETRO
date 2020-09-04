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

#ifndef FLOYD_WARSHALL_H
#define FLOYD_WARSHALL_H

#include <stdlib.h>
#include "graph.h"
#include "list.h"
#include "hash.h"
#include "types.h"

typedef double (*WFunc) (void *vertex1, void *vertex2);

struct _FloydWarshal
{
	Graph   *graph;
	Hash    *nodes;

	HashFunc hash_fun;
	EqualFun match_fun;
	WFunc    weight_fun;

	size_t   size;
	double **dist;
	int    **next;
};

typedef struct _FloydWarshal FloydWarshal;

FloydWarshal * floyd_warshall_new  (const Graph *graph,
		HashFunc hash_fun, EqualFun match_fun,
		WFunc weight_fun);

void  floyd_warshall_free (FloydWarshal *fw);
void  floyd_warshall_run  (FloydWarshal *fw);

double floyd_warshall_dist (FloydWarshal *fw,
		const void *from, const void *to);

#endif /* floyd_warshall.h */
