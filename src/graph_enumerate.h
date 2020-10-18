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

#ifndef GRAPH_ENUMERATE_H
#define GRAPH_ENUMERATE_H

#include "graph_unipath.h"
#include "list.h"
#include "types.h"

enum _GraphLabel
{
	GRAPH_LABEL_NONE,
	GRAPH_LABEL_TIP,
	GRAPH_LABEL_FORK,
	GRAPH_LABEL_CYCLE,
	GRAPH_LABEL_BUBBLE
};

typedef enum _GraphLabel GraphLabel;

struct _GraphEnum
{
	GraphLabel    label;
	int           sense;
	AdjList      *parent;
	AdjList      *child1;
	AdjList      *child2;
};

typedef struct _GraphEnum GraphEnum;

List * graph_enumerate (Unipath *unipath, EqualFun match_fun);

#endif /* graph_enumerate.h */
