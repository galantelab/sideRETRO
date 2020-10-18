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
#include "graph_enumerate.h"

static inline void
__enum_push (List *list, GraphLabel label, int sense, AdjList *parent,
		AdjList *child1, AdjList *child2)
{
	GraphEnum *e = xcalloc (1, sizeof (GraphEnum));

	*e = (GraphEnum) {
		.label  = label,
		.sense  = sense,
		.parent = parent,
		.child1 = child1,
		.child2 = child2
	};

	list_append (list, e);
}

static inline void
__enumerate_adjacent (Unipath *unipath, EqualFun match_fun,
		AdjList *adjlist, List *labels)
{
	AdjList *clr_adjlist1 = NULL;
	AdjList *clr_adjlist2 = NULL;
	ListElmt *cur1 = NULL;
	ListElmt *cur2 = NULL;

	cur1 = list_head (adjlist->adjacent);
	for (; cur1 != NULL; cur1 = list_next (cur1))
		{
			clr_adjlist1 = graph_adjlist (unipath, list_data (cur1));

			if (match_fun (adjlist->vertex, clr_adjlist1->vertex))
				{
					__enum_push (labels, GRAPH_LABEL_CYCLE, 1,
							adjlist, clr_adjlist1, NULL);
					continue;
				}

			cur2 = list_next (cur1);
			for (; cur2 != NULL; cur2 = list_next (cur2))
				{
					clr_adjlist2 = graph_adjlist (unipath, list_data (cur2));

					// Ignore cycle here
					if (match_fun (adjlist->vertex, clr_adjlist2->vertex))
						continue;

					if (match_fun (clr_adjlist1->vertex, clr_adjlist2->vertex))
						{
							__enum_push (labels, GRAPH_LABEL_BUBBLE, 1,
									adjlist, clr_adjlist1, NULL);
						}
					else if (list_size (clr_adjlist1->adjacent) == 0
							&& list_size (clr_adjlist2->adjacent) > 0)
						{
							__enum_push (labels, GRAPH_LABEL_TIP, 1,
									adjlist, clr_adjlist1, clr_adjlist2);
						}
					else if (list_size (clr_adjlist2->adjacent) == 0
							&& list_size (clr_adjlist1->adjacent) > 0)
						{
							__enum_push (labels, GRAPH_LABEL_TIP, 1,
									adjlist, clr_adjlist2, clr_adjlist1);
						}
					else if ((list_size (clr_adjlist1->adjacent) > 0
								&& list_size (clr_adjlist2->adjacent) > 0)
							|| (list_size (clr_adjlist1->adjacent) == 0
								&& list_size (clr_adjlist2->adjacent) == 0))
						{
							__enum_push (labels, GRAPH_LABEL_FORK, 1,
									adjlist, clr_adjlist1, clr_adjlist2);
						}
				}
		}
}

static inline void
__enumerate_parent (Unipath *unipath, EqualFun match_fun,
		AdjList *adjlist, List *labels)
{
	AdjList *clr_adjlist1 = NULL;
	AdjList *clr_adjlist2 = NULL;
	ListElmt *cur1 = NULL;
	ListElmt *cur2 = NULL;

	cur1 = list_head (adjlist->parent);
	for (; cur1 != NULL; cur1 = list_next (cur1))
		{
			clr_adjlist1 = graph_adjlist (unipath, list_data (cur1));

			// Ignore cycle here
			if (match_fun (adjlist->vertex, clr_adjlist1->vertex))
				continue;

			cur2 = list_next (cur1);
			for (; cur2 != NULL; cur2 = list_next (cur2))
				{
					clr_adjlist2 = graph_adjlist (unipath, list_data (cur2));

					// Ignore cycle at all
					if (match_fun (adjlist->vertex, clr_adjlist2->vertex))
						continue;

					// Detect bubbles only in adjacent. Ignore in parent, in order
					// to avoid duplications
					if (match_fun (clr_adjlist1->vertex, clr_adjlist2->vertex))
						{
							continue;
						}
					else if (list_size (clr_adjlist1->parent) == 0
							&& list_size (clr_adjlist2->parent) > 0)
						{
							__enum_push (labels, GRAPH_LABEL_TIP, 0,
									adjlist, clr_adjlist1, clr_adjlist2);
						}
					else if (list_size (clr_adjlist2->parent) == 0
							&& list_size (clr_adjlist1->parent) > 0)
						{
							__enum_push (labels, GRAPH_LABEL_TIP, 0,
									adjlist, clr_adjlist2, clr_adjlist1);
						}
					else if ((list_size (clr_adjlist1->parent) > 0
								&& list_size (clr_adjlist2->parent) > 0)
							|| (list_size (clr_adjlist1->parent) == 0
								&& list_size (clr_adjlist2->parent) == 0))
						{
							__enum_push (labels, GRAPH_LABEL_FORK, 0,
									adjlist, clr_adjlist1, clr_adjlist2);
						}
				}
		}
}

List *
graph_enumerate (Unipath *unipath, EqualFun match_fun)
{
	assert (unipath != NULL && match_fun != NULL);

	AdjList *adjlist = NULL;
	List *labels = NULL;
	GraphIter iter = {};

	labels = list_new (xfree);

	graph_iter_init (&iter, unipath);
	while (graph_iter_next (&iter, &adjlist))
		{
			__enumerate_adjacent (unipath, match_fun, adjlist, labels);
			__enumerate_parent (unipath, match_fun, adjlist, labels);
		}

	return labels;
}
