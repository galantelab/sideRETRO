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
#include "debrujin.h"

struct _DeBrujin
{
	int          k;
	Hash        *k_mers;

	char        *buf;

	Graph       *graph;
};

struct _DeBrujinVetex
{
	int          in_degree;
	int          out_degree;
	int          depth;

	int          index;
	int          lowlink;
	int          on_stack;

	const char  *k_mer_affix;

	ListElmt    *cur_adj;
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
	if (debrujin == NULL)
		return;

	xfree (debrujin->buf);
	hash_free (debrujin->k_mers);
	graph_free (debrujin->graph);

	xfree (debrujin);
}

static DeBrujinVetex *
debrujin_insert_k_mer_affix (DeBrujin *debrujin,
		const char *k_mer_affix)
{
	DeBrujinVetex *vertex = NULL;
	AdjList *adjlist = NULL;

	DeBrujinVetex temp = { .k_mer_affix = k_mer_affix };
	adjlist = graph_adjlist (debrujin->graph, &temp);

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
			pref->depth++;

			// Get non overlap base back
			debrujin->buf[k - 1] = non_overlap;

			// Set suffix k_mer
			suff = debrujin_insert_k_mer_affix (debrujin,
					debrujin->buf + 1);
			suff->depth++;

			// Connect path for this k_mer
			if (graph_ins_edge (debrujin->graph, pref, suff))
				{
					pref->out_degree++;
					suff->in_degree++;
				}
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

static void
debrujin_dfs (DeBrujin *debrujin, AdjList *adjlist, List *tour)
{
	AdjList *clr_adjlist = NULL;
	DeBrujinVetex *clr_vertex = NULL;
	DeBrujinVetex *vertex = NULL;

	vertex = adjlist->vertex;

	while (vertex->cur_adj != NULL)
		{
			// Determine the color of the next adjacent vertex
			clr_vertex = list_data (vertex->cur_adj);

			// Get adjacency list of the clr_vertex
			clr_adjlist = graph_adjlist (debrujin->graph, clr_vertex);

			// Move one vertex deeper
			vertex->cur_adj = list_next (vertex->cur_adj);
			debrujin_dfs (debrujin, clr_adjlist, tour);
		}

	// Sequence grows backwards
	list_prepend (tour, vertex->k_mer_affix);
}

///////////////////////////////////////////
#include "log.h"

static void
strongconnect (DeBrujin *debrujin, AdjList *adjlist,
		int *index_p, List *stack)
{
	AdjList *clr_adjlist = NULL;
	DeBrujinVetex *clr_vertex = NULL;
	DeBrujinVetex *vertex = NULL;
	ListElmt *cur = NULL;

	vertex = adjlist->vertex;
	vertex->index = vertex->lowlink = ++ *index_p;
	vertex->on_stack = 1;

	list_append (stack, adjlist);

	cur = list_head (adjlist->adjacent);
	for (; cur != NULL; cur = list_next (cur))
		{
			clr_vertex = list_data (cur);
			if (clr_vertex->index == 0)
				{
					clr_adjlist = graph_adjlist (debrujin->graph, clr_vertex);
					strongconnect (debrujin, clr_adjlist, index_p, stack);

					if (clr_vertex->lowlink < vertex->lowlink)
						vertex->lowlink = clr_vertex->lowlink;
				}
			else if (clr_vertex->on_stack
					&& clr_vertex->index < vertex->lowlink)
				vertex->lowlink = clr_vertex->index;
		}

	if (vertex->lowlink == vertex->index)
		{
			do
				{
					list_remove (stack, list_tail (stack),
							(void **) &clr_adjlist);

					clr_vertex = clr_adjlist->vertex;
					clr_vertex->on_stack = 0;
					log_info ("[%d] [%d %d]: %s", vertex->index, clr_vertex->index,
							clr_vertex->lowlink, clr_vertex->k_mer_affix);
				}
			while (clr_vertex != vertex);
		}
}

static void
tarjan (DeBrujin *debrujin)
{
	int index = 0;
	List *stack = NULL;
	AdjList *adjlist = NULL;
	DeBrujinVetex *vertex = NULL;
	GraphIter iter;

	stack = list_new (NULL);

	graph_iter_init (&iter, debrujin->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			vertex = adjlist->vertex;
			if (vertex->index == 0)
				strongconnect (debrujin, adjlist, &index, stack);
		}

	log_info (":: stack_size: %ld", list_size (stack));
	list_free (stack);
}

/*void*/
/*debrujin_fix_edges (DeBrujin *debrujin)*/
/*{*/
	/*AdjList *adjlist = NULL;*/
	/*ListElmt *cur = NULL;*/
	/*DeBrujinVetex *vertex = NULL;*/
	/*DeBrujinVetex *clr_vertex = NULL;*/
	/*int diff = 0;*/
	/*int i = 0;*/
	/*GraphIter iter;*/

	/*graph_iter_init (&iter, debrujin->graph);*/
	/*while (graph_iter_next (&iter, &adjlist))*/
		/*{*/
			/*vertex = adjlist->vertex;*/
			/*if (vertex->in_degree > vertex->out_degree)*/
				/*{*/
					/*diff = vertex->in_degree - vertex->out_degree;*/
					/*cur = list_head (adjlist->adjacent);*/
					/*for (; cur != NULL; cur = list_next (cur))*/
						/*{*/
							/*clr_vertex = list_data (cur);*/
							/*[>if ((clr_vertex->out_degree - clr_vertex->in_degree) == diff)<]*/
								/*[>{<]*/
									/*for (i = 0; i < diff; i++)*/
										/*{*/
											/*graph_ins_multi_edge (debrujin->graph,*/
													/*vertex, clr_vertex);*/
											/*vertex->out_degree++;*/
											/*clr_vertex->in_degree++;*/
										/*}*/
									/*[>break;<]*/
								/*[>}<]*/
						/*}*/
				/*}*/
		/*}*/
/*}*/
///////////////////////////////////////////

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

	/*if (!debrujin_has_eulerian_path (debrujin))*/
		/*return NULL;*/

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
