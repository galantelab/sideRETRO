#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "../src/wrapper.h"
#include "../src/array.h"
#include "../src/log.h"
#include "../src/utils.h"
#include "../src/graph.h"
#include "../src/hungarian.h"
#include "../src/debrujin.c"

#define MYSELF "sucateiro"
#define K_MER 55
#define MIN_COV 3

#define MTX_NEW(size1,size2,_type) ({ \
	_type **mtx = NULL; \
	do { \
		size_t i = 0; \
		mtx = xcalloc (size1, sizeof (_type *)); \
		for (i = 0; i < size1; i++) \
			mtx[i] = xcalloc (size2, sizeof (_type)); \
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

void
print_usage (FILE *fp)
{
	fprintf (fp,
			"%s\n"
			"\n"
			"Usage: %s [-h] [-q] [-d] [-k] [-c] FILE\n"
			"\n"
			"Assembler for testing purpose\n"
			"\n"
			"Arguments:\n"
			"    FILE: One read per line\n"
			"\n"
			"Input/Output Options:\n"
			"   -h, --help        Show help options\n"
			"   -q, --quiet       Decrease verbosity\n"
			"   -d, --debug       Increase verbosity\n"
			"                     Pass twice for trace\n"
			"                     verbosity\n"
			"\n"
			"Assembly Options:\n"
			"   -k, --k-mer       K-mer size [%d]\n"
			"   -c, --min-cov     Minimum coverage [%d]\n"
			"\n",
			MYSELF, MYSELF, K_MER, MIN_COV);
}

void
parse_file_k_mers (const char *file, DeBrujin *d)
{
	FILE *fp = NULL;
	char *line = NULL;
	size_t len = 0;
	size_t nread;

	fp = xfopen (file, "r");
	errno = 0;

	while ((nread = getline (&line, &len, fp)) != EOF)
		{
			if (errno)
				log_errno_fatal ("Failed to read '%s'",
						file);

			chomp (line);
			trim (line);

			// Line is empty
			if (*line == '\0')
				continue;

			log_trace ("Insert read %s", line);
			debrujin_insert (d, line);
		}

	if (errno)
		log_errno_fatal ("Failed to read '%s'",
				file);

	xfree (line);
	xfclose (fp);
}

void
print_k_mers_depth (DeBrujin *d)
{
	AdjList *adjlist = NULL;
	DeBrujinVetex *v = NULL;
	DeBrujinVetex *u = NULL;
	ListElmt *cur = NULL;
	GraphIter iter;

	graph_iter_init (&iter, d->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			v = adjlist->vertex;

			cur = list_head (adjlist->adjacent);
			for (; cur != NULL; cur = list_next (cur))
				{
					u = list_data (cur);
					log_debug ("K-MER: %s => %s: %d",
							v->k_mer_affix, u->k_mer_affix, u->depth);
				}
		}
}

void
rem_low_coverage_edges (DeBrujin *d, int min_cov)
{
	DeBrujinVetex *v = NULL;
	DeBrujinVetex *u = NULL;
	DeBrujinVetex *u2 = NULL;
	AdjList *adjlist = NULL;
	AdjList *adjlist2 = NULL;
	ListElmt *cur = NULL;
	List *to_rm = NULL;
	GraphIter iter;

	to_rm = list_new (NULL);

	log_debug ("Find low coverage edges");

	graph_iter_init (&iter, d->graph);
	while (graph_iter_next (&iter, &adjlist))
		{
			v = adjlist->vertex;
			cur = list_head (adjlist->adjacent);
			while (cur != NULL)
				{
					u = list_data (cur);
					cur = list_next (cur);

					if (u->depth < min_cov)
						{
							if (graph_rem_edge (d->graph, v, (void **) &u))
								{
									adjlist2 = graph_adjlist (d->graph, u);
									u2 = adjlist2->vertex;

									v->out_degree--;
									u2->in_degree--;

									if (v->out_degree == 0 && v->in_degree == 0)
										list_append (to_rm, v);

									if (u2->out_degree == 0 && u2->in_degree == 0)
										list_append (to_rm, u2);

									log_trace ("Remove edge between %s <=> %s",
											v->k_mer_affix, u->k_mer_affix);

									xfree (u);
								}
						}
				}
		}

	log_debug ("Remove low coverage vertices");

	for (cur = list_head (to_rm); cur != NULL; cur = list_next (cur))
		{
			v = list_data (cur);

			if (graph_rem_vertex (d->graph, (void **) &v))
				log_trace ("Remove vertex: %s", v->k_mer_affix);

			// LEAK:
			debrujin_vertex_free (v);
		}

	list_free (to_rm);
}

void
route_inspection (DeBrujin *d)
{
	Array *neg = NULL;
	Array *pos = NULL;
	GraphIter iter = {};
	DeBrujinVetex *n = NULL;
	DeBrujinVetex *p = NULL;
	DeBrujinVetex *a = NULL;
	DeBrujinVetex *b = NULL;
	Hungarian *h = NULL;
	AdjList *l = NULL;
	List *path = NULL;
	ListElmt *cur = NULL;
	ListElmt *prev = NULL;
	double **cost = NULL;
	int **assig_res = NULL;
	List ***paths = NULL;
	int i, j;

	if (debrujin_has_eulerian_path (d))
		{
			log_info ("UFA! Graph has eulerian path");
			return;
		}

	log_info ("Graph has no eulerian path. Let's fix");

	log_info ("Get all odd degree vertices");

	neg = array_new (NULL);
	pos = array_new (NULL);

	graph_iter_init (&iter, d->graph);

	while (graph_iter_next (&iter, &l))
		{
			a = l->vertex;
			if (a->out_degree - a->in_degree > 0 && a->in_degree > 0)
				{
					log_debug ("pos: %s: %d", a->k_mer_affix,
							a->out_degree - a->in_degree);
					for (i = 0; i < abs (a->out_degree - a->in_degree); i++)
						array_add (pos, a);
				}
			else if (a->out_degree - a->in_degree < 0 && a->out_degree > 0)
				{
					log_debug ("neg: %s: %d", a->k_mer_affix,
							a->out_degree - a->in_degree);
					for (i = 0; i < abs (a->out_degree - a->in_degree); i++)
						array_add (neg, a);
				}
		}

	// Check for zero values!
	cost = MTX_NEW (array_len (neg), array_len (pos), double);
	paths = MTX_NEW (array_len (neg), array_len (pos), List *);

	log_info ("Build assignment matrix with the shortest paths");

	for (i = 0; i < array_len (neg); i++)
		{
			n = array_get (neg, i);
			for (j = 0; j < array_len (pos); j++)
				{
					p = array_get (pos, j);
					paths[i][j] = debrujin_shortest_path (d, n, p);

					if (paths[i][j] != NULL)
						cost[i][j] = p->dist;
					else
						cost[i][j] = INFINITY;
				}
		}

	log_info ("Solve minimum assignment problem with hungarian method");

	h = hungarian_new (cost, array_len (neg),
			array_len (pos), HUNGARIAN_MODE_MINIMIZE_COST);

	hungarian_solve (h);
	assig_res = hungarian_assignment (h);

	log_info ("Fix edges and make the graph eulerian");

	for (i = 0; i < array_len (neg); i++)
		{
			n = array_get (neg, i);
			for (j = 0; j < array_len (pos); j++)
				{
					p = array_get (pos, j);
					if (assig_res[i][j])
						{
							log_debug ("Find path at %d,%d for %s -> %s",
									i, j, n->k_mer_affix, p->k_mer_affix);

							path = paths[i][j];

							cur = list_head (path);
							prev = cur;

							for (cur = list_next (cur); cur != NULL; cur = list_next (cur))
								{
									a = list_data (prev);
									b = list_data (cur);
									debrujin_insert_multi_edge (d, a, b);
									prev = cur;
								}
						}
				}
		}

	if (!debrujin_has_eulerian_path (d))
		log_fatal ("After all, there is no eulerian path!");

	for (i = 0; i < array_len (neg); i++)
		for (j = 0; j < array_len (pos); j++)
			list_free (paths[i][j]);

	MTX_FREE (cost, array_len (neg));
	MTX_FREE (paths, array_len (neg));

	array_free (pos, 1);
	array_free (neg, 1);
	hungarian_free (h);
}

int
main (int argc, char *argv[])
{
	int k = K_MER;
	int min_cov = MIN_COV;
	int debug = 0;
	int silent = 0;
	const char *file = NULL;

	int option_index = 0;
	int c;

	DeBrujin *d = NULL;
	List *seqs = NULL;
	ListElmt *cur = NULL;

	if (argc == 1)
		{
			print_usage (stdout);
			return EXIT_SUCCESS;
		}

	struct option opt[] =
	{
		{"help",      no_argument,       0, 'h'},
		{"quiet",     no_argument,       0, 'q'},
		{"debug",     no_argument,       0, 'd'},
		{"k-mer",     required_argument, 0, 'k'},
		{"min-cov",   required_argument, 0, 'c'}
	};

	while ((c = getopt_long (argc, argv, "hqdk:c:", opt, &option_index)) >= 0)
		{
			switch (c)
				{
				case 'h':
					{
						print_usage (stdout);
						return EXIT_SUCCESS;
					}
				case 'q':
					{
						silent = 1;
						break;
					}
				case 'd':
					{
						debug ++;
						break;
					}
				case 'k':
					{
						k = atoi (optarg);
						break;
					}
				case 'c':
					{
						min_cov = atoi (optarg);
						break;
					}
				case '?':
				case ':':
					{
						return EXIT_FAILURE;
					}
				}
		}

	if (silent)
		log_set_quiet (1);
	else
		{
			if (debug == 1)
				log_set_level (LOG_DEBUG);
			else if (debug > 1)
				log_set_level (LOG_TRACE);
			else
				log_set_level (LOG_INFO);
		}

	file = argv[optind];

	log_info ("Build debrujin obj and parse file '%s'", file);
	d = debrujin_new (k);
	parse_file_k_mers (file, d);

	log_info ("Remove low coverage k-mers");
	rem_low_coverage_edges (d, min_cov);
	/*print_k_mers_depth (d);*/

	log_info ("Check if graph is eulerian");
	route_inspection (d);

	log_info ("Assemble k-mers");
	seqs = debrujin_assembly (d);

	if (seqs != NULL)
		{
			cur = list_head (seqs);
			for (; cur != NULL; cur = list_next (cur))
				printf ("CONTIG: %s\n", (char *) list_data (cur));
		}

	list_free (seqs);
	debrujin_free (d);

	return EXIT_SUCCESS;
}
