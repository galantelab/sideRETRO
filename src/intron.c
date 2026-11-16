/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br>
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
#include "log.h"
#include "ibitree.h"
#include "intron.h"

static Intron *
intron_new (const char *gene_name, const char *chr, long start,
		long end, int exon_id_left, int exon_id_right)
{
	Intron *intron = xcalloc (1, sizeof (Intron));

	*intron = (Intron) {
		.gene_name     = xstrdup (gene_name),
		.chr           = xstrdup (chr),
		.start         = start,
		.end           = end,
		.exon_id_left  = exon_id_left,
		.exon_id_right = exon_id_right
	};

	return intron;
}

static void
intron_free (Intron *intron)
{
	if (intron == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	xfree (intron->gene_name);
	xfree (intron->chr);
	xfree (intron);
}

static void
intron_tree_insert (IntronTree *intron_tree, int retrocopy_id,
		Intron *intron)
{
	Hash *by_chr = NULL;
	IBiTree *tree = NULL;
	int *id_alloc = NULL;

	// RETROCOPY_ID => CHR => TREE
	by_chr = hash_lookup (intron_tree->idx, &retrocopy_id);

	if (by_chr == NULL)
		{
			by_chr = hash_new_full (str_hash, str_equal, NULL,
					(DestroyNotify) ibitree_free);

			id_alloc = xcalloc (1, sizeof (int));
			*id_alloc = retrocopy_id;

			hash_insert (intron_tree->idx, id_alloc, by_chr);
		}

	tree = hash_lookup (by_chr, intron->chr);

	if (tree == NULL)
		{
			tree = ibitree_new (NULL);
			hash_insert (by_chr, intron->chr, tree);
		}

	ibitree_insert (tree, intron->start, intron->end, intron);
	list_append (intron_tree->introns, intron);

	log_debug ("Retrocopy %d: Intron of '%s' at %s:%li-%li, between exons %d and %d",
			retrocopy_id, intron->gene_name, intron->chr, intron->start,
			intron->end, intron->exon_id_left, intron->exon_id_right);
}

static void
intron_tree_index (IntronTree *intron_tree, sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	Intron *intron = NULL;

	char *gene_name_prev = NULL;
	char *chr_prev = NULL;
	const char *gene_name = NULL;
	const char *chr = NULL;

	int retrocopy_id = 0;
	int retrocopy_id_prev = -1;
	int exon_id = 0;
	int exon_id_left = 0;
	long start = 0;
	long end = 0;
	long max_end = 0;

	/*
	* The exons of the parental genes of each retrocopy,
	* sorted by their start. Walk through them keeping the
	* greatest end so far - the union of the exons: Whenever
	* an exon starts after it, there is a gap, an intron.
	* Its left exon is the first one to reach the greatest
	* end, and its right exon is the longest one starting
	* after the gap
	*/
	const char sql[] =
		"SELECT DISTINCT cm.retrocopy_id, e.gene_name, e.chr,\n"
		"	e.start, e.end, e.id\n"
		"FROM cluster_merging AS cm\n"
		"INNER JOIN cluster AS cl\n"
		"	ON cl.id = cm.cluster_id AND cl.sid = cm.cluster_sid\n"
		"INNER JOIN exon AS e\n"
		"	ON e.gene_name = cl.gene_name\n"
		"ORDER BY cm.retrocopy_id ASC, e.gene_name ASC, e.chr ASC,\n"
		"	e.start ASC, e.end DESC, e.id ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	while (db_step (stmt) == SQLITE_ROW)
		{
			retrocopy_id = db_column_int   (stmt, 0);
			gene_name    = db_column_text  (stmt, 1);
			chr          = db_column_text  (stmt, 2);
			start        = db_column_int64 (stmt, 3);
			end          = db_column_int64 (stmt, 4);
			exon_id      = db_column_int   (stmt, 5);

			// Another gene: Start the union again
			if (retrocopy_id != retrocopy_id_prev
					|| strcmp (gene_name, gene_name_prev)
					|| strcmp (chr, chr_prev))
				{
					xfree (gene_name_prev);
					xfree (chr_prev);

					retrocopy_id_prev = retrocopy_id;
					gene_name_prev = xstrdup (gene_name);
					chr_prev = xstrdup (chr);

					exon_id_left = exon_id;
					max_end = end;

					continue;
				}

			if (start > max_end + 1)
				{
					intron = intron_new (gene_name, chr, max_end + 1,
							start - 1, exon_id_left, exon_id);

					intron_tree_insert (intron_tree, retrocopy_id, intron);
				}

			if (end > max_end)
				{
					exon_id_left = exon_id;
					max_end = end;
				}
		}

	xfree (gene_name_prev);
	xfree (chr_prev);
	db_finalize (stmt);
}

IntronTree *
intron_tree_new (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	IntronTree *intron_tree = xcalloc (1, sizeof (IntronTree));

	intron_tree->idx = hash_new_full (int_hash, int_equal, xfree,
			(DestroyNotify) hash_free);

	intron_tree->introns = list_new ((DestroyNotify) intron_free);

	intron_tree_index (intron_tree, db);

	return intron_tree;
}

void
intron_tree_free (IntronTree *intron_tree)
{
	if (intron_tree == NULL)
		return;

	hash_free (intron_tree->idx);
	list_free (intron_tree->introns);

	xfree (intron_tree);
}

static void
catch_if_inside_intron (IBiTreeLookupData *ldata, void *user_data)
{
	const Intron **intron = user_data;

	// The interval must be entirely inside the intron: An
	// alignment reaching an exon is not an intronic one
	if (ldata->interval_low >= ldata->node_low
			&& ldata->interval_high <= ldata->node_high)
		*intron = ldata->data;
}

const Intron *
intron_tree_lookup (const IntronTree *intron_tree, int retrocopy_id,
		const char *chr, long low, long high)
{
	assert (intron_tree != NULL && chr != NULL && low <= high);

	Hash *by_chr = NULL;
	IBiTree *tree = NULL;
	const Intron *intron = NULL;

	by_chr = hash_lookup (intron_tree->idx, &retrocopy_id);

	if (by_chr == NULL)
		return NULL;

	tree = hash_lookup (by_chr, chr);

	if (tree == NULL)
		return NULL;

	ibitree_lookup (tree, low, high, -1, -1, 0,
			catch_if_inside_intron, &intron);

	return intron;
}

void
intron_dump (sqlite3_stmt *intron_stmt, const Intron *intron,
		int retrocopy_id, int source_id, int reads)
{
	assert (intron_stmt != NULL && intron != NULL);

	log_debug ("Retrocopy %d, source %d: %d reads with the mate at the "
			"intron %s:%li-%li of '%s'", retrocopy_id, source_id, reads,
			intron->chr, intron->start, intron->end, intron->gene_name);

	db_insert_intron_retention (intron_stmt, retrocopy_id, source_id,
			intron->gene_name, intron->chr, intron->start, intron->end,
			intron->exon_id_left, intron->exon_id_right, reads);
}
