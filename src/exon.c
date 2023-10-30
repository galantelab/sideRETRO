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

#include <inttypes.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "gff.h"
#include "ibitree.h"
#include "chr.h"
#include "exon.h"

struct _ExonTreeData
{
	ExonTree *tree;
	int64_t   alignment_id;
};

typedef struct _ExonTreeData ExonTreeData;

ExonTree *
exon_tree_new (sqlite3_stmt *exon_stmt,
		sqlite3_stmt *overlapping_stmt, ChrStd *cs)
{
	assert (exon_stmt != NULL && overlapping_stmt != NULL
			&& cs != NULL);

	ExonTree *exon_tree = xcalloc (1, sizeof (ExonTree));

	exon_tree->exon_stmt = exon_stmt;
	exon_tree->overlapping_stmt = overlapping_stmt;

	exon_tree->cs = cs;

	exon_tree->idx = hash_new (xfree,
			(DestroyNotify) ibitree_free);

	exon_tree->cache = hash_new (xfree, NULL);

	return exon_tree;
}

void
exon_tree_free (ExonTree *exon_tree)
{
	if (exon_tree == NULL)
		return;

	hash_free (exon_tree->idx);
	hash_free (exon_tree->cache);

	xfree (exon_tree);
}

void
exon_tree_index_dump (ExonTree *exon_tree,
		const char *gff_file)
{
	assert (exon_tree != NULL && gff_file != NULL);

	GffFile *gff = gff_open_for_reading (gff_file);
	GffEntry *entry = gff_entry_new ();
	GffFilter *filter = gff_filter_new ();

	gff_filter_insert_feature (filter, "exon");
	gff_filter_insert_hard_attribute (filter,
			"transcript_type", "protein_coding");

	IBiTree *tree = NULL;

	int64_t table_id = 0;
	int64_t *alloc_id = NULL;

	const char *chr_std = NULL;
	const char *gene_name = NULL;
	const char *gene_id = NULL;
	const char *exon_id = NULL;
	const char *exon_id_copy = NULL;
	char strand[1];

	while (gff_read_filtered (gff, entry, filter))
		{
			gene_name = gff_attribute_find (entry, "gene_name");
			gene_id = gff_attribute_find (entry, "gene_id");
			exon_id = gff_attribute_find (entry, "exon_id");

			if (gene_name == NULL || gene_id == NULL || exon_id == NULL)
				{
					log_warn ("Missing gene_name|gene_id|exon_id at line %zu",
						entry->num_line);
					continue;
				}

			if (hash_contains (exon_tree->cache, exon_id))
				{
					continue;
				}
			else
				{
					exon_id_copy = xstrdup (exon_id);
					hash_insert (exon_tree->cache, exon_id_copy,
							exon_id_copy);
				}

			chr_std = chr_std_lookup (exon_tree->cs, entry->seqname);

			log_debug ("Index exon from gene '%s' at %s:%zu-%zu", gene_name,
					chr_std, entry->start, entry->end);

			strand[0] = entry->strand;
			alloc_id = xcalloc (1, sizeof (int64_t));
			* (int64_t *) alloc_id = ++table_id;

			tree = hash_lookup (exon_tree->idx, chr_std);

			if (tree == NULL)
				{
					tree = ibitree_new (xfree);
					hash_insert (exon_tree->idx,
							xstrdup (chr_std), tree);
				}

			ibitree_insert (tree, entry->start, entry->end,
					alloc_id);

			db_insert_exon (exon_tree->exon_stmt,
					table_id, gene_name, chr_std, entry->start,
					entry->end, strand, gene_id, exon_id);
		}

	gff_filter_free (filter);
	gff_entry_free (entry);
	gff_close (gff);
}

static void
dump_if_overlaps_exon (IBiTreeLookupData *ldata,
		void *user_data)
{
	const int64_t *exon_id = ldata->data;
	ExonTreeData *data = user_data;

	log_debug ("Dump overlapping exon [%" PRId64 "] %li-%li with alignment [%" PRId64 "] %li-%li at %li-%li",
			*exon_id, ldata->node_low, ldata->node_high, data->alignment_id,
			ldata->interval_low, ldata->interval_high, ldata->overlap_pos,
			ldata->overlap_pos + ldata->overlap_len - 1);

	db_insert_overlapping (data->tree->overlapping_stmt, *exon_id,
			data->alignment_id, ldata->overlap_pos, ldata->overlap_len);
}

int
exon_tree_lookup_dump (ExonTree *exon_tree, const char *chr,
		long low, long high, float exon_overlap_frac,
		float alignment_overlap_frac, int either,
		int64_t alignment_id)
{
	assert (exon_tree != NULL && chr != NULL);

	int acm = 0;
	IBiTree *tree = NULL;

	tree = hash_lookup (exon_tree->idx, chr);

	if (tree != NULL)
		{
			ExonTreeData data = {exon_tree, alignment_id};
			acm = ibitree_lookup (tree, low, high, exon_overlap_frac,
					alignment_overlap_frac, either, dump_if_overlaps_exon,
					&data);
		}

	return acm;
}
