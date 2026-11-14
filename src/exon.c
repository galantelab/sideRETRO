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
#include "log.h"
#include <stdio.h>
#include "gff.h"
#include "bed.h"
#include "ibitree.h"
#include "chr.h"
#include "exon.h"

struct _ExonTreeData
{
	ExonTree *tree;
	long      alignment_id;
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

static void
exon_tree_insert (ExonTree *exon_tree, long *table_id, const char *seqname,
		long start, long end, char strand, const char *gene_name,
		const char *gene_id, const char *exon_id)
{
	IBiTree *tree = NULL;
	const char *chr_std = NULL;
	char *exon_id_copy = NULL;
	char strand_str[2] = {};
	long *alloc_id = NULL;

	// An exon shared by several transcripts is indexed once
	if (hash_contains (exon_tree->cache, exon_id))
		return;

	exon_id_copy = xstrdup (exon_id);
	hash_insert (exon_tree->cache, exon_id_copy, exon_id_copy);

	chr_std = chr_std_lookup (exon_tree->cs, seqname);

	log_debug ("Index exon '%s' from gene '%s' at %s:%li-%li", exon_id,
			gene_name, chr_std, start, end);

	strand_str[0] = strand;
	alloc_id = xcalloc (1, sizeof (long));
	*alloc_id = ++(*table_id);

	tree = hash_lookup (exon_tree->idx, chr_std);

	if (tree == NULL)
		{
			tree = ibitree_new (xfree);
			hash_insert (exon_tree->idx, xstrdup (chr_std), tree);
		}

	ibitree_insert (tree, start, end, alloc_id);

	db_insert_exon (exon_tree->exon_stmt, *table_id, gene_name, chr_std,
			start, end, strand_str, gene_id, exon_id);
}

static void
exon_tree_index_dump_from_gff (ExonTree *exon_tree, const char *gff_file,
		const GffFilter *filter, const char *gene_name_key, const char *exon_id_key)
{
	GffFile *gff = gff_open_for_reading (gff_file);
	GffEntry *entry = gff_entry_new ();

	const char *gene_name = NULL;
	const char *gene_id = NULL;
	const char *exon_id = NULL;
	char *exon_id_alloc = NULL;
	long table_id = 0;

	while (gff_read_filtered (gff, entry, filter))
		{
			// The gene name groups the exons into the parental
			// gene. Fall back to 'gene_id' if it is missing
			gene_name = gff_attribute_find (entry, gene_name_key);
			gene_id = gff_attribute_find (entry, "gene_id");

			if (gene_name == NULL)
				gene_name = gene_id;

			if (gene_name == NULL)
				{
					log_warn ("Missing '%s' and 'gene_id' at line %zu: Skip the exon",
							gene_name_key, entry->num_line);
					continue;
				}

			if (gene_id == NULL)
				gene_id = gene_name;

			// Without an exon id, use its coordinates
			exon_id = gff_attribute_find (entry, exon_id_key);

			if (exon_id == NULL)
				{
					xasprintf (&exon_id_alloc, "%s:%zu_%zu",
							entry->seqname, entry->start, entry->end);
					exon_id = exon_id_alloc;
				}

			exon_tree_insert (exon_tree, &table_id, entry->seqname,
					entry->start, entry->end, entry->strand,
					gene_name, gene_id, exon_id);

			xfree (exon_id_alloc);
			exon_id_alloc = NULL;
		}

	gff_entry_free (entry);
	gff_close (gff);
}

static void
exon_tree_index_dump_from_bed (ExonTree *exon_tree, const char *bed_file)
{
	BedFile *bed = bed_open_for_reading (bed_file);
	BedEntry *entry = bed_entry_new ();

	char *exon_id = NULL;
	long table_id = 0;
	long start = 0;

	while (bed_read (bed, entry))
		{
			// BED6: chrom, start, end, name (gene), score, strand
			if (entry->num_field < 6)
				log_fatal ("Annotation '%s' must be BED6 (chrom, start, end, gene name, "
						"score, strand): Only %zu fields at line %zu", bed_file,
						entry->num_field, entry->num_line);

			// BED start is 0-based
			start = entry->chrom_start + 1;

			xasprintf (&exon_id, "%s:%li_%li", entry->chrom, start,
					entry->chrom_end);

			exon_tree_insert (exon_tree, &table_id, entry->chrom,
					start, entry->chrom_end, entry->strand,
					entry->name, entry->name, exon_id);

			xfree (exon_id);
			exon_id = NULL;
		}

	bed_entry_free (entry);
	bed_close (bed);
}

void
exon_tree_index_dump (ExonTree *exon_tree, const char *annotation_file,
		const GffFilter *filter, const char *gene_name_key, const char *exon_id_key)
{
	assert (exon_tree != NULL && annotation_file != NULL);

	GffFilter *default_filter = NULL;

	if (bed_looks_like_bed_file (annotation_file))
		{
			log_debug ("Index exons from BED6 '%s'", annotation_file);
			exon_tree_index_dump_from_bed (exon_tree, annotation_file);

			if (hash_size (exon_tree->cache) == 0)
				log_fatal ("No exon found at '%s'", annotation_file);

			return;
		}

	// GENCODE defaults for what is not set
	if (filter == NULL)
		{
			default_filter = gff_filter_new ();
			gff_filter_insert_feature (default_filter, EXON_DEFAULT_FEATURE);
			gff_filter_insert_attribute (default_filter, EXON_DEFAULT_ATTRIBUTE);
			filter = default_filter;
		}

	if (gene_name_key == NULL)
		gene_name_key = EXON_DEFAULT_GENE_NAME;

	if (exon_id_key == NULL)
		exon_id_key = EXON_DEFAULT_EXON_ID;

	log_debug ("Index exons from GTF/GFF3 '%s'", annotation_file);
	exon_tree_index_dump_from_gff (exon_tree, annotation_file,
			filter, gene_name_key, exon_id_key);

	// Nothing can be found without exons: Probably the
	// annotation does not follow the expected attributes
	if (hash_size (exon_tree->cache) == 0)
		log_fatal ("No exon passed the filter at '%s'. Check the feature "
				"(--gff-feature) and the attributes (--gff-attribute): "
				"Ensembl, for instance, uses 'transcript_biotype' instead "
				"of 'transcript_type'", annotation_file);

	log_info ("Indexed %zu exons from '%s'", hash_size (exon_tree->cache),
			annotation_file);

	gff_filter_free (default_filter);
}

static void
dump_if_overlaps_exon (IBiTreeLookupData *ldata,
		void *user_data)
{
	const long *exon_id = ldata->data;
	ExonTreeData *data = user_data;

	log_debug ("Dump overlapping exon [%li] %li-%li with alignment [%li] %li-%li at %li-%li",
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
		long alignment_id)
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
