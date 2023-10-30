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
#include <regex.h>
#include "wrapper.h"
#include "log.h"
#include "gff.h"
#include "bed.h"
#include "ibitree.h"
#include "str.h"
#include "blacklist.h"

struct _BlacklistData
{
	Blacklist *blacklist;
	int64_t    cluster_id;
	int64_t    cluster_sid;
};

typedef struct _BlacklistData BlacklistData;

Blacklist *
blacklist_new (sqlite3_stmt *blacklist_stmt,
		sqlite3_stmt *overlapping_blacklist_stmt,
		ChrStd *cs)
{
	assert (blacklist_stmt != NULL
			&& overlapping_blacklist_stmt != NULL
			&& cs != NULL);

	Blacklist *blacklist = xcalloc (1, sizeof (Blacklist));

	blacklist->blacklist_stmt = blacklist_stmt;
	blacklist->overlapping_blacklist_stmt =
		overlapping_blacklist_stmt;

	blacklist->cs = cs;

	blacklist->idx = hash_new (xfree,
			(DestroyNotify) ibitree_free);

	return blacklist;
}

void
blacklist_free (Blacklist *blacklist)
{
	if (blacklist == NULL)
		return;

	hash_free (blacklist->idx);
	xfree (blacklist);
}

static void
clean_blacklist_tables (sqlite3 *db)
{
	// Delete all values from
	// previous runs
	const char sql[] =
		"DELETE FROM overlapping_blacklist;\n"
		"DELETE FROM blacklist;";

	log_debug ("Clean tables:\n%s", sql);
	db_exec (db, sql);
}

void
blacklist_index_dump_from_gff (Blacklist *blacklist, const char *gff_file, const GffFilter *filter)
{
	assert (blacklist != NULL && gff_file != NULL && filter != NULL);

	log_debug ("Clean blacklist tables");
	clean_blacklist_tables (
			sqlite3_db_handle (blacklist->blacklist_stmt));

	GffFile *gff = gff_open_for_reading (gff_file);
	GffEntry *entry = gff_entry_new ();

	IBiTree *tree = NULL;

	int64_t table_id = 0;
	int64_t *alloc_id = NULL;

	const char *chr_std = NULL;
	const char *gene_name = NULL;

	while (gff_read_filtered (gff, entry, filter))
		{
			gene_name = gff_attribute_find (entry, "gene_name");

			if (gene_name == NULL)
				gene_name = "blacklist";

			chr_std = chr_std_lookup (blacklist->cs, entry->seqname);

			log_debug ("Index blacklist '%s' at %s:%zu-%zu", gene_name,
					chr_std, entry->start, entry->end);

			alloc_id = xcalloc (1, sizeof (int64_t));
			* (int64_t *) alloc_id = ++table_id;

			tree = hash_lookup (blacklist->idx, chr_std);

			if (tree == NULL)
				{
					tree = ibitree_new (xfree);
					hash_insert (blacklist->idx,
							xstrdup (chr_std), tree);
				}

			ibitree_insert (tree, entry->start, entry->end,
					alloc_id);

			db_insert_blacklist (blacklist->blacklist_stmt,
					table_id, gene_name, chr_std, entry->start,
					entry->end);
		}

	gff_entry_free (entry);
	gff_close (gff);
}

void
blacklist_index_dump_from_bed (Blacklist *blacklist, const char *bed_file)
{
	assert (blacklist != NULL && bed_file != NULL);

	log_debug ("Clean blacklist tables");
	clean_blacklist_tables (
			sqlite3_db_handle (blacklist->blacklist_stmt));

	BedFile *bed = bed_open_for_reading (bed_file);
	BedEntry *entry = bed_entry_new ();

	IBiTree *tree = NULL;

	int64_t table_id = 0;
	int64_t *alloc_id = NULL;

	const char *chr_std = NULL;
	const char *name = NULL;

	while (bed_read (bed, entry))
		{
			name = entry->num_field > 3 && entry->name != NULL
				? entry->name
				: "blacklist";

			chr_std = chr_std_lookup (blacklist->cs, entry->chrom);

			log_debug ("Index blacklist '%s' at %s:%zu-%zu", name,
					chr_std, entry->chrom_start, entry->chrom_end);

			alloc_id = xcalloc (1, sizeof (int64_t));
			* (int64_t *) alloc_id = ++table_id;

			tree = hash_lookup (blacklist->idx, chr_std);

			if (tree == NULL)
				{
					tree = ibitree_new (xfree);
					hash_insert (blacklist->idx,
							xstrdup (chr_std), tree);
				}

			ibitree_insert (tree, entry->chrom_start,
					entry->chrom_end, alloc_id);

			db_insert_blacklist (blacklist->blacklist_stmt,
					table_id, name, chr_std, entry->chrom_start,
					entry->chrom_end);
		}

	bed_entry_free (entry);
	bed_close (bed);
}

static void
dump_if_overlaps_blacklist (IBiTreeLookupData *ldata,
		void *user_data)
{
	const int64_t *blacklist_id = ldata->data;
	BlacklistData *data = user_data;

	log_debug ("Dump overlapping blacklist [%" PRId64 "] %li-%li with cluster [%" PRId64 "] %li-%li at %li-%li",
			*blacklist_id, ldata->node_low, ldata->node_high, data->cluster_id,
			ldata->interval_low, ldata->interval_high, ldata->overlap_pos,
			ldata->overlap_pos + ldata->overlap_len - 1);

	db_insert_overlapping_blacklist (data->blacklist->overlapping_blacklist_stmt, *blacklist_id,
			data->cluster_id, data->cluster_sid, ldata->overlap_pos, ldata->overlap_len);
}

int
blacklist_lookup (Blacklist *blacklist, const char *chr,
		long low, long high, long padding, const int64_t cluster_id,
		const int64_t cluster_sid)
{
	assert (blacklist != NULL && chr != NULL
			&& padding >= 0);

	int acm = 0;
	IBiTree *tree = NULL;

	tree = hash_lookup (blacklist->idx, chr);

	if (tree != NULL)
		{
			low -= padding;
			BlacklistData data = {blacklist, cluster_id, cluster_sid};
			acm = ibitree_lookup (tree, low > 0 ? low : 0, high + padding,
					-1, -1, 0, dump_if_overlaps_blacklist, &data);
		}

	return acm;
}
