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
#include <sys/stat.h>
#include "wrapper.h"
#include "utils.h"
#include "log.h"
#include "hash.h"
#include "db_merge.h"

#define NUM_TABLES 5

enum _Tables
{
	BATCH,
	SOURCE,
	EXON,
	ALIGNMENT,
	OVERLAPPING
};

static const char *tables[NUM_TABLES] =
{
	"batch",
	"source",
	"exon",
	"alignment",
	"overlapping"
};

static const char *sel_sql[NUM_TABLES] =
{
	"SELECT * FROM batch",
	"SELECT * FROM source",
	"SELECT * FROM exon",
	"SELECT * FROM alignment",
	"SELECT * FROM overlapping",
};

static const char *max_id_sql[NUM_TABLES] =
{
	"SELECT MAX(id) FROM batch",
	"SELECT MAX(id) FROM source",
	"SELECT MAX(id) FROM exon",
	"SELECT MAX(id) FROM alignment",
	NULL
};

static Hash *ense_h = NULL;
static Hash *exon_id_h = NULL;

static void
merge_batch (sqlite3_stmt *in_stmt, sqlite3_stmt *sel_stmt,
		int max_id[NUM_TABLES])
{
	int id = 0;
	const char *timestamp = NULL;

	db_reset (sel_stmt);

	while (db_step (sel_stmt) == SQLITE_ROW)
		{
			// Catch all values from db2
			id = db_column_int (sel_stmt, 0);
			timestamp = db_column_text (sel_stmt, 1);

			// Insert them into db
			db_insert_batch (in_stmt, id + max_id[BATCH],
					timestamp);
		}
}

static void
merge_souce (sqlite3_stmt *in_stmt, sqlite3_stmt *sel_stmt,
		int max_id[NUM_TABLES])
{
	int id = 0;
	int batch_id = 0;
	const char *path = NULL;

	db_reset (sel_stmt);

	while (db_step (sel_stmt) == SQLITE_ROW)
		{
			// Catch all values from db2
			id = db_column_int (sel_stmt, 0);
			batch_id = db_column_int (sel_stmt, 1);
			path = db_column_text (sel_stmt, 2);

			// Insert them into db
			db_insert_source (in_stmt, id + max_id[SOURCE],
					batch_id + max_id[BATCH], path);
		}
}

static void
merge_exon (sqlite3_stmt *in_stmt, sqlite3_stmt *sel_stmt,
		int max_id[NUM_TABLES])
{
	int id = 0;
	const char *gene_name = NULL;
	const char *chr = NULL;
	long start = 0;
	long end = 0;
	const char *strand = NULL;
	const char *ensg = NULL;
	const char *ense = NULL;

	const int *exon_id = NULL;
	int *id_copy = NULL;
	int *merged_id = NULL;

	db_reset (sel_stmt);

	while (db_step (sel_stmt) == SQLITE_ROW)
		{
			// Catch all values from db2
			id = db_column_int (sel_stmt, 0);
			gene_name = db_column_text (sel_stmt, 1);
			chr = db_column_text (sel_stmt, 2);
			start = db_column_int64 (sel_stmt, 3);
			end = db_column_int64 (sel_stmt, 4);
			strand = db_column_text (sel_stmt, 5);
			ensg = db_column_text (sel_stmt, 6);
			ense = db_column_text (sel_stmt, 7);

			/*
			* exon_id_h maps the exon id at this database to the
			* exon id at the merged database. An exon already
			* merged (the same ense) keeps its id, a new one is
			* shifted by max_id
			*/
			id_copy = xcalloc (1, sizeof (int));
			*id_copy = id;

			exon_id = hash_lookup (ense_h, ense);

			if (exon_id != NULL)
				{
					hash_insert (exon_id_h, id_copy, exon_id);
					continue;
				}

			merged_id = xcalloc (1, sizeof (int));
			*merged_id = id + max_id[EXON];

			hash_insert (ense_h, xstrdup (ense), merged_id);
			hash_insert (exon_id_h, id_copy, merged_id);

			// Insert them into db
			db_insert_exon (in_stmt, *merged_id, gene_name,
					chr, start, end, strand, ensg, ense);
		}
}

static void
merge_alignment (sqlite3_stmt *in_stmt, sqlite3_stmt *sel_stmt,
		int max_id[NUM_TABLES])
{
	db_reset (sel_stmt);

	int id = 0;
	const char *qname = NULL;
	int flag = 0;
	const char *chr = NULL;
	long pos = 0;
	int mapq = 0;
	const char *cigar = NULL;
	int qlen = 0;
	int rlen = 0;
	const char *chr_next = NULL;
	long pos_next = 0;
	int type = 0;
	int source_id = 0;
	const char *seq = NULL;
	const char *qual = NULL;

	while (db_step (sel_stmt) == SQLITE_ROW)
		{
			id = db_column_int (sel_stmt, 0);
			qname = db_column_text (sel_stmt, 1);
			flag = db_column_int (sel_stmt, 2);
			chr = db_column_text (sel_stmt, 3);
			pos = db_column_int64 (sel_stmt, 4);
			mapq = db_column_int (sel_stmt, 5);
			cigar = db_column_text (sel_stmt, 6);
			qlen = db_column_int (sel_stmt, 7);
			rlen = db_column_int (sel_stmt, 8);
			chr_next = db_column_text (sel_stmt, 9);
			pos_next = db_column_int64 (sel_stmt, 10);
			type = db_column_int (sel_stmt, 11);
			source_id = db_column_int (sel_stmt, 12);
			seq = db_column_text (sel_stmt, 13);
			qual = db_column_text (sel_stmt, 14);

			db_insert_alignment (in_stmt, id + max_id[ALIGNMENT], qname,
					flag, chr, pos, mapq, cigar, qlen, rlen, chr_next,
					pos_next, type, source_id + max_id[SOURCE], seq, qual);
		}
}

static void
merge_overlapping (sqlite3_stmt *in_stmt, sqlite3_stmt *sel_stmt,
		int max_id[NUM_TABLES])
{
	db_reset (sel_stmt);

	int exon_id = 0;
	int alignment_id = 0;
	long pos = 0;
	long len = 0;

	int *exon_id_align = NULL;

	while (db_step (sel_stmt) == SQLITE_ROW)
		{
			exon_id = db_column_int (sel_stmt, 0);
			alignment_id = db_column_int (sel_stmt, 1);
			pos = db_column_int64 (sel_stmt, 2);
			len = db_column_int64 (sel_stmt, 3);

			exon_id_align = hash_lookup (exon_id_h, &exon_id);
			assert (exon_id_align != NULL);

			exon_id = *exon_id_align;

			db_insert_overlapping (in_stmt, exon_id,
					alignment_id + max_id[ALIGNMENT],
					pos, len);
		}
}

static void
(*insert[NUM_TABLES]) (sqlite3_stmt *,
		sqlite3_stmt *, int *) =
{
	merge_batch,
	merge_souce,
	merge_exon,
	merge_alignment,
	merge_overlapping
};

static sqlite3_stmt *
(*prepare[NUM_TABLES])(sqlite3 *) =
{
	db_prepare_batch_stmt,
	db_prepare_source_stmt,
	db_prepare_exon_stmt,
	db_prepare_alignment_stmt,
	db_prepare_overlapping_stmt
};

static void
calc_max_id (sqlite3_stmt *max_id_stmt[NUM_TABLES],
		int max_id[NUM_TABLES])
{
	int i = 0;

	for (; i < NUM_TABLES; i++)
		{
			if (max_id_stmt[i] == NULL)
				continue;

			db_reset (max_id_stmt[i]);

			if (db_step (max_id_stmt[i]) == SQLITE_ROW
					&& sqlite3_column_type (max_id_stmt[i], 0) != SQLITE_NULL)
				max_id[i] = db_column_int (max_id_stmt[i], 0);
			else
				max_id[i] = 0;
		}
}

void
exon_id_init (sqlite3 *db)
{
	sqlite3_stmt *sel_stmt = NULL;
	int id = 0;
	int *id_copy = NULL;
	const char *ense = NULL;

	// Prepare select stmt
	sel_stmt = db_prepare (db, sel_sql[EXON]);

	// Init hash ense => exon_id at the merged database.
	// It owns the ids
	ense_h = hash_new_full (str_hash, str_equal, xfree, xfree);

	db_reset (sel_stmt);

	while (db_step (sel_stmt) == SQLITE_ROW)
		{
			// Get database values
			id = db_column_int (sel_stmt, 0);
			ense = db_column_text (sel_stmt, 7);

			id_copy = xcalloc (1, sizeof (int));
			*id_copy = id;

			// Index to the new hash
			hash_insert (ense_h, xstrdup (ense), id_copy);
		}

	db_finalize (sel_stmt);
}

static int
same_file (const char *path1, const char *path2)
{
	struct stat s1;
	struct stat s2;

	// The same file may be reached by different paths
	if (path1 == NULL || path2 == NULL
			|| stat (path1, &s1) < 0 || stat (path2, &s2) < 0)
		return 0;

	return s1.st_dev == s2.st_dev && s1.st_ino == s2.st_ino;
}

static long
count_exons (sqlite3 *db)
{
	sqlite3_stmt *stmt = NULL;
	long n = 0;

	stmt = db_prepare (db, "SELECT COUNT(*) FROM exon");

	if (db_step (stmt) == SQLITE_ROW)
		n = db_column_int64 (stmt, 0);

	db_finalize (stmt);
	return n;
}

static void
index_sources (sqlite3 *db, const char *db_file, Hash *path_h,
		Hash *name_h)
{
	sqlite3_stmt *stmt = NULL;
	const char *path = NULL;
	const char *prev = NULL;
	char *name = NULL;

	stmt = db_prepare (db, "SELECT path FROM source");

	while (db_step (stmt) == SQLITE_ROW)
		{
			path = db_column_text (stmt, 0);

			// The same sample twice: Its reads would be counted
			// twice and it would be twice at the VCF
			prev = hash_lookup (path_h, path);
			if (prev != NULL)
				log_fatal ("The sample '%s' is in both '%s' and '%s': "
						"Merge each sample only once", path, prev, db_file);

			hash_insert (path_h, xstrdup (path), xstrdup (db_file));

			// The VCF names the samples by their file names
			name = path_file (path, 1);
			prev = hash_lookup (name_h, name);

			if (prev != NULL)
				{
					log_warn ("The samples '%s' and '%s' have the same name '%s' "
							"at the VCF", prev, path, name);
					xfree (name);
				}
			else
				hash_insert (name_h, name, xstrdup (path));
		}

	db_finalize (stmt);
}

static void
check_sources (sqlite3 *db, const char *db_file, int argc, char **argv,
		const int *skip)
{
	sqlite3 *db2 = NULL;
	Hash *path_h = NULL;
	Hash *name_h = NULL;
	int i = 0;

	path_h = hash_new_full (str_hash, str_equal, xfree, xfree);
	name_h = hash_new_full (str_hash, str_equal, xfree, xfree);

	// The samples already at the output database (--in-place)
	index_sources (db, db_file, path_h, name_h);

	for (i = 0; i < argc; i++)
		{
			if (skip[i])
				continue;

			db2 = db_connect (argv[i]);
			index_sources (db2, argv[i], path_h, name_h);
			db_close (db2);
		}

	hash_free (path_h);
	hash_free (name_h);
}

/*
 * Merge the databases at argv into db. The same database
 * file is merged only once. Returns the number of merged
 * databases
 */
int
db_merge (sqlite3 *db, int argc, char **argv)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && argc > 0 && argv != NULL);

	sqlite3 *db2 = NULL;

	sqlite3_stmt *in_stmt[NUM_TABLES];
	sqlite3_stmt *sel_stmt[NUM_TABLES];
	sqlite3_stmt *max_id_stmt[NUM_TABLES];

	int max_id[NUM_TABLES] = {};
	int *skip = NULL;
	const char *db_file = NULL;
	long exons_before = 0;
	long exons_db2 = 0;
	long exons_after = 0;
	int merged = 0;
	int i = 0;
	int j = 0;

	// The output database: Empty or --in-place
	db_file = sqlite3_db_filename (db, "main");

	/*
	* The same database file twice - by different paths, or the
	* output database of --in-place again: Merge it only once
	*/
	skip = xcalloc (argc, sizeof (int));

	for (i = 0; i < argc; i++)
		{
			if (same_file (argv[i], db_file))
				{
					log_warn ("'%s' is the output database itself: Skip it", argv[i]);
					skip[i] = 1;
					continue;
				}

			for (j = 0; j < i; j++)
				{
					if (!skip[j] && same_file (argv[i], argv[j]))
						{
							log_warn ("'%s' is the same database as '%s': Merge it only once",
									argv[i], argv[j]);
							skip[i] = 1;
							break;
						}
				}
		}

	// Before writing anything: Each sample only once
	check_sources (db, db_file, argc, argv, skip);

	// Init insert and max_id statement
	for (i = 0; i < NUM_TABLES; i++)
		{
			in_stmt[i] = (*prepare[i])(db);
			max_id_stmt[i] = max_id_sql[i] != NULL
				? db_prepare (db, max_id_sql[i])
				: NULL;
		}

	// Init ense_h and exon_id_h and
	// index the values already present
	// in the database
	exon_id_init (db);

	// Merge all files at argv
	for (i = 0; i < argc; i++)
		{
			if (skip[i])
				continue;

			log_info ("Merge database '%s'", argv[i]);
			db2 = db_connect (argv[i]);

			exons_before = count_exons (db);
			exons_db2 = count_exons (db2);

			// Calculate de max id for all tables
			// to be merged
			calc_max_id (max_id_stmt, max_id);

			// The exon ids are local to each database: Map
			// them to the merged database ids again
			hash_free (exon_id_h);
			exon_id_h = hash_new_full (int_hash, int_equal, xfree, NULL);

			for (j = 0; j < NUM_TABLES; j++)
				{
					log_debug ("Merging table '%s' from database '%s'",
							tables[j], argv[i]);

					// Prepare 'select' - just one time for
					// this database
					sel_stmt[j] = db_prepare (db2, sel_sql[j]);

					// Insert table from db2 to db1 and correct the
					// shift for id with max_id. This way, I can keep
					// the foreign key relation
					(*insert[j]) (in_stmt[j], sel_stmt[j], max_id);

					// Clean this select statement
					db_finalize (sel_stmt[j]);
				}

			/*
			* All databases must come from the same annotation:
			* Otherwise, the reads at the exons missing from one
			* database were not marked as exonic in its samples
			*/
			exons_after = count_exons (db);

			if (exons_before > 0
					&& (exons_after != exons_before || exons_db2 != exons_after))
				log_warn ("The annotation of '%s' (%li exons) differs from the "
						"one of the other databases (%li exons): The samples may "
						"not be comparable. Process all samples with the same "
						"annotation", argv[i], exons_db2, exons_before);

			// Goodbye
			db_close (db2);
			merged++;
		}

	// Cleanup the rest of the mess
	for (i = 0; i < NUM_TABLES; i++)
		{
			if (in_stmt[i] != NULL)
				db_finalize (in_stmt[i]);
			if (max_id_stmt[i] != NULL)
				db_finalize (max_id_stmt[i]);
		}

	hash_free (ense_h);
	hash_free (exon_id_h);
	xfree (skip);

	return merged;
}
