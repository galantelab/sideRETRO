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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/log.h"
#include "../src/hash.h"
#include "../src/ibitree.h"
#include "../src/utils.h"
#include "../src/db.h"
#include "../src/gff.h"
#include "../src/exon.h"

#define TEST_EXON_BUFSIZ 64

struct _TestExonTree
{
	char db_path[TEST_EXON_BUFSIZ];
	char gtf_path[TEST_EXON_BUFSIZ];

	sqlite3 *db;
	sqlite3_stmt *exon_stmt;
	sqlite3_stmt *overlapping_stmt;

	ChrStd *cs;

	ExonTree *exon_tree;
};

typedef struct _TestExonTree TestExonTree;

static const char *gtf =
	"chr1\t.\tgene\t1\t1000\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\";transcript_type \"protein_coding\";\n"
	"chr1\t.\ttranscript\t1\t1000\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\";transcript_type \"protein_coding\";\n"
	"chr1\t.\texon\t1\t100\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; exon_id \"e1\"\n"
	"chr1\t.\texon\t200\t300\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; exon_id \"e2\"\n"
	"chr1\t.\texon\t400\t500\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; exon_id \"e3\"\n"
	"chr1\t.\texon\t600\t700\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; exon_id \"e4\"\n"
	"chr1\t.\texon\t800\t900\t.\t+\t.\t"
	"gene_name \"ponga\"; gene_id \"g1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; exon_id \"e5\"\n";

static const char *gtf_id[] = {"e1", "e2", "e3", "e4", "e5"};

static const int gtf_pos[][2] = {{1, 100}, {200, 300}, {400, 500}, {600, 700}, {800, 900}};

static const int gtf_size = 5;

static sqlite3 *
create_db (char *db_path)
{
	int fd;

	fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
create_gtf (char *gtf_path)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (gtf_path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, gtf, "");

	xfclose (fp);
}

static void
catch_id (IBiTreeLookupData *ldata, void *id)
{
	* (int *) id = * (int *) ldata->data;
}

static sqlite3_stmt *
prepare_exon_search_stmt (sqlite3 *db)
{
	sqlite3_stmt *search_stmt = NULL;
	int rc = 0;

	sqlite3_prepare_v2 (db,
			"SELECT id,start,end,ense FROM exon ORDER BY id ASC",
			-1, &search_stmt, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to prepare search stmt: %s",
				sqlite3_errmsg (db));

	return search_stmt;
}

static void
test_exon_tree_init (TestExonTree *t)
{
	strncpy (t->db_path, "/tmp/ponga.db.XXXXXX",
			TEST_EXON_BUFSIZ - 1);
	strncpy (t->gtf_path, "/tmp/ponga.gtf.XXXXXX",
			TEST_EXON_BUFSIZ - 1);

	t->db_path[TEST_EXON_BUFSIZ - 1] = '\0';
	t->gtf_path[TEST_EXON_BUFSIZ - 1] = '\0';

	t->db = create_db (t->db_path);
	t->exon_stmt = db_prepare_exon_stmt (t->db);
	t->overlapping_stmt = db_prepare_overlapping_stmt (t->db);

	t->cs = chr_std_new ();

	create_gtf (t->gtf_path);

	t->exon_tree = exon_tree_new (t->exon_stmt,
			t->overlapping_stmt, t->cs);
}

static void
test_exon_tree_destroy (TestExonTree *t)
{
	if (t == NULL)
		return;

	db_finalize (t->exon_stmt);
	db_finalize (t->overlapping_stmt);

	db_close (t->db);
	chr_std_free (t->cs);

	exon_tree_free (t->exon_tree);

	// Remove temp files
	xunlink (t->db_path);
	xunlink (t->gtf_path);
}

START_TEST (test_exon_tree_index_dump)
{

	// Init ExonTree struct and create database
	// and gtf files
	TestExonTree t;
	test_exon_tree_init (&t);

	IBiTree *tree = NULL;
	sqlite3_stmt *search_stmt = NULL;

	int tree_id = 0;
	long start = 0;
	long end = 0;
	const char *exon_id = NULL;
	int i = 0;

	/* RUN FOOLS */
	exon_tree_index_dump (t.exon_tree, t.gtf_path, NULL, NULL, NULL);

	ck_assert_int_eq (hash_size (t.exon_tree->idx), 1);
	ck_assert_int_eq (hash_size (t.exon_tree->cache), gtf_size);

	// Test if cache hash contains all exons
	for (i = 0; i < gtf_size; i++)
		ck_assert_int_eq (hash_contains (t.exon_tree->cache,
					gtf_id[i]), 1);

	// Test if idx hash contains chromosome 1
	tree = hash_lookup (t.exon_tree->idx, "chr1");
	ck_assert (tree != NULL);

	// Test if all entries into database and interval tree
	// are the same entries into gtf file
	search_stmt = prepare_exon_search_stmt (t.db);

	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			start = db_column_int64 (search_stmt, 1);
			end = db_column_int64 (search_stmt, 2);
			exon_id = db_column_text (search_stmt, 3);

			ck_assert_int_eq (start, gtf_pos[i][0]);
			ck_assert_int_eq (end, gtf_pos[i][1]);
			ck_assert_str_eq (gtf_id[i], exon_id);

			tree_id = 0;

			// If looking for [start end] into tree,
			// it must find the same id into database
			ibitree_lookup (tree, start, end, -1, -1, 0,
					catch_id, &tree_id);

			ck_assert_int_eq (tree_id, i + 1);
		}

	// Cleanup all the mess
	db_finalize (search_stmt);
	test_exon_tree_destroy (&t);
}
END_TEST

static sqlite3_stmt *
prepare_overlapping_search_stmt (sqlite3 *db)
{
	const char sql[] =
		"SELECT exon_id, alignment_id FROM overlapping ORDER BY exon_id ASC";
	return db_prepare (db, sql);
}

START_TEST (test_exon_tree_lookup_dump)
{

	// Init ExonTree struct and create database
	// and gtf files
	TestExonTree t;
	test_exon_tree_init (&t);

	sqlite3_stmt *search_stmt = NULL;

	int exon_id = 0;
	int alignment_id = 0;

	int i = 0;
	int acm = 0;

	// alignment ids, pos and true 'acm'
	// (number of exons)
	int alignment_size = 5;
	int alignment_acm[] = {1, 1, 0, 2, 0};
	int alignment_ids[] = {1, 2, 3, 4, 5};
	int alignment_pos[][2] = {
		{1,    50},
		{49,   100},
		{110,  190},
		{200,  400},
		{2000, 2100}
	};

	int alignment_id_overlap_exon_id_size = 4;
	int alignment_id_overlap_exon_id[][2] = {
		{1, 1},
		{2, 1},
		{4, 2},
		{4, 3}
	};

	/* RUN FOOLS */
	exon_tree_index_dump (t.exon_tree, t.gtf_path, NULL, NULL, NULL);

	/* ADD IDS: 1, 2, 4, 4 */
	for (i = 0; i < alignment_size; i++)
		{
			acm = exon_tree_lookup_dump (t.exon_tree, "chr1",
					alignment_pos[i][0], alignment_pos[i][1],
					-1, -1, 0, alignment_ids[i]);
			ck_assert_int_eq (acm, alignment_acm[i]);
		}

	search_stmt = prepare_overlapping_search_stmt (t.db);

	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			exon_id = db_column_int (search_stmt, 0);
			ck_assert_int_eq (exon_id,
					alignment_id_overlap_exon_id[i][1]);

			alignment_id = db_column_int (search_stmt, 1);
			ck_assert_int_eq (alignment_id,
					alignment_id_overlap_exon_id[i][0]);
		}

	ck_assert_int_eq (i, alignment_id_overlap_exon_id_size);

	// Cleanup all the mess
	db_finalize (search_stmt);
	test_exon_tree_destroy (&t);
}
END_TEST

static void
index_custom (const char *suffix, const char *body, const GffFilter *filter,
		const char *gene_name_key, const char *exon_id_key,
		const char *expected[][4], int expected_size)
{
	char db_path[] = "/tmp/ponga.db.XXXXXX";
	char path[TEST_EXON_BUFSIZ] = {};
	sqlite3 *db = NULL;
	sqlite3_stmt *exon_stmt = NULL;
	sqlite3_stmt *overlapping_stmt = NULL;
	sqlite3_stmt *search_stmt = NULL;
	ChrStd *cs = NULL;
	ExonTree *exon_tree = NULL;
	FILE *fp = NULL;
	int fd = 0;
	int i = 0;

	snprintf (path, TEST_EXON_BUFSIZ, "/tmp/ponga.XXXXXX%s", suffix);
	fd = mkstemps (path, strlen (suffix));
	ck_assert_int_ge (fd, 0);
	fp = xfdopen (fd, "w");
	xfprintf (fp, "%s", body);
	xfclose (fp);

	db = create_db (db_path);
	exon_stmt = db_prepare_exon_stmt (db);
	overlapping_stmt = db_prepare_overlapping_stmt (db);
	cs = chr_std_new ();
	exon_tree = exon_tree_new (exon_stmt, overlapping_stmt, cs);

	exon_tree_index_dump (exon_tree, path, filter, gene_name_key, exon_id_key);

	search_stmt = db_prepare (db,
			"SELECT gene_name, start || '-' || end, strand, ense FROM exon ORDER BY id ASC");

	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			ck_assert_int_lt (i, expected_size);
			ck_assert_str_eq (db_column_text (search_stmt, 0), expected[i][0]);
			ck_assert_str_eq (db_column_text (search_stmt, 1), expected[i][1]);
			ck_assert_str_eq (db_column_text (search_stmt, 2), expected[i][2]);
			ck_assert_str_eq (db_column_text (search_stmt, 3), expected[i][3]);
		}

	ck_assert_int_eq (i, expected_size);

	db_finalize (search_stmt);
	db_finalize (exon_stmt);
	db_finalize (overlapping_stmt);
	db_close (db);
	chr_std_free (cs);
	exon_tree_free (exon_tree);
	xunlink (db_path);
	xunlink (path);
}

START_TEST (test_exon_tree_index_dump_bed)
{
	// BED is 0-based: The exon id is 'chr:start_end' 1-based.
	// The duplicated exon is indexed once
	const char *bed =
		"chr1\t0\t100\tGENE1\t0\t+\n"
		"chr1\t199\t300\tGENE1\t0\t+\n"
		"chr1\t199\t300\tGENE1\t0\t+\n"
		"chr2\t999\t1100\tGENE2\t0\t-\n";

	const char *expected[][4] = {
		{"GENE1", "1-100",     "+", "chr1:1_100"},
		{"GENE1", "200-300",   "+", "chr1:200_300"},
		{"GENE2", "1000-1100", "-", "chr2:1000_1100"}
	};

	index_custom (".bed", bed, NULL, NULL, NULL, expected, 3);
}
END_TEST

START_TEST (test_exon_tree_index_dump_custom)
{
	// Ensembl like: 'transcript_biotype', no 'gene_name' at the
	// second gene (falls back to 'gene_id') and no exon id at
	// the third exon (falls back to the coordinates)
	const char *gtf_ensembl =
		"chr1\t.\texon\t1\t100\t.\t+\t.\t"
		"gene_id \"G1\"; gene_symbol \"ponga\"; transcript_biotype \"protein_coding\"; ex \"E1\";\n"
		"chr1\t.\texon\t200\t300\t.\t+\t.\t"
		"gene_id \"G1\"; gene_symbol \"ponga\"; transcript_biotype \"lncRNA\"; ex \"E2\";\n"
		"chr1\t.\texon\t400\t500\t.\t-\t.\t"
		"gene_id \"G2\"; transcript_biotype \"protein_coding\"; ex \"E3\";\n"
		"chr1\t.\texon\t600\t700\t.\t-\t.\t"
		"gene_id \"G2\"; gene_symbol \"pongo\"; transcript_biotype \"protein_coding\";\n";

	const char *expected[][4] = {
		{"ponga", "1-100",   "+", "E1"},
		{"G2",    "400-500", "-", "E3"},
		{"pongo", "600-700", "-", "chr1:600_700"}
	};

	GffFilter *filter = gff_filter_new ();
	gff_filter_insert_feature (filter, "exon");
	gff_filter_insert_attribute (filter, "transcript_biotype=protein_coding");

	index_custom (".gtf", gtf_ensembl, filter, "gene_symbol", "ex", expected, 3);

	gff_filter_free (filter);
}
END_TEST

START_TEST (test_exon_tree_index_dump_missing_names)
{
	// No gene name: Use 'gene_id'. No gene name nor 'gene_id':
	// Skip the exon. 'gene_id' missing: Use the gene name
	const char *gtf_missing =
		"chr1\t.\texon\t1\t100\t.\t+\t.\t"
		"gene_id \"G1\"; transcript_type \"protein_coding\"; exon_id \"E1\";\n"
		"chr1\t.\texon\t200\t300\t.\t+\t.\t"
		"transcript_type \"protein_coding\"; exon_id \"E2\";\n"
		"chr1\t.\texon\t400\t500\t.\t+\t.\t"
		"gene_name \"ponga\"; transcript_type \"protein_coding\"; exon_id \"E3\";\n";

	const char *expected[][4] = {
		{"G1",    "1-100",   "+", "E1"},
		{"ponga", "400-500", "+", "E3"}
	};

	index_custom (".gtf", gtf_missing, NULL, NULL, NULL, expected, 2);

	exon_tree_free (NULL);
}
END_TEST

START_TEST (test_exon_tree_index_dump_bed_fatal)
{
	// BED4: No strand
	index_custom (".bed", "chr1\t0\t100\tGENE1\n", NULL, NULL, NULL, NULL, 0);
}
END_TEST

START_TEST (test_exon_tree_index_dump_bed_empty_fatal)
{
	index_custom (".bed", "", NULL, NULL, NULL, NULL, 0);
}
END_TEST

START_TEST (test_exon_tree_index_dump_no_exon_fatal)
{
	// No 'transcript_type=protein_coding'
	index_custom (".gtf",
			"chr1\t.\texon\t1\t100\t.\t+\t.\tgene_name \"ponga\"; exon_id \"E1\";\n",
			NULL, NULL, NULL, NULL, 0);
}
END_TEST

Suite *
make_exon_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("Exon");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_exon_tree_index_dump);
	tcase_add_test (tc_core, test_exon_tree_lookup_dump);
	tcase_add_test (tc_core, test_exon_tree_index_dump_bed);
	tcase_add_test (tc_core, test_exon_tree_index_dump_custom);
	tcase_add_test (tc_core, test_exon_tree_index_dump_missing_names);
	suite_add_tcase (s, tc_core);

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_exit_test (tc_abort, test_exon_tree_index_dump_bed_fatal,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_exon_tree_index_dump_bed_empty_fatal, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_exon_tree_index_dump_no_exon_fatal,   EXIT_FAILURE);
	suite_add_tcase (s, tc_abort);

	return s;
}
