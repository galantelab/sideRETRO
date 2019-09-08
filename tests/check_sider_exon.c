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
	exon_tree_index_dump (t.exon_tree, t.gtf_path);

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
	exon_tree_index_dump (t.exon_tree, t.gtf_path);

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

Suite *
make_exon_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Exon");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_exon_tree_index_dump);
	tcase_add_test (tc_core, test_exon_tree_lookup_dump);
	suite_add_tcase (s, tc_core);

	return s;
}
