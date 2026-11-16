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
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/db.h"
#include "../src/intron.h"

static sqlite3 *
create_db (char *db_path)
{
	int fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
populate_db (sqlite3 *db)
{
	/*
	* Retrocopy 1, from 'geneA' at chr1:
	*   exons 100-200, 150-300 and 180-250 (nested), 301-400
	*   (adjacent), 500-600 and 500-650, 900-1000
	*   union: 100-400, 500-650, 900-1000
	*   introns: 401-499 (exons 4 and 6: The longest starting
	*   at 500) and 651-899 (exons 6 and 7)
	* Retrocopy 2, from the same 'geneA' and from 'geneB',
	* which overlaps it at chr1 and has an exon at chr2 too:
	*   geneB at chr1: 450-470 and 700-720: intron 471-699
	*   geneB at chr2: 100-200: No intron
	* Retrocopy 3, from the single exon 'geneC': No intron
	* 'geneD' has introns, but it is not a parental gene
	*/
	const char sql[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO exon VALUES (1,'geneA','chr1',100,200,'+','ega','eea1');\n"
		"INSERT INTO exon VALUES (2,'geneA','chr1',150,300,'+','ega','eea2');\n"
		"INSERT INTO exon VALUES (3,'geneA','chr1',180,250,'+','ega','eea3');\n"
		"INSERT INTO exon VALUES (4,'geneA','chr1',301,400,'+','ega','eea4');\n"
		"INSERT INTO exon VALUES (5,'geneA','chr1',500,600,'+','ega','eea5');\n"
		"INSERT INTO exon VALUES (6,'geneA','chr1',500,650,'+','ega','eea6');\n"
		"INSERT INTO exon VALUES (7,'geneA','chr1',900,1000,'+','ega','eea7');\n"
		"INSERT INTO exon VALUES (8,'geneB','chr1',450,470,'-','egb','eeb1');\n"
		"INSERT INTO exon VALUES (9,'geneB','chr1',700,720,'-','egb','eeb2');\n"
		"INSERT INTO exon VALUES (10,'geneB','chr2',100,200,'-','egb','eeb3');\n"
		"INSERT INTO exon VALUES (11,'geneC','chr3',100,200,'+','egc','eec1');\n"
		"INSERT INTO exon VALUES (12,'geneD','chr4',100,200,'+','egd','eed1');\n"
		"INSERT INTO exon VALUES (13,'geneD','chr4',300,400,'+','egd','eed2');\n"
		"INSERT INTO source VALUES (1,1,'ponga.bam');\n"
		"INSERT INTO cluster VALUES (1,1,'chr9',1,200,'geneA',31);\n"
		"INSERT INTO cluster VALUES (2,1,'chr9',1000,1200,'geneA',31);\n"
		"INSERT INTO cluster VALUES (3,1,'chr9',1100,1300,'geneB',31);\n"
		"INSERT INTO cluster VALUES (4,1,'chr9',5000,5200,'geneC',31);\n"
		"INSERT INTO cluster_merging VALUES (1,1,1);\n"
		"INSERT INTO cluster_merging VALUES (2,2,1);\n"
		"INSERT INTO cluster_merging VALUES (2,3,1);\n"
		"INSERT INTO cluster_merging VALUES (3,4,1);\n"
		"COMMIT;";

	db_exec (db, sql);
}

static void
check_intron (const Intron *intron, const char *gene_name,
		long start, long end, int exon_id_left, int exon_id_right)
{
	ck_assert_ptr_nonnull (intron);
	ck_assert_str_eq (intron->gene_name, gene_name);
	ck_assert_str_eq (intron->chr, "chr1");
	ck_assert_int_eq (intron->start, start);
	ck_assert_int_eq (intron->end, end);
	ck_assert_int_eq (intron->exon_id_left, exon_id_left);
	ck_assert_int_eq (intron->exon_id_right, exon_id_right);
}

START_TEST (test_intron_tree)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = NULL;
	IntronTree *t = NULL;

	db = create_db (db_file);
	populate_db (db);

	t = intron_tree_new (db);

	// Retrocopy 1: 2 introns; retrocopy 2: 2 of geneA
	// and 1 of geneB; retrocopy 3: none
	ck_assert_int_eq (intron_tree_size (t), 5);

	// Inside each intron, and at its very ends
	check_intron (intron_tree_lookup (t, 1, "chr1", 420, 480),
			"geneA", 401, 499, 4, 6);
	check_intron (intron_tree_lookup (t, 1, "chr1", 401, 499),
			"geneA", 401, 499, 4, 6);
	check_intron (intron_tree_lookup (t, 1, "chr1", 700, 700),
			"geneA", 651, 899, 6, 7);

	// Reaching an exon, inside the exons and outside the gene
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 400, 450));
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 450, 500));
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 350, 380));
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 1, 50));
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 2000, 2010));

	// Nested and adjacent exons make no intron
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 260, 280));

	// Another contig, and a contig of another parental gene
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr2", 420, 480));
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr4", 250, 260));

	// Retrocopy 2 has the introns of both parental genes:
	// 480-490 is intron for both, 660-690 too
	check_intron (intron_tree_lookup (t, 2, "chr1", 410, 440),
			"geneA", 401, 499, 4, 6);
	check_intron (intron_tree_lookup (t, 2, "chr1", 500, 640),
			"geneB", 471, 699, 8, 9);
	ck_assert_ptr_nonnull (intron_tree_lookup (t, 2, "chr1", 480, 490));
	ck_assert_ptr_null (intron_tree_lookup (t, 2, "chr2", 120, 130));

	// The intron of geneB is not for the retrocopy 1
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 500, 640));

	// Single exon gene, and unknown retrocopy
	ck_assert_ptr_null (intron_tree_lookup (t, 3, "chr3", 300, 310));
	ck_assert_ptr_null (intron_tree_lookup (t, 9, "chr1", 420, 480));

	intron_tree_free (t);
	intron_tree_free (NULL);

	db_close (db);
	xunlink (db_file);
}
END_TEST

START_TEST (test_intron_tree_empty)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = NULL;
	IntronTree *t = NULL;

	// No retrocopies at all
	db = create_db (db_file);
	t = intron_tree_new (db);

	ck_assert_int_eq (intron_tree_size (t), 0);
	ck_assert_ptr_null (intron_tree_lookup (t, 1, "chr1", 1, 10));

	intron_tree_free (t);
	db_close (db);
	xunlink (db_file);
}
END_TEST

START_TEST (test_intron_dump)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = NULL;
	sqlite3_stmt *intron_stmt = NULL;
	sqlite3_stmt *search_stmt = NULL;
	IntronTree *t = NULL;

	db = create_db (db_file);
	populate_db (db);

	t = intron_tree_new (db);
	intron_stmt = db_prepare_intron_retention_stmt (db);

	intron_dump (intron_stmt, intron_tree_lookup (t, 1, "chr1", 420, 480), 1, 1, 7);
	intron_dump (intron_stmt, intron_tree_lookup (t, 2, "chr1", 500, 640), 2, 1, 3);

	search_stmt = db_prepare (db,
			"SELECT retrocopy_id, source_id, gene_name, chr, start, end,\n"
			"	exon_id_left, exon_id_right, reads\n"
			"FROM intron_retention\n"
			"ORDER BY retrocopy_id");

	ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
	ck_assert_int_eq (db_column_int (search_stmt, 0), 1);
	ck_assert_int_eq (db_column_int (search_stmt, 1), 1);
	ck_assert_str_eq (db_column_text (search_stmt, 2), "geneA");
	ck_assert_str_eq (db_column_text (search_stmt, 3), "chr1");
	ck_assert_int_eq (db_column_int (search_stmt, 4), 401);
	ck_assert_int_eq (db_column_int (search_stmt, 5), 499);
	ck_assert_int_eq (db_column_int (search_stmt, 6), 4);
	ck_assert_int_eq (db_column_int (search_stmt, 7), 6);
	ck_assert_int_eq (db_column_int (search_stmt, 8), 7);

	ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
	ck_assert_int_eq (db_column_int (search_stmt, 0), 2);
	ck_assert_str_eq (db_column_text (search_stmt, 2), "geneB");
	ck_assert_int_eq (db_column_int (search_stmt, 4), 471);
	ck_assert_int_eq (db_column_int (search_stmt, 5), 699);
	ck_assert_int_eq (db_column_int (search_stmt, 6), 8);
	ck_assert_int_eq (db_column_int (search_stmt, 7), 9);
	ck_assert_int_eq (db_column_int (search_stmt, 8), 3);

	ck_assert_int_eq (db_step (search_stmt), SQLITE_DONE);

	db_finalize (search_stmt);
	db_finalize (intron_stmt);
	intron_tree_free (t);
	db_close (db);
	xunlink (db_file);
}
END_TEST

Suite *
make_intron_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Intron");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_intron_tree);
	tcase_add_test (tc_core, test_intron_tree_empty);
	tcase_add_test (tc_core, test_intron_dump);
	suite_add_tcase (s, tc_core);

	return s;
}
