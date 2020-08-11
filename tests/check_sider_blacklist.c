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
#include "../src/blacklist.h"

#define TEST_BLACKLIST_BUFSIZ 64

struct _TestBlacklist
{
	char db_path[TEST_BLACKLIST_BUFSIZ];
	char file_path[TEST_BLACKLIST_BUFSIZ];

	sqlite3 *db;
	sqlite3_stmt *blacklist_stmt;
	sqlite3_stmt *overlapping_blacklist_stmt;

	ChrStd *cs;

	Blacklist *blacklist;
};

typedef struct _TestBlacklist TestBlacklist;

static const char *gtf =
	"chr1\t.\tgene\t1\t1000\t.\t+\t.\t"
	"gene_name \"ponga1\"; gene_type \"processed_pseudogene\";\n"
	"chr17\t.\tgene\t1\t1000\t.\t+\t.\t"
	"tag \"retrogene\";\n";

static const char *bed =
	"chr1\t1\t1000\tponga1\n"
	"chr17\t1\t1000\n";

static sqlite3 *
create_db (char *db_path)
{
	int fd;

	fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
create_file (char *path, const char *content)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, "%s", content);

	xfclose (fp);
}

static void
test_blacklist_init (TestBlacklist *t, const char *content)
{
	strncpy (t->db_path, "/tmp/ponga.db.XXXXXX",
			TEST_BLACKLIST_BUFSIZ - 1);
	strncpy (t->file_path, "/tmp/ponga.gtf.XXXXXX",
			TEST_BLACKLIST_BUFSIZ - 1);

	t->db_path[TEST_BLACKLIST_BUFSIZ - 1] = '\0';
	t->file_path[TEST_BLACKLIST_BUFSIZ - 1] = '\0';

	t->db = create_db (t->db_path);
	t->blacklist_stmt = db_prepare_blacklist_stmt (t->db);
	t->overlapping_blacklist_stmt = db_prepare_overlapping_blacklist_stmt (t->db);

	t->cs = chr_std_new ();

	create_file (t->file_path, content);

	t->blacklist = blacklist_new (t->blacklist_stmt,
			t->overlapping_blacklist_stmt, t->cs);
}

static void
test_blacklist_destroy (TestBlacklist *t)
{
	if (t == NULL)
		return;

	db_finalize (t->blacklist_stmt);
	db_finalize (t->overlapping_blacklist_stmt);

	db_close (t->db);
	chr_std_free (t->cs);

	blacklist_free (NULL);
	blacklist_free (t->blacklist);

	// Remove temp files
	xunlink (t->db_path);
	xunlink (t->file_path);
}

START_TEST (test_blacklist_index_dump_from_gff)
{
	IBiTree *tree = NULL;
	TestBlacklist t;
	test_blacklist_init (&t, gtf);

	GffFilter *filter = gff_filter_new ();
	gff_filter_insert_feature (filter, "gene");
	gff_filter_insert_soft_attribute (filter, "gene_type", "processed_pseudogene");
	gff_filter_insert_soft_attribute (filter, "tag", "retrogene");

	blacklist_index_dump_from_gff (t.blacklist, t.file_path, filter);
	ck_assert_int_eq (hash_size (t.blacklist->idx), 2);

	// Test if idx hash contains chromosome 2
	tree = hash_lookup (t.blacklist->idx, "chr1");
	ck_assert (tree != NULL);

	tree = hash_lookup (t.blacklist->idx, "chr17");
	ck_assert (tree != NULL);

	gff_filter_free (filter);
	test_blacklist_destroy (&t);
}
END_TEST

START_TEST (test_blacklist_index_dump_from_bed)
{
	IBiTree *tree = NULL;
	TestBlacklist t;
	test_blacklist_init (&t, bed);

	blacklist_index_dump_from_bed (t.blacklist, t.file_path);
	ck_assert_int_eq (hash_size (t.blacklist->idx), 2);

	// Test if idx hash contains chromosome 2
	tree = hash_lookup (t.blacklist->idx, "chr1");
	ck_assert (tree != NULL);

	tree = hash_lookup (t.blacklist->idx, "chr17");
	ck_assert (tree != NULL);

	test_blacklist_destroy (&t);
}
END_TEST

START_TEST (test_blacklist_lookup)
{
	TestBlacklist t;
	test_blacklist_init (&t, bed);

	// alignment ids, pos and true 'acm'
	int alignment_size = 5;
	int alignment_acm[] = {1, 1, 1, 1, 0};
	int alignment_ids[] = {1, 2, 3, 4, 5};
	int alignment_pos[][2] = {
		{1,    50},
		{49,   100},
		{110,  190},
		{200,  400},
		{2000, 2100}
	};

	int i = 0;
	int acm = 0;

	/* RUN FOOLS */
	blacklist_index_dump_from_bed (t.blacklist, t.file_path);

	/* ADD IDS: 1, 2, 4, 4 */
	for (i = 0; i < alignment_size; i++)
		{
			acm = blacklist_lookup (t.blacklist, "chr1",
					alignment_pos[i][0], alignment_pos[i][1],
					0, alignment_ids[i], 1);
			ck_assert_int_eq (acm, alignment_acm[i]);
		}

	test_blacklist_destroy (&t);
}
END_TEST

Suite *
make_blacklist_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Blacklist");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_blacklist_index_dump_from_gff);
	tcase_add_test (tc_core, test_blacklist_index_dump_from_bed);
	tcase_add_test (tc_core, test_blacklist_lookup);
	suite_add_tcase (s, tc_core);

	return s;
}
