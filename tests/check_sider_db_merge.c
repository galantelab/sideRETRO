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

#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/db.h"
#include "../src/db_merge.h"

static void
create_db_sample (char *path, const char *sample)
{
	char *sql = NULL;

	// One sample by database, as process-sample of each one
	xasprintf (&sql,
		"INSERT INTO batch VALUES(1,\"2019-02-31\");\n"
		"INSERT INTO source VALUES(1,1,\"%s\");\n"
		"INSERT INTO exon VALUES(1,\"g1\",\"chr1\",1,200,\"+\",\"ENG0066\",\"ENSE0066\");\n"
		"INSERT INTO alignment VALUES(1,\"r1\",99,\"chr1\",1,20,\"101M\",101,101,\"chr1\",200,0,1,\"*\",\"*\");\n"
		"INSERT INTO overlapping VALUES(1,1,1,101);", sample);

	int fd = xmkstemp (path);
	close (fd);

	sqlite3 *db = db_create (path);
	db_exec (db, sql);
	db_close (db);

	xfree (sql);
}

static void
create_db (char *path)
{
	// A sample named after the database
	create_db_sample (path, path);
}

static long
count_rows (sqlite3 *db, const char *sql)
{
	sqlite3_stmt *stmt = db_prepare (db, sql);
	long n = 0;

	if (db_step (stmt) == SQLITE_ROW)
		n = db_column_int64 (stmt, 0);

	db_finalize (stmt);
	return n;
}

START_TEST (test_db_merge)
{
	int i = 0;
	int num_db = 3;
	char db_path1[] = "/tmp/ponga1.db.XXXXXX";
	char db_path2[] = "/tmp/ponga2.db.XXXXXX";
	char db_path3[] = "/tmp/ponga3.db.XXXXXX";

	char *db_paths[] = {db_path1, db_path2, db_path3};

	for (i = 0; i < num_db; i++)
		create_db (db_paths[i]);

	sqlite3 *db = db_create (":memory:");

	db_merge (db, num_db, db_paths);

	// The exon is shared: Merged once
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM source"), 3);
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM alignment"), 3);
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM exon"), 1);
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM overlapping"), 3);

	db_close (db);

	for (i = 0; i < num_db; i++)
		xunlink (db_paths[i]);
}
END_TEST

START_TEST (test_db_merge_in_line)
{
	int i = 0;
	int num_db = 3;
	char db_path1[] = "/tmp/ponga1.db.XXXXXX";
	char db_path2[] = "/tmp/ponga2.db.XXXXXX";
	char db_path3[] = "/tmp/ponga3.db.XXXXXX";

	char *db_paths[] = {db_path1, db_path2, db_path3};

	for (i = 0; i < num_db; i++)
		create_db (db_paths[i]);

	// Merge  databases to the first
	sqlite3 *db = db_connect (db_path1);

	db_merge (db, num_db - 1, &db_paths[1]);

	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM source"), 3);
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM alignment"), 3);

	db_close (db);

	for (i = 0; i < num_db; i++)
		xunlink (db_paths[i]);
}
END_TEST

START_TEST (test_db_merge_exon_id)
{
	char db_path1[] = "/tmp/ponga1.db.XXXXXX";
	char db_path2[] = "/tmp/ponga2.db.XXXXXX";
	char *db_paths[] = {db_path2};
	sqlite3 *db = NULL;
	sqlite3_stmt *stmt = NULL;
	int i = 0;

	// The second database has the same exon (ENSE0066) with
	// another id: Its alignments must point to the merged one
	const char sql2[] =
		"INSERT INTO batch VALUES(1,\"2019-02-31\");\n"
		"INSERT INTO source VALUES(1,1,\"pongo.bam\");\n"
		"INSERT INTO exon VALUES(1,\"g2\",\"chr2\",1,200,\"+\",\"ENG0099\",\"ENSE0099\");\n"
		"INSERT INTO exon VALUES(2,\"g1\",\"chr1\",1,200,\"+\",\"ENG0066\",\"ENSE0066\");\n"
		"INSERT INTO alignment VALUES(1,\"r2\",99,\"chr1\",1,20,\"101M\",101,101,\"chr1\",200,0,1,\"*\",\"*\");\n"
		"INSERT INTO overlapping VALUES(2,1,1,101);";

	create_db (db_path1);

	int fd = xmkstemp (db_path2);
	close (fd);
	db = db_create (db_path2);
	db_exec (db, sql2);
	db_close (db);

	db = db_connect (db_path1);
	db_merge (db, 1, db_paths);

	stmt = db_prepare (db,
			"SELECT a.qname, e.ense FROM overlapping AS o\n"
			"INNER JOIN exon AS e ON e.id = o.exon_id\n"
			"INNER JOIN alignment AS a ON a.id = o.alignment_id\n"
			"ORDER BY a.qname ASC");

	const char *expected[][2] = {{"r1", "ENSE0066"}, {"r2", "ENSE0066"}};

	for (i = 0; db_step (stmt) == SQLITE_ROW; i++)
		{
			ck_assert_int_lt (i, 2);
			ck_assert_str_eq (db_column_text (stmt, 0), expected[i][0]);
			ck_assert_str_eq (db_column_text (stmt, 1), expected[i][1]);
		}

	ck_assert_int_eq (i, 2);

	db_finalize (stmt);
	db_close (db);
	xunlink (db_path1);
	xunlink (db_path2);
}
END_TEST

START_TEST (test_db_merge_same_file)
{
	char db_path1[] = "/tmp/ponga1.db.XXXXXX";
	char db_path2[] = "/tmp/ponga2.db.XXXXXX";
	char *alias = NULL;
	char *db_paths[3] = {};

	create_db (db_path1);
	create_db (db_path2);

	// The second database again, by another path, and
	// the output database itself: Merged only once
	xasprintf (&alias, "/tmp/../tmp/%s", db_path2 + 5);

	db_paths[0] = db_path2;
	db_paths[1] = alias;
	db_paths[2] = db_path1;

	sqlite3 *db = db_connect (db_path1);
	db_merge (db, 3, db_paths);

	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM source"), 2);
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM alignment"), 2);

	db_close (db);
	xunlink (db_path1);
	xunlink (db_path2);
	xfree (alias);
}
END_TEST

START_TEST (test_db_merge_same_sample)
{
	// The same sample in two databases: Its reads
	// would be counted twice
	char db_path1[] = "/tmp/ponga1.db.XXXXXX";
	char db_path2[] = "/tmp/ponga2.db.XXXXXX";
	char *db_paths[] = {db_path1, db_path2};

	create_db_sample (db_path1, "/data/ponga.bam");
	create_db_sample (db_path2, "/data/ponga.bam");

	sqlite3 *db = db_create (":memory:");
	db_merge (db, 2, db_paths);
}
END_TEST

START_TEST (test_db_merge_warnings)
{
	char db_path1[] = "/tmp/ponga1.db.XXXXXX";
	char db_path2[] = "/tmp/ponga2.db.XXXXXX";
	char *db_paths[] = {db_path1, db_path2};
	sqlite3 *db2 = NULL;

	// The same sample name from different paths, and a
	// database with one more exon (another annotation):
	// Warnings only
	create_db_sample (db_path1, "/data1/ponga.bam");
	create_db_sample (db_path2, "/data2/ponga.bam");

	db2 = db_connect (db_path2);
	db_exec (db2, "INSERT INTO exon VALUES(2,\"g2\",\"chr2\",1,200,\"+\",\"ENG0067\",\"ENSE0067\");");
	db_close (db2);

	sqlite3 *db = db_create (":memory:");
	db_merge (db, 2, db_paths);

	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM source"), 2);
	ck_assert_int_eq (count_rows (db, "SELECT COUNT(*) FROM exon"), 2);

	db_close (db);
	xunlink (db_path1);
	xunlink (db_path2);
}
END_TEST

Suite *
make_db_merge_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("DBMerge");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_db_merge);
	tcase_add_test (tc_core, test_db_merge_in_line);
	tcase_add_test (tc_core, test_db_merge_exon_id);
	tcase_add_test (tc_core, test_db_merge_same_file);
	tcase_add_test (tc_core, test_db_merge_warnings);
	suite_add_tcase (s, tc_core);

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_exit_test (tc_abort, test_db_merge_same_sample, EXIT_FAILURE);
	suite_add_tcase (s, tc_abort);

	return s;
}
