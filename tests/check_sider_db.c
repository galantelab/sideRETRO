#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/db.h"

static void
handle_sigabrt (int sig)
{
	if (sig == SIGABRT)
		exit (1);
}

static sqlite3 *
create_db (char *path)
{
	int fd = xmkstemp (path);
	close (fd);
	return db_open (path,
			SQLITE_OPEN_CREATE|SQLITE_OPEN_READWRITE);
}

START_TEST (test_db_open)
{
	char db_path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = create_db (db_path);
	db_close (db);

	db = db_open (db_path, SQLITE_OPEN_READWRITE);
	ck_assert (db != NULL);
	db_close (db);

	xunlink (db_path);
}
END_TEST

START_TEST (test_db_open_abort)
{
	db_open ("it-must-not-exist", SQLITE_OPEN_READWRITE);
}
END_TEST

START_TEST (test_db_close_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READONLY);
	db_prepare (db, "HEY, CREATE PLEASE TABLE ponga (id INTEGER)");
	db_close (db);
}
END_TEST

START_TEST (test_db_exec)
{
	char db_path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = create_db (db_path);

	db_exec (db, "CREATE TABLE ponga (id INTEGER)");

	db_close (db);
	xunlink (db_path);
}
END_TEST

START_TEST (test_db_exec_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READONLY);
	db_exec (db, "HEY, CREATE PLEASE TABLE ponga (id INTEGER)");
}
END_TEST

START_TEST (test_db_prepare)
{
	char db_path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = create_db (db_path);

	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");

	db_step (stmt);

	db_finalize (stmt);
	db_close (db);
	xunlink (db_path);
}
END_TEST

START_TEST (test_db_prepare_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READONLY);
	db_prepare (db, "HEY, CREATE PLEASE TABLE ponga (id INTEGER)");
}
END_TEST

START_TEST (test_db_step_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");
	db_step (stmt);
	db_step (stmt);
}
END_TEST

START_TEST (test_db_reset_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");
	db_step (stmt);
	sqlite3_step (stmt);
	db_reset (stmt);
}
END_TEST

START_TEST (test_db_bind_int_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");
	db_bind_int (stmt, 10, 66);
}
END_TEST

START_TEST (test_db_bind_int64_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");
	db_bind_int64 (stmt, 10, 666666666666);
}
END_TEST

START_TEST (test_db_bind_double_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");
	db_bind_double (stmt, 10, 66.66);
}
END_TEST

START_TEST (test_db_bind_text_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	sqlite3_stmt *stmt = db_prepare (db,
			"CREATE TABLE ponga (id INTEGER)");
	db_bind_text (stmt, 10, "It must fail");
}
END_TEST

START_TEST (test_db_column_int_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	db_exec (db, "CREATE TABLE ponga (i INTEGER, t TEXT, d REAL)");
	db_exec (db, "INSERT INTO ponga (i,t,d) VALUES (1, \"ponga\", 66.66)");
	sqlite3_stmt *stmt = db_prepare (db,
			"SELECT * FROM ponga");
	db_step (stmt);
	db_column_int (stmt, 1);
}
END_TEST

START_TEST (test_db_column_int64_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	db_exec (db, "CREATE TABLE ponga (i INTEGER, t TEXT, d REAL)");
	db_exec (db, "INSERT INTO ponga (i,t,d) VALUES (1, \"ponga\", 66.66)");
	sqlite3_stmt *stmt = db_prepare (db,
			"SELECT * FROM ponga");
	db_step (stmt);
	db_column_int64 (stmt, 1);
}
END_TEST

START_TEST (test_db_column_double_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	db_exec (db, "CREATE TABLE ponga (i INTEGER, t TEXT, d REAL)");
	db_exec (db, "INSERT INTO ponga (i,t,d) VALUES (1, \"ponga\", 66.66)");
	sqlite3_stmt *stmt = db_prepare (db,
			"SELECT * FROM ponga");
	db_step (stmt);
	db_column_double (stmt, 1);
}
END_TEST

START_TEST (test_db_column_text_abort)
{
	sqlite3 *db = db_open (":memory:", SQLITE_OPEN_READWRITE);
	db_exec (db, "CREATE TABLE ponga (i INTEGER, t TEXT, d REAL)");
	db_exec (db, "INSERT INTO ponga (i,t,d) VALUES (1, \"ponga\", 66.66)");
	sqlite3_stmt *stmt = db_prepare (db,
			"SELECT * FROM ponga");
	db_step (stmt);
	db_column_text (stmt, 0);
}
END_TEST

START_TEST (test_db_schema)
{
	char db_path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = create_db (db_path);
	db_close (db);

	db = db_create (db_path);

	sqlite3_stmt *batch_stmt = db_prepare_batch_stmt (db);
	sqlite3_stmt *source_stmt = db_prepare_source_stmt (db);
	sqlite3_stmt *exon_stmt = db_prepare_exon_stmt (db);
	sqlite3_stmt *alignment_stmt = db_prepare_alignment_stmt (db);
	sqlite3_stmt *overlapping_stmt = db_prepare_overlapping_stmt (db);
	sqlite3_stmt *clustering_stmt = db_prepare_clustering_stmt (db);
	sqlite3_stmt *cluster_stmt = db_prepare_cluster_stmt (db);
	sqlite3_stmt *overlapping_blacklist_stmt =
		db_prepare_overlapping_blacklist_stmt (db);
	sqlite3_stmt *blacklist_stmt = db_prepare_blacklist_stmt (db);
	sqlite3_stmt *cluster_merge_stmt = db_prepare_cluster_merging_stmt (db);
	sqlite3_stmt *retrocopy_stmt = db_prepare_retrocopy_stmt (db);
	sqlite3_stmt *genotype_stmt = db_prepare_genotype_stmt (db);

	db_cache_size (db, DB_DEFAULT_CACHE_SIZE - 1);
	db_begin_transaction (db);

	db_insert_batch (batch_stmt, 1, "2019-02-31");
	db_insert_source (source_stmt, 1, 1, "ponga.bam");
	db_insert_exon (exon_stmt, 1, "PONGA", "chr1", 1, 200, "+",
			"ENSG000666", "ENSE000666");
	db_insert_alignment (alignment_stmt, 1, "run1", 99, "chr1", 1, 20,
			"101M", 101, 101, "chr1", 200, 0, 1);
	db_insert_overlapping (overlapping_stmt, 1, 1, 1, 101);
	db_insert_clustering (clustering_stmt, 1, 1, 1, 0, 1);
	db_insert_cluster (cluster_stmt, 1, 1, "chr1", 1, 101, "PONGA", 1);
	db_insert_blacklist (blacklist_stmt, 1, "blackponga", "chr1", 1, 200);
	db_insert_overlapping_blacklist (overlapping_blacklist_stmt, 1, 1, 1, 1, 101);
	db_insert_cluster_merging (cluster_merge_stmt, 1, 1, 1);
	db_insert_retrocopy (retrocopy_stmt, 1, "chr1", 1, 200, "ponga1/ponga2",
			12, 100, 1, -0.87, 0.00001);
	db_insert_genotype (genotype_stmt, 1, 1, 0);

	db_end_transaction (db);

	db_finalize (batch_stmt);
	db_finalize (source_stmt);
	db_finalize (exon_stmt);
	db_finalize (alignment_stmt);
	db_finalize (overlapping_stmt);
	db_finalize (clustering_stmt);
	db_finalize (cluster_stmt);
	db_finalize (blacklist_stmt);
	db_finalize (overlapping_blacklist_stmt);
	db_finalize (cluster_merge_stmt);
	db_finalize (retrocopy_stmt);
	db_finalize (genotype_stmt);
	db_close (db);
	xunlink (db_path);
}
END_TEST

Suite *
make_db_suite (void)
{
	setup_signal (SIGABRT, handle_sigabrt);

	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("DB");

	/* Core test case */
	tc_core = tcase_create ("Core");

	/* Abort test case */
	tc_abort = tcase_create ("Abort");
	/*tcase_set_tags (tc_abort, "no-valgrind");*/

	tcase_add_test (tc_core, test_db_open);
	tcase_add_test (tc_core, test_db_exec);
	tcase_add_test (tc_core, test_db_prepare);
	tcase_add_test (tc_core, test_db_schema);

	tcase_add_exit_test (tc_abort, test_db_open_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_close_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_exec_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_prepare_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_step_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_reset_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_bind_int_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_bind_int64_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_bind_double_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_bind_text_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_column_int_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_column_int64_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_column_double_abort, 1);
	tcase_add_exit_test (tc_abort, test_db_column_text_abort, 1);

	suite_add_tcase (s, tc_core);
	suite_add_tcase (s, tc_abort);

	return s;
}
