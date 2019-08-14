#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/abnormal.h"
#include "../src/wrapper.h"
#include "../src/utils.h"
#include "../src/db.h"
#include "../src/cluster.h"

static sqlite3 *
create_db (char *db_path)
{
	int fd;

	fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
populate_db (sqlite3_stmt *alignment_stmt, int offset,
		const char *chr, const char **id,
		int (*pos)[2], int size)
{
	int i = 0;

	for (i = 0; i < size; i++)
		{
			db_insert_alignment (alignment_stmt, ++offset, id[i], 66,
					chr, pos[i][0], 60, "100M", pos[i][1], pos[i][1],
					chr, 1, ABNORMAL_EXONIC, 1);
			db_insert_alignment (alignment_stmt, ++offset, id[i], 66,
					chr, pos[i][0], 60, "100M", pos[i][1], pos[i][1],
					chr, 1, 0, 1);
		}
}

sqlite3_stmt *
prepare_query_stmt (sqlite3 *db)
{
	const char sql[] =
		"SELECT cluster_id, alignment_id, label, neighbors\n"
		"FROM clustering ORDER BY alignment_id ASC";
	return db_prepare (db, sql);
}

START_TEST (test_cluster)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";

	sqlite3 *db = NULL;
	sqlite3_stmt *alignment_stmt = NULL;
	sqlite3_stmt *clustering_stmt = NULL;
	sqlite3_stmt *search_stmt = NULL;

	int eps = 500;
	int min_pts = 3;
	int i = 0;
	int j = 0;

	const char *chr1 = "chr1";
	const char *chr2 = "chr2";

	const char *id1[] = {"id1", "id2", "id3",
		"id4", "id5", "id6"};

	const char *id2[] = {"id7", "id8", "id9",
		"id10", "id11", "id12"};

	// Number of positions
	// for each chromosome
	int size = 6;

	int pos[][2] = {
		{1000, 101},
		{1050, 101},
		{1300, 101},
		{2000, 101},
		{2500, 101},
		{2560, 101}
	};

	// True positive
	int true_size_col = 4;

	int true[][4] = {
		{1, 2, 3, 3},
		{1, 4, 3, 3},
		{1, 6, 3, 3},
		{2, 8, 2, 2},
		{2, 10, 3, 3},
		{2, 12, 2, 2},
		{3, 14, 3, 3},
		{3, 16, 3, 3},
		{3, 18, 3, 3},
		{4, 20, 2, 2},
		{4, 22, 3, 3},
		{4, 24, 2, 2}
	};

	db = create_db (db_file);
	alignment_stmt = db_prepare_alignment_stmt (db);
	clustering_stmt = db_prepare_clustering_stmt (db);

	// Populate alignment: To avoid same id constraints,
	// I need to pass an offset for the id: 0 and the
	// size number of fragments (size * 2, cause it's paired-end)
	populate_db (alignment_stmt, 0, chr1, id1, pos, size);
	populate_db (alignment_stmt, size * 2, chr2, id2, pos, size);

	// RUN
	cluster (clustering_stmt, eps, min_pts);

	// Let's get the clustering table values
	search_stmt = prepare_query_stmt (db);

	/* TIME TO TEST */
	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		for (j = 0; j < true_size_col; j++)
			ck_assert_int_eq (db_column_int (search_stmt, j),
					true[i][j]);

	ck_assert_int_eq (i, size * 2);

	db_finalize (alignment_stmt);
	db_finalize (clustering_stmt);
	db_finalize (search_stmt);
	db_close (db);
	xunlink (db_file);
}
END_TEST

Suite *
make_cluster_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Cluster");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_cluster);
	suite_add_tcase (s, tc_core);

	return s;
}
