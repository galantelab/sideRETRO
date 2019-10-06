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
#include "../src/set.h"
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
populate_db (sqlite3 *db)
{
	// Database dump
	static const char schema[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO exon VALUES (1,'gene1','chr11',1,3000,'+','eg1','ee1');\n"
		"INSERT INTO exon VALUES (2,'gene2','chr12',1,3000,'-','eg2','ee2');\n"
		"INSERT INTO alignment VALUES(1,'id1',66,'chr11',1000,60,'100M',101,101,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES(2,'id1',66,'chr1',1000,60,'100M',101,101,'chr1',1,0,1);\n"
		"INSERT INTO alignment VALUES(3,'id2',66,'chr11',1050,60,'100M',101,101,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES(4,'id2',66,'chr1',1050,60,'100M',101,101,'chr1',1,0,1);\n"
		"INSERT INTO alignment VALUES(5,'id3',66,'chr11',1300,60,'100M',101,101,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES(6,'id3',66,'chr1',1300,60,'100M',101,101,'chr1',1,0,1);\n"
		"INSERT INTO alignment VALUES(7,'id4',66,'chr11',2000,60,'100M',101,101,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES(8,'id4',66,'chr1',2000,60,'100M',101,101,'chr1',1,0,1);\n"
		"INSERT INTO alignment VALUES(9,'id5',66,'chr11',2500,60,'100M',101,101,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES(10,'id5',66,'chr1',2500,60,'100M',101,101,'chr1',1,0,1);\n"
		"INSERT INTO alignment VALUES(11,'id6',66,'chr11',2560,60,'100M',101,101,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES(12,'id6',66,'chr1',2560,60,'100M',101,101,'chr1',1,0,1);\n"
		"INSERT INTO alignment VALUES(13,'id7',66,'chr12',1000,60,'100M',101,101,'chr2',1,8,1);\n"
		"INSERT INTO alignment VALUES(14,'id7',66,'chr2',1000,60,'100M',101,101,'chr2',1,0,1);\n"
		"INSERT INTO alignment VALUES(15,'id8',66,'chr12',1050,60,'100M',101,101,'chr2',1,8,1);\n"
		"INSERT INTO alignment VALUES(16,'id8',66,'chr2',1050,60,'100M',101,101,'chr2',1,0,1);\n"
		"INSERT INTO alignment VALUES(17,'id9',66,'chr12',1300,60,'100M',101,101,'chr2',1,8,1);\n"
		"INSERT INTO alignment VALUES(18,'id9',66,'chr2',1300,60,'100M',101,101,'chr2',1,0,1);\n"
		"INSERT INTO alignment VALUES(19,'id10',66,'chr12',2000,60,'100M',101,101,'chr2',1,8,1);\n"
		"INSERT INTO alignment VALUES(20,'id10',66,'chr2',2000,60,'100M',101,101,'chr2',1,0,1);\n"
		"INSERT INTO alignment VALUES(21,'id11',66,'chr12',2500,60,'100M',101,101,'chr2',1,8,1);\n"
		"INSERT INTO alignment VALUES(22,'id11',66,'chr2',2500,60,'100M',101,101,'chr2',1,0,1);\n"
		"INSERT INTO alignment VALUES(23,'id12',66,'chr12',2560,60,'100M',101,101,'chr2',1,8,1);\n"
		"INSERT INTO alignment VALUES(24,'id12',66,'chr2',2560,60,'100M',101,101,'chr2',1,0,1);\n"
		"INSERT INTO overlapping VALUES(1,1,1,100);\n"
		"INSERT INTO overlapping VALUES(1,3,1,100);\n"
		"INSERT INTO overlapping VALUES(1,5,1,100);\n"
		"INSERT INTO overlapping VALUES(1,7,1,100);\n"
		"INSERT INTO overlapping VALUES(1,9,1,100);\n"
		"INSERT INTO overlapping VALUES(1,11,1,100);\n"
		"INSERT INTO overlapping VALUES(2,13,1,100);\n"
		"INSERT INTO overlapping VALUES(2,15,1,100);\n"
		"INSERT INTO overlapping VALUES(2,17,1,100);\n"
		"INSERT INTO overlapping VALUES(2,19,1,100);\n"
		"INSERT INTO overlapping VALUES(2,21,1,100);\n"
		"INSERT INTO overlapping VALUES(2,23,1,100);\n"
		"COMMIT;";

	db_exec (db, schema);
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
	sqlite3_stmt *clustering_stmt = NULL;
	sqlite3_stmt *search_stmt = NULL;

	Set *blacklist_chr = set_new (NULL);

	int eps = 500;
	int min_pts = 3;
	int distance = 10000;
	int i = 0;
	int j = 0;

	// Number of positions
	// for each chromosome
	int size = 6;

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
	clustering_stmt = db_prepare_clustering_stmt (db);

	// Populate database
	populate_db (db);

	// RUN
	cluster (clustering_stmt, eps, min_pts, blacklist_chr, distance);

	// Let's get the clustering table values
	search_stmt = prepare_query_stmt (db);

	/* TIME TO TEST */
	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		for (j = 0; j < true_size_col; j++)
			ck_assert_int_eq (db_column_int (search_stmt, j),
					true[i][j]);

	ck_assert_int_eq (i, size * 2);

	db_finalize (clustering_stmt);
	db_finalize (search_stmt);
	set_free (blacklist_chr);
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
