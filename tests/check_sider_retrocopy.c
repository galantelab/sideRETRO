#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/utils.h"
#include "../src/db.h"
#include "../src/retrocopy.h"

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
		"INSERT INTO exon VALUES (1,'gene1','chr1',1,3000,'+','eg1','ee1');\n"
		"INSERT INTO exon VALUES (2,'gene2_1','chr2',1,3000,'-','eg2','ee2');\n"
		"INSERT INTO exon VALUES (3,'gene2_2','chr2',2000,5000,'-','eg3','ee3');\n"
		"INSERT INTO exon VALUES (4,'gene3_1','chr3',1000,3000,'+','eg4','ee4');\n"
		"INSERT INTO exon VALUES (5,'gene3_2','chr3',5000,8000,'+','eg5','ee5');\n"
		"INSERT INTO exon VALUES (6,'gene4_1','chr4',1000,5000,'+','eg6','ee6');\n"
		"INSERT INTO exon VALUES (7,'gene4_2','chr5',1000,5000,'+','eg7','ee7');\n"
		"INSERT INTO exon VALUES (8,'gene5_1','chr6',1,3000,'-','eg8','ee8');\n"
		"INSERT INTO exon VALUES (9,'gene5_2','chr6',2000,5000,'-','eg9','ee9');\n"
		"INSERT INTO exon VALUES (10,'gene5_3','chr7',10000,13000,'+','eg10','ee10');\n"
		"INSERT INTO exon VALUES (11,'gene5_4','chr7',15000,18000,'+','eg11','ee11');\n"
		"INSERT INTO alignment VALUES (1,'q1',0x800,'chr10',1,20,'100M10S',100,100,'chr10',1,1,1);\n"
		"INSERT INTO alignment VALUES (12,'q1',0x800,'chr10',1,20,'100M10S',100,100,'chr10',1,8,1);\n"
		"INSERT INTO alignment VALUES (2,'q4',97,'chr11',1,20,'110M',100,50,'chr11',1,8,1);\n"
		"INSERT INTO alignment VALUES (3,'q5',97,'chr11',200,20,'110M',100,50,'chr11',1,8,1);\n"
		"INSERT INTO alignment VALUES (4,'q2',0x800,'chr12',250,20,'10H100M',100,100,'chr12',1,8,1);\n"
		"INSERT INTO alignment VALUES (5,'q3',0x800,'chr12',200,20,'100M10S',100,50,'chr1',1,8,1);\n"
		"INSERT INTO alignment VALUES (6,'q6',97,'chr13',1,20,'100M10S',100,50,'chr13',1,8,1);\n"
		"INSERT INTO alignment VALUES (7,'q7',97,'chr13',200,20,'100M10S',100,50,'chr13',1,8,1);\n"
		"INSERT INTO alignment VALUES (8,'q8',97,'chr14',1,20,'100M10S',100,50,'chr14',1,8,1);\n"
		"INSERT INTO alignment VALUES (9,'q9',97,'chr14',200,20,'100M10S',100,50,'chr14',1,8,1);\n"
		"INSERT INTO alignment VALUES (10,'q10',97,'chr14',400,20,'100M10S',100,50,'chr14',1,8,1);\n"
		"INSERT INTO alignment VALUES (11,'q11',97,'chr14',500,20,'100M10S',100,50,'chr14',1,8,1);\n"
		"INSERT INTO clustering VALUES (1,2,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,2,2,3,100);\n"
		"INSERT INTO clustering VALUES (3,2,3,3,100);\n"
		"INSERT INTO clustering VALUES (4,2,4,3,100);\n"
		"INSERT INTO clustering VALUES (5,2,5,3,100);\n"
		"INSERT INTO clustering VALUES (6,2,6,3,100);\n"
		"INSERT INTO clustering VALUES (7,2,7,3,100);\n"
		"INSERT INTO clustering VALUES (8,2,8,3,100);\n"
		"INSERT INTO clustering VALUES (9,2,9,3,100);\n"
		"INSERT INTO clustering VALUES (10,2,10,3,100);\n"
		"INSERT INTO clustering VALUES (11,2,11,3,100);\n"
		"INSERT INTO cluster VALUES (1,2,'chr10',1,300,'gene1',31);\n"
		"INSERT INTO cluster VALUES (2,2,'chr11',1,300,'gene2_1',31);\n"
		"INSERT INTO cluster VALUES (3,2,'chr11',200,500,'gene2_2',31);\n"
		"INSERT INTO cluster VALUES (4,2,'chr12',1,300,'gene3_1',31);\n"
		"INSERT INTO cluster VALUES (5,2,'chr12',200,500,'gene3_2',31);\n"
		"INSERT INTO cluster VALUES (6,2,'chr13',1,300,'gene4_1',31);\n"
		"INSERT INTO cluster VALUES (7,2,'chr13',200,500,'gene4_2',31);\n"
		"INSERT INTO cluster VALUES (8,2,'chr14',1,300,'gene5_1',31);\n"
		"INSERT INTO cluster VALUES (9,2,'chr14',200,500,'gene5_2',31);\n"
		"INSERT INTO cluster VALUES (10,2,'chr14',400,600,'gene5_3',31);\n"
		"INSERT INTO cluster VALUES (11,2,'chr14',500,700,'gene5_4',31);\n"
		"COMMIT;";

	db_exec (db, schema);
}

START_TEST (test_retrocopy)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";

	sqlite3 *db = NULL;
	sqlite3_stmt *cluster_merging_stmt = NULL;
	sqlite3_stmt *retrocopy_stmt = NULL;

	db = create_db (db_file);
	cluster_merging_stmt = db_prepare_cluster_merging_stmt (db);
	retrocopy_stmt = db_prepare_retrocopy_stmt (db);

	// Populate database
	populate_db (db);

	retrocopy (retrocopy_stmt,
			cluster_merging_stmt);

	db_finalize (cluster_merging_stmt);
	db_finalize (retrocopy_stmt);
	db_close (db);

	xunlink (db_file);
}
END_TEST

Suite *
make_retrocopy_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Retrocopy");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_retrocopy);
	suite_add_tcase (s, tc_core);

	return s;
}
