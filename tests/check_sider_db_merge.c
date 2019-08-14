#include "config.h"

#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/db.h"
#include "../src/db_merge.h"

static void
create_db (char *path)
{
	const char sql[] =
		"INSERT INTO batch VALUES(1,\"2019-02-31\");\n"
		"INSERT INTO source VALUES(1,1,\"ponga.bam\");\n"
		"INSERT INTO exon VALUES(1,\"g1\",\"chr1\",1,200,\"+\",\"ENG0066\",\"ENSE0066\");\n"
		"INSERT INTO alignment VALUES(1,\"r1\",99,\"chr1\",1,20,\"101M\",101,101,\"chr1\",200,0,1);\n"
		"INSERT INTO overlapping VALUES(1,1,1,101);";

	int fd = xmkstemp (path);
	close (fd);

	sqlite3 *db = db_create (path);
	db_exec (db, sql);
	db_close (db);
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

	db_close (db);

	for (i = 0; i < num_db; i++)
		xunlink (db_paths[i]);
}
END_TEST

Suite *
make_db_merge_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("DBMerge");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_db_merge);
	suite_add_tcase (s, tc_core);

	return s;
}
