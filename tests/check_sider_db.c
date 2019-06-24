#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/log.h"
#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/db.h"

START_TEST (test_db)
{
	log_set_quiet (1);

	char db_path[] = "/tmp/ponga.db.XXXXXX";
	int fd = xmkstemp (db_path);
	close (fd);

	sqlite3 *db = db_create (db_path);
	ck_assert (db != NULL);

	ck_assert (exists (db_path));

	db_close (db);
	xunlink (db_path);
}
END_TEST

Suite *
make_db_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("DB");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_db);
	suite_add_tcase (s, tc_core);

	return s;
}
