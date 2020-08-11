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
#include "../src/db.h"
#include "../src/dedup.h"

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
	static const char schema[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO alignment VALUES(1,'id1',66,'chr1',1,60,'100M',101,101,'chr2',1,10,1);\n"
		"INSERT INTO alignment VALUES(2,'id1',66,'chr1',1,60,'100M',101,101,'chr2',1,10,1);\n"
		"INSERT INTO alignment VALUES(3,'id2',66,'chr1',1,60,'100M',101,101,'chr2',1,10,1);\n"
		"INSERT INTO alignment VALUES(4,'id3',66,'chr1',1,60,'100M',101,101,'chr3',1,10,1);\n"
		"INSERT INTO alignment VALUES(5,'id4',66,'chr1',2,60,'100M',101,101,'chr2',10,10,1);\n"
		"INSERT INTO alignment VALUES(6,'id2',66,'chr2',1,60,'100M',101,101,'chr1',1,10,1);\n"
		"INSERT INTO alignment VALUES(7,'id3',66,'chr3',1,60,'100M',101,101,'chr1',1,10,1);\n"
		"INSERT INTO alignment VALUES(8,'id4',66,'chr2',10,60,'100M',101,101,'chr1',2,10,1);\n"
		"INSERT INTO alignment VALUES(9,'id5',66,'chr3',1,60,'100M',101,101,'chr3',1000,10,1);\n"
		"INSERT INTO alignment VALUES(10,'id6',66,'chr3',1,60,'100M',101,101,'chr3',1000,10,1);\n"
		"INSERT INTO alignment VALUES(11,'id5',66,'chr3',1000,60,'100M',101,101,'chr3',1,10,1);\n"
		"INSERT INTO alignment VALUES(12,'id6',66,'chr3',1000,60,'100M',101,101,'chr3',1,10,1);\n"
		"INSERT INTO alignment VALUES(13,'id7',66,'chr3',1,60,'100M',101,101,'chr3',1000,10,2);\n"
		"INSERT INTO alignment VALUES(14,'id7',66,'chr3',1000,60,'100M',101,101,'chr3',1,10,2);\n"
		"INSERT INTO alignment VALUES(15,'id8',66,'chr3',1,60,'100M',101,101,'chr3',1000,10,1);\n"
		"INSERT INTO alignment VALUES(16,'id8',66,'chr3',1000,60,'100M',101,101,'chr3',1,10,1);\n"
		"COMMIT;";

	db_exec (db, schema);
}

static sqlite3_stmt *
prepare_query_stmt (sqlite3 *db)
{
	const char sql[] =
		"SELECT DISTINCT qname FROM alignment WHERE type = 0";

	return db_prepare (db, sql);
}

START_TEST (test_dedup)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";

	sqlite3 *db = NULL;
	sqlite3_stmt *stmt = NULL;

	const char *qname = NULL;

	const char *true_positive_qnames[] = {
		"id2", "id6", "id8"
	};

	int i = 0;

	db = create_db (db_file);
	stmt = prepare_query_stmt (db);

	// Populate database
	populate_db (db);

	// Let's dedup
	dedup (db);

	for (i = 0; db_step (stmt) == SQLITE_ROW; i++)
		{
			qname = db_column_text (stmt, 0);
			ck_assert_str_eq (qname, true_positive_qnames[i]);
		}

	db_finalize (stmt);
	db_close (db);

	xunlink (db_file);
}
END_TEST

Suite *
make_dedup_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Dedup");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_dedup);
	suite_add_tcase (s, tc_core);

	return s;
}
