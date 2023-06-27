/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2023 Thiago L. A. Miller <tmiller@mochsl.org.br
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

#include <assert.h>
#include "wrapper.h"
#include "db.h"
#include "kv.h"

struct _KV
{
	const char   *path;
	sqlite3      *db;
	sqlite3_stmt *in_stmt;
	sqlite3_stmt *del_stmt;
	sqlite3_stmt *test_stmt;
	sqlite3_stmt *all_stmt;
	sqlite3_stmt *count_stmt;
};

static void
kv_create_table (sqlite3 *db)
{
	const char sql[] =
		"CREATE TABLE IF NOT EXISTS kv (\n"
		" key TEXT PRIMARY KEY,\n"
		" value BLOB NOT NULL)";

	db_exec (db, sql);
}

static sqlite3 *
kv_create_db (const char *path)
{
	sqlite3 *db = NULL;

	db = db_open (path,
			SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE);

	kv_create_table (db);

	return db;
}

static sqlite3_stmt *
kv_prepare_in_stmt (sqlite3 *db)
{
	const char sql[] =
		"INSERT OR REPLACE INTO kv (key,value)\n"
		" VALUES (?1,?2)";

	return db_prepare (db, sql);
}

static sqlite3_stmt *
kv_prepare_del_stmt (sqlite3 *db)
{
	const char sql[] =
		"DELETE FROM kv\n"
		" WHERE key = ?1";

	return db_prepare (db, sql);
}

static sqlite3_stmt *
kv_prepare_test_stmnt (sqlite3 *db)
{
	const char sql[] =
		"SELECT value\n"
		"FROM kv\n"
		"WHERE key = ?1";

	return db_prepare (db, sql);
}

static sqlite3_stmt *
kv_prepare_all_stmt (sqlite3 * db)
{
	const char sql[] =
		"SELECT key,value\n"
		"FROM kv";

	return db_prepare (db, sql);
}

static sqlite3_stmt *
kv_prepare_count_stmt (sqlite3 * db)
{
	const char sql[] =
		"SELECT COUNT(*)\n"
		"FROM kv";

	return db_prepare (db, sql);
}

KV *
kv_new (const char *path)
{
	assert (path != NULL);

	KV *kv = NULL;

	kv = xcalloc (1, sizeof (KV));

	kv->path       = (const char *) xstrdup (path);
	kv->db         = kv_create_db (path);
	kv->in_stmt    = kv_prepare_in_stmt (kv->db);
	kv->del_stmt   = kv_prepare_del_stmt (kv->db);
	kv->test_stmt  = kv_prepare_test_stmnt (kv->db);
	kv->all_stmt   = kv_prepare_all_stmt (kv->db);
	kv->count_stmt = kv_prepare_count_stmt (kv->db);

	return kv;
}

void
kv_free (KV *kv)
{
	if (kv == NULL)
		return;

	db_finalize (kv->in_stmt);
	db_finalize (kv->del_stmt);
	db_finalize (kv->test_stmt);
	db_finalize (kv->all_stmt);
	db_finalize (kv->count_stmt);

	db_close (kv->db);

	xfree ((char *) kv->path);
	xfree (kv);
}

const char *
kv_path (KV *kv)
{
	assert (kv != NULL);
	return kv->path;
}

int
kv_count (KV *kv)
{
	assert (kv != NULL);

	int count = 0;

	db_reset (kv->count_stmt);

	if (db_step (kv->count_stmt) == SQLITE_ROW)
		count = db_column_int (kv->count_stmt, 0);

	return count;
}

void
kv_insert (KV *kv, const char *key, const void *value, int n)
{
	assert (kv != NULL && key != NULL);

	int default_value = 1;

	if (value == NULL)
		{
			value = &default_value;
			n = sizeof (int);
		}

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (kv->in_stmt)));

	db_reset (kv->in_stmt);
	db_clear_bindings (kv->in_stmt);

	db_bind_text (kv->in_stmt, 1, key);
	db_bind_blob (kv->in_stmt, 2, value, n);

	db_step (kv->in_stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (kv->in_stmt)));
}

void
kv_del_key (KV *kv, const char *key)
{
	assert (kv != NULL && key != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (kv->del_stmt)));

	db_reset (kv->del_stmt);
	db_clear_bindings (kv->del_stmt);

	db_bind_text (kv->del_stmt, 1, key);

	db_step (kv->del_stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (kv->del_stmt)));
}

const void *
kv_get_value (KV *kv, const char *key)
{
	assert (kv != NULL && key != NULL);

	const void *value = NULL;

	db_reset (kv->test_stmt);
	db_clear_bindings (kv->test_stmt);

	db_bind_text (kv->test_stmt, 1, key);

	if (db_step (kv->test_stmt) == SQLITE_ROW)
		value = db_column_blob (kv->test_stmt, 0);

	return value;
}

void
kv_foreach (KV *kv, KVFunc func, void *user_data)
{
	assert (kv != NULL && func != NULL);

	const char *key = NULL;
	const void *value = NULL;

	db_reset (kv->all_stmt);

	while (db_step (kv->all_stmt) == SQLITE_ROW)
		{
			key   = db_column_text (kv->all_stmt, 0);
			value = db_column_blob (kv->all_stmt, 1);

			func (key, value, user_data);
		}
}
