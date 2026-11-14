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

#include <stdlib.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "db.h"

/* Low-level wrapper for sqlite3 interface */

static const char *
db_error_hint (int rc)
{
	// The primary result code is at the lower 8 bits
	switch (rc & 0xff)
		{
		case SQLITE_BUSY:
		case SQLITE_LOCKED:
			return "The database is in use by another process. Wait for it to "
				"finish, or raise the waiting time with the environment variable "
				DB_BUSY_TIMEOUT_ENV " (in milliseconds)";
		case SQLITE_IOERR:
			return "Input/output error: Check the storage. Network filesystems "
				"(NFS, Lustre, GPFS) are prone to it - prefer a local disk for "
				"the database";
		case SQLITE_FULL:
			return "The disk is full: Free some space or use another "
				"output directory";
		case SQLITE_CORRUPT:
		case SQLITE_NOTADB:
			return "The database is damaged or it is not a sideRETRO database: "
				"Generate it again";
		case SQLITE_CANTOPEN:
		case SQLITE_PERM:
		case SQLITE_READONLY:
			return "Check the path and its permissions";
		default:
			return NULL;
		}
}

static void
db_fatal (sqlite3 *db, const char *path, int rc, const char *what,
		const char *detail)
{
	const char *hint = NULL;
	int ext_rc = rc;

	if (db != NULL)
		{
			ext_rc = sqlite3_extended_errcode (db);

			if (path == NULL)
				path = sqlite3_db_filename (db, "main");
		}

	hint = db_error_hint (ext_rc);

	log_fatal ("%s failed at database '%s': %s (SQLite error %d)%s%s%s%s",
			what, path != NULL ? path : "?",
			db != NULL ? sqlite3_errmsg (db) : sqlite3_errstr (rc), ext_rc,
			detail != NULL ? "\n" : "", detail != NULL ? detail : "",
			hint != NULL ? "\n" : "", hint != NULL ? hint : "");
}

int
db_busy_timeout (void)
{
	const char *env = getenv (DB_BUSY_TIMEOUT_ENV);
	char *end = NULL;
	long ms = 0;

	if (env == NULL || *env == '\0')
		return DB_DEFAULT_BUSY_TIMEOUT;

	ms = strtol (env, &end, 10);

	if (*end != '\0' || ms < 0 || ms > 86400000)
		{
			log_warn ("Invalid %s='%s': Use the default of %d ms",
					DB_BUSY_TIMEOUT_ENV, env, DB_DEFAULT_BUSY_TIMEOUT);
			return DB_DEFAULT_BUSY_TIMEOUT;
		}

	return (int) ms;
}

sqlite3 *
db_open (const char *path, int flags)
{
	assert (path != NULL);

	sqlite3 *db = NULL;
	int rc = 0;

	rc = sqlite3_open_v2 (path, &db, flags, NULL);

	if (rc != SQLITE_OK)
		db_fatal (db, path, rc, "Opening", NULL);

	/*
	* A locked database makes SQLite return BUSY at once.
	* Instead, let it wait and retry for a while: Locks
	* are usually brief
	*/
	sqlite3_busy_timeout (db, db_busy_timeout ());

	return db;
}

void
db_close (sqlite3 *db)
{
	assert (db != NULL);

	int rc = sqlite3_close (db);

	if (rc != SQLITE_OK)
		db_fatal (db, NULL, rc, "Closing", NULL);
}

void
db_exec (sqlite3 *db, const char *sql)
{
	assert (db != NULL && sql != NULL);

	int rc = 0;

	rc = sqlite3_exec (db, sql, NULL, NULL, NULL);

	if (rc != SQLITE_OK)
		db_fatal (db, NULL, rc, "Executing", sql);
}

sqlite3_stmt *
db_prepare (sqlite3 *db, const char *sql)
{
	assert (db != NULL && sql != NULL);

	sqlite3_stmt *stmt = NULL;
	int rc = 0;

	rc = sqlite3_prepare_v2 (db, sql, -1, &stmt, NULL);

	if (rc != SQLITE_OK)
		db_fatal (db, NULL, rc, "Preparing", sql);

	return stmt;
}

int
db_step (sqlite3_stmt *stmt)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_step (stmt);

	if (rc != SQLITE_DONE && rc != SQLITE_ROW)
		db_fatal (sqlite3_db_handle (stmt), NULL, rc, "Running", sqlite3_sql (stmt));

	return rc;
}

void
db_finalize (sqlite3_stmt *stmt)
{
	assert (stmt != NULL);

	sqlite3 *db = sqlite3_db_handle (stmt);
	int rc = sqlite3_finalize (stmt);

	if (rc != SQLITE_OK)
		db_fatal (db, NULL, rc, "Finalizing a statement", NULL);
}

void
db_reset (sqlite3_stmt *stmt)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_reset (stmt);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_reset: %s",
				sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

void
db_clear_bindings (sqlite3_stmt *stmt)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_clear_bindings (stmt);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_clear_bindings: %s", // LCOV_EXCL_LINE: it always returns SQLITE_OK
				sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

void
db_bind_int (sqlite3_stmt *stmt,
		int i, int value)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_bind_int (stmt, i, value);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_bind_int '[%d] %d': %s",
				i, value, sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

void
db_bind_int64 (sqlite3_stmt *stmt,
		int i, int64_t value)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_bind_int64 (stmt, i, value);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_bind_int64 at '[%d] %jd': %s",
				i, (intmax_t) value,
				sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

void
db_bind_double (sqlite3_stmt *stmt,
		int i, double value)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_bind_double (stmt, i, value);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_bind_double at '[%d] %f': %s",
				i, value, sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

void
db_bind_text (sqlite3_stmt *stmt,
		int i, const char *value)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_bind_text (stmt, i, value, -1, SQLITE_TRANSIENT);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_bind_text at '[%d] %s': %s",
				i, value, sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

void
db_bind_null (sqlite3_stmt *stmt, int i)
{
	assert (stmt != NULL);

	int rc = 0;

	rc = sqlite3_bind_null (stmt, i);

	if (rc != SQLITE_OK)
		log_fatal ("Failed sqlite3_bind_null at '[%d]': %s",
				i, sqlite3_errmsg (sqlite3_db_handle (stmt)));
}

int
db_column_int (sqlite3_stmt *stmt, int i)
{
	assert (stmt != NULL);

	if (sqlite3_column_type (stmt, i) != SQLITE_INTEGER)
		log_fatal ("Failed sqlite3_column_int at %d: Wrong type", i);

	return sqlite3_column_int (stmt, i);
}

int64_t
db_column_int64 (sqlite3_stmt *stmt, int i)
{
	assert (stmt != NULL);

	if (sqlite3_column_type (stmt, i) != SQLITE_INTEGER)
		log_fatal ("Failed sqlite3_column_int64 at %d: Wrong type", i);

	return sqlite3_column_int64 (stmt, i);
}

double
db_column_double (sqlite3_stmt *stmt, int i)
{
	assert (stmt != NULL);

	if (sqlite3_column_type (stmt, i) != SQLITE_FLOAT)
		log_fatal ("Failed sqlite3_column_double at %d: Wrong type", i);

	return sqlite3_column_double (stmt, i);
}

const char *
db_column_text (sqlite3_stmt *stmt, int i)
{
	assert (stmt != NULL);

	if (sqlite3_column_type (stmt, i) != SQLITE_TEXT)
		log_fatal ("Failed sqlite3_column_text at %d: Wrong type", i);

	return (const char *) sqlite3_column_text (stmt, i);
}

/* db management functions */

static void
db_create_tables (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"DROP TABLE IF EXISTS schema;\n"
		"CREATE TABLE schema (\n"
		"	major_version INTEGER NOT NULL,\n"
		"	minor_version INTEGER NOT NULL);\n"
		"\n"
		"DROP TABLE IF EXISTS batch;\n"
		"CREATE TABLE batch (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	timestamp TEXT NOT NULL);\n"
		"\n"
		"DROP TABLE IF EXISTS source;\n"
		"CREATE TABLE source (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	batch_id INTEGER NOT NULL,\n"
		"	path TEXT NOT NULL,\n"
		"	FOREIGN KEY (batch_id) REFERENCES batch(id));\n"
		"\n"
		"DROP TABLE IF EXISTS exon;\n"
		"CREATE TABLE exon (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	gene_name TEXT NOT NULL,\n"
		"	chr TEXT NOT NULL,\n"
		"	start INTEGER NOT NULL,\n"
		"	end INTEGER NOT NULL,\n"
		"	strand TEXT NOT NULL,\n"
		"	ensg TEXT NOT NULL,\n"
		"	ense TEXT NOT NULL,\n"
		"	UNIQUE (ense));\n"
		"\n"
		"DROP TABLE IF EXISTS alignment;\n"
		"CREATE TABLE alignment (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	qname TEXT NOT NULL,\n"
		"	flag INTEGER NOT NULL,\n"
		"	chr TEXT NOT NULL,\n"
		"	pos INTEGER NOT NULL,\n"
		"	mapq INTEGER NOT NULL,\n"
		"	cigar TEXT NOT NULL,\n"
		"	qlen INTEGER DEFAULT -1,\n"
		"	rlen INTEGER DEFAULT -1,\n"
		"	chr_next TEXT NOT NULL,\n"
		"	pos_next INTEGER NOT NULL,\n"
		"	type INT DEFAULT 0,\n"
		"	source_id INTEGER NOT NULL,\n"
		"	seq TEXT NOT NULL,\n"
		"	qual TEXT NOT NULL,\n"
		"	FOREIGN KEY (source_id) REFERENCES source(id));\n"
		"\n"
		"DROP TABLE IF EXISTS overlapping;\n"
		"CREATE TABLE overlapping (\n"
		"	exon_id INTEGER NOT NULL,\n"
		"	alignment_id INTEGER NOT NULL,\n"
		"	pos INTEGER NOT NULL,\n"
		"	len INTEGER NOT NULL,\n"
		"	FOREIGN KEY (exon_id) REFERENCES exon(id),\n"
		"	FOREIGN KEY (alignment_id) REFERENCES alignment(id),\n"
		"	PRIMARY KEY (exon_id, alignment_id));\n"
		"\n"
		"DROP TABLE IF EXISTS cluster;\n"
		"CREATE TABLE cluster (\n"
		"	id INTEGER NOT NULL,\n"
		"	sid INTEGER NOT NULL,\n"
		"	chr TEXT NOT NULL,\n"
		"	start INTEGER NOT NULL,\n"
		"	end INTEGER NOT NULL,\n"
		"	gene_name TEXT NOT NULL,\n"
		"	filter INTEGER NOT NULL,\n"
		"	PRIMARY KEY (id, sid));\n"
		"\n"
		"DROP TABLE IF EXISTS clustering;\n"
		"CREATE TABLE clustering (\n"
		"	cluster_id INTEGER NOT NULL,\n"
		"	cluster_sid INTEGER NOT NULL,\n"
		"	alignment_id INTEGER NOT NULL,\n"
		"	label INTEGER NOT NULL,\n"
		"	neighbors INTEGER NOT NULL,\n"
		"	FOREIGN KEY (cluster_id, cluster_sid) REFERENCES cluster(id, sid),\n"
		"	FOREIGN KEY (alignment_id) REFERENCES alignment(id),\n"
		"	PRIMARY KEY (cluster_id, cluster_sid, alignment_id));\n"
		"\n"
		"DROP TABLE IF EXISTS blacklist;\n"
		"CREATE TABLE blacklist (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	name TEXT NOT NULL,\n"
		"	chr TEXT NOT NULL,\n"
		"	start INTEGER NOT NULL,\n"
		"	end INTEGER NOT NULL);\n"
		"\n"
		"DROP TABLE IF EXISTS overlapping_blacklist;\n"
		"CREATE TABLE overlapping_blacklist (\n"
		"	blacklist_id INTEGER NOT NULL,\n"
		"	cluster_id INTEGER NOT NULL,\n"
		"	cluster_sid INTEGER NOT NULL,\n"
		"	pos INTEGER NOT NULL,\n"
		"	len INTEGER NOT NULL,\n"
		"	FOREIGN KEY (blacklist_id) REFERENCES blacklist(id),\n"
		"	FOREIGN KEY (cluster_id, cluster_sid) REFERENCES cluster(id, sid),\n"
		"	PRIMARY KEY (blacklist_id, cluster_id, cluster_sid));\n"
		"\n"
		"DROP TABLE IF EXISTS retrocopy;\n"
		"CREATE TABLE retrocopy (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	chr TEXT NOT NULL,\n"
		"	window_start INTEGER NOT NULL,\n"
		"	window_end INTEGER NOT NULL,\n"
		"	parental_gene_name TEXT NOT NULL,\n"
		"	level INTEGER NOT NULL,\n"
		"	insertion_point INTEGER,\n"
		"	insertion_point_type INTEGER,\n"
		"	orientation_rho REAL,\n"
		"	orientation_p_value REAL,\n"
		"	breakpoint_left INTEGER,\n"
		"	breakpoint_left_sr INTEGER NOT NULL DEFAULT 0,\n"
		"	breakpoint_right INTEGER,\n"
		"	breakpoint_right_sr INTEGER NOT NULL DEFAULT 0,\n"
		"	polya_len INTEGER,\n"
		"	polya_lower_bound INTEGER NOT NULL DEFAULT 0,\n"
		"	polya_strand TEXT,\n"
		"	polya_sr INTEGER NOT NULL DEFAULT 0);\n"
		"\n"
		"DROP TABLE IF EXISTS cluster_merging;\n"
		"CREATE TABLE cluster_merging (\n"
		"	retrocopy_id INTEGER NOT NULL,\n"
		"	cluster_id INTEGER NOT NULL,\n"
		"	cluster_sid INTEGER NOT NULL,\n"
		"	FOREIGN KEY (retrocopy_id) REFERENCES retrocopy(id),\n"
		"	FOREIGN KEY (cluster_id, cluster_sid) REFERENCES cluster(id, sid),\n"
		"	PRIMARY KEY (retrocopy_id, cluster_id, cluster_sid));\n"
		"\n"
		"DROP TABLE IF EXISTS junction;\n"
		"CREATE TABLE junction (\n"
		"	retrocopy_id INTEGER NOT NULL,\n"
		"	exon_id1 INTEGER NOT NULL,\n"
		"	exon_id2 INTEGER NOT NULL,\n"
		"	reads INTEGER NOT NULL,\n"
		"	FOREIGN KEY (retrocopy_id) REFERENCES retrocopy(id),\n"
		"	FOREIGN KEY (exon_id1) REFERENCES exon(id),\n"
		"	FOREIGN KEY (exon_id2) REFERENCES exon(id),\n"
		"	PRIMARY KEY (retrocopy_id, exon_id1, exon_id2));\n"
		"\n"
		"DROP TABLE IF EXISTS insertion;\n"
		"CREATE TABLE insertion (\n"
		"	retrocopy_id INTEGER PRIMARY KEY,\n"
		"	five_prime_pos INTEGER,\n"
		"	five_prime_exon_id INTEGER,\n"
		"	five_prime_sr INTEGER NOT NULL DEFAULT 0,\n"
		"	three_prime_pos INTEGER,\n"
		"	three_prime_exon_id INTEGER,\n"
		"	three_prime_sr INTEGER NOT NULL DEFAULT 0,\n"
		"	len INTEGER,\n"
		"	seq TEXT,\n"
		"	left_seq TEXT,\n"
		"	right_seq TEXT,\n"
		"	FOREIGN KEY (retrocopy_id) REFERENCES retrocopy(id),\n"
		"	FOREIGN KEY (five_prime_exon_id) REFERENCES exon(id),\n"
		"	FOREIGN KEY (three_prime_exon_id) REFERENCES exon(id));\n"
		"\n"
		"DROP TABLE IF EXISTS evidence;\n"
		"CREATE TABLE evidence (\n"
		"	retrocopy_id INTEGER NOT NULL,\n"
		"	alignment_id INTEGER NOT NULL,\n"
		"	role INTEGER NOT NULL,\n"
		"	FOREIGN KEY (retrocopy_id) REFERENCES retrocopy(id),\n"
		"	FOREIGN KEY (alignment_id) REFERENCES alignment(id),\n"
		"	PRIMARY KEY (retrocopy_id, alignment_id, role));\n"
		"\n"
		"DROP TABLE IF EXISTS genotype;\n"
		"CREATE TABLE genotype (\n"
		"	source_id INTEGER NOT NULL,\n"
		"	retrocopy_id INTEGER NOT NULL,\n"
		"	reference_depth INTEGER NOT NULL,\n"
		"	alternate_depth INTEGER NOT NULL,\n"
		"	ho_ref_likelihood REAL NOT NULL,\n"
		"	he_likelihood REAL NOT NULL,\n"
		"	ho_alt_likelihood NOT NULL,\n"
		"	FOREIGN KEY (source_id) REFERENCES source(id),\n"
		"	FOREIGN KEY (retrocopy_id) REFERENCES retrocopy(id),\n"
		"	PRIMARY KEY (source_id, retrocopy_id));";

	log_debug ("Database schema:\n%s", sql);
	db_exec (db, sql);
}

static void
db_insert_schema_version (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	sqlite3_stmt *stmt = NULL;

	stmt = db_prepare (db,
		"INSERT INTO schema (major_version,minor_version)\n"
		"VALUES (?1,?2)");

	db_bind_int (stmt, 1, DB_SCHEMA_MAJOR_VERSION);
	db_bind_int (stmt, 2, DB_SCHEMA_MINOR_VERSION);
	db_step (stmt);

	db_finalize (stmt);
}

sqlite3 *
db_create (const char *path)
{
	log_trace ("Inside %s", __func__);
	assert (path != NULL);

	sqlite3 *db = NULL;

	log_debug ("Create and open database '%s'", path);
	db = db_open (path, SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE);

	log_debug ("Create tables into database '%s'", path);
	db_create_tables (db);

	log_debug ("Insert schema version 'v%d.%d' into database '%s'",
			DB_SCHEMA_MAJOR_VERSION, DB_SCHEMA_MINOR_VERSION, path);
	db_insert_schema_version (db);

	return db;
}

void
db_cache_size (sqlite3 *db, size_t size)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && size > 0);

	char sql[64] = {};

	if (size < DB_DEFAULT_CACHE_SIZE)
		log_warn ("Cache size of %zuKiB is smaller than the default of %uKiB",
				size, DB_DEFAULT_CACHE_SIZE);

	xsnprintf (sql, 63, "PRAGMA cache_size = -%zu", size);
	sql[63] = '\0';

	db_exec (db, sql);
}

void
db_begin_transaction (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	// Take the write lock at once: A deferred transaction that
	// becomes a write one may get BUSY without waiting
	db_exec (db, "BEGIN IMMEDIATE TRANSACTION");
}

void
db_end_transaction (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	db_exec (db, "END TRANSACTION");
}

static void
db_check_schema_version (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	sqlite3_stmt *stmt = NULL;
	int major_version = 0;
	int minor_version = 0;

	stmt = db_prepare (db,
			"SELECT major_version,minor_version FROM schema LIMIT 1");

	if (db_step (stmt) == SQLITE_ROW)
		{
			major_version = db_column_int (stmt, 0);
			minor_version = db_column_int (stmt, 1);
		}
	else
		log_fatal ("Failed to request schema version");

	if ((major_version > DB_SCHEMA_MAJOR_VERSION)
			|| (major_version == DB_SCHEMA_MAJOR_VERSION
				&& minor_version > DB_SCHEMA_MINOR_VERSION))
		{
			log_fatal ("Schema version 'v%d.%d' of database '%s' is newer than the version 'v%d.%d' of '%s': Use a newer sideRETRO",
				major_version, minor_version, sqlite3_db_filename (db, "main"),
				DB_SCHEMA_MAJOR_VERSION, DB_SCHEMA_MINOR_VERSION, PACKAGE_STRING);
		}
	else if ((major_version < DB_SCHEMA_MAJOR_VERSION)
			|| (major_version == DB_SCHEMA_MAJOR_VERSION
				&& minor_version < DB_SCHEMA_MINOR_VERSION))
		{
			log_fatal ("Schema version 'v%d.%d' of database '%s' is older than the version 'v%d.%d' of '%s': Generate it again with 'process-sample'",
				major_version, minor_version, sqlite3_db_filename (db, "main"),
				DB_SCHEMA_MAJOR_VERSION, DB_SCHEMA_MINOR_VERSION, PACKAGE_STRING);
		}

	db_finalize (stmt);
}

sqlite3 *
db_connect (const char *path)
{
	log_trace ("Inside %s", __func__);
	assert (path != NULL);

	sqlite3 *db = NULL;

	log_debug ("Connect to database '%s'", path);
	db = db_open (path, SQLITE_OPEN_READWRITE);

	log_debug ("Check database schema version for '%s'", path);
	db_check_schema_version (db);

	return db;
}

sqlite3_stmt *
db_prepare_exon_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO exon (id,gene_name,chr,start,end,strand,ensg,ense)"
		"	VALUES (?1,?2,?3,?4,?5,?6,?7,?8)";

	return db_prepare (db, sql);
}

void
db_insert_exon (sqlite3_stmt *stmt, int id, const char *gene_name,
		const char *chr, long start, long end, const char *strand, const char *ensg,
		const char *ense)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL && gene_name != NULL && chr != NULL
			&& strand != NULL && ensg != NULL && ense != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_text (stmt, 2, gene_name);
	db_bind_text (stmt, 3, chr);
	db_bind_int64 (stmt, 4, start);
	db_bind_int64 (stmt, 5, end);
	db_bind_text (stmt, 6, strand);
	db_bind_text (stmt, 7, ensg);
	db_bind_text (stmt, 8, ense);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_batch_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO batch (id,timestamp) VALUES (?1,?2)";

	return db_prepare (db, sql);
}

void
db_insert_batch (sqlite3_stmt *stmt, int id,
		const char *timestamp)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL && timestamp != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_text (stmt, 2, timestamp);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_source_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO source (id,batch_id,path) VALUES (?1,?2,?3)";

	return db_prepare (db, sql);
}

void
db_insert_source (sqlite3_stmt *stmt, int id,
		int batch_id, const char *path)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL && path != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_int (stmt, 2, batch_id);
	db_bind_text (stmt, 3, path);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_alignment_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO alignment (id,qname,flag,chr,pos,mapq,cigar,qlen,rlen,chr_next,pos_next,type,source_id,seq,qual)\n"
		"VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15)";

	return db_prepare (db, sql);
}

void
db_insert_alignment (sqlite3_stmt *stmt, int id, const char *qname, int flag,
		const char *chr, long pos, int mapq, const char *cigar, int qlen, int rlen,
		const char *chr_next, long pos_next, int type, int source_id,
		const char *seq, const char *qual)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL && qname != NULL && chr != NULL
			&& cigar != NULL && chr_next != NULL
			&& seq != NULL && qual != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_text (stmt, 2, qname);
	db_bind_int (stmt, 3, flag);
	db_bind_text (stmt, 4, chr);
	db_bind_int64 (stmt, 5, pos);
	db_bind_int (stmt, 6, mapq);
	db_bind_text (stmt, 7, cigar);
	db_bind_int (stmt, 8, qlen);
	db_bind_int (stmt, 9, rlen);
	db_bind_text (stmt, 10, chr_next);
	db_bind_int64 (stmt, 11, pos_next);
	db_bind_int (stmt, 12, type);
	db_bind_int (stmt, 13, source_id);
	db_bind_text (stmt, 14, seq);
	db_bind_text (stmt, 15, qual);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_overlapping_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO overlapping (exon_id,alignment_id,pos,len)\n"
		"VALUES (?1,?2,?3,?4)";

	return db_prepare (db, sql);
}

void
db_insert_overlapping (sqlite3_stmt *stmt, int exon_id,
	int alignment_id, long pos, long len)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, exon_id);
	db_bind_int (stmt, 2, alignment_id);
	db_bind_int64 (stmt, 3, pos);
	db_bind_int64 (stmt, 4, len);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_clustering_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO clustering (cluster_id,cluster_sid,alignment_id,label,neighbors)\n"
		"VALUES (?1,?2,?3,?4,?5)";

	return db_prepare (db, sql);
}

void
db_insert_clustering (sqlite3_stmt *stmt, int cluster_id, int cluster_sid,
		int alignment_id, int label, int neighbors)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, cluster_id);
	db_bind_int (stmt, 2, cluster_sid);
	db_bind_int (stmt, 3, alignment_id);
	db_bind_int (stmt, 4, label);
	db_bind_int (stmt, 5, neighbors);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_cluster_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO cluster (id,sid,chr,start,end,gene_name,filter)\n"
		"VALUES (?1,?2,?3,?4,?5,?6,?7)";

	return db_prepare (db, sql);
}

void
db_insert_cluster (sqlite3_stmt *stmt, int id, int sid, const char *chr,
		long start, long end, const char *gene_name, int filter)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_int (stmt, 2, sid);
	db_bind_text (stmt, 3, chr);
	db_bind_int64 (stmt, 4, start);
	db_bind_int64 (stmt, 5, end);
	db_bind_text (stmt, 6, gene_name);
	db_bind_int (stmt, 7, filter);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_blacklist_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO blacklist (id,name,chr,start,end)"
		"	VALUES (?1,?2,?3,?4,?5)";

	return db_prepare (db, sql);
}

void
db_insert_blacklist (sqlite3_stmt *stmt, int id, const char *name,
		const char *chr, long start, long end)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL && name != NULL && chr != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_text (stmt, 2, name);
	db_bind_text (stmt, 3, chr);
	db_bind_int64 (stmt, 4, start);
	db_bind_int64 (stmt, 5, end);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_overlapping_blacklist_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO overlapping_blacklist (blacklist_id,cluster_id,cluster_sid,pos,len)\n"
		"VALUES (?1,?2,?3,?4,?5)";

	return db_prepare (db, sql);
}

void
db_insert_overlapping_blacklist (sqlite3_stmt *stmt, int blacklist_id,
	int cluster_id, int cluster_sid, long pos, long len)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, blacklist_id);
	db_bind_int (stmt, 2, cluster_id);
	db_bind_int (stmt, 3, cluster_sid);
	db_bind_int64 (stmt, 4, pos);
	db_bind_int64 (stmt, 5, len);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_cluster_merging_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO cluster_merging (retrocopy_id,cluster_id,cluster_sid)\n"
		"VALUES (?1,?2,?3)";

	return db_prepare (db, sql);
}

void
db_insert_cluster_merging (sqlite3_stmt *stmt, int retrocopy_id,
		int cluster_id, int cluster_sid)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, retrocopy_id);
	db_bind_int (stmt, 2, cluster_id);
	db_bind_int (stmt, 3, cluster_sid);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_retrocopy_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO retrocopy (id,chr,window_start,window_end,parental_gene_name,level,\n"
		"	insertion_point,insertion_point_type,orientation_rho,orientation_p_value,\n"
		"	breakpoint_left,breakpoint_left_sr,breakpoint_right,breakpoint_right_sr,\n"
		"	polya_len,polya_lower_bound,polya_strand,polya_sr)\n"
		"VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18)";

	return db_prepare (db, sql);
}

void
db_insert_retrocopy (sqlite3_stmt *stmt, int id, const char *chr, long window_start,
		long window_end, const char *parental_gene_name, int level, long insertion_point,
		int insertion_point_type, double orientation_rho, double orientation_p_value,
		long breakpoint_left, int breakpoint_left_sr, long breakpoint_right,
		int breakpoint_right_sr, int polya_len, int polya_lower_bound,
		const char *polya_strand, int polya_sr)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int    (stmt, 1, id);
	db_bind_text   (stmt, 2, chr);
	db_bind_int64  (stmt, 3, window_start);
	db_bind_int64  (stmt, 4, window_end);
	db_bind_text   (stmt, 5, parental_gene_name);
	db_bind_int    (stmt, 6, level);
	db_bind_int64  (stmt, 7, insertion_point);
	db_bind_int    (stmt, 8, insertion_point_type);
	db_bind_double (stmt, 9, orientation_rho);
	db_bind_double (stmt, 10, orientation_p_value);

	// Breakpoints are unknown without split reads
	if (breakpoint_left_sr > 0)
		db_bind_int64 (stmt, 11, breakpoint_left);
	else
		db_bind_null (stmt, 11);

	db_bind_int (stmt, 12, breakpoint_left_sr);

	if (breakpoint_right_sr > 0)
		db_bind_int64 (stmt, 13, breakpoint_right);
	else
		db_bind_null (stmt, 13);

	db_bind_int (stmt, 14, breakpoint_right_sr);

	// Poly-A is unknown without reads crossing its junction
	if (polya_sr > 0)
		{
			db_bind_int (stmt, 15, polya_len);
			db_bind_text (stmt, 17, polya_strand);
		}
	else
		{
			db_bind_null (stmt, 15);
			db_bind_null (stmt, 17);
		}

	db_bind_int (stmt, 16, polya_lower_bound);
	db_bind_int (stmt, 18, polya_sr);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_junction_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO junction (retrocopy_id,exon_id1,exon_id2,reads)\n"
		"VALUES (?1,?2,?3,?4)";

	return db_prepare (db, sql);
}

void
db_insert_junction (sqlite3_stmt *stmt, int retrocopy_id,
		int exon_id1, int exon_id2, int reads)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, retrocopy_id);
	db_bind_int (stmt, 2, exon_id1);
	db_bind_int (stmt, 3, exon_id2);
	db_bind_int (stmt, 4, reads);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_insertion_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO insertion (retrocopy_id,five_prime_pos,five_prime_exon_id,five_prime_sr,\n"
		"	three_prime_pos,three_prime_exon_id,three_prime_sr,len,seq,left_seq,right_seq)\n"
		"VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11)";

	return db_prepare (db, sql);
}

void
db_insert_insertion (sqlite3_stmt *stmt, int retrocopy_id,
		long five_prime_pos, int five_prime_exon_id, int five_prime_sr,
		long three_prime_pos, int three_prime_exon_id, int three_prime_sr,
		int len, const char *seq, const char *left_seq, const char *right_seq)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, retrocopy_id);

	// The ends are unknown without reads matching the exons
	if (five_prime_sr > 0)
		{
			db_bind_int64 (stmt, 2, five_prime_pos);
			db_bind_int (stmt, 3, five_prime_exon_id);
		}
	else
		{
			db_bind_null (stmt, 2);
			db_bind_null (stmt, 3);
		}

	db_bind_int (stmt, 4, five_prime_sr);

	if (three_prime_sr > 0)
		{
			db_bind_int64 (stmt, 5, three_prime_pos);
			db_bind_int (stmt, 6, three_prime_exon_id);
		}
	else
		{
			db_bind_null (stmt, 5);
			db_bind_null (stmt, 6);
		}

	db_bind_int (stmt, 7, three_prime_sr);

	// The predicted sequence needs the 5' end
	if (seq != NULL)
		{
			db_bind_int (stmt, 8, len);
			db_bind_text (stmt, 9, seq);
		}
	else
		{
			db_bind_null (stmt, 8);
			db_bind_null (stmt, 9);
		}

	// Observed sequences at the junctions
	if (left_seq != NULL)
		db_bind_text (stmt, 10, left_seq);
	else
		db_bind_null (stmt, 10);

	if (right_seq != NULL)
		db_bind_text (stmt, 11, right_seq);
	else
		db_bind_null (stmt, 11);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_evidence_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT OR IGNORE INTO evidence (retrocopy_id,alignment_id,role)\n"
		"VALUES (?1,?2,?3)";

	return db_prepare (db, sql);
}

void
db_insert_evidence (sqlite3_stmt *stmt, int retrocopy_id,
		int alignment_id, int role)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, retrocopy_id);
	db_bind_int (stmt, 2, alignment_id);
	db_bind_int (stmt, 3, role);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}

sqlite3_stmt *
db_prepare_genotype_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	const char sql[] =
		"INSERT INTO genotype (source_id,retrocopy_id,reference_depth,alternate_depth,ho_ref_likelihood,he_likelihood,ho_alt_likelihood)\n"
		"VALUES (?1,?2,?3,?4,?5,?6,?7)";

	return db_prepare (db, sql);
}

void
db_insert_genotype (sqlite3_stmt *stmt, int source_id, int retrocopy_id, int reference_depth,
		int alternate_depth, double ho_ref_likelihood, double he_likelihood, double ho_alt_likelihood)
{
	log_trace ("Inside %s", __func__);
	assert (stmt != NULL);

	sqlite3_mutex_enter (sqlite3_db_mutex (sqlite3_db_handle (stmt)));

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, source_id);
	db_bind_int (stmt, 2, retrocopy_id);
	db_bind_int (stmt, 3, reference_depth);
	db_bind_int (stmt, 4, alternate_depth);
	db_bind_double (stmt, 5, ho_ref_likelihood);
	db_bind_double (stmt, 6, he_likelihood);
	db_bind_double (stmt, 7, ho_alt_likelihood);

	db_step (stmt);

	sqlite3_mutex_leave (sqlite3_db_mutex (sqlite3_db_handle (stmt)));
}
