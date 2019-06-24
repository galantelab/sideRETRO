#include <stdlib.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "db.h"

#define DB_DEFAULT_CACHE_SIZE 2000

static void
db_create_tables (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	int rc = 0;
	char *err_msg = NULL;

	const char sql[] =
		"DROP TABLE IF EXISTS gene;\n"
		"CREATE TABLE gene (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	name TEXT NOT NULL,\n"
		"	chr TEXT NOT NULL,\n"
		"	start INTEGER NOT NULL,\n"
		"	end INTEGER NOT NULL,\n"
		"	strand TEXT NOT NULL,\n"
		"	ensembl_id TEXT NOT NULL,\n"
		"	UNIQUE (name, ensembl_id));\n"
		"\n"
		"DROP TABLE IF EXISTS alignment;\n"
		"CREATE TABLE alignment (\n"
		"	id INTEGER PRIMARY KEY,\n"
		"	name TEXT NOT NULL,\n"
		"	flag INTEGER NOT NULL,\n"
		"	chr TEXT NOT NULL,\n"
		"	pos INTEGER NOT NULL,\n"
		"	mapq INTEGER NOT NULL,\n"
		"	cigar TEXT NOT NULL,\n"
		"	qlen INTEGER DEFAULT -1,\n"
		"	rlen INTEGER DEFAULT -1,\n"
		"	chr_next TEXT NOT NULL,\n"
		"	pos_next INTEGER NOT NULL,\n"
		"	type INT DEFAULT 0);\n"
		"\n"
		"DROP TABLE IF EXISTS overlapping;\n"
		"CREATE TABLE overlapping (\n"
		"	gene_id INTEGER NOT NULL,\n"
		"	alignment_id INTEGER NOT NULL,\n"
		"	exons TEXT DEFAULT NULL,\n"
		"	FOREIGN KEY (gene_id) REFERENCES gene(id),\n"
		"	FOREIGN KEY (alignment_id) REFERENCES alignment(id),\n"
		"	PRIMARY KEY (gene_id, alignment_id));";

	log_debug ("Database schema:\n%s", sql);
	rc = sqlite3_exec (db, sql, NULL, NULL, &err_msg);

	if (rc != SQLITE_OK)
		log_fatal ("SQL error: %s", err_msg);
}

sqlite3 *
db_create (const char *path)
{
	log_trace ("Inside %s", __func__);
	assert (path != NULL);

	sqlite3 *db = NULL;
	int rc = 0;

	log_debug ("Create and open database '%s'", path);
	rc = sqlite3_open_v2 (path, &db,
			SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Can't create database '%s': %s",
				sqlite3_errmsg (db));

	log_debug ("Create tables into database '%s'", path);
	db_create_tables (db);

	return db;
}

void
db_cache_size (sqlite3 *db, size_t size)
{
	assert (db != NULL && size > 0);

	int rc = 0;
	int len = 0;
	char *err_msg = NULL;
	char sql[128];

	if (size < DB_DEFAULT_CACHE_SIZE)
		log_warn ("cache size of %dKiB is lesser than the default value of %dKiB",
				size, DB_DEFAULT_CACHE_SIZE);

	len = xsnprintf (sql, 128, "PRAGMA cache_size=-%d", size);
	rc = sqlite3_exec (db, sql, NULL, NULL, &err_msg);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to change default database cache size: %s",
				err_msg);
}

void
db_begin_transaction (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	int rc = 0;
	char *err_msg = NULL;

	rc = sqlite3_exec (db, "BEGIN TRANSACTION",
			NULL, NULL, &err_msg);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to begin database transaction: %s",
				err_msg);
}

void
db_end_transaction (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	int rc = 0;
	char *err_msg = NULL;

	rc = sqlite3_exec (db, "END TRANSACTION",
			NULL, NULL, &err_msg);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to end database transaction: %s",
				err_msg);
}

sqlite3 *
db_connect (const char *path)
{
	log_trace ("Inside %s", __func__);
	assert (path != NULL);

	sqlite3 *db = NULL;
	int rc = 0;

	rc = sqlite3_open_v2 (path, &db,
			SQLITE_OPEN_READWRITE, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Can't connect to database '%s': %s",
				sqlite3_errmsg (db));

	return db;
}

void
db_close (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	int rc = sqlite3_close (db);

	if (rc != SQLITE_OK)
		log_fatal ("Failed close database: %s",
				sqlite3_errmsg (db));
}

void
db_finalize (sqlite3 *db, sqlite3_stmt *stmt)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && stmt != NULL);

	int rc = sqlite3_finalize (stmt);

	if (rc != SQLITE_OK)
		log_fatal ("Failed finalize stmt: %s",
				sqlite3_errmsg (db));
}

sqlite3_stmt *
db_prepare_gene_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	int rc = 0;
	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"INSERT OR IGNORE INTO gene (id,name,chr,start,end,strand,ensembl_id)\n"
		"	VALUES (?1,?2,?3,?4,?5,?6,?7)";

	rc = sqlite3_prepare_v2 (db, sql, -1, &stmt, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to prepare gene stmt: %s",
				sqlite3_errmsg (db));

	return stmt;
}

void
db_insert_gene (sqlite3 *db, sqlite3_stmt *stmt, int id, const char *name,
		const char *chr, long start, long end, const char *strand, const char *ensembl_id)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && stmt != NULL && name != NULL && chr != NULL
			&& strand != NULL && ensembl_id != NULL);

	int rc = 0;

	rc = sqlite3_reset (stmt);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_clear_bindings (stmt);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 1, id);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 2, name, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 3, chr, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int64 (stmt, 4, start);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int64 (stmt, 5, end);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 6, strand, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 7, ensembl_id, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_step (stmt);
	if (rc != SQLITE_DONE)
		log_fatal ("Failed to insert gene data: %s",
				sqlite3_errmsg (db));
}

sqlite3_stmt *
db_prepare_alignment_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	int rc = 0;
	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"INSERT INTO alignment (id,name,flag,chr,pos,mapq,cigar,qlen,rlen,chr_next,pos_next,type)\n"
		"	VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12)";

	rc = sqlite3_prepare_v2 (db, sql, -1, &stmt, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to prepare alignment stmt: %s",
				sqlite3_errmsg (db));

	return stmt;
}

void
db_insert_alignment (sqlite3 *db, sqlite3_stmt *stmt, int id, const char *name,
		int flag, const char *chr, long pos, int mapq, const char *cigar, int qlen,
		int rlen, const char *chr_next, long pos_next, int type)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && stmt != NULL && name != NULL && chr != NULL
			&& cigar != NULL && chr_next != NULL);

	int rc = 0;

	rc = sqlite3_reset (stmt);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_clear_bindings (stmt);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 1, id);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 2, name, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 3, flag);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 4, chr, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int64 (stmt, 5, pos);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 6, mapq);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 7, cigar, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 8, qlen);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 9, rlen);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_text (stmt, 10, chr_next, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int64 (stmt, 11, pos_next);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 12, type);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_step (stmt);
	if (rc != SQLITE_DONE)
		log_fatal ("Failed to insert alignment data: %s",
				sqlite3_errmsg (db));
}

sqlite3_stmt *
db_prepare_overlapping_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	int rc = 0;
	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"INSERT INTO overlapping (gene_id,alignment_id,exons) VALUES (?1,?2,?3)";

	rc = sqlite3_prepare_v2 (db, sql, -1, &stmt, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to prepare overlapping stmt: %s",
				sqlite3_errmsg (db));

	return stmt;
}

void
db_insert_overlapping (sqlite3 *db, sqlite3_stmt *stmt, int gene_id,
		int alignment_id, const char *exons)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && stmt != NULL);

	int rc = 0;

	rc = sqlite3_reset (stmt);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_clear_bindings (stmt);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 1, gene_id);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 2, alignment_id);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	if (exons != NULL)
		rc = sqlite3_bind_text (stmt, 3, exons, -1, SQLITE_TRANSIENT);
	else
		rc = sqlite3_bind_null (stmt, 3);

	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_step (stmt);
	if (rc != SQLITE_DONE)
		log_fatal ("Failed to insert overlapping data: %s",
				sqlite3_errmsg (db));
}
