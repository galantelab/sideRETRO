#ifndef DB_H
#define DB_H

#include <stdlib.h>
#include <stdint.h>
#include <sqlite3.h>

/* Database schema version */
#define DB_SCHEMA_VERSION 0.4

#define DB_DEFAULT_CACHE_SIZE 2000

/* Low-level functions */

sqlite3 *      db_open (const char *path, int flags);
void           db_close (sqlite3 *db);

void           db_exec (sqlite3 *db, const char *sql);
sqlite3_stmt * db_prepare (sqlite3 *db, const char *sql);

int            db_step (sqlite3_stmt *stmt);
void           db_finalize (sqlite3_stmt *stmt);

void           db_reset (sqlite3_stmt *stmt);
void           db_clear_bindings (sqlite3_stmt *stmt);

void           db_bind_int (sqlite3_stmt *stmt, int i, int value);
void           db_bind_int64 (sqlite3_stmt *stmt, int i, int64_t value);
void           db_bind_double (sqlite3_stmt *stmt, int i, double value);
void           db_bind_text (sqlite3_stmt *stmt, int i, const char *value);

int            db_column_int (sqlite3_stmt *stmt, int i);
int64_t        db_column_int64 (sqlite3_stmt *stmt, int i);
double         db_column_double (sqlite3_stmt *stmt, int i);
const char   * db_column_text (sqlite3_stmt *stmt, int i);

/* database interface  */

sqlite3 * db_create            (const char *path);
sqlite3 * db_connect           (const char *path);
void      db_cache_size        (sqlite3 *db, size_t size);
void      db_begin_transaction (sqlite3 *db);
void      db_end_transaction   (sqlite3 *db);

sqlite3_stmt * db_prepare_exon_stmt (sqlite3 *db);
void db_insert_exon (sqlite3_stmt *stmt, int id, const char *gene_name,
		const char *chr, long start, long end, const char *strand, const char *ensg,
		const char *ense);

sqlite3_stmt *db_prepare_batch_stmt (sqlite3 *db);
void db_insert_batch (sqlite3_stmt *stmt, int id,
		const char *timestamp);

sqlite3_stmt * db_prepare_source_stmt (sqlite3 *db);
void db_insert_source (sqlite3_stmt *stmt, int id,
		int batch_id, const char *path);

sqlite3_stmt * db_prepare_alignment_stmt (sqlite3 *db);
void db_insert_alignment (sqlite3_stmt *stmt, int id, const char *name,
		int flag, const char *chr, long pos, int mapq, const char *cigar, int qlen,
		int rlen, const char *chr_next, long pos_next, int type, int source_id);

sqlite3_stmt * db_prepare_overlapping_stmt (sqlite3 *db);
void db_insert_overlapping (sqlite3_stmt *stmt, int exon_id,
	int alignment_id, long pos, long len);

sqlite3_stmt * db_prepare_clustering_stmt (sqlite3 *db);
void db_insert_clustering (sqlite3_stmt *stmt, int id,
	int alignment_id, int label, int neighbors);

sqlite3_stmt * db_prepare_reclustering_stmt (sqlite3 *db);
void db_insert_reclustering (sqlite3_stmt *stmt, int id, int clustering_id,
	int clustering_alignment_id, int label, int neighbors,
	const char *parental_gene_name);

#endif /* db.h */
