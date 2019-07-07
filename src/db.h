#ifndef DB_H
#define DB_H

#include "sqlite3.h"

sqlite3 * db_create  (const char *path);
sqlite3 * db_connect (const char *path);
void      db_close   (sqlite3 *db);

void db_cache_size        (sqlite3 *db, size_t size);
void db_begin_transaction (sqlite3 *db);
void db_end_transaction   (sqlite3 *db);

void db_finalize (sqlite3 *db, sqlite3_stmt *stmt);

sqlite3_stmt * db_prepare_exon_stmt (sqlite3 *db);
void db_insert_exon (sqlite3 *db, sqlite3_stmt *stmt, int id, const char *gene_name,
		const char *chr, long start, long end, const char *strand, const char *ensg,
		const char *ense);

sqlite3_stmt * db_prepare_source_stmt (sqlite3 *db);
void db_insert_source (sqlite3 *db, sqlite3_stmt *stmt, int id,
		const char *name, const char *path);

sqlite3_stmt * db_prepare_alignment_stmt (sqlite3 *db);
void db_insert_alignment (sqlite3 *db, sqlite3_stmt *stmt, int id, const char *name,
		int flag, const char *chr, long pos, int mapq, const char *cigar, int qlen,
		int rlen, const char *chr_next, long pos_next, int type, int source_id);

sqlite3_stmt * db_prepare_overlapping_stmt (sqlite3 *db);
void db_insert_overlapping (sqlite3 *db, sqlite3_stmt *stmt,
		int exon_id, int alignment_id);

#endif /* db.h */
