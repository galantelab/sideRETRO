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

#ifndef DB_H
#define DB_H

#include <stdlib.h>
#include <stdint.h>
#include <sqlite3.h>

/* Database schema version */
#define DB_SCHEMA_MAJOR_VERSION 0
#define DB_SCHEMA_MINOR_VERSION 12

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
void db_insert_exon (sqlite3_stmt *stmt, int64_t id, const char *gene_name,
		const char *chr, int64_t start, int64_t end, const char *strand, const char *ensg,
		const char *ense);

sqlite3_stmt *db_prepare_batch_stmt (sqlite3 *db);
void db_insert_batch (sqlite3_stmt *stmt, int64_t id,
		const char *timestamp);

sqlite3_stmt * db_prepare_source_stmt (sqlite3 *db);
void db_insert_source (sqlite3_stmt *stmt, int64_t id,
		int64_t batch_id, const char *path);

sqlite3_stmt * db_prepare_alignment_stmt (sqlite3 *db);
void db_insert_alignment (sqlite3_stmt *stmt, int64_t id, const char *name,
		int flag, const char *chr, int64_t pos, int mapq, const char *cigar, int qlen,
		int rlen, const char *chr_next, int64_t pos_next, int type, int64_t source_id);

sqlite3_stmt * db_prepare_overlapping_stmt (sqlite3 *db);
void db_insert_overlapping (sqlite3_stmt *stmt, int64_t exon_id,
	int64_t alignment_id, int64_t pos, int64_t len);

sqlite3_stmt * db_prepare_clustering_stmt (sqlite3 *db);
void db_insert_clustering (sqlite3_stmt *stmt, int64_t cluster_id, int64_t cluster_sid,
		int64_t alignment_id, int label, int neighbors);

sqlite3_stmt * db_prepare_cluster_stmt (sqlite3 *db);
void db_insert_cluster (sqlite3_stmt *stmt, int64_t id, int64_t sid, const char *chr,
		int64_t start, int64_t end, const char *gene_name, int filter);

sqlite3_stmt * db_prepare_blacklist_stmt (sqlite3 *db);
void db_insert_blacklist (sqlite3_stmt *stmt, int64_t id, const char *name,
		const char *chr, int64_t start, int64_t end);

sqlite3_stmt * db_prepare_overlapping_blacklist_stmt (sqlite3 *db);
void db_insert_overlapping_blacklist (sqlite3_stmt *stmt, int64_t blacklist_id,
	int64_t cluster_id, int64_t cluster_sid, int64_t pos, int64_t len);

sqlite3_stmt * db_prepare_cluster_merging_stmt (sqlite3 *db);
void db_insert_cluster_merging (sqlite3_stmt *stmt, int64_t retrocopy_id,
		int64_t cluster_id, int64_t cluster_sid);

sqlite3_stmt * db_prepare_retrocopy_stmt (sqlite3 *db);
void db_insert_retrocopy (sqlite3_stmt *stmt, int64_t id, const char *chr, int64_t window_start,
		int64_t window_end, const char *parental_gene_name, int level, int64_t insertion_point,
		int insertion_point_type, double orientation_rho, double orientation_p_value);

sqlite3_stmt * db_prepare_genotype_stmt (sqlite3 *db);
void db_insert_genotype (sqlite3_stmt *stmt, int64_t source_id, int64_t retrocopy_id, int reference_depth,
		int alternate_depth, double ho_ref_likelihood, double he_likelihood, double ho_alt_likelihood);

#endif /* db.h */
