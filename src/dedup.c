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

#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "db.h"
#include "str.h"
#include "abnormal.h"
#include "dedup.h"

#define BLOCK_SIZE 8
#define STR_SIZE   32

struct _DedupData
{
	int64_t   id;
	String   *qname;

	String   *chr;
	long      pos;

	String   *chr_next;
	long      pos_next;

	int64_t   source_id;
};

typedef struct _DedupData DedupData;

static inline void
dedup_data_init (DedupData *data)
{
	data->qname = string_sized_new (STR_SIZE);
	data->chr = string_sized_new (STR_SIZE);
	data->chr_next = string_sized_new (STR_SIZE);
}

static inline void
dedup_data_read (DedupData *data, sqlite3_stmt *stmt)
{
	data->id        = db_column_int64 (stmt, 0);
	data->pos       = db_column_int64 (stmt, 3);
	data->pos_next  = db_column_int64 (stmt, 5);
	data->source_id = db_column_int64 (stmt, 6);

	string_set (data->qname, db_column_text (stmt, 1));
	string_set (data->chr, db_column_text (stmt, 2));
	string_set (data->chr_next, db_column_text (stmt, 4));
}

static inline void
dedup_data_copy (DedupData *to, const DedupData *from)
{
	to->id        = from->id;
	to->pos       = from->pos;
	to->pos_next  = from->pos_next;
	to->source_id = from->source_id;

	string_set (to->qname, from->qname->str);
	string_set (to->chr, from->chr->str);
	string_set (to->chr_next, from->chr_next->str);
}

static inline int
dedup_data_is_dup (const DedupData *data1, const DedupData *data2)
{
	return data1->id != data2->id
		&& !strcmp (data1->chr->str, data2->chr->str)
		&& data1->pos == data2->pos
		&& !strcmp (data1->chr_next->str, data2->chr_next->str)
		&& data1->pos_next == data2->pos_next
		&& data1->source_id == data2->source_id;
}

static inline void
dedup_data_destroy (DedupData *data)
{
	string_free (data->qname, 1);
	string_free (data->chr, 1);
	string_free (data->chr_next, 1);
}

static sqlite3_stmt *
prepare_alignment_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"SELECT id, qname, chr, pos, chr_next, pos_next, source_id\n"
		"FROM alignment\n"
		"ORDER BY source_id ASC,\n"
		"	chr ASC, pos ASC,\n"
		"	chr_next ASC, pos_next ASC,\n"
		"	qname ASC";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static void
create_temp_dup_table (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"CREATE TEMPORARY TABLE dup (\n"
		"	qname TEXT NOT NULL,\n"
		"	source_id INTEGER NOT NULL,\n"
		"	is_primary INTEGER NOT NULL,\n"
		"	PRIMARY KEY (qname, source_id))";

	log_debug ("Create table:\n%s", sql);
	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_temp_dup_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"INSERT OR IGNORE INTO dup (qname,source_id,is_primary)\n"
		"VALUES (?1,?2,?3)";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static inline void
insert_temp_dup (sqlite3_stmt *stmt, const char *qname,
		const int64_t source_id, const int is_primary)
{
	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_text  (stmt, 1, qname);
	db_bind_int64 (stmt, 2, source_id);
	db_bind_int   (stmt, 3, is_primary);

	db_step (stmt);
}

static void
set_dup (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Invalidate duplicated reads
	// without removing it
	const char sql[] =
		"UPDATE alignment\n"
		"SET type = $NONE\n"
		"WHERE (qname, source_id) IN (\n"
		"	SELECT qname, source_id\n"
		"	FROM dup\n"
		"	WHERE is_primary = 0)";

	log_debug ("Set duplicated reads to ABNORMAL_NONE flag:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$NONE"),
			ABNORMAL_NONE);

	// RUN
	db_step (stmt);

	// Clean
	db_finalize (stmt);
}

static void
mark_dup (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *query_stmt = NULL;
	sqlite3_stmt *temp_dup_stmt = NULL;

	// Current entry
	DedupData data = {};

	// Previous entry
	DedupData data_prev = {};

	// Manage the variadic arrays
	int size = 0;
	int alloc = 0;
	int old_alloc = 0;

	// Array of qnames and source_ids
	int64_t *source_ids = NULL;
	String **qnames = NULL;

	int i = 0;

	// Init DedupData
	dedup_data_init (&data);
	dedup_data_init (&data_prev);

	// Prepare query and dup statements
	query_stmt = prepare_alignment_query_stmt (db);
	temp_dup_stmt = prepare_temp_dup_stmt (db);

	while (db_step (query_stmt) == SQLITE_ROW)
		{
			dedup_data_read (&data, query_stmt);

			if (!data_prev.id)
				{
					// First loop
					dedup_data_copy (&data_prev, &data);
				}
			else if (!dedup_data_is_dup (&data, &data_prev))
				{
					if (size > 1)
						{
							// The primary read - Just the first one
							insert_temp_dup (temp_dup_stmt, qnames[0]->str,
									source_ids[0], 1);

							for (i = 1; i < size; i++)
								insert_temp_dup (temp_dup_stmt, qnames[i]->str,
										source_ids[i], 0);
						}

					size = 0;
				}
			else if (!strcmp (data.qname->str, data_prev.qname->str))
				{
					// Ignore weird mates mapping the same genomic coordinate
					continue;
				}

			if (size >= alloc)
				{
					old_alloc = alloc;
					alloc += BLOCK_SIZE;

					source_ids = xrealloc (source_ids,
							sizeof (int64_t) * alloc);

					qnames = xrealloc (qnames,
							sizeof (String *) * alloc);

					for (; old_alloc < alloc; old_alloc++)
						qnames[old_alloc] =
							string_sized_new (STR_SIZE);
				}

			source_ids[size] = data.source_id;

			string_set (qnames[size],
					data.qname->str);

			size++;

			// Update data_prev
			dedup_data_copy (&data_prev, &data);
		}

	for (i = 0; i < alloc; i++)
		string_free (qnames[i], 1);

	xfree (qnames);
	xfree (source_ids);

	dedup_data_destroy (&data);
	dedup_data_destroy (&data_prev);

	db_finalize (query_stmt);
	db_finalize (temp_dup_stmt);
}

void
dedup (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	// Create temp table
	log_debug ("Create temporary table 'dup'");
	create_temp_dup_table (db);

	// Insert duplicated reads to 'dup' table
	log_info ("Mark duplicated reads");
	mark_dup (db);

	// Delete duplicated reads from alignment
	log_info ("Remove duplicated reads");
	set_dup (db);
}
