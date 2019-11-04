#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "db.h"
#include "str.h"
#include "dedup.h"

#define BLOCK_SIZE 8
#define STR_SIZE   32

struct _DedupData
{
	int       id;
	String   *qname;

	String   *chr;
	long      pos;

	String   *chr_next;
	long      pos_next;

	int       source_id;
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
	data->id        = db_column_int   (stmt, 0);
	data->pos       = db_column_int64 (stmt, 3);
	data->pos_next  = db_column_int64 (stmt, 5);
	data->source_id = db_column_int   (stmt, 6);

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
		"	id ASC";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static void
create_temp_dup_table (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"CREATE TEMPORARY TABLE dup (\n"
		"	id INTEGER NOT NULL,\n"
		"	qname TEXT NOT NULL,\n"
		"	source_id INTEGER NOT NULL)";

	log_debug ("Create table:\n%s", sql);
	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_temp_dup_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"INSERT INTO dup (id,qname,source_id)\n"
		"VALUES (?1,?2,?3)";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static inline void
insert_temp_dup (sqlite3_stmt *stmt, const int id,
		const char *qname, const int source_id)
{
	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt, 1, id);
	db_bind_text (stmt, 2, qname);
	db_bind_int (stmt, 3, source_id);

	db_step (stmt);
}

static void
delete_primary_reads_from_dup (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"DELETE FROM dup\n"
		"WHERE (qname, source_id) NOT IN (\n"
		"	SELECT qname, source_id\n"
		"	FROM dup AS d1\n"
		"	INNER JOIN dup AS d2\n"
		"		USING (qname, source_id)\n"
		"	WHERE d1.id < d2.id\n"
		"	GROUP BY d1.id)";

	log_debug ("Delete from table:\n%s", sql);
	db_exec (db, sql);
}

static void
delete_dup (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"DELETE FROM alignment\n"
		"WHERE (qname, source_id) IN (\n"
		"	SELECT qname, source_id\n"
		"	FROM dup)";

	log_debug ("Delete from table:\n%s", sql);
	db_exec (db, sql);
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
	int *source_ids = NULL;
	String **qnames = NULL;

	int id = 0;
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

			// First loop
			if (!data_prev.id)
				dedup_data_copy (&data_prev, &data);

			if (!dedup_data_is_dup (&data, &data_prev))
				{
					if (size > 1)
						{
							id++;
							for (i = 0; i < size; i++)
								insert_temp_dup (temp_dup_stmt,
										id, qnames[i]->str, source_ids[i]);
						}

					dedup_data_copy (&data_prev, &data);

					size = 0;
				}

			if (size >= alloc)
				{
					old_alloc = alloc;
					alloc += BLOCK_SIZE;

					source_ids = xrealloc (source_ids,
							sizeof (int) * alloc);

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

	// Delete primary reads from 'dup' table
	log_info ("Choose primary reads");
	delete_primary_reads_from_dup (db);

	// Delete duplicated reads from alignment
	log_info ("Remove duplicated reads");
	delete_dup (db);
}
