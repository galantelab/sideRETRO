#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "abnormal.h"
#include "dbscan.h"
#include "cluster.h"

struct _Cluster
{
	sqlite3      *db;
	sqlite3_stmt *clustering_stmt;
	int           id;
};

typedef struct _Cluster Cluster;

static sqlite3_stmt *
prepare_query_stmt (sqlite3 *db)
{
	int rc = 0;
	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"SELECT id, chr, pos, pos + rlen - 1\n"
		"FROM alignment\n"
		"WHERE qname IN (\n"
		"	SELECT qname\n"
		"	FROM alignment\n"
		"	WHERE type & ?1\n"
		")\n"
		"AND NOT (type & ?2)\n"
		"ORDER BY chr ASC, pos ASC;";

	log_debug ("Query schema:\n%s", sql);
	rc = sqlite3_prepare_v2 (db, sql, -1, &stmt, NULL);

	if (rc != SQLITE_OK)
		log_fatal ("Failed to prepare query stmt: %s",
				sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 1, ABNORMAL_EXONIC);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	rc = sqlite3_bind_int (stmt, 2, ABNORMAL_EXONIC);
	if (rc != SQLITE_OK) log_fatal ("%s", sqlite3_errmsg (db));

	return stmt;
}

static void
dump_clustering (Point *p, void *user_data)
{
	Cluster *c = user_data;
	int alignment_id = * (int *) p->data;

	log_debug ("Dump cluster%d aid = %d label = %d n = %d",
			p->id + c->id, alignment_id, p->label,
			p->neighbors);

	db_insert_clustering (c->db, c->clustering_stmt,
			p->id + c->id, alignment_id, p->label,
			p->neighbors);
}

void
cluster (sqlite3 *db, sqlite3_stmt *clustering_stmt,
		long eps, int min_pts)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && clustering_stmt != NULL
			&& min_pts > 2);

	DBSCAN *dbscan = NULL;
	sqlite3_stmt *query_stmt = NULL;
	char *chr_prev = NULL;
	const char *chr = NULL;
	int id = 0;
	int *id_alloc = NULL;
	long low = 0;
	long high = 0;
	int acm = 0;

	Cluster c = {
		.db = db,
		.clustering_stmt = clustering_stmt,
		.id = 0
	};

	log_debug ("Prepare query stmt");
	query_stmt = prepare_query_stmt (db);

	log_info ("Clustering abnormal alignments");

	while (sqlite3_step (query_stmt) == SQLITE_ROW)
		{
			id = sqlite3_column_int (query_stmt, 0);
			chr = sqlite3_column_text (query_stmt, 1);
			low = sqlite3_column_int64 (query_stmt, 2);
			high = sqlite3_column_int64 (query_stmt, 3);

			// First loop - Init dbscan and chr_prev
			if (chr_prev == NULL)
				{
					dbscan = dbscan_new (xfree);
					chr_prev = xstrdup (chr);
				}

			// If all alignments at chromosome
			// were catch, then cluster them
			if (strcmp (chr_prev, chr))
				{
					log_debug ("Clustering at '%s'", chr_prev);
					acm = dbscan_cluster (dbscan, eps, min_pts,
							dump_clustering, &c);

					log_debug ("Found %d clusters at %s",
							acm, chr_prev);

					c.id += acm;

					dbscan_free (dbscan);
					dbscan = dbscan_new (xfree);

					xfree (chr_prev);
					chr_prev = xstrdup (chr);
				}

			id_alloc = xcalloc (1, sizeof (int));
			*id_alloc = id;

			// Insert a new point
			dbscan_insert_point (dbscan, low, high, id_alloc);
		}

	// Run the last clustering - or
	// the first, if there is only
	// one cluster
	// Test if there any entry
	if (dbscan != NULL)
		{
			log_debug ("Clustering at '%s'", chr_prev);
			acm = dbscan_cluster (dbscan, eps, min_pts,
					dump_clustering, &c);

			c.id += acm;

			log_debug ("Found %d clusters at %s",
					acm, chr_prev);
		}

	log_info ("Found %d clusters", c.id);

	xfree (chr_prev);
	db_finalize (db, query_stmt);
	dbscan_free (dbscan);
}
