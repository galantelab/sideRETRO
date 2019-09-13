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
	sqlite3_stmt *clustering_stmt;
	int           id;
};

typedef struct _Cluster Cluster;

static void
index_alignment_qname (sqlite3 *db)
{
	const char sql[] =
		"DROP INDEX IF EXISTS alignment_qname_idx;\n"
		"CREATE INDEX alignment_qname_idx\n"
		"	ON alignment(qname)";

	log_debug ("Create index:\n%s", sql);
	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_query_stmt (sqlite3 *db)
{
	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"WITH\n"
		"	alignment_overlaps_exon(id, qname, chr, pos, rlen, type, gene_name) AS (\n"
		"		SELECT a.id, a.qname, a.chr, a.pos, a.rlen, a.type, e.gene_name\n"
		"		FROM alignment AS a\n"
		"		LEFT JOIN overlapping AS o\n"
		"			ON a.id = o.alignment_id\n"
		"		LEFT JOIN exon AS e\n"
		"			ON e.id = o.exon_id\n"
		"	)\n"
		"SELECT DISTINCT aoe1.id, aoe1.chr, aoe1.pos, aoe1.pos + aoe1.rlen - 1\n"
		"FROM alignment_overlaps_exon AS aoe1\n"
		"INNER JOIN alignment_overlaps_exon AS aoe2\n"
		"	USING(qname)\n"
		"WHERE aoe1.id != aoe2.id\n"
		"	AND aoe2.type & ?1\n"
		"	AND ((NOT aoe1.type & ?1)\n"
		"		OR (aoe1.type & ?1 AND aoe1.gene_name IS NOT aoe2.gene_name))\n"
		"ORDER BY aoe1.chr ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);
	db_bind_int (stmt, 1, ABNORMAL_EXONIC);

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

	db_insert_clustering (c->clustering_stmt,
			p->id + c->id, alignment_id, p->label,
			p->neighbors);
}

void
cluster (sqlite3_stmt *clustering_stmt,
		long eps, int min_pts, Set *blacklist_chr)
{
	log_trace ("Inside %s", __func__);
	assert (clustering_stmt != NULL
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
		.clustering_stmt = clustering_stmt,
		.id = 0
	};

	log_info ("Index abnormal alignment qnames");
	index_alignment_qname (sqlite3_db_handle (clustering_stmt));

	log_debug ("Prepare query stmt");
	query_stmt = prepare_query_stmt (
			sqlite3_db_handle (clustering_stmt));

	log_info ("Clustering abnormal alignments");

	while (db_step (query_stmt) == SQLITE_ROW)
		{
			id = db_column_int (query_stmt, 0);
			chr = db_column_text (query_stmt, 1);
			low = db_column_int64 (query_stmt, 2);
			high = db_column_int64 (query_stmt, 3);

			if (set_is_member (blacklist_chr, chr))
				{
					log_debug ("Avoid blacklistied chr: '%s'", chr);
					continue;
				}

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
	db_finalize (query_stmt);
	dbscan_free (dbscan);
}
