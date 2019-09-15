#include "config.h"

#include <assert.h>
#include "log.h"
#include "hash.h"
#include "dbscan.h"
#include "wrapper.h"
#include "recluster.h"

struct _Recluster
{
	sqlite3_stmt *reclustering_stmt;
	int           id;
};

typedef struct _Recluster Recluster;

static sqlite3_stmt *
prepare_query_stmt (sqlite3 *db, const int distance)
{
	sqlite3_stmt *stmt = NULL;

	const char sql[]=
		"SELECT cid, aid, astart, aend, gene_name\n"
		"FROM (\n"
		"	SELECT DISTINCT c.id AS cid,\n"
		"		c.alignment_id AS aid,\n"
		"		a.chr AS achr,\n"
		"		a.pos AS astart,\n"
		"		(a.pos + a.rlen - 1) AS aend,\n"
		"		gene_name,\n"
		"		e.chr AS gchr,\n"
		"		MIN(e.start) AS gstart,\n"
		"		MAX(e.end) AS gend\n"
		"	FROM clustering AS c\n"
		"	INNER JOIN alignment AS a\n"
		"		ON c.alignment_id = a.id\n"
		"	INNER JOIN alignment AS a2\n"
		"		USING (qname)\n"
		"	INNER JOIN overlapping AS o\n"
		"		ON o.alignment_id = a2.id\n"
		"	INNER JOIN exon AS e\n"
		"		ON o.exon_id = e.id\n"
		"	GROUP BY cid, aid, gene_name\n"
		"	ORDER BY cid ASC, gene_name ASC, astart ASC, aend ASC)\n"
		"WHERE achr != gchr\n"
		"	OR NOT (astart <= (gend + ?1)\n"
		"		AND aend >= (gstart - ?1))";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);
	db_bind_int (stmt, 1, distance);

	return stmt;
}

static void
dump_reclustering (Point *p, void *user_data)
{
	Recluster *r = user_data;
	int aid = * (int *) p->data;

	log_debug ("Redump cluster%d aid = %d label = %d n = %d",
			p->id + r->id, aid, p->label, p->neighbors);

	db_insert_reclustering (r->reclustering_stmt,
			p->id + r->id, aid, p->label,
			p->neighbors);
}

static void
cluster_analysis (Recluster *r, Hash *cluster_by_gene,
		const long eps, const int min_pts)
{
	const char *gene_name = NULL;
	DBSCAN *dbscan = NULL;
	int acm = 0;

	HashIter iter;
	hash_iter_init (&iter, cluster_by_gene);

	while (hash_iter_next (&iter, (void **) &gene_name, (void **) &dbscan))
		{
			log_debug ("Cluster for putative parental '%s'", gene_name);
			acm = dbscan_cluster (dbscan, eps, min_pts, dump_reclustering, r);

			log_debug ("Found %d clusters at %s", acm, gene_name);
			r->id += acm;
		}
}

void
recluster (sqlite3_stmt *reclustering_stmt,
		const long eps, const int min_pts,
		const int distance)
{
	log_trace ("Inside %s", __func__);
	assert (reclustering_stmt != NULL
			&& min_pts > 2);

	sqlite3_stmt *query_stmt = NULL;
	Hash *cluster_by_gene = NULL;
	DBSCAN *dbscan = NULL;

	int cid = 0;
	int cid_prev = 0;
	int aid = 0;
	int *aid_alloc = NULL;
	long astart = 0;
	long aend = 0;
	const char *gene_name = NULL;

	Recluster r = {reclustering_stmt, 0};

	log_debug ("Prepare query stmt");
	query_stmt = prepare_query_stmt (
			sqlite3_db_handle (reclustering_stmt),
			distance);

	log_info ("Reclustering clustered alignments");

	while (db_step (query_stmt) == SQLITE_ROW)
		{
			cid = db_column_int (query_stmt, 0);
			aid = db_column_int (query_stmt, 1);
			astart = db_column_int64 (query_stmt, 2);
			aend = db_column_int64 (query_stmt, 3);
			gene_name = db_column_text (query_stmt, 4);

			if (!cid_prev)
				{
					cluster_by_gene = hash_new (xfree, (DestroyNotify) dbscan_free);
					cid_prev = cid;
				}

			if (cid_prev != cid)
				{
					log_debug ("Reclustering cluster '%d'", cid_prev);
					cluster_analysis (&r, cluster_by_gene, eps, min_pts);

					hash_free (cluster_by_gene);
					cluster_by_gene = hash_new (xfree, (DestroyNotify) dbscan_free);

					cid_prev = cid;
				}

			aid_alloc = xcalloc (1, sizeof (int));
			*aid_alloc = aid;

			dbscan = hash_lookup (cluster_by_gene, gene_name);
			if (dbscan == NULL)
				{
					dbscan = dbscan_new (xfree);
					hash_insert (cluster_by_gene, xstrdup (gene_name), dbscan);
				}

			dbscan_insert_point (dbscan, astart, aend, aid_alloc);
		}

	// Run the last reclustering - or
	// the first, if there is only
	// one cluster
	// Test if there any entry
	if (cluster_by_gene != NULL)
		{
			log_debug ("Reclustering cluster '%d'", cid_prev);
			cluster_analysis (&r, cluster_by_gene, eps, min_pts);
		}

	log_info ("Found %d clusters after reclustering", r.id);

	hash_free (cluster_by_gene);
	db_finalize (query_stmt);
}
