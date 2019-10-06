#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "str.h"
#include "abnormal.h"
#include "dbscan.h"
#include "cluster.h"

struct _Cluster
{
	sqlite3_stmt *stmt;
	const char   *gene_name;
	int           id;
};

typedef struct _Cluster Cluster;

static void
index_alignment_qname (sqlite3 *db)
{
	// Index qname for speedup queries
	const char sql[] =
		"DROP INDEX IF EXISTS alignment_qname_idx;\n"
		"CREATE INDEX alignment_qname_idx\n"
		"	ON alignment(qname)";

	log_debug ("Create index:\n%s", sql);
	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_query_stmt (sqlite3 *db, const Set *blacklist_chr, const int distance)
{
	sqlite3_stmt *stmt = NULL;
	String *sql = NULL;
	List *list = NULL;
	ListElmt *cur = NULL;
	int blacklist_size = 0;
	int limit = 0;
	int i = 0;

	blacklist_size = set_size (blacklist_chr);
	limit = sqlite3_limit (db, SQLITE_LIMIT_VARIABLE_NUMBER, -1);

	// Test for set size. Needs to be less than SQLITE_LIMIT_VARIABLE_NUMBER - 2
	if (blacklist_size >= (limit - 2))
		log_fatal ("blacklisted chrs size (%d) is greater than the maximum allowed (%d)",
				blacklist_size, limit);

	/*
	* Where the magic happens!
	* Catch all alignments whose mate overlaps
	* a given exon and filter them by:
	* - blacklisted chromosomes (e.g. chrM)
	* - read cannot be exonic from its own parental
	* - distance from its own parental gene 
	*
	* Sort by chromosome and gene, so the clustering
	* will procide for each gene by its abnormal reads
	* along the chromosomes
	*/
	sql = string_new (
		"WITH\n"
		"	gene_pos (gene_name, chr, start, end) AS (\n"
		"		SELECT gene_name, chr, MIN(start), MAX(end)\n"
		"		FROM exon\n"
		"		GROUP BY gene_name\n"
		"	),\n"
		"	alignment_overlaps_exon(id, qname, chr, pos, rlen, type, gene_name) AS (\n"
		"		SELECT a.id, a.qname, a.chr, a.pos, a.rlen, a.type, e.gene_name\n"
		"		FROM alignment AS a\n"
		"		LEFT JOIN overlapping AS o\n"
		"			ON a.id = o.alignment_id\n"
		"		LEFT JOIN exon AS e\n"
		"			ON e.id = o.exon_id\n"
		"	)\n"
		"SELECT DISTINCT aoe1.id,\n"
		"	aoe1.chr AS achr,\n"
		"	aoe1.pos AS astart,\n"
		"	(aoe1.pos + aoe1.rlen - 1) AS aend,\n"
		"	aoe2.gene_name\n"
		"FROM alignment_overlaps_exon AS aoe1\n"
		"INNER JOIN alignment_overlaps_exon AS aoe2\n"
		"	USING(qname)\n"
		"INNER JOIN gene_pos AS g\n"
		"	ON aoe2.gene_name = g.gene_name\n"
		"WHERE "
	);

	if (blacklist_size)
		{
			string_concat (sql, "aoe1.chr NOT IN (?1");
			for (i = 1; i < blacklist_size; i++)
				string_concat_printf (sql, ",?%d", i + 1);

			string_concat (sql, ")\n\tAND aoe2.chr NOT IN (?1");
			for (i = 1; i < blacklist_size; i++)
				string_concat_printf (sql, ",?%d", i + 1);

			string_concat (sql, ")\n\tAND ");
		}

	string_concat (sql,
			"((achr != g.chr)\n"
			"		OR NOT (astart <= (g.end + $DIST) AND aend >= (g.start - $DIST)))\n"
			"	AND aoe1.id != aoe2.id\n"
			"	AND aoe2.type & $EXONIC\n"
			"	AND ((NOT aoe1.type & $EXONIC)\n"
			"		OR (aoe1.type & $EXONIC AND aoe1.gene_name IS NOT aoe2.gene_name))\n"
			"ORDER BY aoe1.chr ASC, aoe2.gene_name ASC");

	log_debug ("Query schema:\n%s", sql->str);
	stmt = db_prepare (db, sql->str);

	list = set_list (blacklist_chr);
	for (cur = list_head (list), i = 1; cur != NULL; cur = list_next (cur), i++)
		db_bind_text (stmt, i, list_data (cur));

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$DIST"),
			distance);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$EXONIC"),
			ABNORMAL_EXONIC);

	string_free (sql, 1);
	return stmt;
}

static void
dump_clustering (Point *p, void *user_data)
{
	Cluster *c = user_data;
	const int *alignment_id = p->data;

	log_debug ("Dump cluster%d aid = %d label = %d n = %d from = %s",
			p->id + c->id, *alignment_id, p->label,
			p->neighbors, c->gene_name);

	db_insert_clustering (c->stmt, p->id + c->id, *alignment_id,
			p->label, p->neighbors, c->gene_name);
}

static void
clustering (sqlite3_stmt *clustering_stmt, sqlite3_stmt *query_stmt,
		const long eps, const int min_pts)
{
	log_trace ("Inside %s", __func__);

	DBSCAN *dbscan = NULL;

	char *chr_prev = NULL;
	const char *chr = NULL;

	char *gene_name_prev = NULL;
	const char *gene_name = NULL;

	int *aid_alloc = NULL;

	int aid = 0;
	long astart = 0;
	long aend = 0;
	int acm = 0;

	Cluster c = {clustering_stmt, 0};

	while (db_step (query_stmt) == SQLITE_ROW)
		{
			aid       = db_column_int   (query_stmt, 0);
			chr       = db_column_text  (query_stmt, 1);
			astart    = db_column_int64 (query_stmt, 2);
			aend      = db_column_int64 (query_stmt, 3);
			gene_name = db_column_text  (query_stmt, 4);

			// First loop - Init dbscan, chr_prev and gene_name_prev
			if (chr_prev == NULL && gene_name_prev == NULL)
				{
					dbscan = dbscan_new (xfree);
					chr_prev = xstrdup (chr);
					gene_name_prev = xstrdup (gene_name);
				}

			// If all alignments at chromosome and gene_name
			// were cach, then cluster them
			if (strcmp (chr_prev, chr) || strcmp (gene_name_prev, gene_name))
				{
					log_debug ("Clustering by chr at '%s' for '%s'",
							chr_prev, gene_name_prev);

					c.gene_name = gene_name_prev;

					acm = dbscan_cluster (dbscan, eps, min_pts,
							dump_clustering, &c);

					log_debug ("Found %d clusters at '%s' for '%s'",
							acm, chr_prev, gene_name_prev);

					c.id += acm;

					dbscan_free (dbscan);
					dbscan = dbscan_new (xfree);

					xfree (chr_prev);
					xfree (gene_name_prev);

					chr_prev = xstrdup (chr);
					gene_name_prev = xstrdup (gene_name);
				}

			aid_alloc = xcalloc (1, sizeof (int));
			*aid_alloc = aid;

			// Insert a new point
			dbscan_insert_point (dbscan, astart, aend, aid_alloc);
		}

	// Run the last clustering - or
	// the first, if there is only
	// one cluster
	// Test if there any entry
	if (dbscan != NULL)
		{
			log_debug ("Clustering at '%s' for '%s'",
					chr_prev, gene_name_prev);

			c.gene_name = gene_name_prev;

			acm = dbscan_cluster (dbscan, eps, min_pts,
					dump_clustering, &c);

			c.id += acm;

			log_debug ("Found %d clusters at %s for '%s'",
					acm, chr_prev, gene_name_prev);
		}

	log_info ("Found %d clusters", c.id);

	xfree (chr_prev);
	xfree (gene_name_prev);
	dbscan_free (dbscan);
}

static void
insert_cluster_ids (sqlite3 *db)
{
	const char sql[] =
		"INSERT INTO cluster SELECT DISTINCT cluster_id FROM clustering";

	log_debug ("Insert cluster schema:\n%s", sql);
	db_exec (db, sql);
}

void
cluster (sqlite3_stmt *clustering_stmt, const long eps, const int min_pts,
		const Set *blacklist_chr, const int distance)
{
	log_trace ("Inside %s", __func__);
	assert (clustering_stmt != NULL
			&& min_pts > 2
			&& distance >= 0);

	sqlite3 *db = NULL;
	sqlite3_stmt *query_stmt = NULL;

	db = sqlite3_db_handle (clustering_stmt);

	log_debug ("Prepare query stmt");
	query_stmt = prepare_query_stmt (db, blacklist_chr, distance);

	log_info ("Index abnormal alignment qnames");
	index_alignment_qname (db);

	log_info ("Clustering abnormal alignments");
	clustering (clustering_stmt, query_stmt, eps, min_pts);

	log_info ("Fill cluster ids from clustering");
	insert_cluster_ids (db);

	log_info ("Done!");

	db_finalize (query_stmt);
}
