#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "str.h"
#include "hash.h"
#include "abnormal.h"
#include "dbscan.h"
#include "cluster.h"

struct _Cluster
{
	sqlite3_stmt *stmt;
	int           id;
};

typedef struct _Cluster Cluster;

static void
index_alignment_qname (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

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
	log_trace ("Inside %s", __func__);

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
create_temp_clustering_table (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	// Create temporary table as clustering.
	// It is necessary for validating each cluster
	// according to genotype support in reads
	const char sql[] =
		"CREATE TEMPORARY TABLE temp_clustering AS SELECT * FROM clustering";

	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_temp_clustering_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"INSERT INTO temp_clustering (cluster_id,alignment_id,label,neighbors)\n"
		"	VALUES (?1,?2,?3,?4)";

	return db_prepare (db, sql);
}

static sqlite3_stmt *
prepare_temp_clustering_query_stmt (sqlite3 *db,
		const int support)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Filter clustering reads according
	// to the genotype support number
	const char sql[] =
		"WITH\n"
		"	filter (cid, sid) AS (\n"
		"		SELECT cluster_id, source_id\n"
		"		FROM temp_clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON a.id = c.alignment_id\n"
		"		GROUP BY cluster_id, source_id\n"
		"		HAVING COUNT(*) >= ?1\n"
		")\n"
		"SELECT cluster_id, alignment_id, pos, (pos + rlen - 1)\n"
		"FROM temp_clustering AS c\n"
		"INNER JOIN alignment AS a\n"
		"	ON c.alignment_id = a.id\n"
		"INNER JOIN filter AS f\n"
		"	ON f.cid = c.cluster_id\n"
		"WHERE a.source_id = f.sid";

	log_debug ("Temporary clustering query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt, 1, support);
	return stmt;
}

static void
dump_clustering (Point *p, void *user_data)
{
	Cluster *c = user_data;
	const int *alignment_id = p->data;

	log_debug ("Dump cluster%d aid = %d label = %d n = %d",
			p->id + c->id, *alignment_id, p->label,
			p->neighbors);

	db_insert_clustering (c->stmt, p->id + c->id, *alignment_id,
			p->label, p->neighbors);
}

static void
clustering (sqlite3_stmt *clustering_stmt, sqlite3_stmt *query_stmt,
		const long eps, const int min_pts, Hash *cluster_gene)
{
	log_trace ("Inside %s", __func__);

	DBSCAN *dbscan = NULL;

	char *chr_prev = NULL;
	const char *chr = NULL;

	char *gene_name_prev = NULL;
	const char *gene_name = NULL;

	int *aid_alloc = NULL;
	int *cid_alloc = NULL;

	int aid = 0;
	long astart = 0;
	long aend = 0;
	int acm = 0;
	int i = 0;

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

					acm = dbscan_cluster (dbscan, eps, min_pts,
							dump_clustering, &c);

					if (acm)
						{
							log_debug ("Found %d clusters at '%s' for '%s'",
									acm, chr_prev, gene_name_prev);

							// Store cluster_id => gene_name_prev relation
							for (i = 0; i < acm; i++)
								{
									cid_alloc = xcalloc (1, sizeof (int));
									*cid_alloc = ++c.id;
									hash_insert (cluster_gene, cid_alloc,
											xstrdup (gene_name_prev));
								}
						}

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

			acm = dbscan_cluster (dbscan, eps, min_pts,
					dump_clustering, &c);

			if (acm)
				{
					log_debug ("Found %d clusters at '%s' for '%s'",
							acm, chr_prev, gene_name_prev);

					// Store cluster_id => gene_name_prev relation
					for (i = 0; i < acm; i++)
						{
							cid_alloc = xcalloc (1, sizeof (int));
							*cid_alloc = ++c.id;
							hash_insert (cluster_gene, cid_alloc,
									xstrdup (gene_name_prev));
						}
				}
		}

	log_info ("Found %d clusters", c.id);

	xfree (chr_prev);
	xfree (gene_name_prev);
	dbscan_free (dbscan);
}

static void
reclustering (sqlite3_stmt *clustering_stmt, sqlite3_stmt *query_stmt,
		const long eps, const int min_pts, Hash *cluster_gene)
{
	log_trace ("Inside %s", __func__);

	DBSCAN *dbscan = NULL;

	const char *gene_name = NULL;

	int *aid_alloc = NULL;
	int *cid_alloc = NULL;

	int cid = 0;
	int cid_prev = 0;

	int aid = 0;
	long astart = 0;
	long aend = 0;
	int acm = 0;
	int i = 0;

	Cluster c = {clustering_stmt, 0};

	while (db_step (query_stmt) == SQLITE_ROW)
		{
			cid       = db_column_int   (query_stmt, 0);
			aid       = db_column_int   (query_stmt, 1);
			astart    = db_column_int64 (query_stmt, 2);
			aend      = db_column_int64 (query_stmt, 3);

			// First loop - Init dbscan and cid_prev
			if (!cid_prev)
				{
					dbscan = dbscan_new (xfree);
					cid_prev = cid;
				}

			// If all alignment from a given clusters were catch, then
			// recluster them!
			if (cid_prev != cid)
				{
					// Get gene_name from clustering step
					gene_name = hash_lookup (cluster_gene, &cid_prev);

					log_debug ("Reclustering by cluster_id = '%d' from '%s'",
							cid_prev, gene_name);

					acm = dbscan_cluster (dbscan, eps, min_pts,
							dump_clustering, &c);

					if (acm)
						{
							log_debug ("Found %d clusters from cluster_id = '%d' after reclustering",
									acm, cid_prev);

							// Store cluster_id => gene_name relation
							for (i = 0; i < acm; i++)
								{
									cid_alloc = xcalloc (1, sizeof (int));
									*cid_alloc = ++c.id;
									hash_insert (cluster_gene, cid_alloc,
											xstrdup (gene_name));
								}
						}

					dbscan_free (dbscan);
					dbscan = dbscan_new (xfree);

					cid_prev = cid;
				}

			aid_alloc = xcalloc (1, sizeof (int));
			*aid_alloc = aid;

			// Insert a new point
			dbscan_insert_point (dbscan, astart, aend, aid_alloc);
		}

	// Run the last reclustering - or
	// the first, if there is only
	// one cluster
	// Test if there any entry
	if (dbscan != NULL)
		{
			// Get gene_name from clustering step
			gene_name = hash_lookup (cluster_gene, &cid_prev);

			log_debug ("Reclustering by cluster_id = '%d' from '%s'",
					cid_prev, gene_name);

			acm = dbscan_cluster (dbscan, eps, min_pts,
					dump_clustering, &c);

			if (acm)
				{
					log_debug ("Found %d clusters from cluster_id = '%d' after reclustering",
							acm, cid_prev);

					// Store cluster_id => gene_name_prev relation
					for (i = 0; i < acm; i++)
						{
							cid_alloc = xcalloc (1, sizeof (int));
							*cid_alloc = ++c.id;
							hash_insert (cluster_gene, cid_alloc,
									xstrdup (gene_name));
						}
				}
		}

	log_info ("Found %d clusters after reclustering", c.id);
	dbscan_free (dbscan);
}

static sqlite3_stmt *
prepare_clustering_query_stmt (sqlite3 *db)
{
	/*
	* All the clusters are build at end, when the
	* clustering step is done
	*/
	const char sql[] =
		"SELECT cluster_id, chr, MIN(pos), MAX(pos + rlen - 1)\n"
		"FROM clustering AS c\n"
		"INNER JOIN alignment AS a\n"
		"	ON a.id = c.alignment_id\n"
		"GROUP BY cluster_id";

	log_debug ("Clustering query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static void
run_clustering_step (sqlite3 *db, const long eps, const int min_pts,
		const int distance, const Set *blacklist_chr,
		Hash *cluster_gene)
{
	log_trace ("Inside %s", __func__);

	// Query the database for valid abnormal reads
	sqlite3_stmt *query_stmt = NULL;

	// Insert into temporary clustering table
	sqlite3_stmt *clustering_stmt = NULL;

	// Prepare query stmt
	query_stmt = prepare_query_stmt (db, blacklist_chr, distance);

	// Prepare temp clustering
	clustering_stmt = prepare_temp_clustering_stmt (db);

	clustering (clustering_stmt, query_stmt, eps, min_pts,
			cluster_gene);

	db_finalize (query_stmt);
	db_finalize (clustering_stmt);
}

static void
run_reclustering_step (sqlite3 *db, sqlite3_stmt *clustering_stmt,
		const long eps, const int min_pts, const int support,
		Hash *cluster_gene)
{
	log_trace ("Inside %s", __func__);

	// Query temporary clustering table for validation
	sqlite3_stmt *query_stmt = NULL;

	// Prepare temp clustering query stmt
	query_stmt = prepare_temp_clustering_query_stmt (db, support);

	reclustering (clustering_stmt, query_stmt,
			eps, min_pts, cluster_gene);

	db_finalize (query_stmt);
}

static void
dump_cluster (sqlite3 *db, sqlite3_stmt *cluster_stmt,
		Hash *cluster_gene, Blacklist *blacklist)
{
	sqlite3_stmt *query_stmt = NULL;

	int id = 0;
	const char *chr = NULL;
	long start = 0;
	long end = 0;
	const char *gene_name = NULL;
	int acm = 0;

	// Prepare clustering query stmt
	query_stmt = prepare_clustering_query_stmt (db);

	while (sqlite3_step (query_stmt) == SQLITE_ROW)
		{
			id = db_column_int (query_stmt, 0);
			chr = db_column_text (query_stmt, 1);
			start = db_column_int64 (query_stmt, 2);
			end = db_column_int64 (query_stmt, 3);
			gene_name = hash_lookup (cluster_gene, &id);

			acm = blacklist_lookup (blacklist, chr, start, end,
					0, id);

			if (acm)
				log_debug ("Cluster [%d] from %s %s:%li-%li overlaps %d blacklisted regions",
						id, gene_name, chr, start, end, acm);

			log_debug ("Dump cluster %d at %s:%li-%li from %s",
					id, chr, start, end, gene_name);

			db_insert_cluster (cluster_stmt, id, chr,
					start, end, gene_name);
		}

	db_finalize (query_stmt);
}

void
cluster (sqlite3_stmt *cluster_stmt, sqlite3_stmt *clustering_stmt,
		const long eps, const int min_pts, const Set *blacklist_chr,
		const int distance, const int support,
		Blacklist *blacklist)
{
	log_trace ("Inside %s", __func__);
	assert (cluster_stmt != NULL
			&& clustering_stmt != NULL
			&& min_pts > 2
			&& distance >= 0
			&& support >= 0
			&& blacklist != NULL);

	sqlite3 *db = NULL;

	// Hold cluster_id => gene relation
	Hash *cluster_gene = NULL;

	// sqlite3 handle
	db = sqlite3_db_handle (clustering_stmt);

	log_debug ("Create temporary clustering table");
	create_temp_clustering_table (db);

	log_info ("Index abnormal alignment qnames");
	index_alignment_qname (db);

	// cluster => gene relation
	cluster_gene = hash_new_full (int_hash, int_equal,
			xfree, xfree);

	log_info ("Clustering abnormal alignments");
	run_clustering_step (db, eps, min_pts, distance,
			blacklist_chr, cluster_gene);

	log_info ("Filter clusters according to genotype support and recluster them");
	run_reclustering_step (db, clustering_stmt, eps, min_pts,
			support, cluster_gene);

	log_info ("Build clusters from clustering");
	dump_cluster (db, cluster_stmt, cluster_gene, blacklist);

	log_info ("Done!");
	hash_free (cluster_gene);
}
