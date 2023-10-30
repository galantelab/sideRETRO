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
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "str.h"
#include "hash.h"
#include "abnormal.h"
#include "dbscan.h"
#include "cluster.h"

struct _Clustering
{
	Hash         *cluster_h;
	sqlite3_stmt *stmt;
	ClusterFilter filter;
	int64_t       id;
	int64_t       sid;
	int64_t       sub;
};

typedef struct _Clustering Clustering;

static inline Hash *
cluster_filter_new (void)
{
	return hash_new_full (int_hash, int_equal,
			xfree, (DestroyNotify) hash_free);
}

static inline void
cluster_filter_free (Hash *h)
{
	hash_free (h);
}

static inline void
cluster_filter_set (Hash *h, const int64_t id, const int64_t sid,
		const ClusterFilter filter)
{
/*
 * Handle Hash of Hashes schema to keep
 * the relation cluster => subcluster => filter
 */
	Hash *s_h = NULL;
	int *filter_alloc = NULL;

	s_h = hash_lookup (h, &id);

	if (s_h == NULL)
		{
			s_h = hash_new_full (int_hash, int_equal,
					xfree, xfree);

			int64_t *id_alloc = xcalloc (1, sizeof (int64_t));
			*id_alloc = id;

			hash_insert (h, id_alloc, s_h);
		}

	filter_alloc = hash_lookup (s_h, &sid);

	if (filter_alloc == NULL)
		{
			int64_t *sid_alloc = xcalloc (1, sizeof (int64_t));
			*sid_alloc = sid;

			filter_alloc = xcalloc (1,
					sizeof (ClusterFilter));

			hash_insert (s_h, sid_alloc, filter_alloc);
		}

	*filter_alloc |= filter;
}

static inline ClusterFilter *
cluster_filter_get (Hash *h, const int64_t id, const int64_t sid)
{
/*
 * Get id => sid => filter
 *  or NULL
 */
	Hash *s_h = NULL;

	s_h = hash_lookup (h, &id);
	if (s_h == NULL)
		return NULL;

	return hash_lookup (s_h, &sid);
}

static void
clean_clustering_tables (sqlite3 *db)
{
	// Delete all values from
	// previous runs
	const char sql[] =
		"DELETE FROM clustering;\n"
		"DELETE FROM cluster;";

	log_debug ("Clean tables:\n%s", sql);
	db_exec (db, sql);
}

static void
index_alignment_qname (sqlite3 *db)
{
	// Index qname for speedup queries
	const char sql[] =
		"DROP INDEX IF EXISTS alignment_qname_idx;\n"
		"CREATE INDEX alignment_qname_idx\n"
		"	ON alignment(qname,source_id)";

	log_debug ("Create index:\n%s", sql);
	db_exec (db, sql);
}

static void
index_overlapping_alignment (sqlite3 *db)
{
	const char sql[] =
		"DROP INDEX IF EXISTS overlapping_alignment_idx;\n"
		"CREATE INDEX overlapping_alignment_idx\n"
		"	ON overlapping(alignment_id)";

	log_debug ("Create index:\n%s", sql);
	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	/*
	* Where the magic happens!
	* Catch all alignments whose mate overlaps
	* a given exon and sort by chromosome and gene,
	* so the clustering will procide for each gene
	* by its abnormal reads along the chromosomes
	*/
	const char sql[] =
		"WITH\n"
		"	alignment_overlaps_exon(id, qname, source_id, chr, pos, rlen, type, gene_name) AS (\n"
		"		SELECT a.id, a.qname, a.source_id, a.chr, a.pos, a.rlen, a.type, e.gene_name\n"
		"		FROM alignment AS a\n"
		"		LEFT JOIN overlapping AS o\n"
		"			ON a.id = o.alignment_id\n"
		"		LEFT JOIN exon AS e\n"
		"			ON e.id = o.exon_id\n"
		"		WHERE type != $NONE\n"
		"	)\n"
		"SELECT DISTINCT aoe1.id,\n"
		"	aoe1.chr,\n"
		"	aoe1.pos,\n"
		"	CASE\n"
		"		WHEN aoe1.rlen <= 0\n"
		"			THEN (aoe1.pos)\n"
		"		ELSE\n"
		"			(aoe1.pos + aoe1.rlen - 1)\n"
		"	END,\n"
		"	aoe2.gene_name\n"
		"FROM alignment_overlaps_exon AS aoe1\n"
		"INNER JOIN alignment_overlaps_exon AS aoe2\n"
		"	USING (qname, source_id)\n"
		"WHERE aoe1.id != aoe2.id\n"
		"	AND aoe2.type & $EXONIC\n"
		"	AND ((NOT aoe1.type & $EXONIC)\n"
		"		OR (aoe1.type & $EXONIC AND aoe1.gene_name IS NOT aoe2.gene_name))\n"
		"ORDER BY aoe1.chr ASC, aoe2.gene_name ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$NONE"),
			ABNORMAL_NONE);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$EXONIC"),
			ABNORMAL_EXONIC);

	return stmt;
}

static void
dump_clustering (Point *p, void *user_data)
{
	Clustering *c = user_data;
	const int64_t *alignment_id = p->data;

	int64_t id = c->id;
	int64_t sid = c->sid;

	// Clustering or reclustering
	// ('sub'clustering)
	if (c->sub)
		sid += p->id;
	else
		id += p->id;

	log_debug ("Dump cluster [%" PRId64 " %" PRId64 "] aid = %" PRId64 " label = %d n = %d",
			id, sid, *alignment_id, p->label, p->neighbors);

	// Set id => sid => filter hash
	cluster_filter_set (c->cluster_h, id, sid, c->filter);

	db_insert_clustering (c->stmt, id, sid, *alignment_id,
			p->label, p->neighbors);
}

static int
clustering (sqlite3_stmt *clustering_stmt, const long eps,
		const int min_pts, Hash *cluster_h)
{
	log_trace ("Inside %s", __func__);

	DBSCAN *dbscan = NULL;
	sqlite3_stmt *query_stmt = NULL;

	char *chr_prev = NULL;
	const char *chr = NULL;

	char *gene_name_prev = NULL;
	const char *gene_name = NULL;

	int64_t *aid_alloc = NULL;

	int64_t aid = 0;
	long astart = 0;
	long aend = 0;
	int acm = 0;

	// Prepare query stmt
	query_stmt = prepare_query_stmt (
			sqlite3_db_handle (clustering_stmt));

	// STEP 1
	Clustering c = {
		.filter    = CLUSTER_FILTER_NONE,
		.cluster_h = cluster_h,
		.stmt      = clustering_stmt,
		.sub       = 0,
		.id        = 0,
		.sid       = 1
	};

	while (db_step (query_stmt) == SQLITE_ROW)
		{
			aid       = db_column_int64 (query_stmt, 0);
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
							c.id += acm;
							log_debug ("Found %d clusters at '%s' for '%s'",
									acm, chr_prev, gene_name_prev);
						}

					dbscan_free (dbscan);
					dbscan = dbscan_new (xfree);

					xfree (chr_prev);
					xfree (gene_name_prev);

					chr_prev = xstrdup (chr);
					gene_name_prev = xstrdup (gene_name);
				}

			aid_alloc = xcalloc (1, sizeof (int64_t));
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
					c.id += acm;
					log_debug ("Found %d clusters at '%s' for '%s'",
							acm, chr_prev, gene_name_prev);
				}
		}

	xfree (chr_prev);
	xfree (gene_name_prev);
	dbscan_free (dbscan);
	db_finalize (query_stmt);

	return c.id;
}

static sqlite3_stmt *
prepare_filter_support_stmt (sqlite3 *db,
		const int support)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Filter clustering reads according
	// to the genotype support number
	const char sql[] =
		"WITH\n"
		"	filter AS (\n"
		"		SELECT cluster_id, cluster_sid, source_id\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON a.id = c.alignment_id\n"
		"		GROUP BY cluster_id, cluster_sid, source_id\n"
		"		HAVING COUNT(*) >= $SUPPORT\n"
		")\n"
		"SELECT c.cluster_id,\n"
		"	c.cluster_sid,\n"
		"	alignment_id,\n"
		"	pos,\n"
		"	CASE\n"
		"		WHEN rlen <= 0\n"
		"			THEN (pos)\n"
		"		ELSE\n"
		"			(pos + rlen - 1)\n"
		"	END\n"
		"FROM clustering AS c\n"
		"INNER JOIN alignment AS a\n"
		"	ON c.alignment_id = a.id\n"
		"INNER JOIN filter AS f\n"
		"	USING (cluster_id, cluster_sid, source_id)";

	log_debug ("Clustering query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$SUPPORT"),
			support);

	return stmt;
}

static int
reclustering (sqlite3_stmt *clustering_stmt, const int support,
		const long eps, const int min_pts, Hash *cluster_h)
{
	log_trace ("Inside %s", __func__);

	DBSCAN *dbscan = NULL;
	sqlite3_stmt *filter_support_stmt = NULL;

	ClusterFilter *filter_alloc = NULL;
	int64_t *aid_alloc = NULL;

	int64_t cid = 0;
	int64_t cid_prev = 0;

	int64_t sid = 0;
	int64_t sid_prev = 0;

	int64_t aid = 0;
	long astart = 0;
	long aend = 0;
	int acm = 0;

	// Prepare filtering support query stmt
	filter_support_stmt = prepare_filter_support_stmt (
			sqlite3_db_handle (clustering_stmt), support);

	// STEP 2
	Clustering c = {
		.cluster_h = cluster_h,
		.stmt      = clustering_stmt,
		.sub       = 1,
		.id        = 0,
		.sid       = 0
	};

	while (db_step (filter_support_stmt) == SQLITE_ROW)
		{
			cid    = db_column_int64 (filter_support_stmt, 0);
			sid    = db_column_int64 (filter_support_stmt, 1);
			aid    = db_column_int64 (filter_support_stmt, 2);
			astart = db_column_int64 (filter_support_stmt, 3);
			aend   = db_column_int64 (filter_support_stmt, 4);

			// First loop - Init dbscan and cid_prev
			if (!cid_prev)
				{
					dbscan = dbscan_new (xfree);
					cid_prev = cid;
					sid_prev = sid;
				}

			// If all alignment from a given clusters were catch, then
			// recluster them!
			if (cid_prev != cid || sid_prev != sid)
				{
					// Get filter from clustering step
					filter_alloc = cluster_filter_get (cluster_h, cid_prev, sid_prev);
					assert (filter_alloc != NULL);

					log_debug ("Reclustering cluster [%" PRId64 " %" PRId64 "]", cid_prev, sid_prev);

					c.id = cid_prev;
					c.sid = sid_prev;
					c.filter =  *filter_alloc | CLUSTER_FILTER_SUPPORT;

					acm = dbscan_cluster (dbscan, eps, min_pts,
							dump_clustering, &c);

					if (acm)
						log_debug ("Found %d clusters from [%" PRId64 " %" PRId64 "] after reclustering",
								acm, cid_prev, sid_prev);

					dbscan_free (dbscan);
					dbscan = dbscan_new (xfree);

					cid_prev = cid;
					sid_prev = sid;
				}

			aid_alloc = xcalloc (1, sizeof (int64_t));
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
			// Get filter from clustering step
			filter_alloc = cluster_filter_get (cluster_h, cid_prev, sid_prev);
			assert (filter_alloc != NULL);

			log_debug ("Reclustering cluster [%" PRId64 " %" PRId64 "]", cid_prev, sid_prev);

			c.id = cid_prev;
			c.sid = sid_prev;
			c.filter =  *filter_alloc | CLUSTER_FILTER_SUPPORT;

			acm = dbscan_cluster (dbscan, eps, min_pts,
					dump_clustering, &c);

			if (acm)
				log_debug ("Found %d clusters from [%" PRId64 " %" PRId64 "] after reclustering",
						acm, cid_prev, sid_prev);
		}

	dbscan_free (dbscan);
	db_finalize (filter_support_stmt);

	return c.id;
}

static sqlite3_stmt *
prepare_filter_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"WITH\n"
		"	gene (gene_name, chr, start, end) AS (\n"
		"		SELECT gene_name, chr, MIN(start), MAX(end)\n"
		"		FROM exon\n"
		"		GROUP BY gene_name\n"
		"	),\n"
		"	cluster (gene_name, id, sid, chr, start, end) AS (\n"
		"		SELECT e.gene_name, cluster_id, cluster_sid,\n"
		"			a1.chr, MIN(a1.pos), MAX(a1.pos + a1.rlen - 1)\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a1\n"
		"			ON c.alignment_id = a1.id\n"
		"		INNER JOIN alignment AS a2\n"
		"			USING (qname, source_id)\n"
		"		INNER JOIN overlapping AS o\n"
		"			ON a2.id = o.alignment_id\n"
		"		INNER JOIN exon AS e\n"
		"			ON o.exon_id = e.id\n"
		"		GROUP BY cluster_id, cluster_sid\n"
		"	)\n"
		"SELECT c.id, c.sid, c.chr, c.start, c.end,\n"
		"	g.gene_name, g.chr, g.start, g.end\n"
		"FROM cluster AS c\n"
		"INNER JOIN gene AS g\n"
		"	USING (gene_name)";

	log_debug ("Filter query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static int
dump_and_filter_clusters (sqlite3_stmt *cluster_stmt, const int distance,
		const int support, Set *blacklist_chr, Blacklist *blacklist,
		const int padding, Hash *cluster_h)
{
	sqlite3_stmt *filter_query = NULL;
	ClusterFilter *filter = NULL;

	// Cluster info
	int64_t cluster_id = 0;
	int64_t cluster_sid = 0;
	const char *cluster_chr = NULL;
	long cluster_start = 0;
	long cluster_end = 0;

	// Parental info
	const char *gene_name = NULL;
	const char *gene_chr = NULL;
	long gene_start = 0;
	long gene_end = 0;

	int num_clusters = 0;
	int acm = 0;

	// Add the flag CLUSTER_FILTER_SUPPORT
	// if there was not need to reclustering
	const int support_flag =
		support > 1 ? 0 : CLUSTER_FILTER_SUPPORT;

	const int all_filters =
		CLUSTER_FILTER_NONE
		|CLUSTER_FILTER_CHR
		|CLUSTER_FILTER_DIST
		|CLUSTER_FILTER_REGION
		|CLUSTER_FILTER_SUPPORT;

	filter_query = prepare_filter_stmt (
			sqlite3_db_handle (cluster_stmt));

	while (db_step (filter_query) == SQLITE_ROW)
		{
			cluster_id    = db_column_int64 (filter_query, 0);
			cluster_sid   = db_column_int64 (filter_query, 1);
			cluster_chr   = db_column_text  (filter_query, 2);
			cluster_start = db_column_int64 (filter_query, 3);
			cluster_end   = db_column_int64 (filter_query, 4);
			gene_name     = db_column_text  (filter_query, 5);
			gene_chr      = db_column_text  (filter_query, 6);
			gene_start    = db_column_int64 (filter_query, 7);
			gene_end      = db_column_int64 (filter_query, 8);

			filter = cluster_filter_get (cluster_h, cluster_id, cluster_sid);
			assert (filter != NULL);

			// CLUSTER_FILTER_SUPPORT
			*filter |= support_flag;

			// Chromosome filter
			if (!set_is_member (blacklist_chr, cluster_chr)
					&& !set_is_member (blacklist_chr, gene_chr))
				*filter |= CLUSTER_FILTER_CHR;

			// Distance filter
			if (strcmp (cluster_chr, gene_chr)
					|| !(cluster_start <= (gene_end + distance)
						&& cluster_end >= (gene_start - distance)))
				*filter |= CLUSTER_FILTER_DIST;

			// Region filter
			acm = blacklist_lookup (blacklist, cluster_chr, cluster_start,
					cluster_end, padding, cluster_id, cluster_sid);

			if (!acm)
				*filter |= CLUSTER_FILTER_REGION;

			log_debug ("Dump cluster [%" PRId64 " %" PRId64 "] at %s:%li-%li from %s filter %d",
					cluster_id, cluster_sid, cluster_chr, cluster_start,
					cluster_end, gene_name, *filter);

			db_insert_cluster (cluster_stmt, cluster_id, cluster_sid,
					cluster_chr, cluster_start, cluster_end,
					gene_name, *filter);

			if (*filter == all_filters)
				num_clusters++;
		}

	db_finalize (filter_query);
	return num_clusters;
}

int
cluster (sqlite3_stmt *cluster_stmt, sqlite3_stmt *clustering_stmt,
		const long eps, const int min_pts, const int distance,
		const int support, Set *blacklist_chr, Blacklist *blacklist,
		const int padding)
{
	log_trace ("Inside %s", __func__);
	assert (cluster_stmt != NULL
			&& clustering_stmt != NULL
			&& min_pts > 2
			&& distance >= 0
			&& support >= 0
			&& padding >= 0
			&& blacklist != NULL);

	Hash *cluster_h = NULL;
	int num_clusters = 0;

	// cluster_id => sid => filter relation
	cluster_h = cluster_filter_new ();

	log_debug ("Clean clustering tables");
	clean_clustering_tables (
			sqlite3_db_handle (cluster_stmt));

	log_info ("Index abnormal alignment qnames");
	index_alignment_qname (
			sqlite3_db_handle (cluster_stmt));

	log_info ("Index overlapping alignment ids");
	index_overlapping_alignment (
			sqlite3_db_handle (cluster_stmt));

	// First clustering step
	log_info ("Clustering abnormal alignments");
	num_clusters = clustering (clustering_stmt,
			eps, min_pts, cluster_h);

	if (num_clusters)
		log_info ("Found %d clusters", num_clusters);
	else
		goto RET;

	// Next clustering
	if (support > 1)
		{
			log_info ("Filter clusters according to genotype support and recluster them");
			num_clusters = reclustering (clustering_stmt, support,
					eps, min_pts, cluster_h);

			if (num_clusters)
				log_info ("%d clusters left after genotype filtering and reclustering",
						num_clusters);
			else
				goto RET;
		}

	// Finally dump all clusters
	log_info (
			"Build clusters from clustering and filter them "
			"by blacklisted regions, and chromosome, and parental distance");

	num_clusters = dump_and_filter_clusters (cluster_stmt,
			distance, support, blacklist_chr,
			blacklist, padding, cluster_h);

	if (num_clusters)
		log_info ("%d clusters have been passed all controling filters",
				num_clusters);
	else
		goto RET;

RET:
	cluster_filter_free (cluster_h);
	return num_clusters;
}
