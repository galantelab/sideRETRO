#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "array.h"
#include "list.h"
#include "hash.h"
#include "log.h"
#include "db.h"
#include "abnormal.h"
#include "cluster.h"
#include "correlation.h"
#include "retrocopy.h"

#define MAX_DIST 3
#define BLOCK_SIZE 64
#define SEED 17

struct _ClusterEntry
{
	// CLUSTER
	int         cid;
	int         sid;
	const char *cchr;
	long        cstart;
	long        cend;

	// PARENTAL
	const char *gene_name;
	const char *gchr;
	long        gstart;
	long        gend;
	int         dist;
};

typedef struct _ClusterEntry ClusterEntry;

struct _RetrocopyEntry
{
	RetrocopyLevel level;
	double         orientation_rho;
	double         orientation_p_value;
};

typedef struct _RetrocopyEntry RetrocopyEntry;

static ClusterEntry *
cluster_entry_new (int cid, int sid, const char *cchr, long cstart, long cend,
		const char *gene_name, const char *gchr, long gstart, long gend, int dist)
{
	ClusterEntry *c = xcalloc (1, sizeof (ClusterEntry));

	*c = (ClusterEntry) {
		.cid       = cid,
		.sid       = sid,
		.cchr      = xstrdup (cchr),
		.cstart    = cstart,
		.cend      = cend,
		.gene_name = xstrdup (gene_name),
		.gchr      = xstrdup (gchr),
		.gstart    = gstart,
		.gend      = gend,
		.dist      = dist
	};

	return c;
}

static void
cluster_entry_free (ClusterEntry *c)
{
	if (c == NULL)
		return;

	xfree ((void *) c->cchr);
	xfree ((void *) c->gene_name);
	xfree ((void *) c->gchr);

	xfree (c);
}

static void
clean_retrocopy_tables (sqlite3 *db)
{
	// Delete all values from
	// previous runs
	const char sql[] =
		"DELETE FROM cluster_merging;\n"
		"DELETE FROM retrocopy;";

	log_debug ("Clean tables:\n%s", sql);
	db_exec (db, sql);
}

static sqlite3_stmt *
prepare_cluster_query_stmt (sqlite3 *db, const int filter)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"WITH\n"
		"	gene (gene_name, chr, start, end) AS (\n"
		"		SELECT gene_name, chr, MIN(start), MAX(end)\n"
		"		FROM exon\n"
		"		GROUP BY gene_name\n"
		"	),\n"
		"	gene_rank (gene_name, chr, start, end, dist) AS (\n"
		"		SELECT *,\n"
		"			DENSE_RANK() OVER (\n"
		"				PARTITION BY chr\n"
		"				ORDER BY start ASC, end ASC\n"
		"			)\n"
		"		FROM gene\n"
		"	)\n"
		"SELECT c.id, c.sid, c.chr, c.start, c.end,\n"
		"	c.gene_name, g.chr, g.start, g.end, g.dist\n"
		"FROM cluster AS c\n"
		"INNER JOIN gene_rank AS g\n"
		"	USING (gene_name)\n"
		"WHERE c.filter = $FILTER\n"
		"ORDER BY c.chr ASC, c.start ASC, c.end ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$FILTER"),
			filter);

	return stmt;
}

static inline int
overlaps (const char *chr1, const long start1, const long end1,
		const char *chr2, const long start2, const long end2)
{
	if (!strcmp (chr1, chr2) && (start1 <= end2 && end1 >= start2))
		return 1;

	return 0;
}

static inline int
is_near (const char *chr1, const int dist1,
		const char *chr2, const int dist2)
{
	if (!strcmp (chr1, chr2) && abs (dist1 - dist2) <= MAX_DIST)
		return 1;

	return 0;
}

static int
cmp_cluster (const void *p1, const void *p2)
{
	ClusterEntry *c1 = * (ClusterEntry **) p1;
	ClusterEntry *c2 = * (ClusterEntry **) p2;

	int rc = strcmp (c1->gchr, c2->gchr);

	if (!rc)
		{
			if (c1->gstart < c2->gend)
				return -1;
			else if (c1->gstart > c2->gend)
				return 1;
			else
				return 0;
		}

	return rc;
}

static void
dump_cluster_merge (sqlite3_stmt *cluster_merging_stmt,
		const List *to_merge, const int rid)
{
	log_trace ("Inside %s", __func__);

	ClusterEntry *c = NULL;
	ListElmt *cur = NULL;

	for (cur = list_head (to_merge); cur != NULL; cur = list_next (cur))
		{
			c = list_data (cur);

			log_debug ("Dump cluster merging [%d %d] into retrocopy %d",
					c->cid, c->sid, rid);

			db_insert_cluster_merging (cluster_merging_stmt,
					rid, c->cid, c->sid);
		}
}

static inline RetrocopyEntry *
rtc_insert_entry (Hash *h, const int rid)
{
	int *rid_alloc = xcalloc (1, sizeof (int));
	*rid_alloc = rid;

	RetrocopyEntry *e = xcalloc (1, sizeof (RetrocopyEntry));
	hash_insert (h, rid_alloc, e);

	return e;
}

static void
cluster_entry_merge_and_classify (sqlite3_stmt *cluster_merging_stmt,
		Array *stack, Hash *rtc_h, int *rid)
{
	log_trace ("Inside %s", __func__);

	List *to_merge = NULL;

	ClusterEntry *c = NULL;
	ClusterEntry *c_prev = NULL;

	RetrocopyEntry *e = NULL;
	RetrocopyLevel level = 0;
	long gend_prev = 0;

	int i = 0;

	// Sort by parental gene position
	array_sort (stack, (CompareFunc) cmp_cluster);

	c_prev = array_get (stack, 0);
	gend_prev = c_prev->gend;

	to_merge = list_new (NULL);
	list_append (to_merge, c_prev);

	for (i = 1; i < array_len (stack); i++)
		{
			c = array_get (stack, i);

			if (overlaps (c_prev->gchr, c_prev->gstart, gend_prev,
						c->gchr, c->gstart, c->gend))
				{
					list_append (to_merge, c);
					level |= RETROCOPY_OVERLAPPED_PARENTALS;
				}
			else if (is_near (c_prev->gchr, c_prev->dist,
						c->gchr, c->dist))
				{
					list_append (to_merge, c);
					level |= RETROCOPY_NEAR_PARENTALS;
				}
			else
				{
					// Update rid
					(*rid)++;

					dump_cluster_merge (cluster_merging_stmt,
							to_merge, *rid);

					e = rtc_insert_entry (rtc_h, *rid);
					e->level = level | RETROCOPY_HOTSPOT;

					list_free (to_merge);
					to_merge = list_new (NULL);

					list_append (to_merge, c);
					gend_prev = c->gend;

					level = RETROCOPY_HOTSPOT;
				}

			if (c->gend > gend_prev)
				gend_prev = c->gend;

			c_prev = c;
		}

	if (!level)
		level = RETROCOPY_PASS;

	// Update rid
	(*rid)++;

	dump_cluster_merge (cluster_merging_stmt,
			to_merge, *rid);

	e = rtc_insert_entry (rtc_h, *rid);
	e->level = level;

	list_free (to_merge);
}

static inline void
merge_cluster_init (Array **a, char **chr, long *start, long *end,
		const char *chr_v, const long start_v, const long end_v)
{
	array_free (*a, 1);
	*a = array_new ((DestroyNotify) cluster_entry_free);

	xfree (*chr);
	*chr = xstrdup (chr_v);

	*start = start_v;
	*end = end_v;
}

static void
merge_cluster (sqlite3_stmt *cluster_merging_stmt,
		const int filter, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *cluster_query_stmt = NULL;

	Array *stack = NULL;
	ClusterEntry *c = NULL;

	// Retrocopy ids acm
	int rid = 0;

	// CLUSTER
	int cid = 0;
	int sid = 0;
	const char *cchr = NULL;
	long cstart = 0;
	long cend = 0;

	// PARENTAL
	const char *gene_name = NULL;
	const char *gchr = NULL;
	long gstart = 0;
	long gend = 0;
	int dist = 0;

	// Previous
	char *cchr_prev = NULL;
	long cstart_prev = 0;
	long cend_prev = 0;

	// Prepare query
	cluster_query_stmt = prepare_cluster_query_stmt (
			sqlite3_db_handle (cluster_merging_stmt),
			filter);

	while (db_step (cluster_query_stmt) == SQLITE_ROW)
		{
			cid       = db_column_int   (cluster_query_stmt, 0);
			sid       = db_column_int   (cluster_query_stmt, 1);
			cchr      = db_column_text  (cluster_query_stmt, 2);
			cstart    = db_column_int64 (cluster_query_stmt, 3);
			cend      = db_column_int64 (cluster_query_stmt, 4);

			gene_name = db_column_text  (cluster_query_stmt, 5);
			gchr      = db_column_text  (cluster_query_stmt, 6);
			gstart    = db_column_int64 (cluster_query_stmt, 7);
			gend      = db_column_int64 (cluster_query_stmt, 8);
			dist      = db_column_int   (cluster_query_stmt, 9);

			// First loop
			if (stack == NULL)
				merge_cluster_init (&stack, &cchr_prev, &cstart_prev, &cend_prev,
						cchr, cstart, cend);

			// Process stacked clusters
			if (!overlaps (cchr_prev, cstart_prev, cend_prev,
						cchr, cstart, cend))
				{
					cluster_entry_merge_and_classify (cluster_merging_stmt,
							stack, rtc_h, &rid);

					merge_cluster_init (&stack, &cchr_prev, &cstart_prev, &cend_prev,
							cchr, cstart, cend);
				}

			c = cluster_entry_new (cid, sid, cchr, cstart, cend,
					gene_name, gchr, gstart, gend, dist);

			array_add (stack, c);

			/*
			* Increase cluster previous window
			*/
			if (cend > cend_prev)
				cend_prev = cend;
		}

	if (stack != NULL)
		cluster_entry_merge_and_classify (cluster_merging_stmt,
				stack, rtc_h, &rid);

	xfree (cchr_prev);
	array_free (stack, 1);
	db_finalize (cluster_query_stmt);
}

static sqlite3_stmt *
prepare_orientation_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Get all alignment position for read and its
	// mate - side by side
	const char sql[] =
		"WITH\n"
		"	alignment_flag (id, sid, aid, pos) AS (\n"
		"		SELECT c.cluster_id, c.cluster_sid, a.qname, pos\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON c.alignment_id = a.id\n"
		"	),\n"
		"	alignment_overlaps_exon (id, qname, pos) AS (\n"
		"		SELECT id, qname, pos\n"
		"		FROM alignment\n"
		"		WHERE type & $EXONIC\n"
		"	),\n"
		"	gene_flag (id, sid, aid, pos) AS (\n"
		"		SELECT c.cluster_id, c.cluster_sid, a.qname, aoe.pos\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON c.alignment_id = a.id\n"
		"		INNER JOIN alignment_overlaps_exon AS aoe\n"
		"			USING (qname)\n"
		"		WHERE a.id != aoe.id\n"
		"	)\n"
		"SELECT DISTINCT retrocopy_id, a.pos, g.pos\n"
		"FROM cluster_merging AS c\n"
		"INNER JOIN alignment_flag AS a\n"
		"	ON c.cluster_id = a.id\n"
		"		AND c.cluster_sid = a.sid\n"
		"INNER JOIN gene_flag AS g\n"
		"	USING (id, sid, aid)";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$EXONIC"),
			ABNORMAL_EXONIC);

	return stmt;
}

static void
calculate_orientation (sqlite3 *db, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *orientation_stmt = NULL;

	RetrocopyEntry *e = NULL;

	int rid = 0;
	long apos = 0;
	long gpos = 0;

	int rid_prev = 0;

	double *apos_a = NULL;
	double *gpos_a = NULL;
	size_t alloc_a = 0;
	size_t size_a = 0;

	double *work1 = NULL;
	double *work2 = NULL;
	size_t size_work = 0;

	double orientation_rho = 0;
	double orientation_p_value = 0;
	unsigned int seed = SEED;

	orientation_stmt = prepare_orientation_stmt (db);

	while (db_step (orientation_stmt) == SQLITE_ROW)
		{
			rid  = db_column_int   (orientation_stmt, 0);
			apos = db_column_int64 (orientation_stmt, 1);
			gpos = db_column_int64 (orientation_stmt, 2);

			// First loop
			if (!rid_prev)
				rid_prev = rid;

			// Calculate spearman test
			if (rid_prev != rid)
				{
					if ((2 * size_a) > size_work)
						{
							size_work = 2 * size_a;
							work1 = xrealloc (work1, sizeof (double) * size_work);
							work2 = xrealloc (work2, sizeof (double) * size_work);
						}

					// spearman test!
					orientation_rho = spearman (apos_a, gpos_a, size_a, work1);

					// Calculate p-value
					orientation_p_value = spearman_permutation_test (apos_a, gpos_a,
							size_a, work1, work2, &seed, orientation_rho);

					// Set orientation value for the given entry
					e = hash_lookup (rtc_h, &rid_prev);
					assert (e != NULL);

					e->orientation_rho = orientation_rho;
					e->orientation_p_value = orientation_p_value;

					size_a = 0;
					rid_prev = rid;
				}

			if (size_a >= alloc_a)
				{
					alloc_a += BLOCK_SIZE;
					apos_a = xrealloc (apos_a, sizeof (double) * alloc_a);
					gpos_a = xrealloc (gpos_a, sizeof (double) * alloc_a);
				}

			apos_a[size_a] = apos;
			gpos_a[size_a] = gpos;
			size_a++;
		}

	// The last or the first retrocopy
	if (rid_prev)
		{
			if ((2 * size_a) > size_work)
				{
					size_work = 2 * size_a;
					work1 = xrealloc (work1, sizeof (double) * size_work);
					work2 = xrealloc (work2, sizeof (double) * size_work);
				}

			// spearman test!
			orientation_rho = spearman (apos_a, gpos_a, size_a, work1);

			// Calculate p-value
			orientation_p_value = spearman_permutation_test (apos_a, gpos_a,
					size_a, work1, work2, &seed, orientation_rho);

			// Set orientation value for the given entry
			e = hash_lookup (rtc_h, &rid_prev);
			assert (e != NULL);

			e->orientation_rho = orientation_rho;
			e->orientation_p_value = orientation_p_value;
		}

	xfree (apos_a);
	xfree (gpos_a);
	xfree (work1);
	xfree (work2);

	db_finalize (orientation_stmt);
}

static sqlite3_stmt *
prepare_cluster_merging_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Merge clusters and calculate the insetion point
	// ALL IN ONE SHOT!
	const char sql[] =
		"WITH\n"
		"	cluster_mean_point (id, sid, mean) AS (\n"
		"		SELECT id, sid, (start + end) / 2\n"
		"		FROM cluster\n"
		"	),\n"
		"	cluster_cigar_mode (id, sid, pos) AS (\n"
		"		SELECT cluster_id, cluster_sid,\n"
		"			(SELECT CASE\n"
		"				WHEN cigar LIKE '%M%S' OR cigar LIKE '%M%H' THEN\n"
		"					pos + rlen\n"
		"				WHEN cigar LIKE '%S%M' OR cigar LIKE '%H%M' THEN\n"
		"					pos\n"
		"				ELSE\n"
		"					NULL\n"
		"				END AS p\n"
		"			FROM clustering AS c1\n"
		"			INNER JOIN alignment AS a\n"
		"				ON c1.alignment_id = a.id\n"
		"			WHERE c1.cluster_id = c2.cluster_id\n"
		"				AND c1.cluster_sid = c2.cluster_sid\n"
		"				AND flag & 0x800\n"
		"			GROUP BY p\n"
		"			ORDER BY COUNT(*) DESC\n"
		"			LIMIT 1)\n"
		"		FROM (SELECT DISTINCT cluster_id, cluster_sid FROM clustering) AS c2\n"
		"	),\n"
		"	cluster_ip (id, sid, ip, ip_type) AS (\n"
		"		SELECT a.id, a.sid,\n"
		"			CASE\n"
		"				WHEN pos IS NULL THEN\n"
		"					mean\n"
		"				ELSE\n"
		"					pos\n"
		"			END,\n"
		"			CASE\n"
		"				WHEN pos IS NULL THEN\n"
		"					$WINDOW_MEAN\n"
		"				ELSE\n"
		"					$SUPPLEMENTARY_MODE\n"
		"			END\n"
		"		FROM cluster_mean_point AS a\n"
		"		LEFT JOIN cluster_cigar_mode AS b\n"
		"			USING (id, sid)\n"
		"	),\n"
		"	cluster_merge (rid, id, sid, chr, start, end, gene) AS (\n"
		"		SELECT retrocopy_id, c.id, c.sid,\n"
		"			chr, MIN(start), MAX(end),\n"
		"			GROUP_CONCAT(gene_name,'/')\n"
		"		FROM cluster AS c\n"
		"		INNER JOIN cluster_merging AS m\n"
		"			ON c.id = m.cluster_id AND c.sid = m.cluster_sid\n"
		"		GROUP BY retrocopy_id\n"
		"	)\n"
		"SELECT rid, chr, start, end, gene, ip, ip_type\n"
		"FROM cluster_merge AS c\n"
		"INNER JOIN cluster_ip AS i\n"
		"	USING (id, sid)";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$WINDOW_MEAN"),
			RETROCOPY_INSERTION_POINT_WINDOW_MEAN);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$SUPPLEMENTARY_MODE"),
			RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE);

	return stmt;
}

static void
annotate_retrocopy (sqlite3_stmt *retrocopy_stmt, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *cluster_merging_query_stmt = NULL;

	int rid = 0;
	const char *chr = NULL;
	long start = 0;
	long end = 0;
	const char *gene = NULL;
	long ip = 0;

	RetrocopyEntry *e = NULL;
	RetrocopyInsertionPoint ip_type = 0;

	cluster_merging_query_stmt = prepare_cluster_merging_query_stmt (
			sqlite3_db_handle (retrocopy_stmt));

	while (db_step (cluster_merging_query_stmt) == SQLITE_ROW)
		{
			rid     = db_column_int   (cluster_merging_query_stmt, 0);
			chr     = db_column_text  (cluster_merging_query_stmt, 1);
			start   = db_column_int64 (cluster_merging_query_stmt, 2);
			end     = db_column_int64 (cluster_merging_query_stmt, 3);
			gene    = db_column_text  (cluster_merging_query_stmt, 4);
			ip      = db_column_int64 (cluster_merging_query_stmt, 5);
			ip_type = db_column_int   (cluster_merging_query_stmt, 6);

			e = hash_lookup (rtc_h, &rid);
			assert (e != NULL);

			log_debug ("%d %s %li %li %s %d %li %d %.6f %.6f",
					rid, chr, start, end, gene, e->level, ip, ip_type,
					e->orientation_rho, e->orientation_p_value);

			db_insert_retrocopy (retrocopy_stmt, rid, chr, start, end,
					gene, e->level, ip, ip_type, e->orientation_rho,
					e->orientation_p_value);
		}

	db_finalize (cluster_merging_query_stmt);
}

void
retrocopy (sqlite3_stmt *retrocopy_stmt,
		sqlite3_stmt *cluster_merging_stmt)
{
	log_trace ("Inside %s", __func__);
	assert (retrocopy_stmt != NULL
			&& cluster_merging_stmt != NULL);

	// Keep ID => level
	Hash *rtc_h = NULL;

	const int filter =
		CLUSTER_FILTER_NONE
		|CLUSTER_FILTER_CHR
		|CLUSTER_FILTER_DIST
		|CLUSTER_FILTER_REGION
		|CLUSTER_FILTER_SUPPORT;

	rtc_h = hash_new_full (int_hash, int_equal,
			xfree, xfree);

	log_debug ("Clean retrocopy tables");
	clean_retrocopy_tables (
			sqlite3_db_handle (retrocopy_stmt));

	log_info ("Analise and merge clusters into retrocopies");
	merge_cluster (cluster_merging_stmt, filter, rtc_h);

	log_info ("Calculate retrocopies orientation");
	calculate_orientation (
			sqlite3_db_handle (retrocopy_stmt),
			rtc_h);

	log_info ("Annotate retrocopies");
	annotate_retrocopy (retrocopy_stmt, rtc_h);

	hash_free (rtc_h);
}
