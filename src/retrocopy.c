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
#include "cigar.h"
#include "retrocopy.h"

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

	// Breakpoints and the number of
	// split reads supporting them
	long           breakpoint_left;
	int            breakpoint_left_sr;
	long           breakpoint_right;
	int            breakpoint_right_sr;
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
		const char *chr2, const int dist2,
		const int near_gene_dist)
{
	if (!strcmp (chr1, chr2) && abs (dist1 - dist2) <= near_gene_dist)
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
		Array *stack, Hash *rtc_h, const int near_gene_dist, int *rid)
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
					if (strcmp (c_prev->gene_name, c->gene_name))
						level |= RETROCOPY_OVERLAPPED_PARENTALS;
				}
			else if (is_near (c_prev->gchr, c_prev->dist,
						c->gchr, c->dist, near_gene_dist))
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
		const int filter, const int near_gene_dist,
		Hash *rtc_h)
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
							stack, rtc_h, near_gene_dist, &rid);

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
				stack, rtc_h, near_gene_dist, &rid);

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
		"	alignment_flag (id, sid, aid, srcid, pos) AS (\n"
		"		SELECT c.cluster_id, c.cluster_sid, a.qname, a.source_id, a.pos\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON c.alignment_id = a.id\n"
		"	),\n"
		"	alignment_overlaps_exon (id, qname, source_id, pos) AS (\n"
		"		SELECT id, qname, source_id, pos\n"
		"		FROM alignment\n"
		"		WHERE type & $EXONIC\n"
		"	),\n"
		"	gene_flag (id, sid, aid, srcid, pos) AS (\n"
		"		SELECT c.cluster_id, c.cluster_sid, a.qname, a.source_id, aoe.pos\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON c.alignment_id = a.id\n"
		"		INNER JOIN alignment_overlaps_exon AS aoe\n"
		"			USING (qname, source_id)\n"
		"		WHERE a.id != aoe.id\n"
		"	)\n"
		"SELECT DISTINCT retrocopy_id, a.pos, g.pos\n"
		"FROM cluster_merging AS c\n"
		"INNER JOIN alignment_flag AS a\n"
		"	ON c.cluster_id = a.id\n"
		"		AND c.cluster_sid = a.sid\n"
		"INNER JOIN gene_flag AS g\n"
		"	USING (id, sid, aid, srcid)";

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
prepare_breakpoint_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	/*
	* Get all clipped alignments clustered into the
	* retrocopies - primary or supplementary. The reads
	* crossing the poly-A junction usually have no
	* supplementary alignment, so do not require it
	*/
	const char sql[] =
		"SELECT DISTINCT retrocopy_id, a.id, a.pos, a.rlen, a.cigar\n"
		"FROM cluster_merging AS cm\n"
		"INNER JOIN clustering AS c\n"
		"	USING (cluster_id, cluster_sid)\n"
		"INNER JOIN alignment AS a\n"
		"	ON a.id = c.alignment_id\n"
		"WHERE a.cigar LIKE '%S%' OR a.cigar LIKE '%H%'\n"
		"ORDER BY retrocopy_id ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	return stmt;
}

static int
cmp_long (const void *p1, const void *p2)
{
	const long l1 = * (const long *) p1;
	const long l2 = * (const long *) p2;
	return (l1 > l2) - (l1 < l2);
}

static inline void
breakpoint_push (long **a, size_t *size, size_t *alloc,
		const long pos)
{
	if (*size >= *alloc)
		{
			*alloc += BLOCK_SIZE;
			*a = xrealloc (*a, sizeof (long) * *alloc);
		}

	(*a)[(*size)++] = pos;
}

static void
breakpoint_mode (long *a, const size_t size, long *mode, int *support)
{
	long cur = 0;
	int acm = 0;
	size_t i = 0;

	*mode = 0;
	*support = 0;

	if (size == 0)
		return;

	// Sort in order to count repeated positions
	// and keep the first one in case of tie
	qsort (a, size, sizeof (long), cmp_long);

	cur = a[0];
	acm = 1;

	for (i = 1; i <= size; i++)
		{
			if (i < size && a[i] == cur)
				{
					acm++;
					continue;
				}

			if (acm > *support)
				{
					*support = acm;
					*mode = cur;
				}

			if (i < size)
				{
					cur = a[i];
					acm = 1;
				}
		}
}

static inline void
breakpoint_set_entry (Hash *rtc_h, const int rid, long *left_a,
		const size_t left_size, long *right_a, const size_t right_size)
{
	RetrocopyEntry *e = NULL;

	e = hash_lookup (rtc_h, &rid);
	assert (e != NULL);

	breakpoint_mode (left_a, left_size,
			&e->breakpoint_left, &e->breakpoint_left_sr);

	breakpoint_mode (right_a, right_size,
			&e->breakpoint_right, &e->breakpoint_right_sr);

	log_debug ("Retrocopy %d breakpoints: left %li (%d SR) right %li (%d SR)",
			rid, e->breakpoint_left, e->breakpoint_left_sr,
			e->breakpoint_right, e->breakpoint_right_sr);
}

static void
calculate_breakpoint (sqlite3 *db, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *breakpoint_stmt = NULL;

	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	uint32_t clip_left = 0;
	uint32_t clip_right = 0;

	int rid = 0;
	int aid = 0;
	long pos = 0;
	long rlen = 0;
	const char *cigar_str = NULL;

	int rid_prev = 0;

	// Candidates for the left and right breakpoints
	long *left_a = NULL;
	size_t left_size = 0;
	size_t left_alloc = 0;

	long *right_a = NULL;
	size_t right_size = 0;
	size_t right_alloc = 0;

	breakpoint_stmt = prepare_breakpoint_stmt (db);

	while (db_step (breakpoint_stmt) == SQLITE_ROW)
		{
			rid       = db_column_int   (breakpoint_stmt, 0);
			aid       = db_column_int   (breakpoint_stmt, 1);
			pos       = db_column_int64 (breakpoint_stmt, 2);
			rlen      = db_column_int64 (breakpoint_stmt, 3);
			cigar_str = db_column_text  (breakpoint_stmt, 4);

			// First loop
			if (!rid_prev)
				rid_prev = rid;

			// Set breakpoints for the previous retrocopy
			if (rid_prev != rid)
				{
					breakpoint_set_entry (rtc_h, rid_prev, left_a,
							left_size, right_a, right_size);

					left_size = 0;
					right_size = 0;
					rid_prev = rid;
				}

			cigar = cigar_yyparser (cigar_str, &n_cigar);
			if (cigar == NULL)
				{
					log_warn ("Skip alignment %d with invalid cigar '%s'",
							aid, cigar_str);
					continue;
				}

			cigar_clip (cigar, n_cigar, &clip_left, &clip_right);
			xfree (cigar);

			/*
			* A read clipped at its right side comes from
			* the left flank: The breakpoint is the last
			* aligned base. A read clipped at its left side
			* comes from the right flank: The breakpoint is
			* the first aligned base
			*/
			if (clip_right >= RETROCOPY_MIN_CLIP && rlen > 0)
				breakpoint_push (&left_a, &left_size, &left_alloc,
						pos + rlen - 1);

			if (clip_left >= RETROCOPY_MIN_CLIP)
				breakpoint_push (&right_a, &right_size, &right_alloc,
						pos);
		}

	// The last or the first retrocopy
	if (rid_prev)
		breakpoint_set_entry (rtc_h, rid_prev, left_a,
				left_size, right_a, right_size);

	xfree (left_a);
	xfree (right_a);

	db_finalize (breakpoint_stmt);
}

static sqlite3_stmt *
prepare_cluster_merging_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Merge the clusters into the retrocopy window
	const char sql[] =
		"SELECT retrocopy_id, chr, MIN(start), MAX(end),\n"
		"	REPLACE(GROUP_CONCAT(DISTINCT gene_name),',','/')\n"
		"FROM cluster AS c\n"
		"INNER JOIN cluster_merging AS m\n"
		"	ON c.id = m.cluster_id AND c.sid = m.cluster_sid\n"
		"GROUP BY retrocopy_id";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

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

			e = hash_lookup (rtc_h, &rid);
			assert (e != NULL);

			/*
			* The insertion point is the first reference base
			* after the insertion. Prefer the right breakpoint,
			* so the insertion (and its TSD) is left aligned.
			* Without split reads, fall back to the window mean
			*/
			if (e->breakpoint_right_sr > 0)
				{
					ip = e->breakpoint_right;
					ip_type = RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE;
				}
			else if (e->breakpoint_left_sr > 0)
				{
					ip = e->breakpoint_left + 1;
					ip_type = RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE;
				}
			else
				{
					ip = (start + end) / 2;
					ip_type = RETROCOPY_INSERTION_POINT_WINDOW_MEAN;
				}

			log_debug ("%d %s %li %li %s %d %li %d %.6f %.6f",
					rid, chr, start, end, gene, e->level, ip, ip_type,
					e->orientation_rho, e->orientation_p_value);

			db_insert_retrocopy (retrocopy_stmt, rid, chr, start, end,
					gene, e->level, ip, ip_type, e->orientation_rho,
					e->orientation_p_value, e->breakpoint_left,
					e->breakpoint_left_sr, e->breakpoint_right,
					e->breakpoint_right_sr);
		}

	db_finalize (cluster_merging_query_stmt);
}

void
retrocopy (sqlite3_stmt *retrocopy_stmt,
		sqlite3_stmt *cluster_merging_stmt,
		int near_gene_dist)
{
	log_trace ("Inside %s", __func__);
	assert (retrocopy_stmt != NULL
			&& cluster_merging_stmt != NULL
			&& near_gene_dist > 0);

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
	merge_cluster (cluster_merging_stmt, filter, near_gene_dist, rtc_h);

	log_info ("Calculate retrocopies orientation");
	calculate_orientation (
			sqlite3_db_handle (retrocopy_stmt),
			rtc_h);

	log_info ("Calculate retrocopies breakpoints");
	calculate_breakpoint (
			sqlite3_db_handle (retrocopy_stmt),
			rtc_h);

	log_info ("Annotate retrocopies");
	annotate_retrocopy (retrocopy_stmt, rtc_h);

	hash_free (rtc_h);
}
