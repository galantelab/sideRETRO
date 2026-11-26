/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br>
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
#include "fasta.h"
#include "str.h"
#include "utils.h"
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

	// Poly-A tail: Its length, if the length is just
	// a lower bound, the strand of the retrocopy it
	// implies and the number of reads supporting it
	int            polya_len;
	int            polya_lower_bound;
	char           polya_strand;
	int            polya_sr;

	// Overlapped parental genes: The ones left
	// out when one of them is chosen ('g1/g2')
	char          *other_genes;
};

typedef struct _RetrocopyEntry RetrocopyEntry;

static void
rtc_entry_free (RetrocopyEntry *e)
{
	xfree (e->other_genes);
	xfree (e);
}

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
		return; // LCOV_EXCL_LINE: static: never NULL

	xfree ((void *) c->cchr);
	xfree ((void *) c->gene_name);
	xfree ((void *) c->gchr);

	xfree (c);
}

static void
clean_retrocopy_tables (sqlite3 *db)
{
	/*
	* Delete all values from previous runs. The genotypes
	* too: They belong to the retrocopies deleted here. If
	* the genotype step fails afterwards, the old ones must
	* not be left for the new retrocopies
	*/
	const char sql[] =
		"DELETE FROM genotype;\n"
		"DELETE FROM intron_retention;\n"
		"DELETE FROM cluster_merging;\n"
		"DELETE FROM junction;\n"
		"DELETE FROM insertion;\n"
		"DELETE FROM evidence;\n"
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

	if (rc)
		return rc;

	/*
	* Sort by the parental gene position: Start and then end.
	* Comparing the start of one with the end of the other was
	* not a consistent order for overlapping genes - each one
	* came before the other - and qsort needs one
	*/
	if (c1->gstart != c2->gstart)
		return c1->gstart < c2->gstart ? -1 : 1;

	return (c1->gend > c2->gend) - (c1->gend < c2->gend);
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

	/*
	* The position of each read at the insertion site and
	* the position of its mate at the parental gene, side
	* by side. Only the mate: The other piece of a split
	* read (primary and supplementary alignments of the same
	* read, 0x40|0x80) lies at the end of the retrocopy next
	* to the junction, whatever the distance from the read
	* to it - a point against the correlation
	*/
	const char sql[] =
		"SELECT DISTINCT cm.retrocopy_id, a.pos, m.pos\n"
		"FROM cluster_merging AS cm\n"
		"INNER JOIN clustering AS c\n"
		"	ON c.cluster_id = cm.cluster_id\n"
		"		AND c.cluster_sid = cm.cluster_sid\n"
		"INNER JOIN alignment AS a\n"
		"	ON a.id = c.alignment_id\n"
		"INNER JOIN alignment AS m\n"
		"	ON m.qname = a.qname\n"
		"		AND m.source_id = a.source_id\n"
		"WHERE m.id != a.id\n"
		"	AND m.type & $EXONIC\n"
		"	AND (m.flag & 0xC0) != (a.flag & 0xC0)\n"
		"ORDER BY cm.retrocopy_id ASC, a.pos ASC, m.pos ASC";

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

					// Calculate p-value. The same seed for each
					// retrocopy: Its p-value depends only on its
					// own points
					seed = SEED;
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
			seed = SEED;
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
		"SELECT DISTINCT retrocopy_id, a.id, a.pos, a.rlen, a.cigar,\n"
		"	(a.flag & 0x800 OR a.type & $SUPPLEMENTARY)\n"
		"FROM cluster_merging AS cm\n"
		"INNER JOIN clustering AS c\n"
		"	USING (cluster_id, cluster_sid)\n"
		"INNER JOIN alignment AS a\n"
		"	ON a.id = c.alignment_id\n"
		"WHERE a.cigar LIKE '%S%' OR a.cigar LIKE '%H%'\n"
		"ORDER BY retrocopy_id ASC, a.pos ASC, a.rlen ASC, a.cigar ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$SUPPLEMENTARY"),
			ABNORMAL_SUPPLEMENTARY);

	return stmt;
}

/*
 * Breakpoint candidate. A candidate is strong if it
 * comes from a true split read - a fragment with a
 * supplementary alignment - and not just from a
 * clipped read, which may be noise as adapters or
 * quality trimming
 */
struct _Breakpoint
{
	long pos;
	int  strong;
	int  aid;
};

typedef struct _Breakpoint Breakpoint;

static int
cmp_breakpoint (const void *p1, const void *p2)
{
	const Breakpoint *b1 = p1;
	const Breakpoint *b2 = p2;
	return (b1->pos > b2->pos) - (b1->pos < b2->pos);
}

static inline void
breakpoint_push (Breakpoint **a, size_t *size, size_t *alloc,
		const long pos, const int strong, const int aid)
{
	if (*size >= *alloc)
		{
			*alloc += BLOCK_SIZE;
			*a = xrealloc (*a, sizeof (Breakpoint) * *alloc);
		}

	(*a)[*size] = (Breakpoint) {.pos = pos, .strong = strong, .aid = aid};
	(*size)++;
}

static void
breakpoint_mode (Breakpoint *a, const size_t size, long *mode, int *support)
{
	long cur = 0;
	int acm = 0;
	int strong = 0;
	size_t i = 0;

	*mode = 0;
	*support = 0;

	if (size == 0)
		return;

	// Sort in order to count repeated positions
	// and keep the first one in case of tie
	qsort (a, size, sizeof (Breakpoint), cmp_breakpoint);

	cur = a[0].pos;
	acm = 1;
	strong = a[0].strong;

	for (i = 1; i <= size; i++)
		{
			if (i < size && a[i].pos == cur)
				{
					acm++;
					strong |= a[i].strong;
					continue;
				}

			/*
			* Accept a position if at least two clipped reads
			* agree on it, or if a single one is a true split
			* read. Random clips can never make a breakpoint
			* by themselves
			*/
			if (acm > *support && (acm > 1 || strong))
				{
					*support = acm;
					*mode = cur;
				}

			if (i < size)
				{
					cur = a[i].pos;
					acm = 1;
					strong = a[i].strong;
				}
		}
}

static inline void
breakpoint_dump_evidence (sqlite3_stmt *evidence_stmt, const int rid,
		const Breakpoint *a, const size_t size, const long mode,
		const int support, const int role)
{
	size_t i = 0;

	if (support == 0)
		return;

	// The reads supporting the accepted breakpoint
	for (i = 0; i < size; i++)
		if (a[i].pos == mode)
			db_insert_evidence (evidence_stmt, rid, a[i].aid, role);
}

static inline void
breakpoint_set_entry (sqlite3_stmt *evidence_stmt, Hash *rtc_h, const int rid,
		Breakpoint *left_a, const size_t left_size, Breakpoint *right_a,
		const size_t right_size)
{
	RetrocopyEntry *e = NULL;

	e = hash_lookup (rtc_h, &rid);
	assert (e != NULL);

	breakpoint_mode (left_a, left_size,
			&e->breakpoint_left, &e->breakpoint_left_sr);

	breakpoint_mode (right_a, right_size,
			&e->breakpoint_right, &e->breakpoint_right_sr);

	breakpoint_dump_evidence (evidence_stmt, rid, left_a, left_size,
			e->breakpoint_left, e->breakpoint_left_sr,
			RETROCOPY_EVIDENCE_SPLIT_LEFT);

	breakpoint_dump_evidence (evidence_stmt, rid, right_a, right_size,
			e->breakpoint_right, e->breakpoint_right_sr,
			RETROCOPY_EVIDENCE_SPLIT_RIGHT);

	log_debug ("Retrocopy %d breakpoints: left %li (%d SR) right %li (%d SR)",
			rid, e->breakpoint_left, e->breakpoint_left_sr,
			e->breakpoint_right, e->breakpoint_right_sr);
}

static void
calculate_breakpoint (sqlite3 *db, sqlite3_stmt *evidence_stmt, Hash *rtc_h)
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
	int strong = 0;

	int rid_prev = 0;

	// Candidates for the left and right breakpoints
	Breakpoint *left_a = NULL;
	size_t left_size = 0;
	size_t left_alloc = 0;

	Breakpoint *right_a = NULL;
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
			strong    = db_column_int   (breakpoint_stmt, 5);

			// First loop
			if (!rid_prev)
				rid_prev = rid;

			// Set breakpoints for the previous retrocopy
			if (rid_prev != rid)
				{
					breakpoint_set_entry (evidence_stmt, rtc_h, rid_prev,
							left_a, left_size, right_a, right_size);

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
						pos + rlen - 1, strong, aid);

			if (clip_left >= RETROCOPY_MIN_CLIP)
				breakpoint_push (&right_a, &right_size, &right_alloc,
						pos, strong, aid);
		}

	// The last or the first retrocopy
	if (rid_prev)
		breakpoint_set_entry (evidence_stmt, rtc_h, rid_prev,
				left_a, left_size, right_a, right_size);

	xfree (left_a);
	xfree (right_a);

	db_finalize (breakpoint_stmt);
}

static sqlite3_stmt *
prepare_polya_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	/*
	* Get all clipped alignments clustered into the
	* retrocopies, with the sequence of the read. A soft
	* clipped one has it. A hard clipped one - a
	* supplementary alignment, as some aligners write it
	* - has no sequence: It comes from the primary
	* alignment of the same read, reverse complemented if
	* they are at different strands (see clip_read_get).
	* The order is given by the reads - not by their ids,
	* which change with the databases merged: The
	* consensus of the clips keeps the first read in
	* case of tie
	*/
	const char sql[] =
		"SELECT DISTINCT retrocopy_id, a.id, a.pos, a.rlen, a.cigar,\n"
		"	CASE\n"
		"		WHEN a.cigar LIKE '%S%' THEN a.seq\n"
		"		ELSE p.seq\n"
		"	END AS read_seq,\n"
		"	CASE\n"
		"		WHEN a.cigar LIKE '%S%' THEN 0\n"
		"		ELSE (a.flag & 0x10) != (p.flag & 0x10)\n"
		"	END\n"
		"FROM cluster_merging AS cm\n"
		"INNER JOIN clustering AS c\n"
		"	USING (cluster_id, cluster_sid)\n"
		"INNER JOIN alignment AS a\n"
		"	ON a.id = c.alignment_id\n"
		"LEFT JOIN alignment AS p\n"
		"	ON a.cigar NOT LIKE '%S%'\n"
		"		AND p.qname = a.qname\n"
		"		AND p.source_id = a.source_id\n"
		"		AND (p.flag & 0xC0) = (a.flag & 0xC0)\n"
		"		AND (p.flag & 0x900) = 0\n"
		"		AND p.cigar NOT LIKE '%H%'\n"
		"		AND p.seq != '*'\n"
		"WHERE (a.cigar LIKE '%S%' AND a.seq != '*')\n"
		"	OR (a.cigar LIKE '%H%' AND p.id IS NOT NULL)\n"
		"ORDER BY retrocopy_id ASC, a.pos ASC, a.cigar ASC, read_seq ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	return stmt;
}

/*
 * The cigar and the sequence of a clipped alignment from
 * prepare_polya_stmt. A hard clipped one is given as soft
 * clipped: Its 'H' turn into 'S', with the sequence of the
 * primary alignment of the read - reverse complemented if
 * at the other strand. The buffers are reused between the
 * calls. Return 0 if the sequence does not fit the cigar
 */
static int
clip_read_get (sqlite3_stmt *stmt, const char **cigar_str, const char **seq,
		char **cigar_buf, char **seq_buf)
{
	const char *c = NULL;
	long qlen = 0;
	long n = 0;
	int rev = 0;
	int i = 0;

	*cigar_str = db_column_text (stmt, 4);
	*seq = db_column_text (stmt, 5);
	rev = db_column_int (stmt, 6);

	// Soft clipped: As stored
	if (strchr (*cigar_str, 'S') != NULL)
		return 1;

	xfree (*cigar_buf);
	xfree (*seq_buf);

	*cigar_buf = xstrdup (*cigar_str);
	*seq_buf = xstrdup (*seq);

	// The bases of the read along the cigar
	for (c = *cigar_buf, i = 0; c[i] != '\0'; i++)
		{
			if (c[i] >= '0' && c[i] <= '9')
				{
					n = n * 10 + (c[i] - '0');
					continue;
				}

			if (c[i] == 'H')
				(*cigar_buf)[i] = 'S';

			if (strchr ("MIHS=X", c[i]) != NULL)
				qlen += n;

			n = 0;
		}

	if (qlen != (long) strlen (*seq_buf))
		return 0;

	if (rev)
		revcomp (*seq_buf, qlen);

	*cigar_str = *cigar_buf;
	*seq = *seq_buf;

	return 1;
}

/*
 * Poly-A evidence accumulated for one strand
 */
struct _PolyA
{
	int  exact;
	int  lower_bound;
	int  sr;

	// The reads supporting this strand
	int *aids;
	int  aids_alloc;
};

typedef struct _PolyA PolyA;

static inline int
polya_run_at (const char *seq, const int len, const char base,
		const int reverse, const int offset)
{
	int run = 0;
	int i = 0;

	// Count the run of 'base' from the end (reverse)
	// or from the beginning of seq, after offset bases
	for (i = offset; i < len; i++)
		{
			char c = reverse ? seq[len - 1 - i] : seq[i];

			if (c != base && c != (base + 32))
				break;

			run++;
		}

	return run;
}

static inline int
polya_run (const char *seq, const int len, const char base,
		const int reverse, int *offset)
{
	int run = 0;
	int i = 0;

	/*
	* The tail may not begin right at the flank: Allow
	* a few bases before it, but only for a run long
	* enough to be called as the poly-A
	*/
	for (i = 0; i <= RETROCOPY_MAX_POLYA_OFFSET && i < len; i++)
		{
			run = polya_run_at (seq, len, base, reverse, i);
			if (run >= RETROCOPY_MIN_POLYA)
				{
					*offset = i;
					return run;
				}
		}

	// No tail: The run at the flank, if any
	*offset = 0;
	return polya_run_at (seq, len, base, reverse, 0);
}

static inline void
polya_add (PolyA *p, const int run, const int clip_len, const int aid)
{
	if (run < RETROCOPY_MIN_POLYA)
		return;

	if (p->sr >= p->aids_alloc)
		{
			p->aids_alloc += BLOCK_SIZE;
			p->aids = xrealloc (p->aids, sizeof (int) * p->aids_alloc);
		}

	p->aids[p->sr++] = aid;

	// If the whole clip is poly-A, the tail
	// extends beyond the read: Lower bound
	if (run == clip_len)
		{
			if (run > p->lower_bound)
				p->lower_bound = run;
		}
	else if (run > p->exact)
		p->exact = run;
}

static inline void
polya_reset (PolyA *p)
{
	int *aids = p->aids;
	int aids_alloc = p->aids_alloc;

	// Keep the buffer
	*p = (PolyA) {.aids = aids, .aids_alloc = aids_alloc};
}

static inline void
polya_set_entry (sqlite3_stmt *evidence_stmt, Hash *rtc_h, const int rid,
		const PolyA *plus, const PolyA *minus)
{
	RetrocopyEntry *e = NULL;
	const PolyA *p = NULL;
	int i = 0;

	e = hash_lookup (rtc_h, &rid);
	assert (e != NULL);

	if (plus->sr == 0 && minus->sr == 0)
		return;

	// Keep the strand with more reads. Prefer
	// '+' in case of tie
	if (plus->sr >= minus->sr)
		{
			p = plus;
			e->polya_strand = '+';
		}
	else
		{
			p = minus;
			e->polya_strand = '-';
		}

	e->polya_sr = p->sr;

	for (i = 0; i < p->sr; i++)
		db_insert_evidence (evidence_stmt, rid, p->aids[i],
				RETROCOPY_EVIDENCE_POLYA);

	// A read with non poly-A bases before the tail
	// gives its exact length. Otherwise all we know
	// is the longest run
	if (p->exact > 0)
		{
			e->polya_len = p->exact;
			e->polya_lower_bound = 0;
		}
	else
		{
			e->polya_len = p->lower_bound;
			e->polya_lower_bound = 1;
		}

	log_debug ("Retrocopy %d poly-A: %d bases%s at strand '%c' (%d reads)",
			rid, e->polya_len, e->polya_lower_bound ? " (lower bound)" : "",
			e->polya_strand, e->polya_sr);
}

static void
locate_soft_clips (const uint32_t *cigar, const uint32_t n_cigar, const long pos,
		CigarOp *left_op, int *has_left, CigarOp *right_op, int *has_right)
{
	CigarIter iter = {};
	CigarOp op = {};
	int seen_aligned = 0;

	*has_left = 0;
	*has_right = 0;

	// The soft clip at the beginning comes before any
	// aligned operation and the one at the end is followed
	// only by hard clips
	cigar_iter_init (&iter, cigar, n_cigar, pos);

	while (cigar_iter_next (&iter, &op))
		{
			if (op.type == CIGAR_HARD_CLIP)
				continue;

			if (op.type == CIGAR_SOFT_CLIP)
				{
					if (!seen_aligned)
						{
							*left_op = op;
							*has_left = 1;
						}
					else
						{
							*right_op = op;
							*has_right = 1;
						}
				}
			else
				{
					seen_aligned = 1;
					*has_right = 0;
				}
		}
}

static void
calculate_polya (sqlite3 *db, sqlite3_stmt *evidence_stmt, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *polya_stmt = NULL;
	RetrocopyEntry *e = NULL;

	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	CigarOp left_op = {};
	CigarOp right_op = {};
	int has_left = 0;
	int has_right = 0;
	char *clip = NULL;
	int run = 0;
	int offset = 0;

	int rid = 0;
	int aid = 0;
	long pos = 0;
	long rlen = 0;
	const char *cigar_str = NULL;
	const char *seq = NULL;
	char *cigar_buf = NULL;
	char *seq_buf = NULL;

	int rid_prev = 0;

	PolyA plus = {};
	PolyA minus = {};

	polya_stmt = prepare_polya_stmt (db);

	while (db_step (polya_stmt) == SQLITE_ROW)
		{
			rid       = db_column_int   (polya_stmt, 0);
			aid       = db_column_int   (polya_stmt, 1);
			pos       = db_column_int64 (polya_stmt, 2);
			rlen      = db_column_int64 (polya_stmt, 3);

			if (!clip_read_get (polya_stmt, &cigar_str, &seq,
						&cigar_buf, &seq_buf))
				{
					log_warn ("Skip alignment %d: Its read does not fit the cigar", aid);
					continue;
				}

			// First loop
			if (!rid_prev)
				rid_prev = rid;

			// Set poly-A for the previous retrocopy
			if (rid_prev != rid)
				{
					polya_set_entry (evidence_stmt, rtc_h, rid_prev, &plus, &minus);

					polya_reset (&plus);
					polya_reset (&minus);
					rid_prev = rid;
				}

			e = hash_lookup (rtc_h, &rid);
			assert (e != NULL);

			// Without breakpoints there is no junction
			if (e->breakpoint_left_sr == 0 && e->breakpoint_right_sr == 0)
				continue;

			cigar = cigar_yyparser (cigar_str, &n_cigar);
			if (cigar == NULL)
				{
					log_warn ("Skip alignment %d with invalid cigar '%s'",
							aid, cigar_str);
					continue;
				}

			// Locate the soft clips at the beginning and
			// at the end of the read
			locate_soft_clips (cigar, n_cigar, pos,
					&left_op, &has_left, &right_op, &has_right);

			xfree (cigar);

			/*
			* Retrocopy at '+': The poly-A comes right before the
			* right flank. The reads clipped at the left at the right
			* breakpoint end their clips with a run of A.
			* Retrocopy at '-': The poly-A is reverse complemented
			* right after the left flank. The reads clipped at the
			* right at the left breakpoint begin their clips with
			* a run of T
			*/
			if (has_left && e->breakpoint_right_sr > 0
					&& pos == e->breakpoint_right)
				{
					clip = cigar_op_seq (&left_op, seq);
					if (clip != NULL)
						{
							run = polya_run (clip, left_op.qlen, 'A', 1, &offset);
							polya_add (&plus, run, left_op.qlen - offset, aid);
						}
					xfree (clip);
				}

			if (has_right && e->breakpoint_left_sr > 0
					&& (pos + rlen - 1) == e->breakpoint_left)
				{
					clip = cigar_op_seq (&right_op, seq);
					if (clip != NULL)
						{
							run = polya_run (clip, right_op.qlen, 'T', 0, &offset);
							polya_add (&minus, run, right_op.qlen - offset, aid);
						}
					xfree (clip);
				}
		}

	// The last or the first retrocopy
	if (rid_prev)
		polya_set_entry (evidence_stmt, rtc_h, rid_prev, &plus, &minus);

	xfree (plus.aids);
	xfree (minus.aids);
	xfree (cigar_buf);
	xfree (seq_buf);

	db_finalize (polya_stmt);
}

/*
 * A read spanning an exon-exon junction of the retrocopy
 * is aligned to the parental gene in two pieces: One
 * alignment at each exon. So, search the fragments of
 * each retrocopy for reads with two alignments on two
 * different exons - not overlapping each other - of the
 * parental gene. The read is identified by its qname,
 * source and the flag bits 0x40|0x80.
 *
 * Exons of different transcripts overlap, so the same
 * pieces overlap many exon pairs. The junction is the one
 * whose boundaries are met by the pieces: The left piece
 * ends - clipped - at the end of the first exon, and the
 * right piece starts - clipped - at the start of the
 * second exon. A piece may stop up to ?1 bases before the
 * boundary (mismatches near the junction), or go up to ?2
 * bases beyond it, into the intron: The sequence shared by
 * the end of one exon and the intron before the other
 * (microhomology) is aligned at both pieces. The closest
 * pair wins, and the exons with the same boundaries are
 * represented by the lowest ids
 */
#define JUNCTION_QUERY_WITH \
	"WITH\n" \
	"	fragment (rid, gene_name, qname, source_id) AS (\n" \
	"		SELECT DISTINCT cm.retrocopy_id, cl.gene_name,\n" \
	"			a.qname, a.source_id\n" \
	"		FROM cluster_merging AS cm\n" \
	"		INNER JOIN cluster AS cl\n" \
	"			ON cl.id = cm.cluster_id AND cl.sid = cm.cluster_sid\n" \
	"		INNER JOIN clustering AS c\n" \
	"			ON c.cluster_id = cm.cluster_id\n" \
	"				AND c.cluster_sid = cm.cluster_sid\n" \
	"		INNER JOIN alignment AS a\n" \
	"			ON a.id = c.alignment_id\n" \
	"	),\n" \
	"	exonic (qname, source_id, read, aid, apos, aend, cigar,\n" \
	"		eid, gene_name, start, end) AS (\n" \
	"		SELECT a.qname, a.source_id, a.flag & 0xC0, a.id,\n" \
	"			a.pos, a.pos + a.rlen - 1, a.cigar,\n" \
	"			e.id, e.gene_name, e.start, e.end\n" \
	"		FROM alignment AS a\n" \
	"		INNER JOIN overlapping AS o\n" \
	"			ON o.alignment_id = a.id\n" \
	"		INNER JOIN exon AS e\n" \
	"			ON e.id = o.exon_id\n" \
	"	),\n" \
	"	candidate (rid, aid, eid1, eid2, end, start, rank) AS (\n" \
	"		SELECT f.rid, x1.aid, x1.eid, x2.eid, x1.end, x2.start,\n" \
	"			ROW_NUMBER () OVER (\n" \
	"				PARTITION BY f.rid, x1.aid, x2.aid\n" \
	"				ORDER BY ABS (x1.aend - x1.end) + ABS (x2.apos - x2.start),\n" \
	"					x1.eid, x2.eid)\n" \
	"		FROM fragment AS f\n" \
	"		INNER JOIN exonic AS x1\n" \
	"			ON x1.qname = f.qname\n" \
	"				AND x1.source_id = f.source_id\n" \
	"				AND x1.gene_name = f.gene_name\n" \
	"		INNER JOIN exonic AS x2\n" \
	"			ON x2.qname = x1.qname\n" \
	"				AND x2.source_id = x1.source_id\n" \
	"				AND x2.read = x1.read\n" \
	"				AND x2.gene_name = x1.gene_name\n" \
	"				AND x2.aid != x1.aid\n" \
	"		WHERE x1.end < x2.start\n" \
	"			AND x1.aend - x1.end BETWEEN -?1 AND ?2\n" \
	"			AND x2.start - x2.apos BETWEEN -?1 AND ?2\n" \
	"			AND x1.cigar GLOB '*[SH]'\n" \
	"			AND LTRIM (x2.cigar, '0123456789') GLOB '[SH]*'\n" \
	"	),\n" \
	"	matched (rid, aid, eid1, eid2, end, start) AS (\n" \
	"		SELECT rid, aid, eid1, eid2, end, start\n" \
	"		FROM candidate\n" \
	"		WHERE rank = 1\n" \
	"	)\n"

static sqlite3_stmt *
prepare_junction_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Each distinct junction once, with its reads
	const char sql[] =
		JUNCTION_QUERY_WITH
		"SELECT rid, MIN(eid1), MIN(eid2), COUNT(DISTINCT aid)\n"
		"FROM matched\n"
		"GROUP BY rid, end, start\n"
		"ORDER BY rid ASC, end ASC, start ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);
	db_bind_int (stmt, 1, RETROCOPY_JUNCTION_MAX_INNER);
	db_bind_int (stmt, 2, RETROCOPY_JUNCTION_MAX_OUTER);

	return stmt;
}

static void
calculate_junction (sqlite3_stmt *junction_stmt, sqlite3_stmt *evidence_stmt)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *junction_query_stmt = NULL;

	int rid = 0;
	int exon_id1 = 0;
	int exon_id2 = 0;
	int reads = 0;

	junction_query_stmt = prepare_junction_query_stmt (
			sqlite3_db_handle (junction_stmt));

	while (db_step (junction_query_stmt) == SQLITE_ROW)
		{
			rid      = db_column_int (junction_query_stmt, 0);
			exon_id1 = db_column_int (junction_query_stmt, 1);
			exon_id2 = db_column_int (junction_query_stmt, 2);
			reads    = db_column_int (junction_query_stmt, 3);

			log_debug ("Retrocopy %d exon-exon junction %d-%d (%d reads)",
					rid, exon_id1, exon_id2, reads);

			db_insert_junction (junction_stmt, rid,
					exon_id1, exon_id2, reads);
		}

	db_finalize (junction_query_stmt);

	// The reads spanning the junctions
	const char sql[] =
		JUNCTION_QUERY_WITH
		"SELECT DISTINCT rid, aid\n"
		"FROM matched";

	log_debug ("Query schema:\n%s", sql);
	junction_query_stmt = db_prepare (sqlite3_db_handle (junction_stmt), sql);
	db_bind_int (junction_query_stmt, 1, RETROCOPY_JUNCTION_MAX_INNER);
	db_bind_int (junction_query_stmt, 2, RETROCOPY_JUNCTION_MAX_OUTER);

	while (db_step (junction_query_stmt) == SQLITE_ROW)
		db_insert_evidence (evidence_stmt,
				db_column_int (junction_query_stmt, 0),
				db_column_int (junction_query_stmt, 1),
				RETROCOPY_EVIDENCE_JUNCTION);

	db_finalize (junction_query_stmt);
}

/*
 * Exon of the parental gene with its sequence
 * (upper case) from the reference genome
 */
struct _ExonSeq
{
	int          id;
	long         start;
	long         end;
	const char  *seq;
};

typedef struct _ExonSeq ExonSeq;

/*
 * Position of a junction clip into the parental gene.
 * 'sense' is set if the insertion is read 5' => 3'
 * at the forward strand of the reference
 */
struct _Hit
{
	long pos;
	int  exon_id;
	int  sense;
};

typedef struct _Hit Hit;

/*
 * One end of the retrocopy into the parental gene
 */
struct _InsertionEnd
{
	long pos;
	int  exon_id;
	int  sr;
};

typedef struct _InsertionEnd InsertionEnd;

static sqlite3_stmt *
prepare_gene_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Parental genes of each retrocopy
	const char sql[] =
		"SELECT DISTINCT cm.retrocopy_id, cl.gene_name\n"
		"FROM cluster_merging AS cm\n"
		"INNER JOIN cluster AS cl\n"
		"	ON cl.id = cm.cluster_id AND cl.sid = cm.cluster_sid\n"
		"ORDER BY cm.retrocopy_id ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	return stmt;
}

static Hash *
index_parental_gene (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *gene_stmt = NULL;
	Hash *gene_h = NULL;
	char *gene = NULL;
	int *rid_alloc = NULL;
	int rid = 0;

	gene_h = hash_new_full (int_hash, int_equal, xfree, xfree);
	gene_stmt = prepare_gene_query_stmt (db);

	while (db_step (gene_stmt) == SQLITE_ROW)
		{
			rid = db_column_int (gene_stmt, 0);
			gene = hash_lookup (gene_h, &rid);

			// A retrocopy with more than one parental
			// gene is ambiguous: Keep an empty name
			if (gene != NULL)
				{
					gene[0] = '\0';
					continue;
				}

			rid_alloc = xcalloc (1, sizeof (int));
			*rid_alloc = rid;

			hash_insert (gene_h, rid_alloc,
					xstrdup (db_column_text (gene_stmt, 1)));
		}

	db_finalize (gene_stmt);
	return gene_h;
}

static sqlite3_stmt *
prepare_exon_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"SELECT id, chr, start, end, strand\n"
		"FROM exon\n"
		"WHERE gene_name = ?1\n"
		"ORDER BY start ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	return stmt;
}

static sqlite3_stmt *
prepare_exon_region_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	// Exons of all parental genes: The only regions
	// of the genome needed
	const char sql[] =
		"SELECT DISTINCT e.chr, e.start, e.end\n"
		"FROM exon AS e\n"
		"WHERE e.gene_name IN (SELECT DISTINCT gene_name FROM cluster)";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	return stmt;
}

static Hash *
index_exon_regions (sqlite3 *db, const char *fasta_file)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *region_stmt = NULL;
	List *regions = NULL;
	Hash *fidx = NULL;

	regions = list_new ((DestroyNotify) fasta_region_free);
	region_stmt = prepare_exon_region_query_stmt (db);

	while (db_step (region_stmt) == SQLITE_ROW)
		list_append (regions, fasta_region_new (
					db_column_text (region_stmt, 0),
					db_column_int64 (region_stmt, 1),
					db_column_int64 (region_stmt, 2)));

	log_info ("Extract %zu parental exons from '%s'", list_size (regions), fasta_file);
	fidx = fasta_index_regions (fasta_file, regions, NULL);

	list_free (regions);
	db_finalize (region_stmt);

	return fidx;
}

static int
exon_seq_load (sqlite3_stmt *exon_stmt, Hash *fidx, const char *gene_name,
		ExonSeq **exons, char *strand)
{
	const char *seq = NULL;
	const char *chr = NULL;
	ExonSeq *e = NULL;
	char *key = NULL;
	int n_exons = 0;
	int alloc = 0;

	*exons = NULL;
	*strand = '+';

	db_reset (exon_stmt);
	db_clear_bindings (exon_stmt);
	db_bind_text (exon_stmt, 1, gene_name);

	while (db_step (exon_stmt) == SQLITE_ROW)
		{
			if (n_exons >= alloc)
				{
					alloc += BLOCK_SIZE;
					*exons = xrealloc (*exons, sizeof (ExonSeq) * alloc);
				}

			e = &(*exons)[n_exons];

			e->id    = db_column_int   (exon_stmt, 0);
			chr      = db_column_text  (exon_stmt, 1);
			e->start = db_column_int64 (exon_stmt, 2);
			e->end   = db_column_int64 (exon_stmt, 3);
			*strand  = db_column_text  (exon_stmt, 4)[0];

			key = fasta_region_key (chr, e->start, e->end);
			seq = hash_lookup (fidx, key);
			xfree (key);

			if (seq == NULL)
				{
					log_warn ("Exon %d at %s:%li-%li not found at the reference genome",
							e->id, chr, e->start, e->end);
					continue;
				}

			e->seq = seq;
			n_exons++;
		}

	return n_exons;
}

/*
 * Search the seed (and its reverse complement) into
 * the exons. 'side' is 0 for a clip at the left
 * junction, whose first base is the insertion boundary,
 * and 1 for a clip at the right junction, whose last
 * base is the boundary. Return the number of distinct
 * positions found and the first one
 */
static int
search_seed (const ExonSeq *exons, const int n_exons, const char *seed,
		const int seed_len, const int side, const char strand, Hit *hit)
{
	char *fwd = NULL;
	char *rev = NULL;
	const char *found = NULL;
	const char *haystack = NULL;
	const ExonSeq *e = NULL;
	long pos = 0;
	int n_hits = 0;
	int forward = 0;
	int i = 0;
	int k = 0;

	fwd = xcalloc (seed_len + 1, sizeof (char));
	memcpy (fwd, seed, seed_len);

	rev = xcalloc (seed_len + 1, sizeof (char));
	memcpy (rev, seed, seed_len);
	revcomp (rev, seed_len);

	for (i = 0; i < n_exons; i++)
		{
			e = &exons[i];

			for (k = 0; k < 2; k++)
				{
					forward = (k == 0);
					haystack = e->seq;

					while ((found = strstr (haystack, forward ? fwd : rev)) != NULL)
						{
							/*
							* The boundary base of the clip: For a forward
							* match the first base of the seed is its first
							* base at the reference, for a reverse complement
							* match it is the last one
							*/
							if ((side == 0) == forward)
								pos = e->start + (found - e->seq);
							else
								pos = e->start + (found - e->seq) + seed_len - 1;

							// Overlapping exons may give the same position
							if (n_hits == 0 || hit->pos != pos)
								{
									if (n_hits == 0)
										{
											hit->pos = pos;
											hit->exon_id = e->id;
											hit->sense = forward
												? strand == '+'
												: strand == '-';
										}

									n_hits++;
								}

							haystack = found + 1;
						}
				}
		}

	xfree (fwd);
	xfree (rev);

	return n_hits;
}

static int
cmp_hit (const void *p1, const void *p2)
{
	const Hit *h1 = p1;
	const Hit *h2 = p2;
	return (h1->pos > h2->pos) - (h1->pos < h2->pos);
}

static inline void
hit_push (Hit **a, size_t *size, size_t *alloc, const Hit *hit)
{
	if (*size >= *alloc)
		{
			*alloc += BLOCK_SIZE;
			*a = xrealloc (*a, sizeof (Hit) * *alloc);
		}

	(*a)[*size] = *hit;
	(*size)++;
}

static void
hit_mode (Hit *a, const size_t size, InsertionEnd *end)
{
	long cur = 0;
	int exon_id = 0;
	int acm = 0;
	size_t i = 0;

	*end = (InsertionEnd) {};

	if (size == 0)
		return;

	qsort (a, size, sizeof (Hit), cmp_hit);

	cur = a[0].pos;
	exon_id = a[0].exon_id;
	acm = 1;

	for (i = 1; i <= size; i++)
		{
			if (i < size && a[i].pos == cur)
				{
					acm++;
					continue;
				}

			if (acm > end->sr)
				{
					end->sr = acm;
					end->pos = cur;
					end->exon_id = exon_id;
				}

			if (i < size)
				{
					cur = a[i].pos;
					exon_id = a[i].exon_id;
					acm = 1;
				}
		}
}

static char *
predict_insertion (const ExonSeq *exons, const int n_exons, const char strand,
		const long five_pos, const long three_pos, const int polya_len, int *len)
{
	String *ins = NULL;
	char *segment = NULL;
	char *seq = NULL;
	long covered = 0;
	long lo = 0;
	long hi = 0;
	long a = 0;
	long b = 0;
	int i = 0;

	ins = string_sized_new (BUFSIZ);

	lo = strand == '+' ? five_pos : three_pos;
	hi = strand == '+' ? three_pos : five_pos;

	/*
	* Concatenate the exon bases between the ends. Exons
	* are sorted by start, so skip the bases already covered
	* by a previous overlapping exon (alternative exons)
	*/
	for (i = 0; i < n_exons; i++)
		{
			a = exons[i].start;
			b = exons[i].end;

			if (a < lo)
				a = lo;

			if (a <= covered)
				a = covered + 1;

			if (b > hi)
				b = hi;

			if (a > b)
				continue;

			segment = xcalloc (b - a + 2, sizeof (char));
			memcpy (segment, exons[i].seq + (a - exons[i].start), b - a + 1);
			string_concat (ins, segment);
			xfree (segment);

			covered = b;
		}

	// Gene at the reverse strand: The mRNA is the
	// reverse complement of the exons
	if (strand == '-')
		revcomp (ins->str, ins->len);

	// Poly-A tail
	for (i = 0; i < polya_len; i++)
		string_concat (ins, "A");

	*len = ins->len;
	seq = string_free (ins, 0);

	return seq;
}

static char *
clip_consensus (char **clips, const size_t n, const int from_end)
{
	char *cons = NULL;
	const char *clip = NULL;
	size_t max_len = 0;
	size_t len = 0;
	size_t col = 0;
	size_t i = 0;
	int count[5] = {};
	int best = 0;
	int tie = 0;
	int k = 0;
	char c = 0;

	if (n == 0)
		return NULL;

	for (i = 0; i < n; i++)
		{
			len = strlen (clips[i]);
			if (len > max_len)
				max_len = len;
		}

	cons = xcalloc (max_len + 1, sizeof (char));

	/*
	* The clips are anchored at the junction: At their
	* beginning (left junction) or at their end (right
	* junction). Vote base by base from the junction.
	* A tie is unknown: 'N' - it does not depend on the
	* order of the reads
	*/
	for (col = 0; col < max_len; col++)
		{
			memset (count, 0, sizeof (count));

			for (i = 0; i < n; i++)
				{
					clip = clips[i];
					len = strlen (clip);

					if (col >= len)
						continue;

					c = from_end ? clip[len - 1 - col] : clip[col];

					switch (c)
						{
						case 'A': {k = 0; break;}
						case 'C': {k = 1; break;}
						case 'G': {k = 2; break;}
						case 'T': {k = 3; break;}
						default:  {k = 4; break;}
						}

					count[k]++;
				}

			best = 0;
			tie = 0;

			for (k = 1; k < 5; k++)
				{
					if (count[k] > count[best])
						{
							best = k;
							tie = 0;
						}
					else if (count[k] == count[best])
						tie = 1;
				}

			cons[from_end ? max_len - 1 - col : col] = tie ? 'N' : "ACGTN"[best];
		}

	return cons;
}

static void
insertion_set_entry (sqlite3_stmt *insertion_stmt, const int rid,
		const RetrocopyEntry *e, const ExonSeq *exons, const int n_exons,
		const char strand, Hit *left_a, const size_t left_size,
		Hit *right_a, const size_t right_size, char **left_clips,
		const size_t n_left_clips, char **right_clips, const size_t n_right_clips)
{
	InsertionEnd left = {};
	InsertionEnd right = {};
	InsertionEnd *five = NULL;
	InsertionEnd *three = NULL;
	char *seq = NULL;
	char *left_seq = NULL;
	char *right_seq = NULL;
	long three_pos = 0;
	int sense_acm = 0;
	int sense = 0;
	int len = 0;
	size_t i = 0;

	if (n_left_clips == 0 && n_right_clips == 0)
		return;

	// Observed sequences at the junctions
	left_seq = clip_consensus (left_clips, n_left_clips, 0);
	right_seq = clip_consensus (right_clips, n_right_clips, 1);

	/*
	* Is the insertion read 5' => 3' at the forward strand
	* of the reference? The poly-A tells it: Its junction is
	* the 3' end. Without a poly-A, the majority of the hits.
	* They may be misled by an inverted piece at the 5' end
	* of the retrocopy, whose reads hit at the other sense
	*/
	if (e->polya_sr > 0)
		{
			sense = e->polya_strand == '+';
		}
	else
		{
			for (i = 0; i < left_size; i++)
				sense_acm += left_a[i].sense;

			for (i = 0; i < right_size; i++)
				sense_acm += right_a[i].sense;

			sense = 2 * sense_acm >= (int) (left_size + right_size);
		}

	hit_mode (left_a, left_size, &left);
	hit_mode (right_a, right_size, &right);

	// If the insertion is at the sense, the left
	// junction is the 5' end of the retrocopy
	if (sense)
		{
			five = &left;
			three = &right;
		}
	else
		{
			five = &right;
			three = &left;
		}

	// Predict the inserted sequence from the 5' end to
	// the 3' end. If the 3' end was not observed, assume
	// the end of the last exon
	if (five->sr > 0)
		{
			/*
			* The exons are sorted by start: At '-', the last exon
			* of the mRNA is the first one. At '+', look for the
			* greatest end, as an alternative exon may end after
			* the last one
			*/
			three_pos = strand == '+'
				? exons[0].end
				: exons[0].start;

			if (strand == '+')
				for (i = 1; i < (size_t) n_exons; i++)
					if (exons[i].end > three_pos)
						three_pos = exons[i].end;

			if (three->sr > 0)
				three_pos = three->pos;

			seq = predict_insertion (exons, n_exons, strand,
					five->pos, three_pos, e->polya_sr > 0 ? e->polya_len : 0,
					&len);
		}

	log_debug ("Retrocopy %d 5' end %li at exon %d (%d reads), "
			"3' end %li at exon %d (%d reads), predicted length %d",
			rid, five->pos, five->exon_id, five->sr,
			three->pos, three->exon_id, three->sr, len);

	db_insert_insertion (insertion_stmt, rid,
			five->pos, five->exon_id, five->sr,
			three->pos, three->exon_id, three->sr,
			len, seq, left_seq, right_seq);

	xfree (seq);
	xfree (left_seq);
	xfree (right_seq);
}

static inline void
clip_push (char ***a, size_t *size, size_t *alloc, char *clip)
{
	if (*size >= *alloc)
		{
			*alloc += BLOCK_SIZE;
			*a = xrealloc (*a, sizeof (char *) * *alloc);
		}

	(*a)[*size] = clip;
	(*size)++;
}

static inline void
clip_clear (char **a, size_t *size)
{
	size_t i = 0;

	for (i = 0; i < *size; i++)
		xfree (a[i]);

	*size = 0;
}

static void
calculate_insertion (sqlite3 *db, sqlite3_stmt *insertion_stmt,
		const char *fasta_file, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *clip_stmt = NULL;
	sqlite3_stmt *exon_stmt = NULL;
	RetrocopyEntry *e = NULL;

	Hash *fidx = NULL;
	Hash *gene_h = NULL;
	const char *gene = NULL;

	ExonSeq *exons = NULL;
	int n_exons = 0;
	char strand = '+';

	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	CigarOp left_op = {};
	CigarOp right_op = {};
	int has_left = 0;
	int has_right = 0;
	char *clip = NULL;
	int run = 0;
	int offset = 0;
	int seed_len = 0;
	int n_hits = 0;
	Hit hit = {};

	int rid = 0;
	int aid = 0;
	long pos = 0;
	long rlen = 0;
	const char *cigar_str = NULL;
	const char *seq = NULL;
	char *cigar_buf = NULL;
	char *seq_buf = NULL;

	int rid_prev = 0;

	// Hits for the left and right junctions
	Hit *left_a = NULL;
	size_t left_size = 0;
	size_t left_alloc = 0;

	Hit *right_a = NULL;
	size_t right_size = 0;
	size_t right_alloc = 0;

	// Clipped sequences at the left and right junctions
	char **left_clips = NULL;
	size_t n_left_clips = 0;
	size_t left_clips_alloc = 0;

	char **right_clips = NULL;
	size_t n_right_clips = 0;
	size_t right_clips_alloc = 0;

	// Without the reference genome, only the observed
	// sequences at the junctions are kept
	if (fasta_file != NULL)
		{
			fidx = index_exon_regions (db, fasta_file);
			gene_h = index_parental_gene (db);
			exon_stmt = prepare_exon_query_stmt (db);
		}

	// Same soft clipped alignments used for the poly-A
	clip_stmt = prepare_polya_stmt (db);

	while (db_step (clip_stmt) == SQLITE_ROW)
		{
			rid       = db_column_int   (clip_stmt, 0);
			aid       = db_column_int   (clip_stmt, 1);
			pos       = db_column_int64 (clip_stmt, 2);
			rlen      = db_column_int64 (clip_stmt, 3);

			if (!clip_read_get (clip_stmt, &cigar_str, &seq,
						&cigar_buf, &seq_buf))
				continue;

			// New retrocopy: Set the previous one and
			// load the exons of the parental gene
			if (rid_prev != rid)
				{
					if (rid_prev)
						insertion_set_entry (insertion_stmt, rid_prev,
								hash_lookup (rtc_h, &rid_prev), exons, n_exons,
								strand, left_a, left_size, right_a, right_size,
								left_clips, n_left_clips, right_clips, n_right_clips);

					left_size = 0;
					right_size = 0;
					clip_clear (left_clips, &n_left_clips);
					clip_clear (right_clips, &n_right_clips);
					rid_prev = rid;

					xfree (exons);
					exons = NULL;
					n_exons = 0;

					if (fidx != NULL)
						{
							gene = hash_lookup (gene_h, &rid);

							if (gene == NULL || gene[0] == '\0')
								log_debug ("Retrocopy %d has no single parental gene", rid);
							else
								n_exons = exon_seq_load (exon_stmt, fidx, gene,
										&exons, &strand);
						}
				}

			e = hash_lookup (rtc_h, &rid);
			assert (e != NULL);

			// Without breakpoints there is no junction
			if (e->breakpoint_left_sr == 0 && e->breakpoint_right_sr == 0)
				continue;

			cigar = cigar_yyparser (cigar_str, &n_cigar);
			if (cigar == NULL)
				{
					log_warn ("Skip alignment %d with invalid cigar '%s'",
							aid, cigar_str);
					continue;
				}

			locate_soft_clips (cigar, n_cigar, pos,
					&left_op, &has_left, &right_op, &has_right);

			xfree (cigar);

			/*
			* Left junction: The clip begins with the insertion.
			* Keep the clip and, if there are exons, skip a poly-T
			* (the reverse complemented poly-A) and search the seed
			* from the beginning of the clip
			*/
			if (has_right && e->breakpoint_left_sr > 0
					&& (pos + rlen - 1) == e->breakpoint_left)
				{
					clip = cigar_op_seq (&right_op, seq);
					if (clip != NULL && n_exons > 0)
						{
							// Skip the bases before the tail too
							run = polya_run (clip, right_op.qlen, 'T', 0, &offset);
							run = run < RETROCOPY_MIN_POLYA ? 0 : run + offset;

							seed_len = right_op.qlen - run;
							if (seed_len > RETROCOPY_SEED_LEN)
								seed_len = RETROCOPY_SEED_LEN;

							if (seed_len >= RETROCOPY_MIN_SEED)
								{
									n_hits = search_seed (exons, n_exons, clip + run,
											seed_len, 0, strand, &hit);

									if (n_hits == 1)
										hit_push (&left_a, &left_size, &left_alloc, &hit);
									else
										log_debug ("Alignment %d left junction seed has %d hits",
												aid, n_hits);
								}
						}

					if (clip != NULL)
						clip_push (&left_clips, &n_left_clips, &left_clips_alloc, clip);
				}

			/*
			* Right junction: The clip ends with the insertion.
			* Keep the clip and, if there are exons, skip the poly-A
			* and search the seed from the end of the clip
			*/
			if (has_left && e->breakpoint_right_sr > 0
					&& pos == e->breakpoint_right)
				{
					clip = cigar_op_seq (&left_op, seq);
					if (clip != NULL && n_exons > 0)
						{
							run = polya_run (clip, left_op.qlen, 'A', 1, &offset);
							run = run < RETROCOPY_MIN_POLYA ? 0 : run + offset;

							seed_len = left_op.qlen - run;
							if (seed_len > RETROCOPY_SEED_LEN)
								seed_len = RETROCOPY_SEED_LEN;

							if (seed_len >= RETROCOPY_MIN_SEED)
								{
									n_hits = search_seed (exons, n_exons,
											clip + left_op.qlen - run - seed_len,
											seed_len, 1, strand, &hit);

									if (n_hits == 1)
										hit_push (&right_a, &right_size, &right_alloc, &hit);
									else
										log_debug ("Alignment %d right junction seed has %d hits",
												aid, n_hits);
								}
						}

					if (clip != NULL)
						clip_push (&right_clips, &n_right_clips, &right_clips_alloc, clip);
				}
		}

	// The last or the first retrocopy
	if (rid_prev)
		insertion_set_entry (insertion_stmt, rid_prev,
				hash_lookup (rtc_h, &rid_prev), exons, n_exons,
				strand, left_a, left_size, right_a, right_size,
				left_clips, n_left_clips, right_clips, n_right_clips);

	clip_clear (left_clips, &n_left_clips);
	clip_clear (right_clips, &n_right_clips);

	xfree (exons);
	xfree (left_a);
	xfree (right_a);
	xfree (left_clips);
	xfree (right_clips);

	if (exon_stmt != NULL)
		db_finalize (exon_stmt);

	db_finalize (clip_stmt);
	xfree (cigar_buf);
	xfree (seq_buf);

	if (gene_h != NULL)
		hash_free (gene_h);

	if (fidx != NULL)
		hash_free (fidx);
}

/*
 * The candidate parental genes of a retrocopy: Their
 * strand and the number of fragments explained by each
 * one alone - the mate overlaps its exons, but no exon
 * of the other candidates
 */
static sqlite3_stmt *
prepare_candidate_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"WITH\n"
		"	candidate (gene_name) AS (\n"
		"		SELECT DISTINCT cl.gene_name\n"
		"		FROM cluster_merging AS cm\n"
		"		INNER JOIN cluster AS cl\n"
		"			ON cl.id = cm.cluster_id AND cl.sid = cm.cluster_sid\n"
		"		WHERE cm.retrocopy_id = $RID\n"
		"	),\n"
		"	exclusive (gene_name, n) AS (\n"
		"		SELECT gene_name, COUNT(*)\n"
		"		FROM (\n"
		"			SELECT MIN(e.gene_name) AS gene_name\n"
		"			FROM cluster_merging AS cm\n"
		"			INNER JOIN clustering AS c\n"
		"				ON c.cluster_id = cm.cluster_id\n"
		"					AND c.cluster_sid = cm.cluster_sid\n"
		"			INNER JOIN alignment AS a\n"
		"				ON a.id = c.alignment_id\n"
		"			INNER JOIN alignment AS m\n"
		"				ON m.qname = a.qname\n"
		"					AND m.source_id = a.source_id\n"
		"			INNER JOIN overlapping AS o\n"
		"				ON o.alignment_id = m.id\n"
		"			INNER JOIN exon AS e\n"
		"				ON e.id = o.exon_id\n"
		"			WHERE cm.retrocopy_id = $RID\n"
		"				AND m.id != a.id\n"
		"				AND e.gene_name IN (SELECT gene_name FROM candidate)\n"
		"			GROUP BY a.qname, a.source_id\n"
		"			HAVING COUNT(DISTINCT e.gene_name) = 1\n"
		"		)\n"
		"		GROUP BY gene_name\n"
		"	)\n"
		"SELECT c.gene_name,\n"
		"	(SELECT strand FROM exon WHERE gene_name = c.gene_name LIMIT 1),\n"
		"	COALESCE(x.n, 0)\n"
		"FROM candidate AS c\n"
		"LEFT JOIN exclusive AS x\n"
		"	USING (gene_name)\n"
		"ORDER BY c.gene_name ASC";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

/*
 * A retrocopy whose clusters come from overlapped parental
 * genes is ambiguous. But the evidence may point to one of
 * the genes:
 *
 * 1) Fragments explained by one gene alone, when no other
 *    gene has them
 * 2) Strand: The poly-A gives the strand of the retrocopy
 *    at the reference, and the orientation gives it in
 *    relation to the parental gene. For genes at opposite
 *    strands, they agree for only one of them
 *
 * The gene is chosen if at least one signal favours it and
 * none favours another. Then the clusters of the other genes
 * are left out of the retrocopy, so the junctions, the ends
 * and the predicted sequence come from the chosen one
 */
static void
choose_parental (sqlite3 *db, Hash *rtc_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *rid_stmt = NULL;
	sqlite3_stmt *candidate_stmt = NULL;
	sqlite3_stmt *delete_stmt = NULL;

	RetrocopyEntry *e = NULL;
	String *others = NULL;

	char **genes = NULL;
	char *strands = NULL;
	int *exclusive = NULL;
	int alloc = 0;
	int n = 0;

	int rid = 0;
	int by_exclusive = 0;
	int by_strand = 0;
	int chosen = 0;
	int acm = 0;
	int i = 0;
	char implied = 0;

	rid_stmt = db_prepare (db,
			"SELECT DISTINCT retrocopy_id\n"
			"FROM cluster_merging\n"
			"ORDER BY retrocopy_id ASC");

	candidate_stmt = prepare_candidate_stmt (db);

	delete_stmt = db_prepare (db,
			"DELETE FROM cluster_merging\n"
			"WHERE retrocopy_id = ?1\n"
			"	AND EXISTS (\n"
			"		SELECT 1\n"
			"		FROM cluster AS cl\n"
			"		WHERE cl.id = cluster_id\n"
			"			AND cl.sid = cluster_sid\n"
			"			AND cl.gene_name != ?2)");

	while (db_step (rid_stmt) == SQLITE_ROW)
		{
			rid = db_column_int (rid_stmt, 0);

			e = hash_lookup (rtc_h, &rid);
			assert (e != NULL);

			// Only overlapped parental genes - not along
			// with near genes or a hotspot
			if (e->level != RETROCOPY_OVERLAPPED_PARENTALS)
				continue;

			db_reset (candidate_stmt);
			db_clear_bindings (candidate_stmt);
			db_bind_int (candidate_stmt,
					sqlite3_bind_parameter_index (candidate_stmt, "$RID"),
					rid);

			n = 0;

			while (db_step (candidate_stmt) == SQLITE_ROW)
				{
					if (n >= alloc)
						{
							alloc += BLOCK_SIZE;
							genes = xrealloc (genes, sizeof (char *) * alloc);
							strands = xrealloc (strands, sizeof (char) * alloc);
							exclusive = xrealloc (exclusive, sizeof (int) * alloc);
						}

					genes[n] = xstrdup (db_column_text (candidate_stmt, 0));
					strands[n] = db_column_text (candidate_stmt, 1)[0];
					exclusive[n] = db_column_int (candidate_stmt, 2);
					n++;
				}

			// 1) The only gene with exclusive fragments
			by_exclusive = -1;
			acm = 0;

			for (i = 0; i < n; i++)
				{
					if (exclusive[i] > 0)
						acm++;

					if (exclusive[i] >= RETROCOPY_MIN_EXCLUSIVE)
						by_exclusive = i;
				}

			if (acm != 1)
				by_exclusive = -1;

			// 2) The only gene whose strand makes the
			// orientation agree with the poly-A
			by_strand = -1;
			acm = 0;

			if (e->polya_sr > 0 && e->orientation_rho != 0.0
					&& e->orientation_p_value <= RETROCOPY_ORIENTATION_ALPHA)
				{
					for (i = 0; i < n; i++)
						{
							implied = e->orientation_rho >= 0.0
								? strands[i]
								: strands[i] == '+' ? '-' : '+';

							if (implied == e->polya_strand)
								{
									by_strand = i;
									acm++;
								}
						}
				}

			if (acm != 1)
				by_strand = -1;

			// At least one signal, and no disagreement
			if (by_exclusive < 0 && by_strand < 0)
				chosen = -1;
			else if (by_exclusive >= 0 && by_strand >= 0
					&& by_exclusive != by_strand)
				chosen = -1;
			else
				chosen = by_exclusive >= 0 ? by_exclusive : by_strand;

			if (chosen >= 0)
				{
					others = string_new (NULL);

					for (i = 0; i < n; i++)
						if (i != chosen)
							string_concat_printf (others, "%s%s",
									others->len ? "/" : "", genes[i]);

					log_debug ("Retrocopy %d: Parental gene '%s' chosen over '%s' "
							"(exclusive fragments: %d, by strand: %s)",
							rid, genes[chosen], others->str, exclusive[chosen],
							by_strand >= 0 ? "yes" : "no");

					e->level |= RETROCOPY_PARENTAL_CHOSEN;
					e->other_genes = string_free (others, 0);

					db_reset (delete_stmt);
					db_clear_bindings (delete_stmt);
					db_bind_int (delete_stmt, 1, rid);
					db_bind_text (delete_stmt, 2, genes[chosen]);
					db_step (delete_stmt);
				}
			else
				log_debug ("Retrocopy %d: No parental gene chosen among %d "
						"overlapped ones", rid, n);

			for (i = 0; i < n; i++)
				xfree (genes[i]);
		}

	xfree (genes);
	xfree (strands);
	xfree (exclusive);

	db_finalize (rid_stmt);
	db_finalize (candidate_stmt);
	db_finalize (delete_stmt);
}

/*
 * A retrocopy with a single breakpoint: One flank of a
 * deletion at the target site, if another retrocopy of
 * the same parental gene has the other flank
 */
struct _TsiteSide
{
	int   rid;
	char *chr;
	char *gene;
	long  breakpoint;
	int   left;

	// The only compatible retrocopy at the other side,
	// and how many there are
	int   partner;
	int   n_partners;

	// Fragments of each sample (source id)
	int  *sources;
	int  *fragments;
	int   n_sources;
	int   alloc;
};

typedef struct _TsiteSide TsiteSide;

static void
tsite_side_free (TsiteSide *t)
{
	xfree (t->chr);
	xfree (t->gene);
	xfree (t->sources);
	xfree (t->fragments);
	xfree (t);
}

static int
tsite_side_fragments (const TsiteSide *t, const int source_id)
{
	int i = 0;

	for (i = 0; i < t->n_sources; i++)
		if (t->sources[i] == source_id)
			return t->fragments[i];

	return 0;
}

/*
 * The samples must agree: At least one of them with
 * fragments at both flanks, and none with many fragments
 * at one flank and nothing at the other
 */
static int
tsite_side_samples_agree (const TsiteSide *a, const TsiteSide *b)
{
	int shared = 0;
	int n = 0;
	int i = 0;

	for (i = 0; i < a->n_sources; i++)
		{
			n = tsite_side_fragments (b, a->sources[i]);

			if (n > 0)
				shared = 1;
			else if (a->fragments[i] >= RETROCOPY_TSITE_MIN_CONFLICT)
				return 0;
		}

	for (i = 0; i < b->n_sources; i++)
		if (b->fragments[i] >= RETROCOPY_TSITE_MIN_CONFLICT
				&& tsite_side_fragments (a, b->sources[i]) == 0)
			return 0;

	return shared;
}

static int
tsite_side_compatible (const TsiteSide *l, const TsiteSide *r,
		Hash *rtc_h, const long max_tsd_del)
{
	const RetrocopyEntry *el = NULL;
	const RetrocopyEntry *er = NULL;

	// The left flank comes first, and the deletion
	// has at least one base
	if (strcmp (l->chr, r->chr) || strcmp (l->gene, r->gene)
			|| l->breakpoint + 1 >= r->breakpoint
			|| r->breakpoint - l->breakpoint - 1 > max_tsd_del)
		return 0;

	el = hash_lookup (rtc_h, &l->rid);
	er = hash_lookup (rtc_h, &r->rid);

	// The poly-A is found at the left flank for '-' and
	// at the right one for '+': One retrocopy per strand
	if (el->polya_sr > 0 && er->polya_sr > 0
			&& el->polya_strand != er->polya_strand)
		return 0;

	// The orientation of the reads must agree
	if (el->orientation_p_value <= RETROCOPY_ORIENTATION_ALPHA
			&& er->orientation_p_value <= RETROCOPY_ORIENTATION_ALPHA
			&& (el->orientation_rho >= 0.0) != (er->orientation_rho >= 0.0))
		return 0;

	return tsite_side_samples_agree (l, r);
}

static void
tsite_side_pair (TsiteSide *t, const int rid)
{
	t->partner = rid;
	t->n_partners++;
}

static int
merge_tsite_deletion (sqlite3 *db, Hash *rtc_h, const long max_tsd_del)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	sqlite3_stmt *update_merging_stmt = NULL;
	sqlite3_stmt *update_evidence_stmt = NULL;

	RetrocopyEntry *e = NULL;
	RetrocopyEntry *er = NULL;
	TsiteSide *t = NULL;
	TsiteSide *l = NULL;
	TsiteSide *r = NULL;

	// The one-sided retrocopies: ID => side
	Hash *side_h = NULL;
	TsiteSide **sides = NULL;
	int n_sides = 0;
	int alloc = 0;

	int rid = 0;
	int merged = 0;
	int i = 0;
	int j = 0;

	side_h = hash_new_full (int_hash, int_equal, NULL,
			(DestroyNotify) tsite_side_free);

	/*
	* A single parental gene - not overlapped, near genes
	* or a hotspot - and a single breakpoint
	*/
	stmt = db_prepare (db,
			"SELECT m.retrocopy_id, MIN(c.chr), MIN(c.gene_name),\n"
			"	COUNT(DISTINCT c.gene_name)\n"
			"FROM cluster_merging AS m\n"
			"INNER JOIN cluster AS c\n"
			"	ON c.id = m.cluster_id AND c.sid = m.cluster_sid\n"
			"GROUP BY m.retrocopy_id\n"
			"ORDER BY m.retrocopy_id ASC");

	while (db_step (stmt) == SQLITE_ROW)
		{
			rid = db_column_int (stmt, 0);

			e = hash_lookup (rtc_h, &rid);
			assert (e != NULL);

			if (e->level != RETROCOPY_PASS
					|| db_column_int (stmt, 3) != 1
					|| (e->breakpoint_left_sr > 0) == (e->breakpoint_right_sr > 0))
				continue;

			t = xcalloc (1, sizeof (TsiteSide));

			t->rid        = rid;
			t->chr        = xstrdup (db_column_text (stmt, 1));
			t->gene       = xstrdup (db_column_text (stmt, 2));
			t->left       = e->breakpoint_left_sr > 0;
			t->breakpoint = t->left ? e->breakpoint_left : e->breakpoint_right;

			hash_insert (side_h, &t->rid, t);

			if (n_sides >= alloc)
				{
					alloc += BLOCK_SIZE;
					sides = xrealloc (sides, sizeof (TsiteSide *) * alloc);
				}

			sides[n_sides++] = t;
		}

	db_finalize (stmt);

	if (n_sides < 2)
		goto Exit;

	// The fragments of each sample
	stmt = db_prepare (db,
			"SELECT m.retrocopy_id, a.source_id, COUNT(DISTINCT a.qname)\n"
			"FROM cluster_merging AS m\n"
			"INNER JOIN clustering AS c\n"
			"	USING (cluster_id, cluster_sid)\n"
			"INNER JOIN alignment AS a\n"
			"	ON a.id = c.alignment_id\n"
			"GROUP BY m.retrocopy_id, a.source_id");

	while (db_step (stmt) == SQLITE_ROW)
		{
			rid = db_column_int (stmt, 0);

			t = hash_lookup (side_h, &rid);
			if (t == NULL)
				continue;

			if (t->n_sources >= t->alloc)
				{
					t->alloc += BLOCK_SIZE;
					t->sources = xrealloc (t->sources, sizeof (int) * t->alloc);
					t->fragments = xrealloc (t->fragments, sizeof (int) * t->alloc);
				}

			t->sources[t->n_sources] = db_column_int (stmt, 1);
			t->fragments[t->n_sources] = db_column_int (stmt, 2);
			t->n_sources++;
		}

	db_finalize (stmt);

	// Pair the left flanks with the right ones
	for (i = 0; i < n_sides; i++)
		{
			if (!sides[i]->left)
				continue;

			for (j = 0; j < n_sides; j++)
				{
					if (sides[j]->left
							|| !tsite_side_compatible (sides[i], sides[j],
								rtc_h, max_tsd_del))
						continue;

					tsite_side_pair (sides[i], sides[j]->rid);
					tsite_side_pair (sides[j], sides[i]->rid);
				}
		}

	update_merging_stmt = db_prepare (db,
			"UPDATE cluster_merging SET retrocopy_id = ?1\n"
			"WHERE retrocopy_id = ?2");

	update_evidence_stmt = db_prepare (db,
			"UPDATE evidence SET retrocopy_id = ?1\n"
			"WHERE retrocopy_id = ?2");

	/*
	* Merge only when both flanks have a single compatible
	* retrocopy at the other side: With more than one, the
	* deletion cannot be told
	*/
	for (i = 0; i < n_sides; i++)
		{
			l = sides[i];

			if (!l->left || l->n_partners != 1)
				continue;

			r = hash_lookup (side_h, &l->partner);
			assert (r != NULL);

			if (r->n_partners != 1)
				continue;

			log_debug ("Retrocopies %d and %d: Flanks of a deletion of %li "
					"bases at the target site", l->rid, r->rid,
					r->breakpoint - l->breakpoint - 1);

			db_reset (update_merging_stmt);
			db_clear_bindings (update_merging_stmt);
			db_bind_int (update_merging_stmt, 1, l->rid);
			db_bind_int (update_merging_stmt, 2, r->rid);
			db_step (update_merging_stmt);

			db_reset (update_evidence_stmt);
			db_clear_bindings (update_evidence_stmt);
			db_bind_int (update_evidence_stmt, 1, l->rid);
			db_bind_int (update_evidence_stmt, 2, r->rid);
			db_step (update_evidence_stmt);

			e = hash_lookup (rtc_h, &l->rid);
			er = hash_lookup (rtc_h, &r->rid);

			e->breakpoint_right = er->breakpoint_right;
			e->breakpoint_right_sr = er->breakpoint_right_sr;

			if (e->polya_sr == 0 && er->polya_sr > 0)
				{
					e->polya_len = er->polya_len;
					e->polya_lower_bound = er->polya_lower_bound;
					e->polya_strand = er->polya_strand;
					e->polya_sr = er->polya_sr;
				}

			e->level |= RETROCOPY_TSITE_MERGED;
			hash_remove (rtc_h, &r->rid);

			merged++;
		}

	db_finalize (update_merging_stmt);
	db_finalize (update_evidence_stmt);

Exit:
	xfree (sides);
	hash_free (side_h);

	return merged;
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
	char *genes = NULL;
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

			// A parental gene chosen among overlapped
			// ones: 'chosen/others'
			if (e->other_genes != NULL)
				{
					xasprintf (&genes, "%s/%s", gene, e->other_genes);
					gene = genes;
				}

			/*
			* The insertion point is the first reference base
			* after the insertion. Prefer the right breakpoint,
			* so the insertion (and its TSD) is left aligned.
			* With a deletion at the target site, it is the first
			* deleted base: The left-most breakpoint, as VCF asks
			* for a structural variant. Without split reads, fall
			* back to the window mean
			*/
			if (e->breakpoint_left_sr > 0 && e->breakpoint_right_sr > 0
					&& e->breakpoint_left + 1 < e->breakpoint_right)
				{
					ip = e->breakpoint_left + 1;
					ip_type = RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE;
				}
			else if (e->breakpoint_right_sr > 0)
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
					e->breakpoint_right_sr, e->polya_len,
					e->polya_lower_bound,
					e->polya_strand == '-' ? "-" : "+",
					e->polya_sr);

			xfree (genes);
			genes = NULL;
		}

	db_finalize (cluster_merging_query_stmt);
}

void
retrocopy (sqlite3_stmt *retrocopy_stmt,
		sqlite3_stmt *cluster_merging_stmt,
		sqlite3_stmt *junction_stmt,
		sqlite3_stmt *insertion_stmt,
		sqlite3_stmt *evidence_stmt,
		const char *fasta_file,
		int near_gene_dist,
		long max_tsd_del)
{
	log_trace ("Inside %s", __func__);
	assert (retrocopy_stmt != NULL
			&& cluster_merging_stmt != NULL
			&& junction_stmt != NULL
			&& insertion_stmt != NULL
			&& evidence_stmt != NULL
			&& near_gene_dist > 0
			&& max_tsd_del >= 0);

	// Keep ID => level
	Hash *rtc_h = NULL;

	const int filter =
		CLUSTER_FILTER_NONE
		|CLUSTER_FILTER_CHR
		|CLUSTER_FILTER_DIST
		|CLUSTER_FILTER_REGION
		|CLUSTER_FILTER_SUPPORT;

	rtc_h = hash_new_full (int_hash, int_equal,
			xfree, (DestroyNotify) rtc_entry_free);

	log_debug ("Clean retrocopy tables");
	clean_retrocopy_tables (
			sqlite3_db_handle (retrocopy_stmt));

	log_info ("Merge the clusters into retrocopies");
	merge_cluster (cluster_merging_stmt, filter, near_gene_dist, rtc_h);

	log_info ("Calculate the orientation of the retrocopies");
	calculate_orientation (
			sqlite3_db_handle (retrocopy_stmt),
			rtc_h);

	log_info ("Find the breakpoints and target site duplications");
	calculate_breakpoint (
			sqlite3_db_handle (retrocopy_stmt),
			evidence_stmt, rtc_h);

	log_info ("Find the poly-A tails");
	calculate_polya (
			sqlite3_db_handle (retrocopy_stmt),
			evidence_stmt, rtc_h);

	log_info ("Choose the parental gene among overlapped ones");
	choose_parental (
			sqlite3_db_handle (retrocopy_stmt),
			rtc_h);

	/*
	* A deletion at the target site longer than the cluster
	* window splits the insertion into two retrocopies, one
	* at each flank, with a single breakpoint each. Merge
	* them, and recalculate the orientation of the merged
	* ones with the reads of both flanks
	*/
	log_info ("Merge the flanks of deletions at the target site");
	if (merge_tsite_deletion (sqlite3_db_handle (retrocopy_stmt),
				rtc_h, max_tsd_del) > 0)
		calculate_orientation (
				sqlite3_db_handle (retrocopy_stmt),
				rtc_h);

	log_info ("Annotate the insertion points");
	annotate_retrocopy (retrocopy_stmt, rtc_h);

	log_info ("Search for reads spanning exon-exon junctions");
	calculate_junction (junction_stmt, evidence_stmt);

	if (fasta_file != NULL)
		log_info ("Search for the retrocopy ends in the parental genes");
	else
		log_info ("No reference genome (--reference-file): Skip the retrocopy ends and predicted sequences");

	calculate_insertion (
			sqlite3_db_handle (retrocopy_stmt),
			insertion_stmt, fasta_file, rtc_h);

	hash_free (rtc_h);
}
