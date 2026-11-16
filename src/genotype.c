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

#include <htslib/sam.h>
#include <math.h>
#include <assert.h>
#include "wrapper.h"
#include "db.h"
#include "log.h"
#include "chr.h"
#include "ibitree.h"
#include "hash.h"
#include "array.h"
#include "list.h"
#include "thpool.h"
#include "abnormal.h"
#include "intron.h"
#include "genotype.h"

struct _Region
{
	const char *chr;

	long        window_start;
	long        window_end;

	long        insertion_point;
};

typedef struct _Region Region;

struct _Genotype
{
	Region     *region;

	int         retrocopy_id;
	int         source_id;

	int         ploidy;

	Array      *abnormal_scores;
	Array      *normal_scores;

	// INTRON => number of reads with
	// the mate inside it
	Hash       *introns;
};

typedef struct _Genotype Genotype;

struct _ZygosityData
{
	sqlite3_stmt *stmt;
	sqlite3_stmt *intron_stmt;

	// Introns of the parental genes. It is shared
	// among the threads, only for reading
	IntronTree   *intron_tree;

	// TID => standardized contig name
	const char  **tid_chr;
	int           n_targets;

	List         *genotype;
	ChrStd       *cs;

	int           phred_quality;

	const char   *path;
};

typedef struct _ZygosityData ZygosityData;

struct _CrossWindowLinear
{
	bam1_t       *align;
	ZygosityData *zd;
	int           phred_quality;
};

typedef struct _CrossWindowLinear CrossWindowLinear;

static Region *
region_new (const char *chr, const long window_start,
		const long window_end, const long insertion_point)
{
	Region *r = xcalloc (1, sizeof (Region));

	*r = (Region) {
		.chr             = xstrdup (chr),
		.window_start    = window_start,
		.window_end      = window_end,
		.insertion_point = insertion_point
	};

	return r;
}

static void
region_free (Region *r)
{
	if (r == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	xfree ((void *) r->chr);
	xfree (r);
}

static Genotype *
genotype_new (const int retrocopy_id, const int source_id,
		const Region *r, int ploidy)
{
	Genotype *g = xcalloc (1, sizeof (Genotype));

	*g = (Genotype) {
		.region          = (Region *) r,
		.ploidy          = ploidy,
		.retrocopy_id    = retrocopy_id,
		.source_id       = source_id,
		.abnormal_scores = array_new (xfree),
		.normal_scores   = array_new (xfree),
		.introns         = hash_new_full (direct_hash, direct_equal,
				NULL, xfree)
	};

	return g;
}

static void
genotype_free (Genotype *g)
{
	if (g == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	array_free (g->abnormal_scores, 1);
	array_free (g->normal_scores, 1);
	hash_free (g->introns);

	xfree (g);
}

static ZygosityData *
zygosity_data_new (void)
{
	ZygosityData *zd = xcalloc (1, sizeof (ZygosityData));
	zd->genotype = list_new ((DestroyNotify) genotype_free);

	return zd;
}

static void
zygosity_data_free (ZygosityData *zd)
{
	if (zd == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	list_free (zd->genotype);
	xfree ((char *) zd->path);
	xfree (zd);
}

static void
clean_genotype_table (sqlite3 *db)
{
	// Delete all values from
	// previous runs
	const char sql[] =
		"DELETE FROM genotype;\n"
		"DELETE FROM intron_retention";

	log_debug ("Clean tables:\n%s", sql);
	db_exec (db, sql);
}

static Hash *
genotype_index_retrocopy (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	Hash *ri = NULL;

	Region *r = NULL;
	int *id = NULL;

	// Database search
	int retrocopy_id = 0;
	const char *chr = NULL;
	long window_start = 0;
	long window_end = 0;
	long insertion_point = 0;

	const char sql[] =
		"SELECT id, chr, window_start, window_end, insertion_point\n"
		"FROM retrocopy";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	ri = hash_new_full (int_hash, int_equal, xfree,
			(DestroyNotify) region_free);

	while (db_step (stmt) == SQLITE_ROW)
		{
			retrocopy_id    = db_column_int   (stmt, 0);
			chr             = db_column_text  (stmt, 1);
			window_start    = db_column_int64 (stmt, 2);
			window_end      = db_column_int64 (stmt, 3);
			insertion_point = db_column_int64 (stmt, 4);

			r = region_new (chr, window_start, window_end,
					insertion_point);

			id = xcalloc (1, sizeof (int));
			*id = retrocopy_id;

			log_debug ("Index retrocopy region [%d] %s:%li-%li in %li",
					retrocopy_id, chr, window_start, window_end,
					insertion_point);

			hash_insert (ri, id, r);
		}

	db_finalize (stmt);

	return ri;
}

static sqlite3_stmt *
prepare_alignment_score_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"SELECT a.mapq\n"
		"FROM retrocopy AS r\n"
		"INNER JOIN cluster_merging AS cm\n"
		"	ON r.id = cm.retrocopy_id\n"
		"INNER JOIN clustering AS c\n"
		"	USING (cluster_id, cluster_sid)\n"
		"INNER JOIN alignment AS a\n"
		"	ON c.alignment_id = a.id\n"
		"WHERE r.id = $RID\n"
		"	AND a.source_id = $SID";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static inline double
dephred_score (const int mapq)
{
	/*
	* MAPQ 0 would give an error probability of 1 and
	* log10 (0) for the reads supporting the genotype:
	* Take it as the least informative score
	*/
	return pow (10.0, -1.0 * (mapq > 0 ? mapq : 1) / 10);
}

static double
likelihood_HE (const int len, const int ploidy)
{
	// log10 ((1 / ploidy) ^ len), without the power: It
	// underflows to zero (log10 => -inf) at a deep site
	return len * log10 (1.0 / ploidy);
}

static double
likelihood_HO (const Array *a1, const Array *a2,
		const int ploidy)
{
	double l = 0.0;
	const int *q = NULL;
	int i = 0;

	for (i = 0; i < array_len (a1); i++)
		{
			q = array_get (a1, i);
			l += log10 (1.0 * ploidy * dephred_score (*q));
		}

	for (i = 0; i < array_len (a2); i++)
		{
			q = array_get (a2, i);
			l += log10 (ploidy * (1.0 - dephred_score (*q)));
		}

	return likelihood_HE (array_len (a1) +
			array_len (a2), ploidy) + l;
}

static void
genotype_get_abnormal_scores (sqlite3_stmt *stmt, const int retrocopy_id,
		const int source_id, Array *a)
{
	int mapq = 0;
	int *q = NULL;

	db_reset (stmt);
	db_clear_bindings (stmt);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$RID"),
			retrocopy_id);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$SID"),
			source_id);

	while (db_step (stmt) == SQLITE_ROW)
		{
			mapq = db_column_int (stmt, 0);

			q = xcalloc (1, sizeof (int));
			*q = mapq;

			array_add (a, q);
		}
}

static Hash *
genotype_index_zygosity_data (sqlite3 *db, Hash *retrocopy_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	sqlite3_stmt *mapq_stmt = NULL;

	Hash *zi = NULL;
	ZygosityData *zd = NULL;
	Genotype *g = NULL;
	Region *r = NULL;
	List *gl = NULL;

	HashIter iter = {};

	const int *retrocopy_id = NULL;
	int *source_id_copy = NULL;

	int source_id = 0;
	const char *path = NULL;

	const char sql[] =
		"SELECT id, path\n"
		"FROM source";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	// Prepare mapq query for normal_scores
	mapq_stmt = prepare_alignment_score_query_stmt (db);

	// SOURCE_ID => @ZYGOSITYDATA. Not by path: Two sources
	// may come from the same file (e.g. merged databases)
	zi = hash_new_full (int_hash, int_equal, xfree,
			(DestroyNotify) zygosity_data_free);

	while (db_step (stmt) == SQLITE_ROW)
		{
			source_id = db_column_int (stmt, 0);
			path = db_column_text (stmt, 1);

			log_debug ("Index source path [%d] %s",
					source_id, path);

			source_id_copy = xcalloc (1, sizeof (int));
			*source_id_copy = source_id;

			zd = zygosity_data_new ();
			zd->path = xstrdup (path);

			hash_insert (zi, source_id_copy, zd);

			gl = zd->genotype;
			hash_iter_init (&iter, retrocopy_h);

			// All paths have a list of all retrocopies
			while (hash_iter_next (&iter, (void **) &retrocopy_id, (void **) &r))
				{
					g = genotype_new (*retrocopy_id, source_id, r, 2);

					genotype_get_abnormal_scores (mapq_stmt, *retrocopy_id,
							source_id, g->abnormal_scores);

					list_append (gl, g);
				}
		}

	db_finalize (stmt);
	db_finalize (mapq_stmt);

	return zi;
}

static Hash *
chr_std2tid (ChrStd *cs, const bam_hdr_t *hdr, const char **tid_chr)
{
	log_trace ("Inside %s", __func__);

	Hash *chr_tid = NULL;
	const char *chr_std = NULL;
	int *tid_copy = NULL;
	int tid = 0;

	chr_tid = hash_new_full (str_hash, str_equal,
			NULL, xfree);

	for (; tid < hdr->n_targets; tid++)
		{
			chr_std = chr_std_lookup (cs, hdr->target_name[tid]);
			tid_chr[tid] = chr_std;

			tid_copy = xcalloc (1, sizeof (int));
			*tid_copy = tid;

			hash_insert (chr_tid, chr_std, tid_copy);
		}

	return chr_tid;
}

static inline void
calculate_align_start_end (const bam1_t *align, long *start, long *end)
{
	uint32_t *cigar = bam_get_cigar (align);
	int rlen = bam_cigar2rlen (align->core.n_cigar, cigar);

	*start = align->core.pos + 1;

	*end = rlen < 1
		? *start
		: *start + rlen - 1;
}

static inline int
cross_insertion_point (const bam1_t *align, const long insertion_point,
		const int phred_quality)
{
	/*
	* Only test with proper aligned reads
	* - paired-end
	* - proper-pair
	* - mapped
	* - mate mapped
	* - not a duplication
	* - not a supplementary
	* - phred quality
	*/
	if (!(align->core.flag & 0x1)
			|| !(align->core.flag & 0x2)
			|| (align->core.flag & 0x4)
			|| (align->core.flag & 0x8)
			|| (align->core.flag & 0x400)
			|| (align->core.flag & 0x800)
			|| (align->core.qual < phred_quality))
		{
			return 0;
		}

	long start = 0;
	long end = 0;
	int half_decil = 0;

	calculate_align_start_end (align, &start, &end);
	half_decil = (end - start + 1) / 20;

	// The read needs to cover at least the window comprising
	// the insertion_point +/- read half_decil
	return (insertion_point - half_decil) >= start
		&& (insertion_point + half_decil) <= end;
}

static Hash *
index_region (const List *genotype, Hash *chr_tid)
{
	log_trace ("Inside %s", __func__);

	Hash *ir = NULL;

	IBiTree *tree = NULL;
	ListElmt *cur = NULL;

	Genotype *g = NULL;
	Region *r = NULL;
	const int *tid = NULL;

	// Init TID => TREE
	ir = hash_new_full (int_hash, int_equal,
			NULL, (DestroyNotify) ibitree_free);

	cur = list_head (genotype);
	for (; cur != NULL; cur = list_next (cur))
		{
			g = list_data (cur);
			r = g->region;

			// Get standardized tid
			tid = hash_lookup (chr_tid, r->chr);
			if (tid == NULL)
				{
					log_warn ("Contig '%s' of retrocopy %d not found at the alignment header",
							r->chr, g->retrocopy_id);
					continue;
				}

			tree = hash_lookup (ir, tid);
			if (tree == NULL)
				{
					tree = ibitree_new (NULL);
					hash_insert (ir, tid, tree);
				}

			ibitree_insert (tree, r->window_start,
					r->window_end, g);
		}

	return ir;
}

static void
dump_genotype (sqlite3_stmt *stmt, const Genotype *g)
{
	double ho_ref, he, ho_alt;
	ho_ref = he = ho_alt = 0.0;

	ho_ref = likelihood_HO (g->abnormal_scores, g->normal_scores, g->ploidy);
	ho_alt = likelihood_HO (g->normal_scores, g->abnormal_scores, g->ploidy);
	he = likelihood_HE (array_len (g->normal_scores) +
			array_len (g->abnormal_scores), g->ploidy);

	log_debug ("retrocopy [%d %d] %.2f,%.2f,%.2f",
			g->retrocopy_id, g->source_id, ho_ref, he, ho_alt);

	db_insert_genotype (stmt, g->source_id, g->retrocopy_id, array_len (g->normal_scores),
			array_len (g->abnormal_scores), ho_ref, he, ho_alt);
}

static void
dump_introns (sqlite3_stmt *stmt, const Genotype *g)
{
	HashIter iter = {};
	const Intron *intron = NULL;
	const int *reads = NULL;

	hash_iter_init (&iter, g->introns);

	while (hash_iter_next (&iter, (void **) &intron, (void **) &reads))
		intron_dump (stmt, intron, g->retrocopy_id, g->source_id, *reads);
}

/*
 * Length of the mate at the reference, from the SAM tag
 * 'MC' (mate CIGAR). Return 0 if the tag is absent, or it
 * cannot be read
 */
static long
mate_rlen (const bam1_t *align)
{
	const uint8_t *aux = NULL;
	const char *mc = NULL;
	long rlen = 0;
	long n = 0;

	aux = bam_aux_get (align, "MC");
	if (aux == NULL)
		return 0;

	mc = bam_aux2Z (aux);
	if (mc == NULL)
		return 0;

	for (; *mc != '\0'; mc++)
		{
			if (*mc >= '0' && *mc <= '9')
				{
					n = n * 10 + (*mc - '0');
					continue;
				}

			// Operations consuming the reference
			if (*mc == 'M' || *mc == 'D' || *mc == 'N'
					|| *mc == '=' || *mc == 'X')
				rlen += n;

			n = 0;
		}

	return rlen;
}

/*
 * Intron retention: An abnormal read at the window of the
 * retrocopy whose mate is inside an intron of the parental
 * gene. The alignment carries the position of its mate, so
 * the mate itself is not read. The mapping quality of the
 * mate (tag 'MQ') and its end (tag 'MC') are used when the
 * aligner wrote them: Otherwise, the mate is taken with the
 * quality and the length of this alignment
 */
static void
count_intron_mate (const ZygosityData *zd, Genotype *g, const bam1_t *align)
{
	const Intron *intron = NULL;
	const uint8_t *aux = NULL;
	const char *chr = NULL;
	int *reads = NULL;
	long start = 0;
	long rlen = 0;

	// The same quality control of the abnormal reads. The
	// secondary and supplementary alignments have the same
	// mate as the primary one: Do not count it twice
	if ((align->core.flag & (0x100 | 0x800))
			|| !abnormal_quality_control (align, zd->phred_quality,
				ABNORMAL_DEFAULT_MAX_BASE_FREQ))
		return;

	if (align->core.mtid < 0 || align->core.mtid >= zd->n_targets)
		return; // LCOV_EXCL_LINE: A mapped mate has its contig

	aux = bam_aux_get (align, "MQ");
	if (aux != NULL && bam_aux2i (aux) < zd->phred_quality)
		return;

	rlen = mate_rlen (align);

	if (rlen < 1)
		rlen = bam_cigar2rlen (align->core.n_cigar, bam_get_cigar (align));

	if (rlen < 1)
		rlen = 1; // LCOV_EXCL_LINE: A mapped read consumes the reference

	chr = zd->tid_chr[align->core.mtid];
	start = align->core.mpos + 1;

	intron = intron_tree_lookup (zd->intron_tree, g->retrocopy_id,
			chr, start, start + rlen - 1);

	if (intron == NULL)
		return;

	reads = hash_lookup (g->introns, intron);

	if (reads == NULL)
		{
			reads = xcalloc (1, sizeof (int));
			hash_insert (g->introns, intron, reads);
		}

	(*reads)++;
}

static void
cross_window (IBiTreeLookupData *ldata, void *user_data)
{
	const CrossWindowLinear *c = user_data;
	Genotype *g = ldata->data;
	int *q = NULL;

	if (cross_insertion_point (c->align, g->region->insertion_point, c->phred_quality))
		{
			q = xcalloc (1, sizeof (int));
			*q = c->align->core.qual;
			array_add (g->normal_scores, q);
		}

	count_intron_mate (c->zd, g, c->align);
}

static void
zygosity_linear_search (samFile *fp, bam_hdr_t *hdr, bam1_t *align,
		ZygosityData *zd, Hash *chr_tid)
{
	log_trace ("Inside %s", __func__);

	// Indexed regions: TID => TREE
	Hash *ir = NULL;
	IBiTree *tree = NULL;

	Genotype *g = NULL;
	ListElmt *cur = NULL;

	CrossWindowLinear c = {.zd = zd, .phred_quality = zd->phred_quality};

	long start = 0;
	long end = 0;

	int rc = 0;

	// Index all retrocopies into a tree by chr
	ir = index_region (zd->genotype, chr_tid);

	while ((rc = sam_read1 (fp, hdr, align)) >= 0)
		{
			tree = hash_lookup (ir, &align->core.tid);
			if (tree == NULL)
				continue;

			calculate_align_start_end (align, &start, &end);
			c.align = align;

			ibitree_lookup (tree, start, end,
					-1, -1, 0, cross_window, &c);
		}

	if (rc < -1)
		log_fatal ("Failed to read sam alignment");

	cur = list_head (zd->genotype);
	for (; cur != NULL; cur = list_next (cur))
		{
			g = list_data (cur);
			dump_genotype (zd->stmt, g);
			dump_introns (zd->intron_stmt, g);
		}

	hash_free (ir);
}

static void
zygosity_indexed_search (samFile *fp, bam_hdr_t *hdr, bam1_t *align,
		hts_idx_t *idx, ZygosityData *zd, Hash *chr_tid)
{
	log_trace ("Inside %s", __func__);

	hts_itr_t *itr = NULL;

	Genotype *g = NULL;
	Region *r = NULL;

	ListElmt *cur = NULL;

	const int *tid = NULL;
	int *q = NULL;
	int rc = 0;

	cur = list_head (zd->genotype);
	for (; cur != NULL; cur = list_next (cur))
		{
			g = list_data (cur);
			r = g->region;

			// Get standardized tid
			tid = hash_lookup (chr_tid, r->chr);

			// Contig not at the alignment header: No reads
			// for the reference allele, as at the linear search
			if (tid == NULL)
				{
					log_warn ("Contig '%s' of retrocopy %d not found at the header of '%s'",
							r->chr, g->retrocopy_id, zd->path);
					dump_genotype (zd->stmt, g);
					continue;
				}

			// Query position CHR:START-END
			// START is 0-based
			// END is 1-based
			itr = sam_itr_queryi (idx, *tid,
					r->window_start - 1, r->window_end);

			if (itr == NULL)
				log_fatal ("Failed to look for '%s:%li-%li' at '%s'",
						r->chr, r->window_start, r->window_end, zd->path);

			while ((rc = sam_itr_next (fp, itr, align)) >= 0)
				{
					if (cross_insertion_point (align, r->insertion_point, zd->phred_quality))
						{
							q = xcalloc (1, sizeof (int));
							*q = align->core.qual;
							array_add (g->normal_scores, q);
						}

					count_intron_mate (zd, g, align);
				}

			if (rc < -1)
				log_fatal ("Failed to read sam alignment");

			dump_genotype (zd->stmt, g);
			dump_introns (zd->intron_stmt, g);

			sam_itr_destroy (itr);
		}
}

static void
zygosity (ZygosityData *zd)
{
	log_trace ("Inside %s", __func__);

	samFile *fp = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *align = NULL;

	hts_idx_t *idx = NULL;

	Hash *chr_tid = NULL;

	// Open BAM file for reading
	fp = sam_open (zd->path, "rb");
	if (fp == NULL)
		log_errno_fatal ("Failed to open '%s' for reading", zd->path);

	// If all OK until now, read header
	hdr = sam_hdr_read (fp);
	if (hdr == NULL)
		log_fatal ("Failed to read sam header");

	// And allocate BAM align entry
	align = bam_init1 ();
	if (align == NULL)
		log_errno_fatal ("Failed to create bam_init1");

	// Get standardized tid, and the
	// standardized contig of each tid
	zd->n_targets = hdr->n_targets;
	zd->tid_chr = xcalloc (hdr->n_targets > 0 ? hdr->n_targets : 1,
			sizeof (char *));

	chr_tid = chr_std2tid (zd->cs, hdr, zd->tid_chr);

	// Look for the index
	idx = sam_index_load (fp, zd->path);

	if (idx == NULL)
		{
			log_warn ("No index found for '%s': Read the whole file (slower)", zd->path);
			zygosity_linear_search (fp, hdr, align, zd, chr_tid);
		}
	else
		{
			log_debug ("Use the index of '%s'", zd->path);
			zygosity_indexed_search (fp, hdr, align, idx, zd, chr_tid);
		}

	if (sam_close (fp) < 0)
		log_errno_fatal ("Failed to close '%s'", zd->path);

	hash_free (chr_tid);
	xfree (zd->tid_chr);
	hts_idx_destroy (idx);
	bam_hdr_destroy (hdr);
	bam_destroy1 (align);
}

void
genotype (sqlite3_stmt *genotype_stmt, sqlite3_stmt *intron_stmt,
		int threads, int phred_quality)
{
	log_trace ("Inside %s", __func__);
	assert (genotype_stmt != NULL
			&& intron_stmt != NULL
			&& threads > 0
			&& phred_quality >= 0);

	sqlite3 *db = NULL;
	threadpool thpool = NULL;

	ChrStd *cs = NULL;
	IntronTree *intron_tree = NULL;
	Hash *retrocopy_h, *zygosity_h;
	retrocopy_h = zygosity_h = NULL;

	HashIter itr = {};
	ZygosityData *zd = NULL;

	// Get DB handle
	db = sqlite3_db_handle (genotype_stmt);

	// Necessary for get the right chr tid
	cs = chr_std_new ();

	// Alloc n threads into the pool
	thpool = thpool_init (threads);

	log_debug ("Clean genotype table");
	clean_genotype_table (db);

	log_debug ("Index the introns of the parental genes");
	intron_tree = intron_tree_new (db);

	log_info ("Genotype the retrocopies in each sample");

	// RETROCOPY_ID => REGION
	retrocopy_h = genotype_index_retrocopy (db);

	log_debug ("Index the retrocopies by alignment file");

	// SOURCE_ID => @ZYGOSITYDATA
	zygosity_h = genotype_index_zygosity_data (db, retrocopy_h);

	// Iterator through sources
	hash_iter_init (&itr, zygosity_h);

	while (hash_iter_next (&itr, NULL, (void **) &zd))
		{
			log_debug ("Look for retrocopies zygosity of file '%s'", zd->path);

			// Don't forget to give chromosome standardization and
			// genotype statement
			zd->cs = cs;
			zd->stmt = genotype_stmt;
			zd->intron_stmt = intron_stmt;
			zd->intron_tree = intron_tree;
			zd->phred_quality = phred_quality;

			// Let's rock!
			thpool_add_work (thpool, (void *) zygosity, (void *) zd);
		}

	// Wait all threads to return
	thpool_wait (thpool);

	// Clean up
	thpool_destroy (thpool);
	chr_std_free (cs);
	intron_tree_free (intron_tree);
	hash_free (retrocopy_h);
	hash_free (zygosity_h);
}
