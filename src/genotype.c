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
};

typedef struct _Genotype Genotype;

struct _ZygosityData
{
	sqlite3_stmt *stmt;

	List         *genotype;
	ChrStd       *cs;

	int           phred_quality;

	const char   *path;
};

typedef struct _ZygosityData ZygosityData;

struct _CrossWindowLinear
{
	bam1_t       *align;
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
		return;

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
		.normal_scores   = array_new (xfree)
	};

	return g;
}

static void
genotype_free (Genotype *g)
{
	if (g == NULL)
		return;

	array_free (g->abnormal_scores, 1);
	array_free (g->normal_scores, 1);

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
		return;

	list_free (zd->genotype);
	xfree (zd);
}

static void
clean_genotype_table (sqlite3 *db)
{
	// Delete all values from
	// previous runs
	const char sql[] =
		"DELETE FROM genotype";

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
	return pow (10.0, -1.0 * mapq / 10);
}

static double
likelihood_HE (const int len, const int ploidy)
{
	return log10 (pow (1.0 / ploidy, len));
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
	const char *path_copy = NULL;

	int source_id = 0;
	const char *path = NULL;

	const char sql[] =
		"SELECT id, path\n"
		"FROM source";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	// Prepare mapq query for normal_scores
	mapq_stmt = prepare_alignment_score_query_stmt (db);

	zi = hash_new_full (str_hash, str_equal, xfree,
			(DestroyNotify) zygosity_data_free);

	while (db_step (stmt) == SQLITE_ROW)
		{
			source_id = db_column_int (stmt, 0);
			path = db_column_text (stmt, 1);

			log_debug ("Index source path [%d] %s",
					source_id, path);

			zd = hash_lookup (zi, path);
			assert (zd == NULL);

			path_copy = xstrdup (path);

			zd = zygosity_data_new ();
			zd->path = path_copy;

			hash_insert (zi, path_copy, zd);

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
chr_std2tid (ChrStd *cs, const bam_hdr_t *hdr)
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

	calculate_align_start_end (align, &start, &end);

	return insertion_point >= start && insertion_point <= end;
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
					log_warn ("No %s contig from retrocopy [%d] found in BAM header",
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

	CrossWindowLinear c = {.phred_quality = zd->phred_quality};

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

			// Query position CHR:START-END
			// START is 0-based
			// END is 1-based
			itr = sam_itr_queryi (idx, *tid,
					r->window_start - 1, r->window_end);

			if (itr == NULL)
				log_fatal ("Failed to look for '%s:%li-%li' at %s\n",
						r->chr, r->window_start, r->window_end, zd->path);

			while ((rc = sam_itr_next (fp, itr, align)) >= 0)
				{
					if (cross_insertion_point (align, r->insertion_point, zd->phred_quality))
						{
							q = xcalloc (1, sizeof (int));
							*q = align->core.qual;
							array_add (g->normal_scores, q);
						}
				}

			if (rc < -1)
				log_fatal ("Failed to read sam alignment");

			dump_genotype (zd->stmt, g);

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

	// Get standardized tid
	chr_tid = chr_std2tid (zd->cs, hdr);

	// Look for the index
	idx = sam_index_load (fp, zd->path);

	if (idx == NULL)
		{
			log_warn ("Failed to open BAM INDEX for '%s'. Make a linear search", zd->path);
			zygosity_linear_search (fp, hdr, align, zd, chr_tid);
		}
	else
		{
			log_info ("Open BAM INDEX for '%s'. Make an indexed search", zd->path);
			zygosity_indexed_search (fp, hdr, align, idx, zd, chr_tid);
		}

	if (sam_close (fp) < 0)
		log_errno_fatal ("Failed to close '%s'", zd->path);

	hash_free (chr_tid);
	hts_idx_destroy (idx);
	bam_hdr_destroy (hdr);
	bam_destroy1 (align);
}

void
genotype (sqlite3_stmt *genotype_stmt, int threads, int phred_quality)
{
	log_trace ("Inside %s", __func__);
	assert (genotype_stmt != NULL
			&& threads > 0
			&& phred_quality >= 0);

	sqlite3 *db = NULL;
	threadpool thpool = NULL;

	ChrStd *cs = NULL;
	Hash *retrocopy_h, *zygosity_h;
	retrocopy_h = zygosity_h = NULL;

	HashIter itr = {};
	const char *path = NULL;
	ZygosityData *zd = NULL;

	// Get DB handle
	db = sqlite3_db_handle (genotype_stmt);

	// Necessary for get the right chr tid
	cs = chr_std_new ();

	// Alloc n threads into the pool
	thpool = thpool_init (threads);

	log_debug ("Clean genotype table");
	clean_genotype_table (db);

	log_info ("Index all retrocopies");

	// RETROCOPY_ID => REGION
	retrocopy_h = genotype_index_retrocopy (db);

	log_info ("Index all BAM path => retrocopy relationship");

	// PATH => @ZYGOSITYDATA
	zygosity_h = genotype_index_zygosity_data (db, retrocopy_h);

	// Iterator through paths
	hash_iter_init (&itr, zygosity_h);

	while (hash_iter_next (&itr, (void **) &path, (void **) &zd))
		{
			log_debug ("Look for retrocopies zygosity of file '%s'", path);

			// Don't forget to give chromosome standardization and
			// genotype statement
			zd->cs = cs;
			zd->stmt = genotype_stmt;
			zd->phred_quality = phred_quality;

			// Let's rock!
			thpool_add_work (thpool, (void *) zygosity, (void *) zd);
		}

	// Wait all threads to return
	thpool_wait (thpool);

	// Clean up
	thpool_destroy (thpool);
	chr_std_free (cs);
	hash_free (retrocopy_h);
	hash_free (zygosity_h);
}
