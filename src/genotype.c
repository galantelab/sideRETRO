#include "config.h"

#include <htslib/sam.h>
#include <assert.h>
#include "wrapper.h"
#include "db.h"
#include "log.h"
#include "chr.h"
#include "ibitree.h"
#include "hash.h"
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

	int         acm;
	Zygosity    z;
};

typedef struct _Genotype Genotype;

struct _ZygosityData
{
	sqlite3_stmt *stmt;

	List         *genotype;
	ChrStd       *cs;

	int           crossing_reads;
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

static ZygosityData *
zygosity_data_new (void)
{
	ZygosityData *zd = xcalloc (1, sizeof (ZygosityData));
	zd->genotype = list_new (xfree);

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

static sqlite3_stmt *
prepare_genotype_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"WITH\n"
		"	genotype (id, sid, source_id, path) AS (\n"
		"		SELECT cluster_id, cluster_sid, source_id, path\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON c.alignment_id = a.id\n"
		"		INNER JOIN source AS s\n"
		"			ON a.source_id = s.id\n"
		"	)\n"
		"SELECT DISTINCT retrocopy_id, source_id,\n"
		"	chr, window_start, window_end,\n"
		"	insertion_point, path\n"
		"FROM retrocopy AS r\n"
		"INNER JOIN cluster_merging AS c\n"
		"	ON r.id = c.retrocopy_id\n"
		"INNER JOIN genotype AS g\n"
		"	ON c.cluster_id = g.id\n"
		"		AND c.cluster_sid = g.sid";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static void
genotype_search_db (sqlite3 *db, Hash *retrocopy_h, Hash *path_h)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *genotype_query_stmt = NULL;

	Region *r = NULL;
	Genotype *g = NULL;
	ZygosityData *zd = NULL;

	int *retrocopy_id_alloc = NULL;
	char *path_copy = NULL;

	// Database search
	int retrocopy_id = 0;
	int source_id = 0;
	const char *chr = NULL;
	long window_start = 0;
	long window_end = 0;
	long insertion_point = 0;
	const char *path = NULL;

	// Prepare query
	genotype_query_stmt = prepare_genotype_query_stmt (db);

	while (db_step (genotype_query_stmt) == SQLITE_ROW)
		{
			retrocopy_id    = db_column_int   (genotype_query_stmt, 0);
			source_id       = db_column_int   (genotype_query_stmt, 1);
			chr             = db_column_text  (genotype_query_stmt, 2);
			window_start    = db_column_int64 (genotype_query_stmt, 3);
			window_end      = db_column_int64 (genotype_query_stmt, 4);
			insertion_point = db_column_int64 (genotype_query_stmt, 5);
			path            = db_column_text  (genotype_query_stmt, 6);

			log_debug ("Index retrocopy [%d] %s:%li-%li into path %s",
					retrocopy_id, chr, window_start, window_end, path);

			r = hash_lookup (retrocopy_h, &retrocopy_id);
			if (r == NULL)
				{
					r = region_new (chr, window_start,
							window_end, insertion_point);

					retrocopy_id_alloc = xcalloc (1, sizeof (int));
					*retrocopy_id_alloc = retrocopy_id;

					hash_insert (retrocopy_h, retrocopy_id_alloc, r);
				}

			zd = hash_lookup (path_h, path);
			if (zd == NULL)
				{
					path_copy = xstrdup (path);

					zd = zygosity_data_new ();
					zd->path = path_copy;

					hash_insert (path_h, path_copy, zd);
				}

			g = xcalloc (1, sizeof (Genotype));

			*g = (Genotype) {
				.region       = r,
				.retrocopy_id = retrocopy_id,
				.source_id    = source_id,
				.acm          = 0,
				.z            = HOMOZYGOUS
			};

			list_append (zd->genotype, g);
		}

	// Clean up
	db_finalize (genotype_query_stmt);
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

static inline void
dump_genotype (sqlite3_stmt *stmt, const Genotype *g)
{
	switch (g->z)
		{
		case HETEROZYGOUS:
			{
				log_debug ("retrocopy [%d %d] is heterozygous", g->retrocopy_id, g->source_id);
				db_insert_genotype (stmt, g->source_id, g->retrocopy_id, g->acm, 1);
				break;
			}
		case HOMOZYGOUS:
			{
				log_debug ("retrocopy [%d %d] is homozygous", g->retrocopy_id, g->source_id);
				db_insert_genotype (stmt, g->source_id, g->retrocopy_id, g->acm, 0);
				break;
			}
		}
}

static void
cross_window (IBiTreeLookupData *ldata, void *user_data)
{
	const CrossWindowLinear *c = user_data;
	Genotype *g = ldata->data;

	if (cross_insertion_point (c->align, g->region->insertion_point, c->phred_quality))
		g->acm++;
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

			if (g->acm >= zd->crossing_reads)
				g->z = HETEROZYGOUS;

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
				if (cross_insertion_point (align, r->insertion_point, zd->phred_quality))
					g->acm++;

			if (rc < -1)
				log_fatal ("Failed to read sam alignment");

			if (g->acm >= zd->crossing_reads)
				g->z = HETEROZYGOUS;

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
genotype (sqlite3_stmt *genotype_stmt, int threads,
		int crossing_reads, int phred_quality)
{
	log_trace ("Inside %s", __func__);
	assert (genotype_stmt != NULL && threads > 0
			&& crossing_reads > 0
			&& phred_quality >= 0);

	threadpool thpool = NULL;

	ChrStd *cs = NULL;
	Hash *retrocopy_h, *path_h;
	retrocopy_h = path_h = NULL;

	HashIter itr = {};
	const char *path = NULL;
	ZygosityData *zd = NULL;

	// Necessary for get the right chr tid
	cs = chr_std_new ();

	// Alloc n threads into the pool
	thpool = thpool_init (threads);

	// RETROCOPY_ID => REGION
	retrocopy_h = hash_new_full (int_hash, int_equal,
			xfree, (DestroyNotify) region_free);

	// PATH => @ZYGOSITYDATA
	path_h = hash_new_full (str_hash, str_equal,
			xfree, (DestroyNotify) zygosity_data_free);

	log_debug ("Clean genotype table");
	clean_genotype_table (sqlite3_db_handle (genotype_stmt));

	log_info ("Index all BAM path => retrocopy relationship");
	genotype_search_db (sqlite3_db_handle (genotype_stmt),
			retrocopy_h, path_h);

	// Iterator through paths
	hash_iter_init (&itr, path_h);

	while (hash_iter_next (&itr, (void **) &path, (void **) &zd))
		{
			log_debug ("Look for retrocopies zygosity of file '%s'", path);

			// Don't forget to give chromosome standardization,
			// genotype statement and crossing_reads
			zd->cs = cs;
			zd->stmt = genotype_stmt;
			zd->crossing_reads = crossing_reads;
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
	hash_free (path_h);
}
