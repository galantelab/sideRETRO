#include "config.h"

#include <htslib/sam.h>
#include <assert.h>
#include "db.h"
#include "log.h"
#include "genotype.h"

#define MAX_CROSS 10

static sqlite3_stmt *
prepare_genotype_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"WITH\n"
		"	genotype (id, sid, source_id, path) AS (\n"
		"		SELECT DISTINCT cluster_id, cluster_sid, source_id, path\n"
		"		FROM clustering AS c\n"
		"		INNER JOIN alignment AS a\n"
		"			ON c.alignment_id = a.id\n"
		"		INNER JOIN source AS s\n"
		"			ON a.source_id = s.id\n"
		"	)\n"
		"SELECT retrocopy_id, source_id,\n"
		"	chr || ':' || window_start || '-' || window_end,\n"
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

static int
cross_insertion_point (const bam1_t *align, const long insertion_point)
{
	/*
	* Only test with proper aligned reads
	* - paired-end
	* - proper-pair
	* - mapped
	* - mate mapped
	* - not a duplication
	*/
	if (!(align->core.flag & 0x1)
			|| !(align->core.flag & 0x2)
			|| (align->core.flag & 0x4)
			|| (align->core.flag & 0x8)
			|| (align->core.flag & 0x400))
		{
			return 0;
		}

	uint32_t *cigar = bam_get_cigar (align);
	int rlen = bam_cigar2rlen (align->core.n_cigar, cigar);

	long end = rlen <= 0
		? align->core.pos
		: align->core.pos + rlen - 1;

	return insertion_point >= align->core.pos && insertion_point <= end;
}

static Zygosity
zygosity (const char *pos, const long insertion_point, const char *path)
{
	samFile *fp = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *align = NULL;

	hts_idx_t *idx = NULL;
	hts_itr_t *itr = NULL;

	Zygosity z = HOMOZYGOUS;
	int acm = 0;
	int rc = 0;

	// Open BAM file for reading
	fp = sam_open (path, "rb");
	if (fp == NULL)
		log_errno_fatal ("Failed to open '%s' for reading", path);

	// If all OK until now, read header
	hdr = sam_hdr_read (fp);
	if (hdr == NULL)
		log_fatal ("Failed to read sam header");

	// And allocate BAM align entry
	align = bam_init1 ();
	if (align == NULL)
		log_errno_fatal ("Failed to create bam_init1");

	// Look for the index
	idx = sam_index_load (fp, path);
	if (idx == NULL)
		log_fatal ("Failed to open BAM INDEX for '%s'", path);

	// Query position CHR:START-END
	itr = sam_itr_querys (idx, hdr, pos);
	if (itr == NULL)
		log_fatal ("Failed to look for '%s' at %s\n", pos, path);

	while ((rc = sam_itr_next (fp, itr, align)) >= 0)
		{
			if (cross_insertion_point (align, insertion_point))
				acm++;

			if (acm > MAX_CROSS)
				{
					z = HETEROZYGOUS;
					break;
				}
		}

	if (rc < -1)
		log_fatal ("Failed to read sam alignment");

	if (sam_close (fp) < 0)
		log_errno_fatal ("Failed to close '%s'", path);

	sam_itr_destroy (itr);
	bam_hdr_destroy (hdr);
	bam_destroy1 (align);

	return z;
}

void
genotype (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL);

	sqlite3_stmt *genotype_query_stmt = NULL;
	Zygosity z = 0;

	int retrocopy_id = 0;
	int source_id = 0;
	const char *pos = NULL;
	long insertion_point = 0;
	const char *path = NULL;

	genotype_query_stmt = prepare_genotype_query_stmt (db);

	while (db_step (genotype_query_stmt) == SQLITE_ROW)
		{
			retrocopy_id    = db_column_int   (genotype_query_stmt, 0);
			source_id       = db_column_int   (genotype_query_stmt, 1);
			pos             = db_column_text  (genotype_query_stmt, 2);
			insertion_point = db_column_int64 (genotype_query_stmt, 3);
			path            = db_column_text  (genotype_query_stmt, 4);

			log_debug ("Look for zygosity of retrocopy [%d %d], position %s at %s",
					retrocopy_id, source_id, pos, path);

			z = zygosity (pos, insertion_point, path);
			switch (z)
				{
				case HETEROZYGOUS:
					{
						log_debug ("retrocopy [%d %d] is heterozygous", retrocopy_id, source_id);
						break;
					}
				case HOMOZYGOUS:
					{
						log_debug ("retrocopy [%d %d] is homozygous", retrocopy_id, source_id);
						break;
					}
				}
		}

	db_finalize (genotype_query_stmt);
}
