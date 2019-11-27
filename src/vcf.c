#include "config.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "list.h"
#include "array.h"
#include "hash.h"
#include "utils.h"
#include "log.h"
#include "str.h"
#include "retrocopy.h"
#include "fasta.h"
#include "vcf.h"

#define VCF_VERSION "VCFv4.2"

struct _VCFHeader
{
	int         id;
	const char *name;
};

typedef struct _VCFHeader VCFHeader;

struct _VCFBody
{
	// Retrocopy ID
	int         id;

	// Contig
	const char *chr;

	// Confidence interval for
	// IMPRECISE
	long        window_start;
	long        window_end;

	// Parental gene
	const char *parental_gene_name;
	const char *parental_strand;
	int         level;

	// Pos
	long        insertion_point;
	int         insertion_point_type;

	// Orientation/Strand
	double      orientation_rho;
	double      orientation_p_value;

	// Depth
	int         acm;
	int         sr_acm;

	// Host and Near genes
	const char *exonic;
	const char *intragenic;
	const char *near;
};

typedef struct _VCFBody VCFBody;

struct _VCFGenotype
{
	int heterozygous;
	int reference_depth;
	int alternative_depth;
};

typedef struct _VCFGenotype VCFGenotype;

static VCFHeader *
vcf_header_new (const int id, const char *name)
{
	VCFHeader *h = xcalloc (1, sizeof (VCFHeader));

	*h = (VCFHeader) {
		.id   = id, .name = name
	};

	return h;
}

static void
vcf_header_free (VCFHeader *h)
{
	if (h == NULL)
		return;

	xfree ((void *) h->name);
	xfree (h);
}

static void
fasta_string_free (String *s)
{
	string_free (s, 1);
}

static Hash *
vcf_index_fasta (const char *fasta_file)
{
	log_trace ("Inside %s", __func__);

	Hash *idx = NULL;
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;

	fasta = fasta_open_for_reading (fasta_file);
	entry = fasta_entry_new ();

	idx = hash_new_full (str_hash, str_equal, xfree,
			(DestroyNotify) fasta_string_free);

	while (fasta_read (fasta, entry))
		hash_insert (idx, xstrdup (entry->contig->str),
				string_new (entry->sequence->str));

	if (!hash_size (idx))
		log_fatal ("FASTA file '%s' has no entries",
				fasta_file);

	fasta_entry_free (entry);
	fasta_close (fasta);

	return idx;
}

static List *
vcf_get_header_line (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	VCFHeader *h = NULL;
	List *hl = NULL;

	int id = 0;
	const char *path = NULL;
	char *basename = NULL;

	const char sql[] =
		"SELECT id, path\n"
		"FROM source\n"
		"ORDER BY id ASC";

	hl = list_new ((DestroyNotify) vcf_header_free);

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	while (db_step (stmt) == SQLITE_ROW)
		{
			id   = db_column_int  (stmt, 0);
			path = db_column_text (stmt, 1);

			basename = path_file (path, 1);

			h = vcf_header_new (id, basename);
			list_append (hl, h);
		}

	db_finalize (stmt);
	return hl;
}

static void
vcf_print_header (const List *hl, Hash *fidx,
		FILE *fp, VCFOption *opt)
{
	log_trace ("Inside %s", __func__);

	const VCFHeader *h = NULL;
	const ListElmt *cur = NULL;

	Array *contigs = NULL;
	const char *contig = NULL;
	const String *seq = NULL;
	int i = 0;

	char timestamp[32] = {};
	time_t t = 0;
	struct tm *lt = NULL;

	// Get current local datetome as timestamp
	t = time (NULL);
	lt = localtime (&t);
	timestamp[strftime (timestamp, sizeof (timestamp),
			"%Y-%m-%d %H:%M:%S", lt)] = '\0';

	xfprintf (fp,
		"##fileformat=%s\n"
		"##fileDate=%s\n"
		"##source=%sv%s\n",
		VCF_VERSION, timestamp, PACKAGE_NAME, PACKAGE_VERSION);

	if (fidx != NULL && opt->fasta_file != NULL)
		{
			xfprintf (fp,
				"##reference=file://%s\n",
				opt->fasta_file);

			contigs = hash_get_keys_as_array (fidx);
			array_sort (contigs, cmpstringp);

			for (i = 0; i < contigs->len; i++)
				{
					contig = array_get (contigs, i);
					seq = hash_lookup (fidx, contig);

					assert (seq != NULL);

					xfprintf (fp,
						"##contig=<ID=%s,length=%zu>\n",
						contig, seq->len);
				}
		}

	// IMPRECISE = No one SR at breakpoint
	xfprintf (fp,
		"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth of segment containing breakpoint\">\n"
		"##INFO=<ID=EXONIC,Number=1,Type=String,Description=\"Exon IDs separated by '/' for intragenic retrocopy\">\n"
		"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
		"##INFO=<ID=INTRONIC,Number=1,Type=String,Description=\"Intron IDs separated by '/' for intragenic retrocopy\">\n"
		"##INFO=<ID=NEAR,Number=1,Type=String,Description=\"Near Gene IDs separated by '/' for intergenic retrocopy\">\n"
		"##INFO=<ID=ORHO,Number=1,Type=Float,Description=\"Spearman's rho to detect the polarity\">\n"
		"##INFO=<ID=PG,Number=1,Type=String,Description=\"Parental Gene IDs separated by '/'\">\n"
		"##INFO=<ID=PGTYPE,Number=1,Type=String,Description=\"Provides information about parental gene:"
		" 1 = Single parental gene; 2 = Overlapped parental genes; 4 = Near parental genes;"
		" 8 = Hotspot - Multiple parental genes with retrocopy at the same segment\">\n"
		"##INFO=<ID=POLARITY,Number=1,Type=Character,Description=\"Mobile element polarity (+/-)\">\n"
		"##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Total number of SRs at the estimated breakpoint for this site\">\n"
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
		"##ALT=<ID=INS:ME:RTC,Description=\"Insertion of a Retrocopy\">\n"
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth of segment containing breakpoint\">\n"
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

	for (cur = list_head (hl); cur != NULL; cur = list_next (cur))
		{
			h = list_data (cur);
			xfprintf (fp, "\t%s", h->name);
		}

	xfprintf (fp, "\n");

	array_free (contigs, 1);
}

static sqlite3_stmt *
prepare_retrocopy_query_stmt (sqlite3 *db, const long near_gene_dist)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"WITH\n"
		"	gene (gene_name, strand, chr, start, end) AS (\n"
		"		SELECT gene_name, strand, chr, MIN(start), MAX(end)\n"
		"		FROM exon\n"
		"		GROUP BY gene_name\n"
		"	),\n"
		"	intragenic (rid, host) AS (\n"
		"		SELECT DISTINCT r.id, g.gene_name\n"
		"		FROM retrocopy AS r\n"
		"		CROSS JOIN gene AS g\n"
		"		WHERE r.chr = g.chr\n"
		"			AND r.insertion_point BETWEEN g.start AND g.end\n"
		"	),\n"
		"	intragenic_g (rid, host) AS (\n"
		"		SELECT rid, GROUP_CONCAT(host, '/')\n"
		"		FROM intragenic\n"
		"		GROUP BY rid\n"
		"	),\n"
		"	exonic (rid, host) AS (\n"
		"		SELECT DISTINCT r.id, e.gene_name\n"
		"		FROM retrocopy AS r\n"
		"		CROSS JOIN exon AS e\n"
		"		WHERE r.chr = e.chr\n"
		"			AND r.insertion_point BETWEEN e.start AND e.end\n"
		"	),\n"
		"	exonic_g (rid, host) AS (\n"
		"		SELECT rid, GROUP_CONCAT(host, '/')\n"
		"		FROM exonic\n"
		"		GROUP BY rid\n"
		"	),\n"
		"	near (rid, gene) AS (\n"
		"		SELECT DISTINCT r.id, g.gene_name\n"
		"		FROM retrocopy AS r\n"
		"		CROSS JOIN gene AS g\n"
		"		WHERE r.chr = g.chr\n"
		"			AND ((r.insertion_point BETWEEN g.start - 1\n"
		"					AND g.start - $DIST)\n"
		"				OR (r.insertion_point BETWEEN g.end + 1\n"
		"					AND g.end + $DIST))\n"
		"	),\n"
		"	near_g (rid, gene) AS (\n"
		"		SELECT rid, GROUP_CONCAT(gene, '/')\n"
		"		FROM near\n"
		"		GROUP BY rid\n"
		"	),\n"
		"	genotype (retrocopy_id, acm) AS (\n"
		"		SELECT retrocopy_id, COUNT(*)\n"
		"		FROM (\n"
		"			SELECT DISTINCT retrocopy_id, source_id, alignment_id\n"
		"			FROM retrocopy AS r\n"
		"			INNER JOIN cluster_merging AS cm\n"
		"				ON r.id = cm.retrocopy_id\n"
		"			INNER JOIN clustering AS c\n"
		"				USING (cluster_id, cluster_sid)\n"
		"			INNER JOIN alignment AS a\n"
		"				ON a.id = c.alignment_id\n"
		"		)\n"
		"		GROUP BY retrocopy_id\n"
		"	),\n"
		"	genotype_sr (retrocopy_id, sr_acm) AS (\n"
		"		SELECT retrocopy_id, COUNT(*)\n"
		"		FROM (\n"
		"			SELECT DISTINCT retrocopy_id, source_id, alignment_id\n"
		"			FROM retrocopy AS r\n"
		"			INNER JOIN cluster_merging AS cm\n"
		"				ON r.id = cm.retrocopy_id\n"
		"			INNER JOIN clustering AS c\n"
		"				USING (cluster_id, cluster_sid)\n"
		"			INNER JOIN alignment AS a\n"
		"				ON a.id = c.alignment_id\n"
		"			WHERE a.flag & 0x800\n"
		"				AND (\n"
		"					((cigar LIKE '%M%S' OR cigar LIKE '%M%H')\n"
		"							AND (a.pos + a.rlen) = insertion_point)\n"
		"						OR ((cigar LIKE '%S%M' OR cigar LIKE '%H%M')\n"
		"							AND a.pos = insertion_point)\n"
		"				)\n"
		"		)\n"
		"		GROUP BY retrocopy_id\n"
		"	)\n"
		"SELECT r.id, r.chr, window_start, window_end,\n"
		"	parental_gene_name,\n"
		"	CASE\n"
		"		WHEN strand IS NOT NULL\n"
		"			THEN strand\n"
		"		ELSE '?'\n"
		"	END,\n"
		"	level,\n"
		"	insertion_point, insertion_point_type,\n"
		"	CASE\n"
		"		WHEN orientation_rho IS NOT NULL\n"
		"			THEN orientation_rho\n"
		"		ELSE 0.00\n"
		"	END,\n"
		"	orientation_p_value,\n"
		"	acm,\n"
		"	CASE\n"
		"		WHEN sr_acm IS NOT NULL\n"
		"			THEN sr_acm\n"
		"		ELSE 0\n"
		"	END,\n"
		"	CASE\n"
		"		WHEN e.host IS NOT NULL\n"
		"			THEN e.host\n"
		"		ELSE '?'\n"
		"	END,\n"
		"	CASE\n"
		"		WHEN i.host IS NOT NULL\n"
		"			THEN i.host\n"
		"		ELSE '?'\n"
		"	END,\n"
		"	CASE\n"
		"		WHEN n.gene IS NOT NULL\n"
		"			THEN n.gene\n"
		"		ELSE '?'\n"
		"	END\n"
		"FROM retrocopy AS r\n"
		"LEFT JOIN gene AS g\n"
		"	ON r.parental_gene_name = g.gene_name\n"
		"INNER JOIN genotype AS gn\n"
		"	ON r.id = gn.retrocopy_id\n"
		"LEFT JOIN genotype_sr AS gn_sr\n"
		"	ON r.id = gn_sr.retrocopy_id\n"
		"LEFT JOIN near_g AS n\n"
		"	ON r.id = n.rid\n"
		"LEFT JOIN intragenic_g AS i\n"
		"	ON r.id = i.rid\n"
		"LEFT JOIN exonic_g AS e\n"
		"	ON r.id = e.rid\n"
		"ORDER BY r.chr ASC, insertion_point ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int64 (stmt,
			sqlite3_bind_parameter_index (stmt, "$DIST"),
			near_gene_dist);

	return stmt;
}

static sqlite3_stmt *
prepare_genotype_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"WITH\n"
		"	genotype_acm (retrocopy_id, source_id, acm) AS (\n"
		"		SELECT retrocopy_id, source_id, COUNT(*)\n"
		"		FROM (\n"
		"			SELECT DISTINCT retrocopy_id, source_id, alignment_id\n"
		"			FROM retrocopy AS r\n"
		"			INNER JOIN cluster_merging AS cm\n"
		"				ON r.id = cm.retrocopy_id\n"
		"			INNER JOIN clustering AS c\n"
		"				USING (cluster_id, cluster_sid)\n"
		"			INNER JOIN alignment AS a\n"
		"				ON a.id = c.alignment_id\n"
		"		)\n"
		"		GROUP BY retrocopy_id, source_id\n"
		"	)\n"
		"SELECT source_id, heterozygous, reference_depth, acm\n"
		"FROM genotype AS g\n"
		"INNER JOIN genotype_acm AS ga\n"
		"	USING (retrocopy_id, source_id)\n"
		"WHERE retrocopy_id = ?1";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static Hash *
human_haploid_chr (void)
{
	Hash *hc = NULL;
	char *chr = NULL;

	hc = hash_new_full (str_hash, str_equal,
			xfree, NULL);

	chr = xstrdup ("chrY");
	hash_insert (hc, chr, chr);

	chr = xstrdup ("chrM");
	hash_insert (hc, chr, chr);

	return hc;
}

static Hash *
vcf_index_genotype (sqlite3_stmt *stmt, const int retrocopy_id)
{
	log_trace ("Inside %s", __func__);

	Hash *gi = NULL;
	VCFGenotype *g = NULL;

	int source_id = 0;
	int heterozygous = 0;
	int reference_depth = 0;
	int alternative_depth = 0;

	int *source_id_alloc = NULL;

	db_reset (stmt);
	db_clear_bindings (stmt);
	db_bind_int (stmt, 1, retrocopy_id);

	gi = hash_new_full (int_hash, int_equal, xfree, xfree);

	while (db_step (stmt) == SQLITE_ROW)
		{
			source_id         = db_column_int (stmt, 0);
			heterozygous      = db_column_int (stmt, 1);
			reference_depth   = db_column_int (stmt, 2);
			alternative_depth = db_column_int (stmt, 3);

			source_id_alloc = xcalloc (1, sizeof (int));
			*source_id_alloc = source_id;

			g = xcalloc (1, sizeof (VCFGenotype));
			*g = (VCFGenotype) {
				.heterozygous      = heterozygous,
				.reference_depth   = reference_depth,
				.alternative_depth = alternative_depth
			};

			hash_insert (gi, source_id_alloc, g);
		}

	return gi;
}

static void
vcf_get_body_line (sqlite3_stmt *stmt, VCFBody *b)
{
	b->id                   = db_column_int    (stmt, 0);
	b->chr                  = db_column_text   (stmt, 1);
	b->window_start         = db_column_int64  (stmt, 2);
	b->window_end           = db_column_int64  (stmt, 3);
	b->parental_gene_name   = db_column_text   (stmt, 4);
	b->parental_strand      = db_column_text   (stmt, 5);
	b->level                = db_column_int    (stmt, 6);
	b->insertion_point      = db_column_int64  (stmt, 7);
	b->insertion_point_type = db_column_int    (stmt, 8);
	b->orientation_rho      = db_column_double (stmt, 9);
	b->orientation_p_value  = db_column_double (stmt, 10);
	b->acm                  = db_column_int    (stmt, 11);
	b->sr_acm               = db_column_int    (stmt, 12);
	b->exonic               = db_column_text   (stmt, 13);
	b->intragenic           = db_column_text   (stmt, 14);
	b->near                 = db_column_text   (stmt, 15);
}

static void
vcf_print_body (sqlite3 *db, const List *hl, Hash *fidx, FILE *fp, VCFOption *opt)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *retrocopy_stmt = NULL;
	sqlite3_stmt *genotype_stmt = NULL;

	ListElmt *cur = NULL;
	VCFHeader *h = NULL;
	Hash *gi = NULL;
	Hash *hc = NULL;

	VCFGenotype *g = NULL;
	VCFBody b = {};

	String *seq = NULL;
	long pos = 0;
	char base = 0;

	// List of haploid chromosomes
	hc = human_haploid_chr ();

	// Prepare retrocopy query
	retrocopy_stmt = prepare_retrocopy_query_stmt (db,
			opt->near_gene_dist);

	// Prepare genotype query
	genotype_stmt = prepare_genotype_query_stmt (db);

	while (db_step (retrocopy_stmt) == SQLITE_ROW)
		{
			vcf_get_body_line (retrocopy_stmt, &b);
			gi = vcf_index_genotype (genotype_stmt, b.id);

			// Pos is 1 position before the insertion
			pos = b.insertion_point == 1
				? b.insertion_point
				: b.insertion_point - 1;

			if (fidx != NULL)
				{
					seq = hash_lookup (fidx, b.chr);

					if (seq == NULL)
						{
							log_warn ("Contig %s has no reference at '%s'",
									b.chr, opt->fasta_file);
							base = 'N';
						}
					else if (seq->len < pos)
						{
							log_warn ("Retrocopy at %s:%li is outside contig range %li",
									b.chr, b.insertion_point, seq->len);
							base = 'N';
						}
					else
						base = seq->str[pos - 1];
				}
			else
				base = 'N';

			xfprintf (fp,
				"%s\t%li\t.\t%c\t<INS:ME:RTC>\t.\tPASS\tSVTYPE=INS",
				b.chr, pos, base >= 'a' && base <= 'z' ? base - 32 : base);

			// Imprecise retrocopies
			if (b.insertion_point_type == RETROCOPY_INSERTION_POINT_WINDOW_MEAN)
				{
					xfprintf (fp,
						";IMPRECISE;CIPOS=%li,%li",
						b.window_start - b.insertion_point,
						b.window_end - b.insertion_point);
				}

			// Polarity
			if (b.level == RETROCOPY_PASS
					&& b.orientation_p_value <= opt->orientation_error)
				{
					xfprintf (fp,
						";ORHO=%f;POLARITY=%c",
						b.orientation_rho,
						b.orientation_rho >= 0.0
							? !strcmp (b.parental_strand, "+")
								? '+'
								: '-'
							: !strcmp (b.parental_strand, "+")
								? '-'
								: '+');
				}

			// Parental gene
			xfprintf (fp,
				";PG=%s;PGTYPE=%d",
				b.parental_gene_name, b.level);

			// Print HOST and NEAR genes
			if (strcmp (b.intragenic, "?"))
				{
					// Exonic
					if (strcmp (b.exonic, "?"))
						xfprintf (fp, ";EXONIC=%s", b.exonic);

					// If not exonic, but intragenic - or
					// exonic different than intragenic, then
					// it is intronic
					if (strcmp (b.exonic, b.intragenic))
						xfprintf (fp, ";INTRONIC=%s", b.intragenic);
				}
			// If not intragenic, it may be near
			// some gene
			else if (strcmp (b.near, "?"))
				xfprintf (fp, ";NEAR=%s", b.near);

			// Depth
			xfprintf (fp, ";DP=%d", b.acm);

			// Splitted reads for precise retrocopies
			if (b.insertion_point_type
					== RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE)
				{
					xfprintf (fp, ";SR=%d", b.sr_acm);
				}

			// FORMAT
			xfprintf (fp, "\tGT:DP");

			// Print genotype
			for (cur = list_head (hl); cur != NULL; cur = list_next (cur))
				{
					h = list_data (cur);

					g = hash_lookup (gi, &h->id);

					// Is not haploid?
					if (hash_lookup (hc, b.chr) == NULL)
						{
							// Print diploid chromosomes
							if (g == NULL)
								xfprintf (fp, "\t0/0:0");
							else
								{
									xfprintf (fp,
										"\t%s:%d",
										g->heterozygous
											? "0/1"
											: "1/1",
										g->alternative_depth);
								}
						}
					else
						{
							// Print haploid chromosomes
							if (g == NULL)
								xfprintf (fp, "\t0:0");
							else
								{
									xfprintf (fp,
										"\t%s:%d",
										"1", g->alternative_depth);
								}
						}
				}

			xfprintf (fp, "\n");

			hash_free (gi);
			gi = NULL;
		}

	hash_free (hc);
	db_finalize (retrocopy_stmt);
	db_finalize (genotype_stmt);
}

void
vcf (sqlite3 *db, const char *output_file, VCFOption *opt)
{
	assert (db != NULL && opt != NULL);
	log_trace ("Inside %s", __func__);

	FILE *fp = NULL;
	List *hl = NULL;
	Hash *fidx = NULL;

	log_info ("Create VCF file '%s'", output_file);
	fp = xfopen (output_file, "w");

	if (opt->fasta_file != NULL)
		{
			log_info ("Index genome '%s'", opt->fasta_file);
			fidx = vcf_index_fasta (opt->fasta_file);
		}
	else
		log_warn (
			"With no reference genome, "
			"it is not possible to determine the VCF's REF field");

	log_info ("Get VCF header line");
	hl = vcf_get_header_line (db);

	log_info ("Write VCF header");
	vcf_print_header (hl, fidx, fp, opt);

	log_info ("Write VCF body");
	vcf_print_body (db, hl, fidx, fp, opt);

	xfclose (fp);
	list_free (hl);
	hash_free (fidx);
}
