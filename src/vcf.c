#include "config.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "list.h"
#include "hash.h"
#include "log.h"
#include "retrocopy.h"
#include "vcf.h"

#define VCF_VERSION "VCFv4.2"
#define ALPHA_ERROR 0.05

struct _Header
{
	int         id;
	const char *name;
};

typedef struct _Header Header;

struct _Genotype
{
	int heterozygous;
	int reference_depth;
	int alternative_depth;
};

typedef struct _Genotype Genotype;

static Header *
header_new (const int id, const char *name)
{
	Header *h = xcalloc (1, sizeof (Header));

	*h = (Header) {
		.id   = id,
		.name = name
	};

	return h;
}

static void
header_free (Header *h)
{
	if (h == NULL)
		return;

	xfree ((void *) h->name);
	xfree (h);
}

static List *
vcf_get_header_line (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	Header *h = NULL;
	List *hl = NULL;

	int id = 0;
	const char *path = NULL;
	char *basename = NULL;

	const char sql[] =
		"SELECT id, path\n"
		"FROM source\n"
		"ORDER BY id ASC";

	hl = list_new ((DestroyNotify) header_free);

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	while (db_step (stmt) == SQLITE_ROW)
		{
			id   = db_column_int  (stmt, 0);
			path = db_column_text (stmt, 1);

			basename = path_file (path, 1);

			h = header_new (id, basename);
			list_append (hl, h);
		}

	db_finalize (stmt);
	return hl;
}

static void
vcf_print_header (const List *hl, FILE *fp)
{
	log_trace ("Inside %s", __func__);

	const Header *h = NULL;
	const ListElmt *cur = NULL;

	char timestamp[32] = {};
	time_t t = 0;
	struct tm *lt = NULL;

	// Get current local datetome as timestamp
	t = time (NULL);
	lt = localtime (&t);
	timestamp[strftime (timestamp, sizeof (timestamp),
			"%Y-%m-%d %H:%M:%S", lt)] = '\0';

	// IMPRECISE = No one SR at breakpoint
	fprintf (fp,
		"##fileformat=%s\n"
		"##fileDate=%s\n"
		"##source=%sv%s\n"
		"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth of segment containing breakpoint\">\n"
		"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
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
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
		VCF_VERSION, timestamp, PACKAGE_NAME, PACKAGE_VERSION);

	for (cur = list_head (hl); cur != NULL; cur = list_next (cur))
		{
			h = list_data (cur);
			fprintf (fp, "\t%s", h->name);
		}

	fprintf (fp, "\n");
}

static sqlite3_stmt *
prepare_retrocopy_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	const char sql[] =
		"WITH\n"
		"	gene (gene_name, strand) AS (\n"
		"		SELECT DISTINCT gene_name, strand\n"
		"		FROM exon\n"
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
		"					((cigar LIKE '%M%S' OR cigar LIKE '%M%H') AND (a.pos + a.rlen) = insertion_point)\n"
		"						OR ((cigar LIKE '%S%M' OR cigar LIKE '%H%M') AND a.pos = insertion_point)\n"
		"				)\n"
		"		)\n"
		"		GROUP BY retrocopy_id\n"
		"	)\n"
		"SELECT r.id, chr, window_start, window_end,\n"
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
		"	END\n"
		"FROM retrocopy AS r\n"
		"LEFT JOIN gene AS g\n"
		"	ON r.parental_gene_name = g.gene_name\n"
		"INNER JOIN genotype AS gn\n"
		"	ON r.id = gn.retrocopy_id\n"
		"LEFT JOIN genotype_sr AS gn_sr\n"
		"	ON r.id = gn_sr.retrocopy_id";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
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
index_genotype (sqlite3_stmt *stmt, const int retrocopy_id)
{
	log_trace ("Inside %s", __func__);

	Hash *gi = NULL;
	Genotype *g = NULL;

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

			g = xcalloc (1, sizeof (Genotype));
			*g = (Genotype) {
				.heterozygous      = heterozygous,
				.reference_depth   = reference_depth,
				.alternative_depth = alternative_depth
			};

			hash_insert (gi, source_id_alloc, g);
		}

	return gi;
}

static void
vcf_print_body (sqlite3 *db, const List *hl, FILE *fp)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *retrocopy_stmt = NULL;
	sqlite3_stmt *genotype_stmt = NULL;

	ListElmt *cur = NULL;
	Header *h = NULL;
	Hash *gi = NULL;
	Hash *hc = NULL;

	Genotype *g = NULL;

	int id = 0;
	const char *chr = NULL;
	long window_start = 0;
	long window_end = 0;
	const char *parental_gene_name = NULL;
	const char *parental_strand = NULL;
	int level = 0;
	long insertion_point = 0;
	int insertion_point_type = 0;
	double orientation_rho = 0;
	double orientation_p_value = 0;
	int acm = 0;
	int sr_acm = 0;

	// List of haploid chromosomes
	hc = human_haploid_chr ();

	retrocopy_stmt = prepare_retrocopy_query_stmt (db);
	genotype_stmt = prepare_genotype_query_stmt (db);

	while (db_step (retrocopy_stmt) == SQLITE_ROW)
		{
			id                   = db_column_int    (retrocopy_stmt, 0);
			chr                  = db_column_text   (retrocopy_stmt, 1);
			window_start         = db_column_int64  (retrocopy_stmt, 2);
			window_end           = db_column_int64  (retrocopy_stmt, 3);
			parental_gene_name   = db_column_text   (retrocopy_stmt, 4);
			parental_strand      = db_column_text   (retrocopy_stmt, 5);
			level                = db_column_int    (retrocopy_stmt, 6);
			insertion_point      = db_column_int64  (retrocopy_stmt, 7);
			insertion_point_type = db_column_int    (retrocopy_stmt, 8);
			orientation_rho      = db_column_double (retrocopy_stmt, 9);
			orientation_p_value  = db_column_double (retrocopy_stmt, 10);
			acm                  = db_column_int    (retrocopy_stmt, 11);
			sr_acm               = db_column_int    (retrocopy_stmt, 12);

			gi = index_genotype (genotype_stmt, id);

			fprintf (fp,
				"%s\t%li\t.\tN\t<INS:ME:RTC>\t.\tPASS\tSVTYPE=INS",
				chr, insertion_point == 1 ? insertion_point : insertion_point - 1);

			// Imprecise retrocopies
			if (insertion_point_type == RETROCOPY_INSERTION_POINT_WINDOW_MEAN)
				{
					fprintf (fp,
						";IMPRECISE;CIPOS=%li,%li",
						window_start - insertion_point,
						window_end - insertion_point);
				}

			// Polarity
			if (level == RETROCOPY_PASS
					&& orientation_p_value <= ALPHA_ERROR)
				{
					fprintf (fp,
						";ORHO=%f;POLARITY=%c",
						orientation_rho,
						orientation_rho >= 0.0
							? !strcmp (parental_strand, "+")
								? '+'
								: '-'
							: !strcmp (parental_strand, "+")
								? '-'
								: '+');
				}

			fprintf (fp,
				";PG=%s;PGTYPE=%d;DP=%d",
				parental_gene_name, level, acm);

			// Splitted reads for precise retrocopies
			if (insertion_point_type
					== RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE)
				{
					fprintf (fp, ";SR=%d", sr_acm);
				}

			// FORMAT
			fprintf (fp, "\tGT:DP");

			// Print genotype
			for (cur = list_head (hl); cur != NULL; cur = list_next (cur))
				{
					h = list_data (cur);

					g = hash_lookup (gi, &h->id);

					// Is not haploid?
					if (hash_lookup (hc, chr) == NULL)
						{
							// Print diploid chromosomes
							if (g == NULL)
								fprintf (fp, "\t0/0:0");
							else
								{
									fprintf (fp,
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
								fprintf (fp, "\t0:0");
							else
								{
									fprintf (fp,
										"\t%s:%d",
										"1", g->alternative_depth);
								}
						}
				}

			fprintf (fp, "\n");

			hash_free (gi);
			gi = NULL;
		}

	hash_free (hc);
	db_finalize (retrocopy_stmt);
	db_finalize (genotype_stmt);
}

void
vcf (sqlite3 *db, const char *fasta_file, const char *output_file)
{
	log_trace ("Inside %s", __func__);

	FILE *fp = NULL;
	List *hl = NULL;

	fp = xfopen (output_file, "w");
	hl = vcf_get_header_line (db);

	vcf_print_header (hl, fp);
	vcf_print_body (db, hl, fp);

	xfclose (fp);
	list_free (hl);
}
