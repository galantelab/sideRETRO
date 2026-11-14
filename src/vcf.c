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

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "list.h"
#include "array.h"
#include "hash.h"
#include "log.h"
#include "str.h"
#include "retrocopy.h"
#include "fasta.h"
#include "vcf.h"

#define VCF_VERSION "VCFv4.2"

// Do not print the TSD sequence if its length is
// greater than this value - the breakpoints are
// probably inconsistent
#define VCF_MAX_TSD_LEN 100

// Retrocopy ID at the VCF, rnames and FASTA files
#define VCF_ID_FMT "RTC%06d"

// FASTA line width
#define VCF_FASTA_WIDTH 60

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

	// Breakpoints
	long        breakpoint_left;
	int         breakpoint_left_sr;
	long        breakpoint_right;
	int         breakpoint_right_sr;

	// Poly-A
	int         polya_len;
	int         polya_lower_bound;
	const char *polya_strand;
	int         polya_sr;

	// Exons of the parental genes and
	// exon-exon junctions
	const char *pg_exons;
	int         eej;
	const char *eej_list;

	// Retrocopy ends into the parental gene,
	// predicted sequence and observed junctions
	long        five_prime_pos;
	const char *five_prime_exon;
	int         five_prime_sr;
	long        three_prime_pos;
	const char *three_prime_exon;
	int         three_prime_sr;
	int         ins_len;
	const char *ins_seq;
	const char *left_seq;
	const char *right_seq;

	// Host and Near genes
	const char *exonic;
	const char *intragenic;
	const char *near;
};

typedef struct _VCFBody VCFBody;

struct _VCFGenotype
{
	int    reference_depth;
	int    alternate_depth;
	double ho_ref_likelihood;
	double he_likelihood;
	double ho_alt_likelihood;
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
vcf_print_header (const List *hl, Hash *contig_len,
		FILE *fp, VCFOption *opt)
{
	log_trace ("Inside %s", __func__);

	const VCFHeader *h = NULL;
	const ListElmt *cur = NULL;

	Array *contigs = NULL;
	const char *contig = NULL;
	const long *len = NULL;
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

	if (contig_len != NULL && opt->fasta_file != NULL)
		{
			xfprintf (fp,
				"##reference=file://%s\n",
				opt->fasta_file);

			contigs = hash_get_keys_as_array (contig_len);
			array_sort (contigs, cmpstringp);

			for (i = 0; i < contigs->len; i++)
				{
					contig = array_get (contigs, i);
					len = hash_lookup (contig_len, contig);

					assert (len != NULL);

					xfprintf (fp,
						"##contig=<ID=%s,length=%li>\n",
						contig, *len);
				}
		}

	// IMPRECISE = No one SR at breakpoint
	xfprintf (fp,
		"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth of segment containing breakpoint\">\n"
		"##INFO=<ID=EEJ,Number=1,Type=Integer,Description=\"Number of reads spanning exon-exon junctions of the parental gene:"
		" Evidence of a processed transcript\">\n"
		"##INFO=<ID=EEJLIST,Number=1,Type=String,Description=\"Exon-exon junctions spanned by reads as 'exon_id1-exon_id2' separated by '/'\">\n"
		"##INFO=<ID=EXONIC,Number=1,Type=String,Description=\"Exon IDs separated by '/' for intragenic retrocopy\">\n"
		"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
		"##INFO=<ID=INTRONIC,Number=1,Type=String,Description=\"Intron IDs separated by '/' for intragenic retrocopy\">\n"
		"##INFO=<ID=LEFT_SVINSSEQ,Number=1,Type=String,Description=\"Observed sequence of the insertion at the left junction:"
		" Consensus of the reads clipped at LBKP\">\n"
		"##INFO=<ID=LBKP,Number=1,Type=Integer,Description=\"Left breakpoint: Last reference base of the left flank, from split reads clipped at the right\">\n"
		"##INFO=<ID=LSR,Number=1,Type=Integer,Description=\"Number of split reads supporting LBKP\">\n"
		"##INFO=<ID=NEAR,Number=1,Type=String,Description=\"Near Gene IDs separated by '/' for intergenic retrocopy\">\n"
		"##INFO=<ID=ORHO,Number=1,Type=Float,Description=\"Spearman's rho to detect the polarity\">\n"
		"##INFO=<ID=PG,Number=1,Type=String,Description=\"Parental Gene IDs separated by '/'\">\n"
		"##INFO=<ID=PGEXONS,Number=1,Type=String,Description=\"Number of exons of each parental gene as 'gene:n' separated by '/'."
		" A single exon gene has no exon-exon junction\">\n"
		"##INFO=<ID=PG5P,Number=1,Type=String,Description=\"5' end of the retrocopy into the parental gene as 'exon_id:position:reads',"
		" where position is at the reference genome\">\n"
		"##INFO=<ID=PG3P,Number=1,Type=String,Description=\"3' end of the retrocopy into the parental gene, before the poly-A,"
		" as 'exon_id:position:reads'\">\n"
		"##INFO=<ID=PGTYPE,Number=1,Type=String,Description=\"Provides information about parental gene:"
		" 1 = Single parental gene; 2 = Overlapped parental genes; 4 = Near parental genes;"
		" 8 = Hotspot - Multiple parental genes with retrocopy at the same segment\">\n"
		"##INFO=<ID=POLARITY,Number=1,Type=Character,Description=\"Mobile element polarity (+/-):"
		" From the poly-A junction when POLYA is present, otherwise from ORHO\">\n"
		"##INFO=<ID=POLYA,Number=1,Type=Integer,Description=\"Poly-A tail length observed in the reads crossing its junction\">\n"
		"##INFO=<ID=POLYAMIN,Number=0,Type=Flag,Description=\"POLYA is a lower bound: The poly-A extends beyond the longest clipped read\">\n"
		"##INFO=<ID=POLYASR,Number=1,Type=Integer,Description=\"Number of reads supporting POLYA\">\n"
		"##INFO=<ID=RBKP,Number=1,Type=Integer,Description=\"Right breakpoint: First reference base of the right flank, from split reads clipped at the left\">\n"
		"##INFO=<ID=RIGHT_SVINSSEQ,Number=1,Type=String,Description=\"Observed sequence of the insertion at the right junction:"
		" Consensus of the reads clipped at RBKP\">\n"
		"##INFO=<ID=RSR,Number=1,Type=Integer,Description=\"Number of split reads supporting RBKP\">\n"
		"##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Total number of SRs at the estimated breakpoints for this site\">\n"
		"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Predicted length of the inserted retrocopy: The reference exon bases"
		" between its ends plus the poly-A. The 3' end is assumed at the end of the last exon when not observed\">\n"
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
		"##INFO=<ID=TSD,Number=1,Type=String,Description=\"Target site duplication sequence: Reference bases from RBKP to LBKP\">\n"
		"##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description=\"Target site duplication length (LBKP - RBKP + 1):"
		" Positive for duplication, zero for a clean insertion and negative for target site deletion\">\n"
		"##ALT=<ID=INS:ME:RTC,Description=\"Insertion of a Retrocopy\">\n"
		"##FORMAT=<ID=DP2,Number=2,Type=Integer,Description=\"Read depth at this position for this sample separated by comma:"
		" The former is reference depth and the last is alternate depth\">\n"
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihood\">\n"
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
		"	gene_exons (gene_name, n) AS (\n"
		"		SELECT gene_name, COUNT(*)\n"
		"		FROM exon\n"
		"		GROUP BY gene_name\n"
		"	),\n"
		"	pg_exons (rid, txt) AS (\n"
		"		SELECT rid, GROUP_CONCAT(gene_name || ':' || n, '/')\n"
		"		FROM (\n"
		"			SELECT DISTINCT cm.retrocopy_id AS rid, cl.gene_name\n"
		"			FROM cluster_merging AS cm\n"
		"			INNER JOIN cluster AS cl\n"
		"				ON cl.id = cm.cluster_id AND cl.sid = cm.cluster_sid\n"
		"			ORDER BY cm.retrocopy_id, cl.gene_name\n"
		"		)\n"
		"		INNER JOIN gene_exons\n"
		"			USING (gene_name)\n"
		"		GROUP BY rid\n"
		"	),\n"
		"	junction_g (rid, eej, txt) AS (\n"
		"		SELECT rid, SUM(reads), GROUP_CONCAT(ense1 || '-' || ense2, '/')\n"
		"		FROM (\n"
		"			SELECT j.retrocopy_id AS rid, j.reads, e1.ense AS ense1, e2.ense AS ense2\n"
		"			FROM junction AS j\n"
		"			INNER JOIN exon AS e1\n"
		"				ON e1.id = j.exon_id1\n"
		"			INNER JOIN exon AS e2\n"
		"				ON e2.id = j.exon_id2\n"
		"			ORDER BY j.retrocopy_id, e1.start, e2.start\n"
		"		)\n"
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
		"	)\n"
		"SELECT r.id, r.chr, window_start, window_end,\n"
		"	parental_gene_name,\n"
		"	CASE\n"
		"		WHEN g.strand IS NOT NULL\n"
		"			THEN g.strand\n"
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
		"	breakpoint_left_sr + breakpoint_right_sr,\n"
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
		"	END,\n"
		"	COALESCE(breakpoint_left, 0), breakpoint_left_sr,\n"
		"	COALESCE(breakpoint_right, 0), breakpoint_right_sr,\n"
		"	COALESCE(polya_len, 0), polya_lower_bound,\n"
		"	COALESCE(polya_strand, '?'), polya_sr,\n"
		"	COALESCE(p.txt, '?'), COALESCE(j.eej, 0), COALESCE(j.txt, '?'),\n"
		"	COALESCE(ins.five_prime_pos, 0), COALESCE(e5.ense, '?'),\n"
		"	COALESCE(ins.five_prime_sr, 0),\n"
		"	COALESCE(ins.three_prime_pos, 0), COALESCE(e3.ense, '?'),\n"
		"	COALESCE(ins.three_prime_sr, 0), COALESCE(ins.len, 0),\n"
		"	COALESCE(ins.seq, '?'), COALESCE(ins.left_seq, '?'),\n"
		"	COALESCE(ins.right_seq, '?')\n"
		"FROM retrocopy AS r\n"
		"LEFT JOIN pg_exons AS p\n"
		"	ON r.id = p.rid\n"
		"LEFT JOIN junction_g AS j\n"
		"	ON r.id = j.rid\n"
		"LEFT JOIN insertion AS ins\n"
		"	ON r.id = ins.retrocopy_id\n"
		"LEFT JOIN exon AS e5\n"
		"	ON e5.id = ins.five_prime_exon_id\n"
		"LEFT JOIN exon AS e3\n"
		"	ON e3.id = ins.three_prime_exon_id\n"
		"LEFT JOIN gene AS g\n"
		"	ON r.parental_gene_name = g.gene_name\n"
		"INNER JOIN genotype AS gn\n"
		"	ON r.id = gn.retrocopy_id\n"
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
		"SELECT source_id,\n"
		"	reference_depth,\n"
		"	alternate_depth,\n"
		"	ho_ref_likelihood,\n"
		"	he_likelihood,\n"
		"	ho_alt_likelihood\n"
		"FROM genotype\n"
		"WHERE retrocopy_id = ?1";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static Hash *
vcf_index_genotype (sqlite3_stmt *stmt, const int retrocopy_id)
{
	log_trace ("Inside %s", __func__);

	Hash *gi = NULL;
	VCFGenotype *g = NULL;

	int source_id = 0;
	int reference_depth = 0;
	int alternate_depth = 0;
	double ho_ref_likelihood = 0.0;
	double he_likelihood = 0.0;
	double ho_alt_likelihood = 0.0;

	int *source_id_alloc = NULL;

	db_reset (stmt);
	db_clear_bindings (stmt);
	db_bind_int (stmt, 1, retrocopy_id);

	gi = hash_new_full (int_hash, int_equal, xfree, xfree);

	while (db_step (stmt) == SQLITE_ROW)
		{
			source_id          = db_column_int    (stmt, 0);
			reference_depth    = db_column_int    (stmt, 1);
			alternate_depth    = db_column_int    (stmt, 2);
			ho_ref_likelihood  = db_column_double (stmt, 3);
			he_likelihood      = db_column_double (stmt, 4);
			ho_alt_likelihood  = db_column_double (stmt, 5);

			source_id_alloc = xcalloc (1, sizeof (int));
			*source_id_alloc = source_id;

			g = xcalloc (1, sizeof (VCFGenotype));
			*g = (VCFGenotype) {
				.reference_depth   = reference_depth,
				.alternate_depth   = alternate_depth,
				.ho_ref_likelihood = ho_ref_likelihood,
				.he_likelihood     = he_likelihood,
				.ho_alt_likelihood = ho_alt_likelihood
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
	b->breakpoint_left      = db_column_int64  (stmt, 16);
	b->breakpoint_left_sr   = db_column_int    (stmt, 17);
	b->breakpoint_right     = db_column_int64  (stmt, 18);
	b->breakpoint_right_sr  = db_column_int    (stmt, 19);
	b->polya_len            = db_column_int    (stmt, 20);
	b->polya_lower_bound    = db_column_int    (stmt, 21);
	b->polya_strand         = db_column_text   (stmt, 22);
	b->polya_sr             = db_column_int    (stmt, 23);
	b->pg_exons             = db_column_text   (stmt, 24);
	b->eej                  = db_column_int    (stmt, 25);
	b->eej_list             = db_column_text   (stmt, 26);
	b->five_prime_pos       = db_column_int64  (stmt, 27);
	b->five_prime_exon      = db_column_text   (stmt, 28);
	b->five_prime_sr        = db_column_int    (stmt, 29);
	b->three_prime_pos      = db_column_int64  (stmt, 30);
	b->three_prime_exon     = db_column_text   (stmt, 31);
	b->three_prime_sr       = db_column_int    (stmt, 32);
	b->ins_len              = db_column_int    (stmt, 33);
	b->ins_seq              = db_column_text   (stmt, 34);
	b->left_seq             = db_column_text   (stmt, 35);
	b->right_seq            = db_column_text   (stmt, 36);
}

static const char *
genotype_likelihood (const double ho_ref, const double he, const double ho_alt)
{
	return
		((ho_ref > he || fequal (ho_ref, he)) && (ho_ref > ho_alt))
		|| ((ho_ref > ho_alt || fequal (ho_ref, ho_alt)) && (ho_ref > he))
		? "0/0"
		: (!fequal (ho_ref, 0) && fequal (ho_ref, he) && fequal (ho_ref, ho_alt))
			|| (he > ho_ref && (he > ho_alt || fequal (he, ho_alt)))
			? "0/1"
			: ho_alt > ho_ref && ho_alt > he
				? "1/1"
				: "./.";
}

static Hash *
vcf_index_regions (sqlite3 *db, const char *fasta_file, Hash **contig_len)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;
	List *regions = NULL;
	Hash *fidx = NULL;

	const char *chr = NULL;
	long ip = 0;
	long left = 0;
	long right = 0;
	long pos = 0;

	// The only regions needed: The base before the
	// insertion (REF) and the TSD
	const char sql[] =
		"SELECT chr, insertion_point,\n"
		"	COALESCE(breakpoint_left, 0), breakpoint_left_sr,\n"
		"	COALESCE(breakpoint_right, 0), breakpoint_right_sr\n"
		"FROM retrocopy";

	regions = list_new ((DestroyNotify) fasta_region_free);

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	while (db_step (stmt) == SQLITE_ROW)
		{
			chr   = db_column_text  (stmt, 0);
			ip    = db_column_int64 (stmt, 1);
			left  = db_column_int64 (stmt, 2);
			right = db_column_int64 (stmt, 4);

			pos = ip == 1 ? ip : ip - 1;
			list_append (regions, fasta_region_new (chr, pos, pos));

			if (db_column_int (stmt, 3) > 0 && db_column_int (stmt, 5) > 0
					&& left >= right && (left - right + 1) <= VCF_MAX_TSD_LEN)
				list_append (regions, fasta_region_new (chr, right, left));
		}

	log_info ("Extract %zu regions from '%s'", list_size (regions), fasta_file);
	fidx = fasta_index_regions (fasta_file, regions, contig_len);

	list_free (regions);
	db_finalize (stmt);

	return fidx;
}

static const char *
vcf_region_lookup (Hash *fidx, const char *chr, long start, long end)
{
	const char *seq = NULL;
	char *key = NULL;

	if (fidx == NULL)
		return NULL;

	key = fasta_region_key (chr, start, end);
	seq = hash_lookup (fidx, key);
	xfree (key);

	return seq;
}

static sqlite3_stmt *
prepare_rnames_query_stmt (sqlite3 *db)
{
	log_trace ("Inside %s", __func__);

	// Reads used to detect the retrocopy by source, with
	// their roles in the evidence: The roles are flags,
	// so their distinct sum is a bitwise OR
	const char sql[] =
		"WITH\n"
		"	role (source_id, qname, role) AS (\n"
		"		SELECT a.source_id, a.qname, SUM(DISTINCT e.role)\n"
		"		FROM evidence AS e\n"
		"		INNER JOIN alignment AS a\n"
		"			ON a.id = e.alignment_id\n"
		"		WHERE e.retrocopy_id = ?1\n"
		"		GROUP BY a.source_id, a.qname\n"
		"	),\n"
		"	read (source_id, qname) AS (\n"
		"		SELECT DISTINCT a.source_id, a.qname\n"
		"		FROM cluster_merging AS cm\n"
		"		INNER JOIN clustering AS c\n"
		"			ON c.cluster_id = cm.cluster_id\n"
		"				AND c.cluster_sid = cm.cluster_sid\n"
		"		INNER JOIN alignment AS a\n"
		"			ON a.id = c.alignment_id\n"
		"		WHERE cm.retrocopy_id = ?1\n"
		"	)\n"
		"SELECT r.source_id, r.qname, COALESCE(o.role, 0)\n"
		"FROM read AS r\n"
		"LEFT JOIN role AS o\n"
		"	USING (source_id, qname)\n"
		"ORDER BY r.source_id ASC, r.qname ASC";

	log_debug ("Query schema:\n%s", sql);
	return db_prepare (db, sql);
}

static const char *
vcf_evidence_roles (const int role, char *buf, const size_t size)
{
	buf[0] = '\0';

	// No role: A discordant read
	if (role == 0)
		{
			snprintf (buf, size, "discordant");
			return buf;
		}

	snprintf (buf, size, "%s%s%s%s",
			role & RETROCOPY_EVIDENCE_SPLIT_LEFT  ? ",split_left"  : "",
			role & RETROCOPY_EVIDENCE_SPLIT_RIGHT ? ",split_right" : "",
			role & RETROCOPY_EVIDENCE_POLYA       ? ",polya"       : "",
			role & RETROCOPY_EVIDENCE_JUNCTION    ? ",junction"    : "");

	// Skip the leading comma
	memmove (buf, buf + 1, strlen (buf));
	return buf;
}

static void
vcf_print_rnames (sqlite3_stmt *stmt, const List *hl, FILE *fp,
		const int retrocopy_id)
{
	const VCFHeader *h = NULL;
	const ListElmt *cur = NULL;
	int source_id = 0;
	int source_id_prev = 0;

	db_reset (stmt);
	db_clear_bindings (stmt);
	db_bind_int (stmt, 1, retrocopy_id);

	char roles[64] = {};

	/*
	* One read per line: retrocopy id, sample, read name and
	* its roles, separated by tab. A read name (QNAME) never
	* has tab, so it can be written as it is - e.g. Illumina
	* names with colons
	*/
	while (db_step (stmt) == SQLITE_ROW)
		{
			source_id = db_column_int (stmt, 0);

			// Look for the sample name only when it changes
			if (source_id != source_id_prev)
				{
					for (cur = list_head (hl); cur != NULL; cur = list_next (cur))
						{
							h = list_data (cur);
							if (h->id == source_id)
								break;
						}

					assert (cur != NULL);
					source_id_prev = source_id;
				}

			xfprintf (fp, VCF_ID_FMT "\t%s\t%s\t%s\n", retrocopy_id,
					h->name, db_column_text (stmt, 1),
					vcf_evidence_roles (db_column_int (stmt, 2),
						roles, sizeof (roles)));
		}
}

static void
vcf_print_insertion (FILE *fp, const VCFBody *b, const long pos)
{
	const char *seq = b->ins_seq;
	size_t len = strlen (seq);
	size_t i = 0;

	xfprintf (fp, ">" VCF_ID_FMT " %s:%li PG=%s SVLEN=%d predicted\n",
			b->id, b->chr, pos, b->parental_gene_name, b->ins_len);

	for (i = 0; i < len; i += VCF_FASTA_WIDTH)
		xfprintf (fp, "%.*s\n", VCF_FASTA_WIDTH, seq + i);
}

static void
vcf_print_body (sqlite3 *db, const List *hl, Hash *fidx, FILE *fp,
		FILE *rnames_fp, FILE *ins_fp, VCFOption *opt)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *retrocopy_stmt = NULL;
	sqlite3_stmt *genotype_stmt = NULL;

	ListElmt *cur = NULL;
	VCFHeader *h = NULL;
	Hash *gi = NULL;

	VCFGenotype *g = NULL;
	VCFBody b = {};

	long pos = 0;
	char base = 0;

	long tsd_len = 0;
	char polarity = 0;

	sqlite3_stmt *rnames_stmt = NULL;
	const char *ref = NULL;

	// Prepare retrocopy query
	retrocopy_stmt = prepare_retrocopy_query_stmt (db,
			opt->near_gene_dist);

	// Prepare genotype query
	genotype_stmt = prepare_genotype_query_stmt (db);

	// Prepare read names query
	rnames_stmt = prepare_rnames_query_stmt (db);

	while (db_step (retrocopy_stmt) == SQLITE_ROW)
		{
			vcf_get_body_line (retrocopy_stmt, &b);
			gi = vcf_index_genotype (genotype_stmt, b.id);

			// Pos is 1 position before the insertion
			pos = b.insertion_point == 1
				? b.insertion_point
				: b.insertion_point - 1;

			base = 'N';

			if (fidx != NULL)
				{
					ref = vcf_region_lookup (fidx, b.chr, pos, pos);

					if (ref == NULL)
						log_warn ("Retrocopy at %s:%li not found at the reference '%s': REF is 'N'",
								b.chr, b.insertion_point, opt->fasta_file);
					else
						base = ref[0];
				}

			xfprintf (fp,
				"%s\t%li\t" VCF_ID_FMT "\t%c\t<INS:ME:RTC>\t.\tPASS\tSVTYPE=INS",
				b.chr, pos, b.id, base);

			// Imprecise retrocopies
			if (b.insertion_point_type == RETROCOPY_INSERTION_POINT_WINDOW_MEAN)
				{
					xfprintf (fp,
						";IMPRECISE;CIPOS=%li,%li",
						b.window_start - b.insertion_point,
						b.window_end - b.insertion_point);
				}

			// Polarity: The poly-A junction is direct evidence.
			// Otherwise use the orientation correlation
			polarity = '?';

			if (b.level == RETROCOPY_PASS
					&& b.orientation_p_value <= opt->orientation_error)
				{
					xfprintf (fp, ";ORHO=%f", b.orientation_rho);

					polarity = b.orientation_rho >= 0.0
						? !strcmp (b.parental_strand, "+")
							? '+'
							: '-'
						: !strcmp (b.parental_strand, "+")
							? '-'
							: '+';
				}

			if (b.polya_sr > 0)
				polarity = b.polya_strand[0];

			if (polarity != '?')
				xfprintf (fp, ";POLARITY=%c", polarity);

			// Parental gene
			xfprintf (fp,
				";PG=%s;PGTYPE=%d",
				b.parental_gene_name, b.level);

			if (strcmp (b.pg_exons, "?"))
				xfprintf (fp, ";PGEXONS=%s", b.pg_exons);

			// Exon-exon junctions
			if (b.eej > 0)
				xfprintf (fp, ";EEJ=%d;EEJLIST=%s", b.eej, b.eej_list);

			// Retrocopy ends into the parental gene
			if (b.five_prime_sr > 0)
				xfprintf (fp, ";PG5P=%s:%li:%d", b.five_prime_exon,
						b.five_prime_pos, b.five_prime_sr);

			if (b.three_prime_sr > 0)
				xfprintf (fp, ";PG3P=%s:%li:%d", b.three_prime_exon,
						b.three_prime_pos, b.three_prime_sr);

			if (b.ins_len > 0)
				xfprintf (fp, ";SVLEN=%d", b.ins_len);

			// Observed sequences at the junctions
			if (strcmp (b.left_seq, "?"))
				xfprintf (fp, ";LEFT_SVINSSEQ=%s", b.left_seq);

			if (strcmp (b.right_seq, "?"))
				xfprintf (fp, ";RIGHT_SVINSSEQ=%s", b.right_seq);

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

			// Breakpoints
			if (b.breakpoint_left_sr > 0)
				xfprintf (fp, ";LBKP=%li;LSR=%d",
						b.breakpoint_left, b.breakpoint_left_sr);

			if (b.breakpoint_right_sr > 0)
				xfprintf (fp, ";RBKP=%li;RSR=%d",
						b.breakpoint_right, b.breakpoint_right_sr);

			// Target site duplication: The reference
			// bases between the two breakpoints
			if (b.breakpoint_left_sr > 0 && b.breakpoint_right_sr > 0)
				{
					tsd_len = b.breakpoint_left - b.breakpoint_right + 1;
					xfprintf (fp, ";TSDLEN=%li", tsd_len);

					if (tsd_len > 0 && tsd_len <= VCF_MAX_TSD_LEN)
						{
							ref = vcf_region_lookup (fidx, b.chr,
									b.breakpoint_right, b.breakpoint_left);

							if (ref != NULL)
								xfprintf (fp, ";TSD=%s", ref);
						}
				}

			// Poly-A
			if (b.polya_sr > 0)
				{
					xfprintf (fp, ";POLYA=%d", b.polya_len);

					if (b.polya_lower_bound)
						xfprintf (fp, ";POLYAMIN");

					xfprintf (fp, ";POLYASR=%d", b.polya_sr);
				}

			// FORMAT
			xfprintf (fp, "\tGT:DP2:GL");

			// Print genotype
			for (cur = list_head (hl); cur != NULL; cur = list_next (cur))
				{
					h = list_data (cur);

					g = hash_lookup (gi, &h->id);
					assert (g != NULL);

					// Print diploid chromosomes
					xfprintf (fp,
						"\t%s:%d,%d:%g,%g,%g",
						genotype_likelihood (g->ho_ref_likelihood,
							g->he_likelihood, g->ho_alt_likelihood),
						g->reference_depth,
						g->alternate_depth,
						g->ho_ref_likelihood,
						g->he_likelihood,
						g->ho_alt_likelihood);
				}

			xfprintf (fp, "\n");

			// Reads used to detect this retrocopy
			if (rnames_fp != NULL)
				vcf_print_rnames (rnames_stmt, hl, rnames_fp, b.id);

			// Predicted inserted sequence
			if (ins_fp != NULL && strcmp (b.ins_seq, "?"))
				{
					vcf_print_insertion (ins_fp, &b, pos);
					opt->num_insertions++;
				}

			opt->num_records++;

			hash_free (gi);
			gi = NULL;
		}

	db_finalize (retrocopy_stmt);
	db_finalize (genotype_stmt);
	db_finalize (rnames_stmt);
}

void
vcf (sqlite3 *db, const char *output_file, VCFOption *opt)
{
	assert (db != NULL && opt != NULL);
	log_trace ("Inside %s", __func__);

	FILE *fp = NULL;
	FILE *rnames_fp = NULL;
	FILE *ins_fp = NULL;
	List *hl = NULL;
	Hash *fidx = NULL;
	Hash *contig_len = NULL;

	log_info ("Write the VCF '%s'", output_file);
	fp = xfopen (output_file, "w");

	if (opt->rnames_file != NULL)
		{
			log_debug ("Create the read names file '%s'", opt->rnames_file);
			rnames_fp = xfopen (opt->rnames_file, "w");
			xfprintf (rnames_fp, "#retrocopy_id\tsample\tread_name\trole\n");
		}

	if (opt->insertion_file != NULL)
		{
			log_debug ("Create the predicted insertions file '%s'", opt->insertion_file);
			ins_fp = xfopen (opt->insertion_file, "w");
		}

	if (opt->fasta_file != NULL)
		fidx = vcf_index_regions (db, opt->fasta_file, &contig_len);
	else
		log_warn ("No reference genome (--reference-file): REF is 'N' and TSD is not reported");

	log_debug ("Get the VCF header line");
	hl = vcf_get_header_line (db);

	log_debug ("Write the VCF header");
	vcf_print_header (hl, contig_len, fp, opt);

	log_debug ("Write the VCF body");
	vcf_print_body (db, hl, fidx, fp, rnames_fp, ins_fp, opt);

	xfclose (fp);

	if (rnames_fp != NULL)
		xfclose (rnames_fp);

	if (ins_fp != NULL)
		xfclose (ins_fp);

	list_free (hl);

	if (fidx != NULL)
		hash_free (fidx);

	if (contig_len != NULL)
		hash_free (contig_len);
}
