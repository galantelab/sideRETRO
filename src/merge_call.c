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

#include <time.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "db.h"
#include "io.h"
#include "log.h"
#include "logger.h"
#include "utils.h"
#include "wrapper.h"
#include "list.h"
#include "set.h"
#include "gff.h"
#include "str.h"
#include "chr.h"
#include "array.h"
#include "blacklist.h"
#include "cluster.h"
#include "db_merge.h"
#include "retrocopy.h"
#include "genotype.h"
#include "summary.h"
#include "merge_call.h"

#define DEFAULT_CACHE_SIZE            200000 /* 200MiB */
#define DEFAULT_PREFIX                "out"
#define DEFAULT_OUTPUT_DIR            "."
#define DEFAULT_LOG_SILENT            0
#define DEFAULT_LOG_LEVEL             LOG_INFO
#define DEFAULT_IN_PLACE              0
#define DEFAULT_EPS                   300
#define DEFAULT_MIN_PTS               10
#define DEFAULT_BLACKLIST_CHR         "chrM"
#define DEFAULT_PARENTAL_DISTANCE     1000000
#define DEFAULT_SUPPORT               1
#define DEFAULT_BLACKLIST_PADDING     0
#define DEFAULT_GFF_FEATURE           "gene"
#define DEFAULT_GFF_ATTRIBUTE         "gene_type=processed_pseudogene|tag=retrogene"
#define DEFAULT_GFF_GENE_NAME         "gene_name"
#define DEFAULT_NEAR_GENE_RANK        3
#define DEFAULT_THREADS               1
#define DEFAULT_PHRED_QUALITY         8

struct _MergeCall
{
	// Mandatory
	Array       *db_files;

	// I/O
	const char  *input_file;
	const char  *output_dir;
	const char  *output_file;
	const char  *prefix;
	int          in_place;
	int          num_merged;

	// Log
	Logger      *logger;
	const char  *log_file;
	int          log_level;
	int          silent;

	// SQLite3
	int          cache_size;

	// Clustering
	int          epsilon;
	int          min_pts;

	// Filtering & Annotation
	const char  *blacklist_region;
	Set         *blacklist_chr;
	GffFilter   *filter;
	List        *attributes;
	const char  *gene_name_key;
	ChrStd      *cs;
	int          padding;
	int          parental_dist;
	int          support;
	int          near_gene_rank;
	const char  *fasta_file;

	// Genotyping
	int          threads;
	int          phred_quality;
};

typedef struct _MergeCall MergeCall;

static long
query_count (sqlite3 *db, const char *sql)
{
	sqlite3_stmt *stmt = NULL;
	long n = 0;

	stmt = db_prepare (db, sql);

	if (db_step (stmt) == SQLITE_ROW)
		n = db_column_int64 (stmt, 0);

	db_finalize (stmt);
	return n;
}

static void
print_summary (const MergeCall *mc, sqlite3 *db, const char *db_file,
		int num_clusters, time_t start)
{
	Summary *s = NULL;
	sqlite3_stmt *stmt = NULL;
	char *name = NULL;
	char b1[SUMMARY_NUM_LEN] = {};
	char b2[SUMMARY_NUM_LEN] = {};
	long n = 0;
	long total = 0;

	s = summary_new ("merge-call summary", start);

	summary_section (s, "Input");

	// --in-place: The output database is an input
	// database too, but it is not at the list
	summary_line (s, "Databases", "%d",
			mc->num_merged + (mc->in_place ? 1 : 0));
	summary_line (s, "Samples", "%s",
			summary_count (query_count (db, "SELECT COUNT(*) FROM source"), b1));
	summary_line (s, "Stored alignments", "%s",
			summary_count (query_count (db, "SELECT COUNT(*) FROM alignment"), b1));
	summary_line (s, "Clusters", "%s (after all filters)",
			summary_count (num_clusters, b1));

	total = query_count (db, "SELECT COUNT(*) FROM retrocopy");

	summary_section (s, "Retrocopies");
	summary_line (s, "Total", "%s", summary_count (total, b1));

#define SUMMARY_RTC(key, sql) \
	do { \
		n = query_count (db, sql); \
		summary_line (s, key, "%s (%s)", summary_count (n, b1), \
				summary_percent (n, total, b2)); \
	} while (0)

	SUMMARY_RTC ("Single parental gene",
			"SELECT COUNT(*) FROM retrocopy WHERE level = 1");
	SUMMARY_RTC ("Precise",
			"SELECT COUNT(*) FROM retrocopy WHERE insertion_point_type = 2");
	SUMMARY_RTC ("Both breakpoints",
			"SELECT COUNT(*) FROM retrocopy\n"
			"WHERE breakpoint_left_sr > 0 AND breakpoint_right_sr > 0");
	SUMMARY_RTC ("Target site duplication",
			"SELECT COUNT(*) FROM retrocopy\n"
			"WHERE breakpoint_left_sr > 0 AND breakpoint_right_sr > 0\n"
			"	AND breakpoint_left >= breakpoint_right");
	SUMMARY_RTC ("Poly-A tail",
			"SELECT COUNT(*) FROM retrocopy WHERE polya_sr > 0");
	SUMMARY_RTC ("Exon-exon junctions",
			"SELECT COUNT(DISTINCT retrocopy_id) FROM junction");

	if (mc->fasta_file != NULL)
		{
			SUMMARY_RTC ("Predicted sequence",
					"SELECT COUNT(*) FROM insertion WHERE seq IS NOT NULL");
		}
	else
		summary_line (s, "Predicted sequence", "- (pass --reference-file)");

#undef SUMMARY_RTC

	// Retrocopies with reads supporting the insertion by sample
	summary_section (s, "Retrocopies by sample (with alternate reads)");

	stmt = db_prepare (db,
			"SELECT s.path, COUNT(g.retrocopy_id)\n"
			"FROM source AS s\n"
			"LEFT JOIN genotype AS g\n"
			"	ON g.source_id = s.id AND g.alternate_depth > 0\n"
			"GROUP BY s.id\n"
			"ORDER BY s.id ASC");

	while (db_step (stmt) == SQLITE_ROW)
		{
			name = path_file (db_column_text (stmt, 0), 1);
			summary_line (s, name, "%s",
					summary_count (db_column_int64 (stmt, 1), b1));
			xfree (name);
		}

	db_finalize (stmt);

	summary_section (s, "Output");
	summary_line (s, "Database", "%s", db_file);
	summary_line (s, "Next step", "%s make-vcf%s%s %s", PACKAGE,
			mc->fasta_file != NULL ? " -r " : "",
			mc->fasta_file != NULL ? mc->fasta_file : "",
			db_file);

	summary_log (s);
	summary_free (s);
}

static void
run (MergeCall *mc)
{
	log_trace ("Inside %s", __func__);

	sqlite3 *db = NULL;

	sqlite3_stmt *cluster_stmt = NULL;
	sqlite3_stmt *clustering_stmt = NULL;
	sqlite3_stmt *blacklist_stmt = NULL;
	sqlite3_stmt *overlapping_blacklist_stmt = NULL;
	sqlite3_stmt *retrocopy_stmt = NULL;
	sqlite3_stmt *junction_stmt = NULL;
	sqlite3_stmt *insertion_stmt = NULL;
	sqlite3_stmt *evidence_stmt = NULL;
	sqlite3_stmt *cluster_merging_stmt = NULL;
	sqlite3_stmt *genotype_stmt = NULL;

	Blacklist *blacklist = NULL;

	int num_clusters = 0;
	char *db_file = NULL;
	time_t start = time (NULL);

	// Merge all databases
	if (mc->output_file != NULL)
		{
			// Use first database
			db_file = xstrdup (mc->output_file);

			log_info ("Open the database '%s'", db_file);

			// Connect to database
			db = db_connect (db_file);
		}
	else
		{
			// Assemble database output filename
			xasprintf_concat (&db_file, "%s/%s.db",
					mc->output_dir, mc->prefix);

			log_debug ("Create output dir '%s'", mc->output_dir);
			mkdir_p (mc->output_dir);

			log_info ("Create the database '%s'", db_file);

			// Create a new database
			db = db_create (db_file);
		}

	// Increase the cache size
	db_cache_size (db, mc->cache_size);

	// Create cluster statement
	cluster_stmt = db_prepare_cluster_stmt (db);

	// Create clustering statement
	clustering_stmt = db_prepare_clustering_stmt (db);

	// Create blacklist statement
	blacklist_stmt = db_prepare_blacklist_stmt (db);

	// Create overlapping blacklist statement
	overlapping_blacklist_stmt =
		db_prepare_overlapping_blacklist_stmt (db);

	// Create retrocopy statement
	retrocopy_stmt = db_prepare_retrocopy_stmt (db);

	log_debug ("Prepare junction statement");
	junction_stmt = db_prepare_junction_stmt (db);

	log_debug ("Prepare insertion statement");
	insertion_stmt = db_prepare_insertion_stmt (db);

	log_debug ("Prepare evidence statement");
	evidence_stmt = db_prepare_evidence_stmt (db);

	// Create cluster merging statement
	cluster_merging_stmt =
		db_prepare_cluster_merging_stmt (db);

	// Create genotype statement
	genotype_stmt = db_prepare_genotype_stmt (db);

	// Index blacklisted regions
	blacklist = blacklist_new (blacklist_stmt,
			overlapping_blacklist_stmt, mc->cs);

	// If the user passed blacklisted regions
	if (mc->blacklist_region != NULL)
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			if (gff_looks_like_gff_file (mc->blacklist_region))
				{
					log_info ("Index the blacklisted regions of '%s'",
							mc->blacklist_region);
					blacklist_index_dump_from_gff (blacklist,
							mc->blacklist_region, mc->filter, mc->gene_name_key);
				}
			else
				{
					log_info ("Index the blacklisted regions of '%s'",
							mc->blacklist_region);
					blacklist_index_dump_from_bed (blacklist,
							mc->blacklist_region);
				}

			// Commit
			db_end_transaction (db);
		}

	// If there are files to merge with ...
	if (array_len (mc->db_files))
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			// Time to merge them all!
			log_info ("Merge the databases into '%s'", db_file);
			mc->num_merged = db_merge (db, array_len (mc->db_files),
					(char **) array_data (mc->db_files));

			// Commit
			db_end_transaction (db);
		}

	// Begin transaction to speed up
	db_begin_transaction (db);

	// Clustering
	log_info ("Run the clustering step");
	num_clusters = cluster (cluster_stmt, clustering_stmt,
			mc->epsilon, mc->min_pts, mc->parental_dist, mc->support,
			mc->blacklist_chr, blacklist, mc->padding);

	// Commit
	db_end_transaction (db);

	if (!num_clusters)
		{
			log_warn ("No cluster found: There is no retrocopy to annotate");
		}
	else
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			// Filtering and Annotation
			log_info ("Run the annotation step");
			retrocopy (retrocopy_stmt, cluster_merging_stmt,
					junction_stmt, insertion_stmt, evidence_stmt, mc->fasta_file,
					mc->near_gene_rank);

			// Commit
			db_end_transaction (db);

			// Begin transaction to speed up
			db_begin_transaction (db);

			// Genotyping
			log_info ("Run the genotyping step");
			genotype (genotype_stmt, mc->threads, mc->phred_quality);

			// Commit
			db_end_transaction (db);

			print_summary (mc, db, db_file, num_clusters, start);
		}

	// Cleanup
	xfree (db_file);
	db_finalize (cluster_stmt);
	db_finalize (clustering_stmt);
	db_finalize (blacklist_stmt);
	db_finalize (overlapping_blacklist_stmt);
	db_finalize (retrocopy_stmt);
	db_finalize (junction_stmt);
	db_finalize (insertion_stmt);
	db_finalize (evidence_stmt);
	db_finalize (cluster_merging_stmt);
	db_finalize (genotype_stmt);
	db_close (db);

	blacklist_free (blacklist);
}

static void
print_usage (FILE *fp)
{
	int pkg_len = strlen (PACKAGE);
	fprintf (fp,
		"%s\n"
		"\n"
		"Usage: %s merge-call [-h] [-q] [-d] [-l FILE] [-o DIR] [-p STR] [-I]\n"
		"       %*c            [-c INT] [-e INT] [-m INT] [-b STR] [-B FILE]\n"
		"       %*c            [-P INT] [-T STR] [-A KEY=VALUE] [-N STR]\n"
		"       %*c            [-x INT] [-g INT] [-n INT] [-r FILE]\n"
		"       %*c            [-t INT] [-Q INT] [-i FILE] <FILE> ...\n"
		"\n"
		"Discover and annotate retrocopy insertions\n"
		"\n"
		"Examples:\n"
		"   $ sider mc -l mc.log -I -r hg38.fa in.db\n"
		"   $ sider mc -e 500 -m 20 -p in1_plus_in2 in1.db in2.db\n"
		"   $ sider mc -b chrY -B gencode.gff3 -t 10 -Q 20 -i list.txt\n"
		"\n"
		"Output:\n"
		"   The annotated retrocopies, at the first input database (with '-I')\n"
		"   or at a new one with all the inputs merged. Next, run 'make-vcf'\n"
		"\n"
		"Arguments:\n"
		"   One or more SQLite3 databases generated by 'process-sample'\n"
		"\n"
		"Mandatory Options:\n"
		"   -i, --input-file           File with a newline separated list of SQLite3\n"
		"                              databases. Not mandatory if the databases are\n"
		"                              passed as arguments. If both are given, their\n"
		"                              union is used\n"
		"\n"
		"Input/Output Options:\n"
		"   -h, --help                 Show help options\n"
		"   -q, --quiet                Decrease verbosity to error messages only,\n"
		"                              or suppress terminal output at all if\n"
		"                              'log-file' is passed\n"
		"       --silent               Same as '--quiet'\n"
		"   -d, --debug                Increase verbosity to debug level\n"
		"   -l, --log-file             Print log messages to a file\n"
		"   -o, --output-dir           Output directory. It is created if it does not\n"
		"                              exist [default:\"%s\"]\n"
		"   -p, --prefix               Prefix of the output files [default:\"%s\"]\n"
		"   -I, --in-place             Merge all databases into the first one, instead\n"
		"                              of creating a new database\n"
		"\n"
		"SQLite3 Options:\n"
		"   -c, --cache-size           SQLite3 cache size in KiB [default:\"%d\"]\n"
		"\n"
		"Clustering Options:\n"
		"   -e, --epsilon              DBSCAN: Maximum distance between two alignments\n"
		"                              inside a cluster [default:\"%d\"]\n"
		"   -m, --min-pts              DBSCAN: Minimum number of alignments required to\n"
		"                              form a dense region [default:\"%d\"]\n"
		"\n"
		"Blacklist Options:\n"
		"   -b, --blacklist-chr        Do not cluster from or to this chromosome. It may\n"
		"                              be passed multiple times [default:\"%s\"]\n"
		"   -B, --blacklist-region     Regions where clusters are discarded, in GTF/GFF3\n"
		"                              or BED format. A GTF/GFF3 file is filtered by\n"
		"                              '-T' and '-A'\n"
		"   -P, --blacklist-padding    Extend the blacklisted regions by N bases at both\n"
		"                              sides [default:\"%d\"]\n"
		"   -T, --gff-feature          The 'feature' (3rd column) of the blacklisted\n"
		"                              regions [default:\"%s\"]\n"
		"   -A, --gff-attribute        An 'attribute' (9th column) to filter in, as\n"
		"                              key=value. Alternatives separated by '|' must\n"
		"                              match any (OR): 'k=v1|v2' or 'k1=v1|k2=v2'. If\n"
		"                              passed multiple times, all of them must match\n"
		"                              (AND). Values match exactly\n"
		"                              [default:\"%s\"]\n"
		"   -N, --gff-gene-name        The attribute naming the blacklisted regions\n"
		"                              [default:\"%s\"]\n"
		"\n"
		"Annotation Options:\n"
		"   -x, --parental-distance    Minimum distance between a cluster and its\n"
		"                              parental gene. Closer clusters are discarded\n"
		"                              [default:\"%d\"]\n"
		"   -g, --genotype-support     Minimum number of reads from a single sample\n"
		"                              (SAM/BAM/CRAM) to keep a cluster [default:\"%d\"]\n"
		"   -n, --near-gene-rank       Maximum distance, in number of genes, between\n"
		"                              the parental genes of overlapping clusters to\n"
		"                              merge them as near parental genes\n"
		"                              [default:\"%d\"]\n"
		"   -r, --reference-file       FASTA file of the reference genome. With it,\n"
		"                              the retrocopy ends are searched in the parental\n"
		"                              exons and its sequence is predicted\n"
		"\n"
		"Genotyping Options:\n"
		"   -t, --threads              Number of threads [default:\"%d\"]\n"
		"   -Q, --phred-quality        Minimum mapping quality (MAPQ) of the reads\n"
		"                              supporting the reference allele [default:\"%d\"]\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ', pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_CACHE_SIZE, DEFAULT_EPS, DEFAULT_MIN_PTS,
		DEFAULT_BLACKLIST_CHR, DEFAULT_BLACKLIST_PADDING, DEFAULT_GFF_FEATURE, DEFAULT_GFF_ATTRIBUTE,
		DEFAULT_GFF_GENE_NAME, DEFAULT_PARENTAL_DISTANCE, DEFAULT_SUPPORT, DEFAULT_NEAR_GENE_RANK,
		DEFAULT_THREADS, DEFAULT_PHRED_QUALITY);
}

static void
print_try_help (FILE *fp)
{
	fprintf (fp, "Try '%s merge-call --help' for more information\n",
			PACKAGE);
}

static void
merge_call_init (MergeCall *mc)
{
	*mc = (MergeCall) {
		.db_files         = array_new (xfree),
		.input_file       = NULL,
		.output_dir       = DEFAULT_OUTPUT_DIR,
		.output_file      = NULL,
		.prefix           = DEFAULT_PREFIX,
		.in_place         = DEFAULT_IN_PLACE,
		.logger           = NULL,
		.log_file         = NULL,
		.log_level        = DEFAULT_LOG_LEVEL,
		.silent           = DEFAULT_LOG_SILENT,
		.cache_size       = DEFAULT_CACHE_SIZE,
		.epsilon          = DEFAULT_EPS,
		.min_pts          = DEFAULT_MIN_PTS,
		.blacklist_region = NULL,
		.blacklist_chr    = set_new_full (str_hash, str_equal, NULL),
		.filter           = gff_filter_new (),
		.attributes       = list_new (xfree),
		.gene_name_key    = DEFAULT_GFF_GENE_NAME,
		.cs               = chr_std_new (),
		.padding          = DEFAULT_BLACKLIST_PADDING,
		.parental_dist    = DEFAULT_PARENTAL_DISTANCE,
		.support          = DEFAULT_SUPPORT,
		.near_gene_rank   = DEFAULT_NEAR_GENE_RANK,
		.fasta_file       = NULL,
		.threads          = DEFAULT_THREADS,
		.phred_quality    = DEFAULT_PHRED_QUALITY
	};
}

static void
merge_call_destroy (MergeCall *mc)
{
	if (mc == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	array_free (mc->db_files, 1);
	set_free (mc->blacklist_chr);
	gff_filter_free (mc->filter);
	list_free (mc->attributes);
	chr_std_free (mc->cs);
	logger_free (mc->logger);
	xfree ((void *) mc->output_file);

	memset (mc, 0, sizeof (MergeCall));
}

static int
merge_call_validate (MergeCall *mc)
{
	int rc = EXIT_SUCCESS;
	int i = 0;

	/*Validate arguments and mandatory options*/

	// If no one file was passed, throw an error
	if (array_len (mc->db_files) == 0)
		{
			fprintf (stderr, "%s: Missing SQLite3 databases\n", PACKAGE);
			print_try_help (stderr);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Test if all alignment files exist
	for (i = 0; i < array_len (mc->db_files); i++)
		{
			const char *db_file = array_get (mc->db_files, i);
			if (!exists (db_file))
				{
					fprintf (stderr, "%s: SQLite3 database '%s': No such file\n", PACKAGE, db_file);
					rc = EXIT_FAILURE; goto Exit;
				}
		}

	/*Validate options*/

	// Validate blacklist_region file
	if (mc->blacklist_region != NULL && !exists (mc->blacklist_region))
		{
			fprintf (stderr, "%s: --blacklist-region '%s': No such file\n", PACKAGE,
					mc->blacklist_region);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate cache_size >= DEFAULT_CACHE_SIZE
	if (mc->cache_size < DEFAULT_CACHE_SIZE)
		{
			fprintf (stderr, "%s: --cache-size must be greater or equal to %uKiB\n",
					PACKAGE, DEFAULT_CACHE_SIZE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->epsilon < 0)
		{
			fprintf (stderr, "%s: --epsilon must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->min_pts < 3)
		{
			fprintf (stderr, "%s: --min-pts must be greater than 2\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->parental_dist < 0)
		{
			fprintf (stderr, "%s: --parental-distance must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->support < 0)
		{
			fprintf (stderr, "%s: --genotype-support must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->padding < 0)
		{
			fprintf (stderr, "%s: --blacklist-padding must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->near_gene_rank < 1)
		{
			fprintf (stderr, "%s: --near-gene-rank must be greater or equal to 1\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->fasta_file != NULL && !exists (mc->fasta_file))
		{
			fprintf (stderr, "%s: --reference-file '%s': No such file\n", PACKAGE,
					mc->fasta_file);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->phred_quality < 0)
		{
			fprintf (stderr, "%s: --phred-quality must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->threads < 1)
		{
			fprintf (stderr, "%s: --threads must be greater or equal to 1\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*Final settings*/

	// Add default blacklisted chr if none
	// has been passed
	if (set_size (mc->blacklist_chr) == 0)
		set_insert (mc->blacklist_chr,
				chr_std_lookup (mc->cs, DEFAULT_BLACKLIST_CHR));

	// If gff_feature was not set
	if (mc->filter->feature == NULL)
		gff_filter_insert_feature (mc->filter, DEFAULT_GFF_FEATURE);

	// If gff_attribute was not set
	if (gff_filter_attribute_size (mc->filter) == 0)
		{
			list_append (mc->attributes, xstrdup (DEFAULT_GFF_ATTRIBUTE));
			gff_filter_insert_attribute (mc->filter, DEFAULT_GFF_ATTRIBUTE);
		}

	// Copy the name of possible output file, if the
	// user chose --in-place
	if (mc->in_place)
		{
			// Copy first file to use it as the final
			// database
			mc->output_file = xstrdup (array_get (mc->db_files, 0));
			i = 0;

			// Remove all repetitive output_file from list
			while (array_find_with_equal_fun (mc->db_files,
						mc->output_file, equalstring, &i))
				array_remove_index (mc->db_files, i);
		}

	// Avoid to include repetitive files
	array_uniq (mc->db_files, cmpstringp);

	// If it's silent and no log file
	// was passed, then set log_level
	// to LOG_ERROR - At least print
	// errors
	if (mc->silent && mc->log_file == NULL)
		{
			mc->silent = 0;
			mc->log_level = LOG_ERROR;
		}

	mc->logger = logger_new (mc->log_file,
			mc->log_level, mc->silent, 1);

Exit:
	return rc;
}

static void
merge_call_print (const MergeCall *mc)
{
	String *msg = NULL;
	ListElmt *cur = NULL;
	int i = 0;

	msg = string_sized_new (BUFSIZ);

	string_concat_printf (msg, ">> Merge Call step <<\n"
		"\n"
		"#\n"
		"# %s\n"
		"#\n"
		"\n"
		"## Command line parsing with default values\n"
		"\n"
		"# Input SQLite3 databases\n"
		"$ cat my-inputfile.txt\n",
		PACKAGE_STRING);

	if (mc->output_file != NULL)
		string_concat_printf (msg, "%s\n", mc->output_file);

	for (i = 0; i < array_len (mc->db_files); i++)
		string_concat_printf (msg, "%s\n",
				(char *) array_get (mc->db_files, i));

	string_concat_printf (msg,
		"\n"
		"# Run %s\n"
		"$ %s merge-call\n"
		"  --input-file='my-inputfile.txt' \\\n"
		"  --output-dir='%s' \\\n"
		"  --prefix='%s' \\\n",
		PACKAGE, PACKAGE, mc->output_dir, mc->prefix);

	if (mc->log_file != NULL)
		string_concat_printf (msg, "  --log-file='%s' \\\n",
				mc->log_file);

	if (mc->log_level <= LOG_DEBUG)
		string_concat_printf (msg, "  --debug \\\n");

	if (mc->silent)
		string_concat_printf (msg, "  --silent \\\n");

	if (mc->in_place)
		string_concat_printf (msg, "  --in-place \\\n");

	string_concat_printf (msg,
		"  --cache-size=%d \\\n"
		"  --epsilon=%d \\\n"
		"  --min-pts=%d \\\n",
		mc->cache_size, mc->epsilon, mc->min_pts);

	cur = list_head (set_list (mc->blacklist_chr));
	for (; cur != NULL; cur = list_next (cur))
		string_concat_printf (msg, "  --blacklist-chr='%s' \\\n",
				(char *) list_data (cur));

	if (mc->blacklist_region != NULL)
		{
			string_concat_printf (msg,
				"  --blacklist-region='%s' \\\n"
				"  --blacklist-padding=%d \\\n",
				mc->blacklist_region, mc->padding);

			if (gff_looks_like_gff_file (mc->blacklist_region))
				{
					string_concat_printf (msg, "  --gff-feature='%s' \\\n",
							mc->filter->feature);

					cur = list_head (mc->attributes);
					for (; cur != NULL; cur = list_next (cur))
						string_concat_printf (msg, "  --gff-attribute='%s' \\\n",
								(char *) list_data (cur));

					string_concat_printf (msg, "  --gff-gene-name='%s' \\\n",
							mc->gene_name_key);
				}
		}

	if (mc->fasta_file != NULL)
		string_concat_printf (msg,
			"  --reference-file=%s \\\n",
			mc->fasta_file);

	string_concat_printf (msg,
		"  --parental-distance=%d \\\n"
		"  --genotype-support=%d \\\n"
		"  --near-gene-rank=%d \\\n"
		"  --threads=%d \\\n"
		"  --phred-quality=%d\n",
		mc->parental_dist, mc->support, mc->near_gene_rank,
		mc->threads, mc->phred_quality);

	log_info ("%s", msg->str);
	string_free (msg, 1);
}

int
parse_merge_call_command_opt (int argc, char **argv)
{
	assert (argc > 1 && argv != NULL && *argv != NULL);

	// No options or arguments
	// Print usage
	if (argc == 2)
		{
			print_usage (stdout);
			return EXIT_SUCCESS;
		}

	struct option opt[] =
	{
		{"help",               no_argument,       0, 'h'},
		{"quiet",              no_argument,       0, 'q'},
		{"silent",             no_argument,       0, 'q'},
		{"debug",              no_argument,       0, 'd'},
		{"in-place",           no_argument,       0, 'I'},
		{"log-file",           required_argument, 0, 'l'},
		{"output-dir",         required_argument, 0, 'o'},
		{"prefix",             required_argument, 0, 'p'},
		{"cache-size",         required_argument, 0, 'c'},
		{"input-file",         required_argument, 0, 'i'},
		{"epsilon",            required_argument, 0, 'e'},
		{"min-pts",            required_argument, 0, 'm'},
		{"parental-distance",  required_argument, 0, 'x'},
		{"genotype-support",   required_argument, 0, 'g'},
		{"blacklist-chr",      required_argument, 0, 'b'},
		{"blacklist-region",   required_argument, 0, 'B'},
		{"blacklist-padding",  required_argument, 0, 'P'},
		{"gff-feature",        required_argument, 0, 'T'},
		{"gff-attribute",      required_argument, 0, 'A'},
		{"gff-gene-name",      required_argument, 0, 'N'},
		{"near-gene-distance", required_argument, 0, 'n'},
		{"reference-file",     required_argument, 0, 'r'},
		{"phred-quality",      required_argument, 0, 'Q'},
		{"threads",            required_argument, 0, 't'},
		{0,                    0,                 0,  0 }
	};

	// Init variables to default values
	MergeCall mc = {};
	merge_call_init (&mc);

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c, i;

	while ((c = getopt_long (argc, argv, "hqdIl:o:p:c:e:m:b:B:P:T:A:N:x:g:n:r:Q:t:i:", opt, &option_index)) >= 0)
		{
			switch (c)
				{
				case 'h':
					{
						print_usage (stdout);
						goto Exit;
						break;
					}
				case 'q':
					{
						mc.silent = 1;
						break;
					}
				case 'd':
					{
						mc.log_level = LOG_DEBUG;
						break;
					}
				case 'I':
					{
						mc.in_place = 1;
						break;
					}
				case 'l':
					{
						mc.log_file = optarg;
						break;
					}
				case 'o':
					{
						mc.output_dir = optarg;
						break;
					}
				case 'p':
					{
						mc.prefix = optarg;
						break;
					}
				case 'c':
					{
						mc.cache_size = atoi (optarg);
						break;
					}
				case 'e':
					{
						mc.epsilon = atoi (optarg);
						break;
					}
				case 'm':
					{
						mc.min_pts = atoi (optarg);
						break;
					}
				case 'b':
					{
						set_insert (mc.blacklist_chr,
								chr_std_lookup (mc.cs, optarg));
						break;
					}
				case 'x':
					{
						mc.parental_dist = atoi (optarg);
						break;
					}
				case 'g':
					{
						mc.support = atoi (optarg);
						break;
					}
				case 'n':
					{
						mc.near_gene_rank = atoi (optarg);
						break;
					}
				case 'r':
					{
						mc.fasta_file = optarg;
						break;
					}
				case 'Q':
					{
						mc.phred_quality = atoi (optarg);
						break;
					}
				case 't':
					{
						mc.threads = atoi (optarg);
						break;
					}
				case 'B':
					{
						mc.blacklist_region = optarg;
						break;
					}
				case 'P':
					{
						mc.padding = atoi (optarg);
						break;
					}
				case 'T':
					{
						gff_filter_insert_feature (mc.filter, optarg);
						break;
					}
				case 'A':
					{
						if (!gff_filter_insert_attribute (mc.filter, optarg))
							{
								fprintf (stderr,
										"%s: --gff-attribute '%s': Expected "
										"KEY=VALUE[|VALUE|KEY=VALUE...]\n", PACKAGE, optarg);
								print_try_help (stderr);
								rc = EXIT_FAILURE; goto Exit;
							}

						list_append (mc.attributes, xstrdup (optarg));
						break;
					}
				case 'N':
					{
						mc.gene_name_key = optarg;
						break;
					}
				case 'i':
					{
						mc.input_file = optarg;
						break;
					}
				case '?':
				case ':':
					{
						print_try_help (stderr);
						rc = EXIT_FAILURE; goto Exit;
						break;
					}
				}
		}

	// Catch all alignment files passed
	// as argument
	for (i = optind + 1; i < argc; i++)
		array_add (mc.db_files, xstrdup (argv[i]));

	// Catch all alignment files passed into
	// --input-file
	if (mc.input_file != NULL)
		read_file_lines (mc.db_files, mc.input_file);

	// Validate and init logger
	rc = merge_call_validate (&mc);

	// If no error
	if (rc == EXIT_SUCCESS)
		{
			// RUN FOOLS
			merge_call_print (&mc);
			run (&mc);
		}

Exit:
	merge_call_destroy (&mc);
	return rc;
}
