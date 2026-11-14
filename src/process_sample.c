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
#include <time.h>
#include <getopt.h>
#include <assert.h>
#include "wrapper.h"
#include "array.h"
#include "utils.h"
#include "chr.h"
#include "db.h"
#include "log.h"
#include "logger.h"
#include "io.h"
#include "str.h"
#include "thpool.h"
#include "exon.h"
#include "summary.h"
#include "list.h"
#include "bed.h"
#include "gff.h"
#include "abnormal.h"
#include "dedup.h"
#include "process_sample.h"

#define DEFAULT_MAX_DISTANCE    10000
#define DEFAULT_CACHE_SIZE      200000 /* 200MiB */
#define DEFAULT_THREADS         1
#define DEFAULT_SORTED          0
#define DEFAULT_DEDUPLICATE     0
#define DEFAULT_EXON_FRAC       1e-09
#define DEFAULT_ALIGNMENT_FRAC  1e-09
#define DEFAULT_EITHER          0
#define DEFAULT_RECIPROCAL      0
#define DEFAULT_PREFIX          "out"
#define DEFAULT_OUTPUT_DIR      "."
#define DEFAULT_LOG_SILENT      0
#define DEFAULT_LOG_LEVEL       LOG_INFO
#define DEFAULT_PHRED_QUALITY   8
#define DEFAULT_MAX_BASE_FREQ   0.75

struct _ProcessSample
{
	// Mandatory
	Array       *sam_files;
	const char  *gff_file;

	// Annotation
	GffFilter   *filter;
	List        *attributes;
	const char  *gene_name_key;
	const char  *exon_id_key;

	// I/O
	const char  *input_file;
	const char  *output_dir;
	const char  *prefix;

	// Log
	Logger      *logger;
	const char  *log_file;
	int          log_level;
	int          silent;

	// SQLite3
	int          cache_size;

	// Read Quality
	float        max_base_freq;
	int          phred_quality;
	int          deduplicate;

	// Processing
	int          threads;
	int          sorted;
	int          max_distance;
	float        exon_frac;
	float        alignment_frac;
	int          alignment_frac_set;
	int          reciprocal;
	int          either;
};

typedef struct _ProcessSample ProcessSample;

static long
count_alignments (sqlite3 *db)
{
	sqlite3_stmt *stmt = NULL;
	long n = 0;

	stmt = db_prepare (db, "SELECT COUNT(*) FROM alignment");

	if (db_step (stmt) == SQLITE_ROW)
		n = db_column_int64 (stmt, 0);

	db_finalize (stmt);
	return n;
}

static void
print_summary (const ProcessSample *ps, const char *db_file, size_t num_exons,
		const AbnormalStat *stats, int num_files, long num_alignments,
		long num_abnormal, long num_exonic, long num_dedup, time_t start)
{
	Summary *s = NULL;
	char *name = NULL;
	char b1[SUMMARY_NUM_LEN] = {};
	char b2[SUMMARY_NUM_LEN] = {};
	char b3[SUMMARY_NUM_LEN] = {};
	char b4[SUMMARY_NUM_LEN] = {};
	int i = 0;

	s = summary_new ("process-sample summary", start);

	summary_section (s, "Annotation");
	summary_line (s, "File", "%s", ps->gff_file);
	summary_line (s, "Exons indexed", "%s", summary_count (num_exons, b1));

	summary_section (s, "Alignment files");
	summary_text (s, "%-24s %15s %13s %13s", "", "Alignments", "Abnormal", "Exonic");

	for (i = 0; i < num_files; i++)
		{
			name = path_file (array_get (ps->sam_files, i), 0);

			summary_text (s, "%-24.24s %15s %13s %13s (%s)", name,
					summary_count (stats[i].alignments, b1),
					summary_count (stats[i].abnormal, b2),
					summary_count (stats[i].exonic, b3),
					summary_percent (stats[i].exonic, stats[i].abnormal, b4));

			xfree (name);
		}

	if (num_files > 1)
		summary_text (s, "%-24s %15s %13s %13s (%s)", "Total",
				summary_count (num_alignments, b1),
				summary_count (num_abnormal, b2),
				summary_count (num_exonic, b3),
				summary_percent (num_exonic, num_abnormal, b4));

	if (num_dedup >= 0)
		summary_line (s, "After deduplication", "%s abnormal alignments (%s removed)",
				summary_count (num_dedup, b1),
				summary_count (num_abnormal - num_dedup, b2));

	summary_section (s, "Output");
	summary_line (s, "Database", "%s", db_file);
	summary_line (s, "Next step", "%s merge-call %s", PACKAGE, db_file);

	summary_log (s);
	summary_free (s);

	if (num_abnormal == 0)
		log_warn ("No abnormal alignment was found: No retrocopy can be called");
}

static void
run (ProcessSample *ps)
{
	log_trace ("Inside %s", __func__);

	sqlite3 *db = NULL;
	sqlite3_stmt *batch_stmt = NULL;
	sqlite3_stmt *source_stmt = NULL;
	sqlite3_stmt *alignment_stmt = NULL;
	sqlite3_stmt *exon_stmt = NULL;
	sqlite3_stmt *overlapping_stmt = NULL;

	ExonTree *exon_tree = NULL;
	ChrStd *cs = NULL;

	threadpool thpool = NULL;

	const int batch_id = 1;
	char timestamp[32] = {};
	time_t t = 0;
	struct tm *lt = NULL;

	const char *sam_file = NULL;
	const int num_files = array_len (ps->sam_files);
	int pool_threads = 0;
	int hts_threads = 0;
	char *db_file = NULL;
	int i = 0;

	time_t start = time (NULL);
	AbnormalStat *stats = NULL;
	long num_alignments = 0;
	long num_abnormal = 0;
	long num_exonic = 0;
	long num_dedup = -1;

	// Assemble database output filename
	xasprintf_concat (&db_file, "%s/%s.db",
			ps->output_dir, ps->prefix);

	log_debug ("Create output dir '%s'", ps->output_dir);
	mkdir_p (ps->output_dir);

	// Create and connect to database
	log_info ("Create the database '%s'", db_file);
	db = db_create (db_file);
	batch_stmt = db_prepare_batch_stmt (db);
	source_stmt = db_prepare_source_stmt (db);
	alignment_stmt = db_prepare_alignment_stmt (db);
	exon_stmt = db_prepare_exon_stmt (db);
	overlapping_stmt = db_prepare_overlapping_stmt (db);

	// Increase the cache size
	db_cache_size (db, ps->cache_size);

	// Begin transaction to speed up
	db_begin_transaction (db);

	// Get current local datetime as timestamp
	// for batch table
	t = time (NULL);
	lt = localtime (&t);
	timestamp[strftime (timestamp, sizeof (timestamp),
			"%Y-%m-%d %H:%M:%S", lt)] = '\0';

	// Dump batch entry
	log_debug ("Dump batch entry %d => %s", batch_id, timestamp);
	db_insert_batch (batch_stmt, batch_id, timestamp);

	// Get chromosome standardization
	cs = chr_std_new ();

	// Index protein coding genes into the database
	// and its exons into an intervalar tree by
	// chromosome
	log_info ("Index the exons of the annotation '%s'", ps->gff_file);
	exon_tree = exon_tree_new (exon_stmt, overlapping_stmt, cs);
	exon_tree_index_dump (exon_tree, ps->gff_file, ps->filter,
			ps->gene_name_key, ps->exon_id_key);

	// Each file is processed by one thread. If there are
	// more threads than files, distribute the spare ones
	// among the files to decompress BAM/CRAM
	pool_threads = ps->threads < num_files ? ps->threads : num_files;
	hts_threads = (ps->threads - pool_threads) / num_files;

	log_debug ("Create thread pool with %d threads", pool_threads);
	thpool = thpool_init (pool_threads);

	if (hts_threads > 0)
		log_info ("Use %d extra threads per file to decompress BAM/CRAM",
				hts_threads);

	stats = xcalloc (num_files, sizeof (AbnormalStat));
	AbnormalArg *ab_args = xcalloc (num_files, sizeof (AbnormalArg));

	for (i = 0; i < num_files; i++)
		{
			sam_file = array_get (ps->sam_files, i);

			// Init abnormal_filter arg
			ab_args[i] = (AbnormalArg)
			{
				.tid              = i + 1,
				.inc_step         = num_files,
				.sam_file         = sam_file,
				.either           = ps->either,
				.exon_frac        = ps->exon_frac,
				.alignment_frac   = ps->alignment_frac,
				.exon_tree        = exon_tree,
				.cs               = cs,
				.alignment_stmt   = alignment_stmt,
				.queryname_sorted = ps->sorted,
				.max_distance     = ps->max_distance,
				.phred_quality    = ps->phred_quality,
				.max_base_freq    = ps->max_base_freq,
				.hts_threads      = hts_threads,
				.stat             = &stats[i]
			};

			log_debug ("Dump source entry '%s'", sam_file);
			db_insert_source (source_stmt, i + 1, batch_id, sam_file);

			log_debug ("Queue '%s' to the abnormal filter", sam_file);
			thpool_add_work (thpool, (void *) abnormal_filter,
					(void *) &ab_args[i]);
		}

	// Wait all threads to return
	thpool_wait (thpool);

	// Commit database
	db_end_transaction (db);

	if (ps->deduplicate)
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			// Run deduplication step
			log_info ("Remove duplicated reads from '%s'", db_file);
			dedup (db);

			// Commit database
			db_end_transaction (db);

			num_dedup = count_alignments (db);
		}

	for (i = 0; i < num_files; i++)
		{
			num_alignments += stats[i].alignments;
			num_abnormal += stats[i].abnormal;
			num_exonic += stats[i].exonic;
		}

	print_summary (ps, db_file, hash_size (exon_tree->cache), stats,
			num_files, num_alignments, num_abnormal, num_exonic,
			num_dedup, start);

	// Time to clean
	db_finalize (exon_stmt);
	db_finalize (batch_stmt);
	db_finalize (source_stmt);
	db_finalize (alignment_stmt);
	db_finalize (overlapping_stmt);
	db_close (db);

	xfree (db_file);
	xfree (ab_args);
	xfree (stats);

	chr_std_free (cs);
	exon_tree_free (exon_tree);

	thpool_destroy (thpool);
}

static void
print_usage (FILE *fp)
{
	int pkg_len = strlen (PACKAGE);
	fprintf (fp,
		"%s\n"
		"\n"
		"Usage: %s process-sample [-h] [-q] [-d] [-l FILE] [-o DIR] [-p STR]\n"
		"       %*c                [-T STR] [-A KEY=VALUE] [-N STR] [-E STR]\n"
		"       %*c                [-c INT] [-Q INT] [-M FLOAT] [-D] [-s]\n"
		"       %*c                [-t INT] [-m INT] [-f FLOAT] [-F FLOAT | -r]\n"
		"       %*c                [-e] [-i FILE] -a FILE <FILE> ...\n"
		"\n"
		"Extract the abnormal alignments that may come from retrocopy insertions\n"
		"\n"
		"Examples:\n"
		"   $ sider ps -l ps.log -a gencode.gff3.gz -o result in.bam\n"
		"   $ sider ps -t 3 -a gencode.gtf in1.bam in2.sam in3.cram\n"
		"   $ sider ps -a ensembl.gtf -A transcript_biotype=protein_coding in.bam\n"
		"   $ sider ps -t 5 -m 15000 -Q 20 -F 0.9 -a exons.bed -i list.txt\n"
		"\n"
		"Output:\n"
		"   A SQLite3 database to be processed by 'merge-call'\n"
		"\n"
		"Arguments:\n"
		"   One or more alignment files in SAM/BAM/CRAM format\n"
		"\n"
		"Mandatory Options:\n"
		"   -a, --annotation-file   Exons of the protein coding genes, in GTF/GFF3\n"
		"                           or BED6 format. See 'Annotation Options'\n"
		"   -i, --input-file        File with a newline separated list of alignment\n"
		"                           files in SAM/BAM/CRAM format. Not mandatory if\n"
		"                           the files are passed as arguments. If both are\n"
		"                           given, their union is used\n"
		"\n"
		"Input/Output Options:\n"
		"   -h, --help              Show help options\n"
		"   -q, --quiet             Decrease verbosity to error messages only,\n"
		"                           or suppress terminal output at all if\n"
		"                           'log-file' is passed\n"
		"       --silent            Same as '--quiet'\n"
		"   -d, --debug             Increase verbosity to debug level\n"
		"   -l, --log-file          Print log messages to a file\n"
		"   -o, --output-dir        Output directory. It is created if it does not\n"
		"                           exist [default:\"%s\"]\n"
		"   -p, --prefix            Prefix of the output files [default:\"%s\"]\n"
		"\n"
		"Annotation Options:\n"
		"   A BED6 file (name ending in '.bed' or '.bed.gz') has one exon per\n"
		"   line, the gene name at the 4th column and the strand at the 6th one.\n"
		"   A GTF/GFF3 file is filtered by the options below\n"
		"   -T, --gff-feature       The 'feature' (3rd column) of the exons\n"
		"                           [default:\"%s\"]\n"
		"   -A, --gff-attribute     An 'attribute' (9th column) to filter in, as\n"
		"                           key=value. Alternatives separated by '|' must\n"
		"                           match any (OR): 'k=v1|v2' or 'k1=v1|k2=v2'. If\n"
		"                           passed multiple times, all of them must match\n"
		"                           (AND). Values match exactly\n"
		"                           [default:\"%s\"]\n"
		"   -N, --gff-gene-name     The attribute with the gene name. If missing,\n"
		"                           'gene_id' is used [default:\"%s\"]\n"
		"   -E, --gff-exon-id       The attribute with the exon id. If missing, the\n"
		"                           exon is named 'chr:start_end' [default:\"%s\"]\n"
		"\n"
		"SQLite3 Options:\n"
		"   -c, --cache-size        SQLite3 cache size in KiB [default:\"%d\"]\n"
		"\n"
		"Read Quality Options:\n"
		"   -Q, --phred-quality     Minimum mapping quality (MAPQ) of the reads\n"
		"                           [default:\"%d\"]\n"
		"   -M, --max-base-freq     Maximum frequency of a single base in a read.\n"
		"                           Reads above it are discarded as low complexity\n"
		"                           sequences, such as poly-A [default:\"%.2f\"]\n"
		"   -D, --deduplicate       Remove duplicated reads: Fragments whose reads\n"
		"                           and mates share the same 5' positions\n"
		"\n"
		"Processing Options:\n"
		"   -s, --sorted            Assume the reads are grouped by read name, even\n"
		"                           without the header tag 'SO:queryname'. Grouped\n"
		"                           files are read once, the others twice\n"
		"   -t, --threads           Number of threads. Each file uses one thread and\n"
		"                           the spare ones decompress BAM/CRAM [default:\"%d\"]\n"
		"   -m, --max-distance      Maximum distance between the reads of a pair.\n"
		"                           Farther pairs are abnormal [default:\"%d\"]\n"
		"   -f, --exon-frac         Minimum overlap between an alignment and an\n"
		"                           exon, as a fraction of the exon\n"
		"                           [default:\"%g\"; 1 base]\n"
		"   -F, --alignment-frac    Minimum overlap between an alignment and an\n"
		"                           exon, as a fraction of the alignment\n"
		"                           [default:\"%g\"; 1 base]\n"
		"   -e, --either            Either '-f' OR '-F' must be satisfied. Without\n"
		"                           '-e', both must be\n"
		"   -r, --reciprocal        The overlap must be reciprocal: '-F' is set to\n"
		"                           the value of '-f'\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ', pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, EXON_DEFAULT_FEATURE, EXON_DEFAULT_ATTRIBUTE,
		EXON_DEFAULT_GENE_NAME, EXON_DEFAULT_EXON_ID, DEFAULT_CACHE_SIZE, DEFAULT_PHRED_QUALITY,
		DEFAULT_MAX_BASE_FREQ, DEFAULT_THREADS, DEFAULT_MAX_DISTANCE,
		DEFAULT_EXON_FRAC, DEFAULT_ALIGNMENT_FRAC);
}

static void
print_try_help (FILE *fp)
{
	fprintf (fp, "Try '%s process-sample --help' for more information\n",
			PACKAGE);
}

static void
process_sample_init (ProcessSample *ps)
{
	*ps = (ProcessSample) {
		.sam_files          = array_new (xfree),
		.gff_file           = NULL,
		.filter             = gff_filter_new (),
		.attributes         = list_new (xfree),
		.gene_name_key      = EXON_DEFAULT_GENE_NAME,
		.exon_id_key        = EXON_DEFAULT_EXON_ID,
		.input_file         = NULL,
		.output_dir         = DEFAULT_OUTPUT_DIR,
		.prefix             = DEFAULT_PREFIX,
		.logger             = NULL,
		.log_file           = NULL,
		.log_level          = DEFAULT_LOG_LEVEL,
		.silent             = DEFAULT_LOG_SILENT,
		.cache_size         = DEFAULT_CACHE_SIZE,
		.max_base_freq      = DEFAULT_MAX_BASE_FREQ,
		.phred_quality      = DEFAULT_PHRED_QUALITY,
		.deduplicate        = DEFAULT_DEDUPLICATE,
		.threads            = DEFAULT_THREADS,
		.sorted             = DEFAULT_SORTED,
		.max_distance       = DEFAULT_MAX_DISTANCE,
		.exon_frac          = DEFAULT_EXON_FRAC,
		.alignment_frac     = DEFAULT_ALIGNMENT_FRAC,
		.alignment_frac_set = 0,
		.reciprocal         = DEFAULT_RECIPROCAL,
		.either             = DEFAULT_EITHER
	};
}

static void
process_sample_destroy (ProcessSample *ps)
{
	if (ps == NULL)
		return;

	array_free (ps->sam_files, 1);
	gff_filter_free (ps->filter);
	list_free (ps->attributes);
	logger_free (ps->logger);

	memset (ps, 0, sizeof (ProcessSample));
}

static int
process_sample_validate (ProcessSample *ps)
{
	int rc = EXIT_SUCCESS;
	int i = 0;

	/*
	* Validate arguments and mandatory options
	*/

	// If no one file was passed, throw an error
	if (array_len (ps->sam_files) == 0)
		{
			fprintf (stderr, "%s: Missing alignment files (SAM/BAM/CRAM)\n", PACKAGE);
			print_try_help (stderr);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Test if all alignment files exist
	for (i = 0; i < array_len (ps->sam_files); i++)
		{
			const char *sam_file = array_get (ps->sam_files, i);
			if (!exists (sam_file))
				{
					fprintf (stderr, "%s: alignment file '%s': No such file\n", PACKAGE, sam_file);
					rc = EXIT_FAILURE; goto Exit;
				}
		}

	// If no annotation-file was passed, throw an error
	if (ps->gff_file == NULL)
		{
			fprintf (stderr, "%s: Missing annotation file (GTF/GFF3)\n", PACKAGE);
			print_try_help (stderr);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Test if annotation file exists
	if (!exists (ps->gff_file))
		{
			fprintf (stderr, "%s: annotation file '%s': No such file\n", PACKAGE, ps->gff_file);
			rc = EXIT_FAILURE; goto Exit;
		}

	// GENCODE defaults for the annotation filter
	if (ps->filter->feature == NULL)
		gff_filter_insert_feature (ps->filter, EXON_DEFAULT_FEATURE);

	if (gff_filter_attribute_size (ps->filter) == 0)
		{
			list_append (ps->attributes, xstrdup (EXON_DEFAULT_ATTRIBUTE));
			gff_filter_insert_attribute (ps->filter, EXON_DEFAULT_ATTRIBUTE);
		}

	/*
	* Validate options
	*/

	// If --reciprocal, set alignment_frac to exon_frac
	if (ps->reciprocal)
		{
			if (ps->alignment_frac_set)
				{
					fprintf (stderr, "%s: --reciprocal must be set solely with '-f'\n", PACKAGE);
					rc = EXIT_FAILURE; goto Exit;
				}
			ps->alignment_frac = ps->exon_frac;
		}

	// Validate exon_frac: ]0, 1]
	if (ps->exon_frac > 1 || ps->exon_frac <= 0)
		{
			fprintf (stderr, "%s: --exon-frac must be in the interval of ]0, 1]\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate alignment_frac ]0, 1]
	if (ps->alignment_frac > 1 || ps->alignment_frac <= 0)
		{
			fprintf (stderr, "%s: --alignment-frac must be in the interval of ]0, 1]\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate threads > 0
	if (ps->threads < 1)
		{
			fprintf (stderr, "%s: --threads must be greater or equal to 1\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate cache_size >= DEFAULT_CACHE_SIZE
	if (ps->cache_size < DEFAULT_CACHE_SIZE)
		{
			fprintf (stderr, "%s: --cache-size must be greater or equal to %uKiB\n",
					PACKAGE, DEFAULT_CACHE_SIZE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate max_distance >= 0
	if (ps->max_distance < 0)
		{
			fprintf (stderr, "%s: --max-distance must be a positive value\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (ps->max_base_freq <= 0.25 || ps->max_base_freq > 1)
		{
			fprintf (stderr, "%s: --max-base-frac must be in the range of ]0.25, 1]\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate phred_quality >= 0
	if (ps->phred_quality < 0)
		{
			fprintf (stderr, "%s: --phred-quality must be a positive value\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*
	* Final settings
	*/

	// Avoid to include repetitive files
	array_uniq (ps->sam_files, cmpstringp);

	// If it's silent and no log file
	// was passed, then set log_level
	// to LOG_ERROR - At least print
	// errors
	if (ps->silent && ps->log_file == NULL)
		{
			ps->silent = 0;
			ps->log_level = LOG_ERROR;
		}

	ps->logger = logger_new (ps->log_file, ps->log_level,
			ps->silent, 1);

Exit:
	return rc;
}

static void
process_sample_print (const ProcessSample *ps)
{
	String *msg = string_sized_new (BUFSIZ);
	int i = 0;

	string_concat_printf (msg, ">> Process Sample step <<\n"
		"\n"
		"#\n"
		"# %s\n"
		"#\n"
		"\n"
		"## Command line parsing with default values\n"
		"\n"
		"# Input SAM/BAM/CRAM files\n"
		"$ cat my-inputfile.txt\n",
		PACKAGE_STRING);

	for (i = 0; i < array_len (ps->sam_files); i++)
		string_concat_printf (msg, "%s\n",
				(char *) array_get (ps->sam_files, i));

	string_concat_printf (msg,
		"\n"
		"# Run %s\n"
		"$ %s process-sample\n"
		"  --input-file='my-inputfile.txt' \\\n"
		"  --annotation-file='%s' \\\n"
		"  --output-dir='%s' \\\n"
		"  --prefix='%s' \\\n",
		PACKAGE, PACKAGE, ps->gff_file, ps->output_dir, ps->prefix);

	if (!bed_looks_like_bed_file (ps->gff_file))
		{
			ListElmt *cur = NULL;

			string_concat_printf (msg, "  --gff-feature='%s' \\\n",
					ps->filter->feature);

			for (cur = list_head (ps->attributes); cur != NULL; cur = list_next (cur))
				string_concat_printf (msg, "  --gff-attribute='%s' \\\n",
						(char *) list_data (cur));

			string_concat_printf (msg,
					"  --gff-gene-name='%s' \\\n"
					"  --gff-exon-id='%s' \\\n",
					ps->gene_name_key, ps->exon_id_key);
		}

	if (ps->log_file != NULL)
		string_concat_printf (msg, "  --log-file='%s' \\\n", ps->log_file);

	if (ps->log_level <= LOG_DEBUG)
		string_concat_printf (msg, "  --debug \\\n");

	if (ps->silent)
		string_concat_printf (msg, "  --silent \\\n");

	string_concat_printf (msg,
		"  --cache-size=%d \\\n"
		"  --max-base-freq=%.2f \\\n"
		"  --phred-quality=%d \\\n",
		ps->cache_size, ps->max_base_freq,  ps->phred_quality);

	if (ps->deduplicate)
		string_concat_printf (msg, "  --deduplicate \\\n");

	if (ps->sorted)
		string_concat_printf (msg, "  --sorted \\\n");

	string_concat_printf (msg,
		"  --threads=%d \\\n"
		"  --max-distance=%d \\\n"
		"  --exon-frac=%e \\\n"
		"  --alignment-frac=%e",
		ps->threads, ps->max_distance, ps->exon_frac,
		ps->alignment_frac);

	if (ps->either)
		{
			string_concat_printf (msg, " \\\n");
			string_concat_printf (msg, "  --either");
		}

	string_concat_printf (msg, "\n");

	log_info ("%s", msg->str);
	string_free (msg, 1);
}

int
parse_process_sample_command_opt (int argc, char **argv)
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
		{"help",            no_argument,       0, 'h'},
		{"quiet",           no_argument,       0, 'q'},
		{"silent",          no_argument,       0, 'q'},
		{"debug",           no_argument,       0, 'd'},
		{"log-file",        required_argument, 0, 'l'},
		{"annotation-file", required_argument, 0, 'a'},
		{"output-dir",      required_argument, 0, 'o'},
		{"prefix",          required_argument, 0, 'p'},
		{"threads",         required_argument, 0, 't'},
		{"phred-quality",   required_argument, 0, 'Q'},
		{"max-distance",    required_argument, 0, 'm'},
		{"max-base-freq",   required_argument, 0, 'M'},
		{"cache-size",      required_argument, 0, 'c'},
		{"sorted",          no_argument,       0, 's'},
		{"deduplicate",     no_argument,       0, 'D'},
		{"exon-frac",       required_argument, 0, 'f'},
		{"alignment-frac",  required_argument, 0, 'F'},
		{"either",          no_argument,       0, 'e'},
		{"reciprocal",      no_argument,       0, 'r'},
		{"input-file",      required_argument, 0, 'i'},
		{"gff-feature",     required_argument, 0, 'T'},
		{"gff-attribute",   required_argument, 0, 'A'},
		{"gff-gene-name",   required_argument, 0, 'N'},
		{"gff-exon-id",     required_argument, 0, 'E'},
		{0,                 0,                 0,  0 }
	};

	// Init variables to default values
	ProcessSample ps = {};
	process_sample_init (&ps);

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c, i;

	while ((c = getopt_long (argc, argv, "hqdsDl:a:o:p:t:m:M:c:Q:f:F:eri:T:A:N:E:", opt, &option_index)) >= 0)
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
						ps.silent = 1;
						break;
					}
				case 'd':
					{
						ps.log_level = LOG_DEBUG;
						break;
					}
				case 'l':
					{
						ps.log_file = optarg;
						break;
					}
				case 'a':
					{
						ps.gff_file = optarg;
						break;
					}
				case 'o':
					{
						ps.output_dir = optarg;
						break;
					}
				case 'p':
					{
						ps.prefix = optarg;
						break;
					}
				case 't':
					{
						ps.threads = atoi (optarg);
						break;
					}
				case 'c':
					{
						ps.cache_size = atoi (optarg);
						break;
					}
				case 'Q':
					{
						ps.phred_quality = atoi (optarg);
						break;
					}
				case 'm':
					{
						ps.max_distance = atoi (optarg);
						break;
					}
				case 'M':
					{
						ps.max_base_freq = atof (optarg);
						break;
					}
				case 'D':
					{
						ps.deduplicate = 1;
						break;
					}
				case 's':
					{
						ps.sorted = 1;
						break;
					}
				case 'f':
					{
						ps.exon_frac = atof (optarg);
						break;
					}
				case 'F':
					{
						ps.alignment_frac = atof (optarg);
						ps.alignment_frac_set = 1;
						break;
					}
				case 'e':
					{
						ps.either = 1;
						break;
					}
				case 'r':
					{
						ps.reciprocal = 1;
						break;
					}
				case 'T':
					{
						gff_filter_insert_feature (ps.filter, optarg);
						break;
					}
				case 'A':
					{
						if (!gff_filter_insert_attribute (ps.filter, optarg))
							{
								fprintf (stderr,
										"%s: --gff-attribute '%s': Expected "
										"KEY=VALUE[|VALUE|KEY=VALUE...]\n", PACKAGE, optarg);
								print_try_help (stderr);
								rc = EXIT_FAILURE; goto Exit;
							}

						list_append (ps.attributes, xstrdup (optarg));
						break;
					}
				case 'N':
					{
						ps.gene_name_key = optarg;
						break;
					}
				case 'E':
					{
						ps.exon_id_key = optarg;
						break;
					}
				case 'i':
					{
						ps.input_file = optarg;
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
		array_add (ps.sam_files, xstrdup (argv[i]));

	// Catch all alignment files passed into
	// --input-file
	if (ps.input_file != NULL)
		read_file_lines (ps.sam_files, ps.input_file);

	// Validate and init logger
	rc = process_sample_validate (&ps);

	// If no error
	if (rc == EXIT_SUCCESS)
		{
			/*
			* RUN FOOLS
			*/
			process_sample_print (&ps);
			run (&ps);
		}

Exit:
	process_sample_destroy (&ps);
	return rc;
}
