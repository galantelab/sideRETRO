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
	char *db_file = NULL;
	int i = 0;

	// Assemble database output filename
	xasprintf_concat (&db_file, "%s/%s.db",
			ps->output_dir, ps->prefix);

	log_info ("Create output dir '%s'", ps->output_dir);
	mkdir_p (ps->output_dir);

	// Create and connect to database
	log_info ("Create and connect to database '%s'", db_file);
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
	log_info ("Index annotation file '%s'", ps->gff_file);
	exon_tree = exon_tree_new (exon_stmt, overlapping_stmt, cs);
	exon_tree_index_dump (exon_tree, ps->gff_file);

	log_info ("Create thread pool");
	thpool = thpool_init (ps->threads);

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
				.max_base_freq    = ps->max_base_freq
			};

			log_debug ("Dump source entry '%s'", sam_file);
			db_insert_source (source_stmt, i + 1, batch_id, sam_file);

			log_info ("Run abnormal filter for '%s'", sam_file);
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
			log_info ("Run deduplication step for '%s'", db_file);
			dedup (db);

			// Commit database
			db_end_transaction (db);
		}

	log_info ("Process Sample at '%s' is finished. "
		"Run merge-call command to discover somatic retrocopies",
		db_file);

	// Time to clean
	db_finalize (exon_stmt);
	db_finalize (batch_stmt);
	db_finalize (source_stmt);
	db_finalize (alignment_stmt);
	db_finalize (overlapping_stmt);
	db_close (db);

	xfree (db_file);
	xfree (ab_args);

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
		"Usage: %s process-sample [-h] [-q] [-d] [-s] [-l FILE] [-o DIR]\n"
		"       %*c                [-p STR] [-t INT] [-c INT] [-Q INT]\n"
		"       %*c                [-m INT] [-f FLOAT] [-F FLOAT | -r]\n"
		"       %*c                [-D] [-M FLOAT] [-e] [-i FILE]\n"
		"       %*c                -a FILE <FILE> ...\n"
		"\n"
		"Extract alignments related to event of retrocopy\n"
		"\n"
		"Examples:\n"
		"   $ sider ps -l ps.log -a gencode.gff3.gz -o result in.bam\n"
		"   $ sider ps -t 3 -a gencode.gtf in1.bam in2.sam in3.bam\n"
		"   $ sider ps -t 5 -m 15000 -Q 20 -F 0.9 -a exon.gtf -i list.txt\n"
		"\n"
		"Output:\n"
		"   A SQLite3 database that can be processed at 'merge-call' step\n"
		"\n"
		"Arguments:\n"
		"   One or more alignment file in SAM/BAM/CRAM format\n"
		"\n"
		"Mandatory Options:\n"
		"   -a, --annotation-file   Gene annotation on the reference genome\n"
		"                           in GTF/GFF3 format\n"
		"   -i, --input-file        File containing a newline separated list of\n"
		"                           alignment files in SAM/BAM/CRAM format.\n"
		"                           This option is not mandatory if one or more\n"
		"                           SAM/BAM/CRAM files are passed as argument.\n"
		"                           If 'input-file' and arguments are set\n"
		"                           concomitantly, then the union of all alignment\n"
		"                           files is used\n"
		"\n"
		"Input/Output Options:\n"
		"   -h, --help              Show help options\n"
		"   -q, --quiet             Decrease verbosity to error messages only\n"
		"                           or suppress terminal outputs at all if\n"
		"                           'log-file' is passed\n"
		"       --silent            Same as '--quiet'\n"
		"   -d, --debug             Increase verbosity to debug level\n"
		"   -l, --log-file          Print log messages to a file\n"
		"   -o, --output-dir        Output directory. Create the directory if it does\n"
		"                           not exist [default:\"%s\"]\n"
		"   -p, --prefix            Prefix output files [default:\"%s\"]\n"
		"\n"
		"SQLite3 Options:\n"
		"   -c, --cache-size        Set SQLite3 cache size in KiB [default:\"%d\"]\n"
		"\n"
		"Read Quality Options:\n"
		"   -Q, --phred-quality     Minimum mapping quality of the reads required\n"
		"                           [default:\"%d\"]\n"
		"   -M, --max-base-freq     Maximum base frequency fraction allowed\n"
		"                           [default:\"%.2f\"]\n"
		"   -D, --deduplicate       Remove duplicated reads. Reads are considered\n"
		"                           duplicates when they share the 5 prime positions\n"
		"                           of both reads and read-pairs\n"
		"\n"
		"Processing Options:\n"
		"   -s, --sorted            Assume all reads are grouped by queryname, even if\n"
		"                           there is no SAM/BAM/CRAM header tag 'SO:queryname'\n"
		"   -t, --threads           Number of threads [default:\"%d\"]\n"
		"   -m, --max-distance      Maximum distance allowed between paired-end reads\n"
		"                           [default:\"%d\"]\n"
		"   -f, --exon-frac         Minimum overlap required as a fraction of exon\n"
		"                           [default:\"%.0e\"; 1 base]\n"
		"   -F, --alignment-frac    Minimum overlap required as a fraction of\n"
		"                           alignment [default:\"%.0e\"; 1 base]\n"
		"   -e, --either            The minimum fraction must be satisfied for at least\n"
		"                           exon OR alignment. Without '-e', both fractions would\n"
		"                           have to be satisfied\n"
		"   -r, --reciprocal        The fraction overlap must be reciprocal for exon and\n"
		"                           alignment. If '-f' is 0.5, then '-F' will be set to\n"
		"                           0.5 as well\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ', pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_CACHE_SIZE, DEFAULT_PHRED_QUALITY,
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
		{0,                 0,                 0,  0 }
	};

	// Init variables to default values
	ProcessSample ps = {};
	process_sample_init (&ps);

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c, i;

	while ((c = getopt_long (argc, argv, "hqdsDl:a:o:p:t:m:M:c:Q:f:F:eri:", opt, &option_index)) >= 0)
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
