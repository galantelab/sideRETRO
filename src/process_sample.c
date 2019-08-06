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
#include "thpool.h"
#include "exon.h"
#include "abnormal.h"
#include "process_sample.h"

#define DEFAULT_MAX_DISTANCE    10000
#define DEFAULT_CACHE_SIZE      DB_DEFAULT_CACHE_SIZE
#define DEFAULT_THREADS         1
#define DEFAULT_SORTED          0
#define DEFAULT_EXON_FRAC       1e-09
#define DEFAULT_ALIGNMENT_FRAC  1e-09
#define DEFAULT_EITHER          0
#define DEFAULT_RECIPROCAL      0
#define DEFAULT_PREFIX          "out"
#define DEFAULT_OUTPUT_DIR      "."
#define DEFAULT_LOG_SILENT      0
#define DEFAULT_LOG_LEVEL       LOG_INFO
#define DEFAULT_LOG_FILE        NULL
#define DEFAULT_GFF_FILE        NULL
#define DEFAULT_INPUT_FILE      NULL

static void
process_sample (const char *output_dir, const char *prefix,
		const Array *sam_files, const char *gff_file,
		int threads, int cache_size, int sorted,
		int max_distance, float exon_frac,
		float alignment_frac, int either)
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
	const int num_files = array_len (sam_files);
	char *db_file = NULL;
	int i = 0;

	// Assemble database output filename
	xasprintf_concat (&db_file, "%s/%s.db", output_dir, prefix);

	log_info (">>> Process Sample step <<<");

	log_info ("Create output dir '%s'", output_dir);
	mkdir_p (output_dir);

	// Create and connect to database
	log_info ("Create and connect to database '%s'", db_file);
	db = db_create (db_file);
	batch_stmt = db_prepare_batch_stmt (db);
	source_stmt = db_prepare_source_stmt (db);
	alignment_stmt = db_prepare_alignment_stmt (db);
	exon_stmt = db_prepare_exon_stmt (db);
	overlapping_stmt = db_prepare_overlapping_stmt (db);

	// Increase the cache size
	db_cache_size (db, cache_size);

	// Begin transaction to speed up
	db_begin_transaction (db);

	// Get current local datetome as timestamp
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
	log_info ("Index annotation file '%s'", gff_file);
	exon_tree = exon_tree_new (exon_stmt, overlapping_stmt, cs);
	exon_tree_index_dump (exon_tree, gff_file);

	log_info ("Create thread pool");
	thpool = thpool_init (threads);

	AbnormalArg *ab_args = xcalloc (num_files, sizeof (AbnormalArg));

	for (i = 0; i < array_len (sam_files); i++)
		{
			sam_file = array_get (sam_files, i);

			// Init abnormal_filter arg
			ab_args[i] = (AbnormalArg)
			{
				.tid              = i + 1,
				.num_threads      = num_files,
				.sam_file         = sam_file,
				.either           = either,
				.exon_frac        = exon_frac,
				.alignment_frac   = alignment_frac,
				.exon_tree        = exon_tree,
				.cs               = cs,
				.alignment_stmt   = alignment_stmt,
				.queryname_sorted = sorted,
				.max_distance     = max_distance
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
		"       %*c                [-p STR] [-t INT] [-c INT] [-m INT]\n"
		"       %*c                [-f FLOAT] [-F FLOAT | -r] [-e]\n"
		"       %*c                [-i FILE] -a FILE <FILE> ...\n"
		"\n"
		"Arguments:\n"
		"   One or more alignment files SAM/BAM\n"
		"\n"
		"Mandatory Options:\n"
		"   -a, --annotation-file   GTF/GFF3 annotation file\n"
		"   -i, --input-file        FILE containing a list of aligned files SAM/BAM\n"
		"                           to be processed. They must be separated by newline.\n"
		"                           This option is not manditory if one or more SAM/BAM\n"
		"                           files are passed as argument. If 'input-file' and\n"
		"                           arguments are set concomitantly, then the union of\n"
		"                           all alignment files is used\n"
		"\n"
		"Options:\n"
		"   -h, --help              Show help options\n"
		"   -q, --quiet             Decrease verbosity to error messages only\n"
		"                           or supress terminal outputs at all if\n"
		"                           'log-file' is passed\n"
		"       --silent            Same as '--quiet'\n"
		"   -d, --debug             Increase verbosity to debug level\n"
		"   -l, --log-file          Print log messages to a FILE\n"
		"   -o, --output-dir        Output directory. Create the directory if it does\n"
		"                           not exist [default:\"%s\"]\n"
		"   -p, --prefix            Prefix output files [default:\"%s\"]\n"
		"   -t, --threads           Number of threads [default:\"%d\"]\n"
		"   -c, --cache-size        Set SQLite3 cache size in KiB [default: \"%d\"]\n"
		"   -s, --sorted            Assume all reads are grouped by queryname, even if\n"
		"                           there is no SAM/BAM header tag 'SO:queryname'\n"
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
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_THREADS, DEFAULT_CACHE_SIZE,
		DEFAULT_MAX_DISTANCE, DEFAULT_EXON_FRAC, DEFAULT_ALIGNMENT_FRAC);
}

static void
print_try_help (FILE *fp)
{
	fprintf (fp, "Try '%s process-sample --help' for more information\n",
			PACKAGE);
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
		{"max-distance",    required_argument, 0, 'm'},
		{"cache-size",      required_argument, 0, 'c'},
		{"sorted",          no_argument,       0, 's'},
		{"exon-frac",       required_argument, 0, 'f'},
		{"alignment-frac",  required_argument, 0, 'F'},
		{"either",          no_argument,       0, 'e'},
		{"reciprocal",      no_argument,       0, 'r'},
		{"input-file",      required_argument, 0, 'i'},
		{0,                 0,                 0,  0 }
	};

	// Init variables to default values
	int         silent         = DEFAULT_LOG_SILENT;
	int         log_level      = DEFAULT_LOG_LEVEL;
	int         threads        = DEFAULT_THREADS;
	int         sorted         = DEFAULT_SORTED;
	int         max_distance   = DEFAULT_MAX_DISTANCE;
	int         cache_size     = DEFAULT_CACHE_SIZE;
	int         either         = DEFAULT_EITHER;
	int         reciprocal     = DEFAULT_RECIPROCAL;
	float       exon_frac      = DEFAULT_EXON_FRAC;
	float       alignment_frac = DEFAULT_ALIGNMENT_FRAC;
	const char *output_dir     = DEFAULT_OUTPUT_DIR;
	const char *prefix         = DEFAULT_PREFIX;
	const char *log_file       = DEFAULT_LOG_FILE;
	const char *gff_file       = DEFAULT_GFF_FILE;
	const char *input_file     = DEFAULT_INPUT_FILE;

	Array *sam_files = array_new (xfree);
	Logger *logger = NULL;

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int alignment_frac_set = 0;
	int c, i;

	while ((c = getopt_long (argc, argv, "hqdsl:a:o:p:t:m:c:f:F:eri:", opt, &option_index)) >= 0)
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
						silent = 1;
						break;
					}
				case 'd':
					{
						log_level = LOG_DEBUG;
						break;
					}
				case 'l':
					{
						log_file = optarg;
						break;
					}
				case 'a':
					{
						gff_file = optarg;
						break;
					}
				case 'o':
					{
						output_dir = optarg;
						break;
					}
				case 'p':
					{
						prefix = optarg;
						break;
					}
				case 't':
					{
						threads = atoi (optarg);
						break;
					}
				case 'c':
					{
						cache_size = atoi (optarg);
						break;
					}
				case 'm':
					{
						max_distance = atoi (optarg);
						break;
					}
				case 's':
					{
						sorted = 1;
						break;
					}
				case 'f':
					{
						exon_frac = atof (optarg);
						break;
					}
				case 'F':
					{
						alignment_frac = atof (optarg);
						alignment_frac_set = 1;
						break;
					}
				case 'e':
					{
						either = 1;
						break;
					}
				case 'r':
					{
						reciprocal = 1;
						break;
					}
				case 'i':
					{
						input_file = optarg;
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
		array_add (sam_files, xstrdup (argv[i]));

	// Catch all alignment files passed into
	// --input-file
	if (input_file != NULL)
		read_file_lines (sam_files, input_file);

	/*
	* Validate arguments and mandatory options
	*/

	// If no one file was passed, throw an error
	if (array_len (sam_files) == 0)
		{
			fprintf (stderr, "%s: Missing alignment files (SAM/BAM)\n", PACKAGE);
			print_try_help (stderr);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Test if all alignment files axist
	for (i = 0; i < array_len (sam_files); i++)
		{
			const char *sam_file = array_get (sam_files, i);
			if (!exists (sam_file))
				{
					fprintf (stderr, "%s: alignment file '%s': No such file\n", PACKAGE, sam_file);
					rc = EXIT_FAILURE; goto Exit;
				}
		}

	// If no annotation-file was passed, throw an error
	if (gff_file == NULL)
		{
			fprintf (stderr, "%s: Missing annotation file (GTF/GFF3)\n", PACKAGE);
			print_try_help (stderr);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Test if annotation file exists
	if (!exists (gff_file))
		{
			fprintf (stderr, "%s: annotation file '%s': No such file\n", PACKAGE, gff_file);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*
	* Validate options
	*/

	// If --reciprocal, set alignment_frac to exon_frac
	if (reciprocal)
		{
			if (alignment_frac_set)
				{
					fprintf (stderr, "%s: --reciprocal must be set solely with '-f'\n", PACKAGE);
					rc = EXIT_FAILURE; goto Exit;
				}
			alignment_frac = exon_frac;
		}

	// Validate exon_frac: ]0, 1]
	if (exon_frac > 1 || exon_frac <= 0)
		{
			fprintf (stderr, "%s: --exon-frac must be in the interval of ]0, 1]\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate alignment_frac ]0, 1]
	if (alignment_frac > 1 || alignment_frac <= 0)
		{
			fprintf (stderr, "%s: --alignment-frac must be in the interval of ]0, 1]\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate threads > 0
	if (threads < 1)
		{
			fprintf (stderr, "%s: --threads must be greater or equal to 1\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate cache_size >= DEFAULT_CACHE_SIZE
	if (cache_size < DEFAULT_CACHE_SIZE)
		{
			fprintf (stderr, "%s: --cache-size must be greater or equal to %uKiB\n",
					PACKAGE, DEFAULT_CACHE_SIZE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate max_distance >= 0
	if (max_distance < 0)
		{
			fprintf (stderr, "%s: --max-distance must be a positive value\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*
	* Final settings
	*/

	// Avoid to include repetitive files
	array_uniq (sam_files, cmpstringp);

	// If it's silent and no log file
	// was passsed, then set log_level
	// to LOG_ERROR - At least print
	// errors
	if (silent && log_file == NULL)
		{
			silent = 0;
			log_level = LOG_ERROR;
		}

	logger = logger_new (log_file, log_level, silent, 1);

	/*
	* RUN FOOLS
	*/
	process_sample (output_dir, prefix, sam_files, gff_file, threads,
			cache_size, sorted, max_distance, exon_frac, alignment_frac,
			either);

Exit:
	logger_free (logger);
	array_free (sam_files, 1);
	return rc;
}
