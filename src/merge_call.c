#include "config.h"

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
#include "set.h"
#include "chr.h"
#include "array.h"
#include "blacklist.h"
#include "cluster.h"
#include "db_merge.h"
#include "merge_call.h"

#define DEFAULT_CACHE_SIZE            DB_DEFAULT_CACHE_SIZE
#define DEFAULT_PREFIX                "out"
#define DEFAULT_OUTPUT_DIR            "."
#define DEFAULT_LOG_SILENT            0
#define DEFAULT_LOG_LEVEL             LOG_INFO
#define DEFAULT_IN_PLACE              0
#define DEFAULT_EPS                   300
#define DEFAULT_MIN_PTS               10
#define DEFAULT_LOG_FILE              NULL
#define DEFAULT_INPUT_FILE            NULL
#define DEFAULT_BLACKLIST_CHR         "chrM"
#define DEFAULT_DISTANCE              10000
#define DEFAULT_SUPPORT               1
#define DEFAULT_BLACKLIST_REGION      NULL
#define DEFAULT_GFF_FEATURE           "gene"
#define DEFAULT_GFF_ATTRIBUTE         "gene_type"
#define DEFAULT_GFF_ATTRIBUTE_VALUE   "processed_pseudogene"

static void
merge_call (const char *output_dir, const char *prefix, Array *db_files,
		const char *output_file, int cache_size, int epsilon, int min_pts,
		Set *blacklist_chr, int distance, int support, const char *blacklist_region,
		const char *gff_feature, const char *gff_attribute,
		const Set *gff_attribute_values)
{
	log_trace ("Inside %s", __func__);

	sqlite3 *db = NULL;

	sqlite3_stmt *cluster_stmt = NULL;
	sqlite3_stmt *clustering_stmt = NULL;
	sqlite3_stmt *blacklist_stmt = NULL;
	sqlite3_stmt *overlapping_blacklist_stmt = NULL;

	Blacklist *blacklist = NULL;
	ChrStd *cs = NULL;

	char *db_file = NULL;

	log_info (">>> Merge Call step <<<");

	log_info ("Create output dir '%s'", output_dir);
	mkdir_p (output_dir);

	// Merge all databases
	if (output_file != NULL)
		{
			// Use first database
			db_file = xstrdup (output_file);

			log_info ("Connect to database '%s'", db_file);

			// Connect to database
			db = db_connect (db_file);
		}
	else
		{
			// Assemble database output filename
			xasprintf_concat (&db_file, "%s/%s.db",
					output_dir, prefix);

			log_info ("Create and connect to database '%s'", db_file);

			// Create a new database
			db = db_create (db_file);
		}

	// Increase the cache size
	db_cache_size (db, cache_size);

	// Create cluster statement
	cluster_stmt = db_prepare_cluster_stmt (db);

	// Create clustering statement
	clustering_stmt = db_prepare_clustering_stmt (db);

	// Create blacklist statement
	blacklist_stmt = db_prepare_blacklist_stmt (db);

	// Create overlapping blacklist statement
	overlapping_blacklist_stmt =
		db_prepare_overlapping_blacklist_stmt (db);

	// Get chromosome standardization
	cs = chr_std_new ();

	// Index blacklisted regions
	blacklist = blacklist_new (blacklist_stmt,
			overlapping_blacklist_stmt, cs);

	// If the user passed blacklisted regions
	if (blacklist_region != NULL)
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			log_info ("Index blacklisted regions from file '%s'",
					blacklist_region);

			blacklist_index_dump (blacklist, blacklist_region, gff_feature,
					gff_attribute, gff_attribute_values);

			// Commit
			db_end_transaction (db);
		}

	// If there are files to merge with ...
	if (array_len (db_files))
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			// Time to merge them all!
			log_info ("Merge all files with '%s'", db_file);
			db_merge (db, array_len (db_files),
					(char **) array_data (db_files));

			// Commit
			db_end_transaction (db);
		}

	// Begin transaction to speed up
	db_begin_transaction (db);

	// RUN
	log_info ("Run clustering step for '%s'", db_file);
	cluster (cluster_stmt, clustering_stmt, epsilon, min_pts,
			distance, support, blacklist_chr, blacklist);

	// Commit
	db_end_transaction (db);

	// Cleanup
	xfree (db_file);
	db_finalize (cluster_stmt);
	db_finalize (clustering_stmt);
	db_finalize (blacklist_stmt);
	db_finalize (overlapping_blacklist_stmt);
	db_close (db);

	chr_std_free (cs);
	blacklist_free (blacklist);
}

static void
print_usage (FILE *fp)
{
	int pkg_len = strlen (PACKAGE);
	fprintf (fp,
		"%s\n"
		"\n"
		"Usage: %s merge-call [-h] [-q] [-d] [-l FILE] [-o DIR] [-p STR]\n"
		"       %*c            [-c INT] [-I] [-e INT] [-m INT] [-b STR]\n"
		"       %*c            [-a FILE] [[-F STR] [-A KEY=VALUE]]\n"
		"       %*c            [-x INT] [-g INT] [-i FILE]\n"
		"       %*c            <FILE> ...\n"
		"\n"
		"Arguments:\n"
		"   One or more SQLite3 databases generated in the 'process-sample' step\n"
		"\n"
		"Mandatory Options:\n"
		"   -i, --input-file           FILE containing a list of SQLite3 databases\n"
		"                              to be processed. They must be separated by newline.\n"
		"                              This option is not manditory if one or more SQLite3\n"
		"                              databases are passed as argument. If 'input-file' and\n"
		"                              arguments are set concomitantly, then the union of\n"
		"                              all files is used\n"
		"\n"
		"Options:\n"
		"   -h, --help                 Show help options\n"
		"   -q, --quiet                Decrease verbosity to error messages only\n"
		"                              or supress terminal outputs at all if\n"
		"                              'log-file' is passed\n"
		"       --silent               Same as '--quiet'\n"
		"   -d, --debug                Increase verbosity to debug level\n"
		"   -l, --log-file             Print log messages to a FILE\n"
		"   -o, --output-dir           Output directory. Create the directory if it does\n"
		"                              not exist [default:\"%s\"]\n"
		"   -p, --prefix               Prefix output files [default:\"%s\"]\n"
		"   -I, --in-place             Merge all databases with the first one of the list,\n"
		"                              instead of creating a new file\n"
		"   -c, --cache-size           Set SQLite3 cache size in KiB [default:\"%d\"]\n"
		"   -e, --epsilon              DBSCAN: Maximum distance between two alignments\n"
		"                              inside a cluster [default:\"%d\"]\n"
		"   -m, --min-pts              DBSCAN: Minimum number of points required to form a\n"
		"                              dense region [default:\"%d\"]\n"
		"   -b, --blacklist-chr        Avoid clustering from and to this chromosome. This\n"
		"                              option may be passed multiple times [default:\"%s\"]\n"
		"   -a, --blacklist-region     GTF/GFF3/BED blacklisted regions. If the file is in\n"
		"                              GTF/GFF3 format, the user may indicate the 'feature'\n"
		"                              (third column), the 'attribute' (ninth column) and\n"
		"                              its values\n"
		"   -F, --gff-feature          The value of 'feature' (third column) for GTF/GFF3\n"
		"                              '--blacklist-region' file [default:\"%s\"]\n"
		"   -A, --gff-attribute        The 'attribute' (ninth column) for GTF/GFF3\n"
		"                              '--blacklist-region' file. It may be passed in the\n"
		"                              format key=value, where 'value' is comma separated\n"
		"                              values (e.g. gene_type=lincRNA,pseudogene,sRNA).\n"
		"                              Each value will match as regex, so 'pseudogene'\n"
		"                              can capture IG_C_pseudogene, IG_V_pseudogene etc\n"
		"                              [default:\"%s=%s\"]\n"
		"   -x, --parental-distance    Minimum distance allowed between a cluster and\n"
		"                              its putative parental gene [default:\"%d\"]\n"
		"   -g, --genotype-support     Minimum number of reads comming from a given source\n"
		"                              (BAM) within a cluster [default:\"%d\"]\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ', pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_CACHE_SIZE, DEFAULT_EPS,
		DEFAULT_MIN_PTS, DEFAULT_BLACKLIST_CHR, DEFAULT_GFF_FEATURE, DEFAULT_GFF_ATTRIBUTE,
		DEFAULT_GFF_ATTRIBUTE_VALUE, DEFAULT_DISTANCE, DEFAULT_SUPPORT);
}

static void
print_try_help (FILE *fp)
{
	fprintf (fp, "Try '%s merge-call --help' for more information\n",
			PACKAGE);
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
		{"help",              no_argument,       0, 'h'},
		{"quiet",             no_argument,       0, 'q'},
		{"silent",            no_argument,       0, 'q'},
		{"debug",             no_argument,       0, 'd'},
		{"in-place",          no_argument,       0, 'I'},
		{"log-file",          required_argument, 0, 'l'},
		{"output-dir",        required_argument, 0, 'o'},
		{"prefix",            required_argument, 0, 'p'},
		{"cache-size",        required_argument, 0, 'c'},
		{"input-file",        required_argument, 0, 'i'},
		{"epsilon",           required_argument, 0, 'e'},
		{"min-pts",           required_argument, 0, 'm'},
		{"parental-distance", required_argument, 0, 'x'},
		{"genotype-support",  required_argument, 0, 'g'},
		{"blacklist-chr",     required_argument, 0, 'b'},
		{"blacklist-region",  required_argument, 0, 'a'},
		{"gff-feature",       required_argument, 0, 'F'},
		{"gff-attribute",     required_argument, 0, 'A'},
		{0,                   0,                 0,  0 }
	};

	// Init variables to default values
	int         silent           = DEFAULT_LOG_SILENT;
	int         log_level        = DEFAULT_LOG_LEVEL;
	int         cache_size       = DEFAULT_CACHE_SIZE;
	int         in_place         = DEFAULT_IN_PLACE;
	int         epsilon          = DEFAULT_EPS;
	int         min_pts          = DEFAULT_MIN_PTS;
	int         distance         = DEFAULT_DISTANCE;
	int         support          = DEFAULT_SUPPORT;
	const char *blacklist_region = DEFAULT_BLACKLIST_REGION;
	const char *gff_feature      = DEFAULT_GFF_FEATURE;
	const char *gff_attribute    = DEFAULT_GFF_ATTRIBUTE;
	const char *output_dir       = DEFAULT_OUTPUT_DIR;
	const char *prefix           = DEFAULT_PREFIX;
	const char *log_file         = DEFAULT_LOG_FILE;
	const char *input_file       = DEFAULT_INPUT_FILE;


	char *output_file = NULL;
	int index_ = 0;

	ChrStd *cs = chr_std_new ();
	Set *blacklist_chr = set_new_full (str_hash, str_equal, NULL);
	Array *db_files = array_new (xfree);
	Set *gff_attribute_values = NULL;
	Logger *logger = NULL;

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c, i;

	while ((c = getopt_long (argc, argv, "hqdIl:o:p:c:e:m:b:a:F:A:x:g:i:", opt, &option_index)) >= 0)
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
				case 'I':
					{
						in_place = 1;
						break;
					}
				case 'l':
					{
						log_file = optarg;
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
				case 'c':
					{
						cache_size = atoi (optarg);
						break;
					}
				case 'e':
					{
						epsilon = atoi (optarg);
						break;
					}
				case 'm':
					{
						min_pts = atoi (optarg);
						break;
					}
				case 'b':
					{
						set_insert (blacklist_chr,
								chr_std_lookup (cs, optarg));
						break;
					}
				case 'x':
					{
						distance = atoi (optarg);
						break;
					}
				case 'g':
					{
						support = atoi (optarg);
						break;
					}
				case 'a':
					{
						blacklist_region = optarg;
						break;
					}
				case 'F':
					{
						gff_feature = optarg;
						break;
					}
				case 'A':
					{
						char *scratch1 = NULL;
						char *scratch2 = NULL;
						char *token = NULL;
						char *subtoken = NULL;

						token = strtok_r (optarg, "=", &scratch1);
						gff_attribute = token == NULL
							? DEFAULT_GFF_ATTRIBUTE
							: token;

						token = strtok_r (NULL, "=", &scratch1);
						if (token != NULL)
							{
								set_free (gff_attribute_values);
								gff_attribute_values = set_new_full (str_hash,
										str_equal, NULL);

								subtoken = strtok_r (token, ",", &scratch2);

								while (subtoken != NULL)
									{
										set_insert (gff_attribute_values, subtoken);
										subtoken = strtok_r (NULL, ",", &scratch2);
									}
							}

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
		array_add (db_files, xstrdup (argv[i]));

	// Catch all alignment files passed into
	// --input-file
	if (input_file != NULL)
		read_file_lines (db_files, input_file);

	/*Validate arguments and mandatory options*/

	// If no one file was passed, throw an error
	if (array_len (db_files) == 0)
		{
			fprintf (stderr, "%s: Missing SQLite3 databases\n", PACKAGE);
			print_try_help (stderr);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Test if all alignment files axist
	for (i = 0; i < array_len (db_files); i++)
		{
			const char *db_file = array_get (db_files, i);
			if (!exists (db_file))
				{
					fprintf (stderr, "%s: SQLite3 database '%s': No such file\n", PACKAGE, db_file);
					rc = EXIT_FAILURE; goto Exit;
				}
		}

	/*Validate options*/

	// Validate blacklist_region file
	if (blacklist_region != NULL && !exists (blacklist_region))
		fprintf (stderr, "%s: --blacklist-region '%s': No such file\n", PACKAGE,
				blacklist_region);

	// Validate cache_size >= DEFAULT_CACHE_SIZE
	if (cache_size < DEFAULT_CACHE_SIZE)
		{
			fprintf (stderr, "%s: --cache-size must be greater or equal to %uKiB\n",
					PACKAGE, DEFAULT_CACHE_SIZE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (epsilon < 0)
		{
			fprintf (stderr, "%s: --epsilon must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (min_pts < 3)
		{
			fprintf (stderr, "%s: --min_pts must be greater than 2\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (distance < 0)
		{
			fprintf (stderr, "%s: --parental-distance must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (support < 0)
		{
			fprintf (stderr, "%s: --genotype-support must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*Final settings*/

	// Add default blacklisted chr if none
	// has been passed
	if (set_size (blacklist_chr) == 0)
		set_insert (blacklist_chr,
				chr_std_lookup (cs, DEFAULT_BLACKLIST_CHR));

	// If gff_attribute_values was not set
	if (gff_attribute_values == NULL)
		gff_attribute_values = set_new_full (str_hash,
				str_equal, NULL);

	// Or it exists and has no values
	if (set_size (gff_attribute_values) == 0)
		set_insert (gff_attribute_values,
				DEFAULT_GFF_ATTRIBUTE_VALUE);

	// Copy the name of possible output file, if the
	// user chose --in-place
	if (in_place)
		{
			// Copy first file to use it as the final
			// database
			output_file = xstrdup (array_get (db_files, 0));

			// Remove all repetitive output_file from list
			while (array_find_with_equal_fun (db_files,
						output_file, equalstring, &index_))
				array_remove_index (db_files, index_);
		}

	// Avoid to include repetitive files
	array_uniq (db_files, cmpstringp);

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

	// RUN FOOLS
	merge_call (output_dir, prefix, db_files, output_file,
			cache_size, epsilon, min_pts, blacklist_chr,
			distance, support, blacklist_region, gff_feature,
			gff_attribute, gff_attribute_values);

Exit:
	logger_free (logger);
	chr_std_free (cs);
	set_free (blacklist_chr);
	set_free (gff_attribute_values);
	array_free (db_files, 1);
	xfree (output_file);
	return rc;
}
