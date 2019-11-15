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
#include "gff.h"
#include "chr.h"
#include "array.h"
#include "blacklist.h"
#include "cluster.h"
#include "db_merge.h"
#include "retrocopy.h"
#include "genotype.h"
#include "merge_call.h"

#define DEFAULT_CACHE_SIZE            DB_DEFAULT_CACHE_SIZE
#define DEFAULT_PREFIX                "out"
#define DEFAULT_OUTPUT_DIR            "."
#define DEFAULT_LOG_SILENT            0
#define DEFAULT_LOG_LEVEL             LOG_INFO
#define DEFAULT_IN_PLACE              0
#define DEFAULT_EPS                   300
#define DEFAULT_MIN_PTS               10
#define DEFAULT_BLACKLIST_CHR         "chrM"
#define DEFAULT_PARENTAL_DISTANCE     10000
#define DEFAULT_SUPPORT               1
#define DEFAULT_BLACKLIST_PADDING     0
#define DEFAULT_GFF_FEATURE           "gene"
#define DEFAULT_GFF_ATTRIBUTE1        "gene_type"
#define DEFAULT_GFF_ATTRIBUTE_VALUE1  "processed_pseudogene"
#define DEFAULT_GFF_ATTRIBUTE2        "tag"
#define DEFAULT_GFF_ATTRIBUTE_VALUE2  "retrogene"
#define DEFAULT_NEAR_GENE_DISTANCE    3
#define DEFAULT_THREADS               1
#define DEFAULT_CROSSING_READS        10

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
	ChrStd      *cs;
	int          padding;
	int          parental_dist;
	int          support;
	int          near_gene_dist;

	// Genotyping
	int          crossing_reads;
	int          threads;
};

typedef struct _MergeCall MergeCall;

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
	sqlite3_stmt *cluster_merging_stmt = NULL;
	sqlite3_stmt *genotype_stmt = NULL;

	Blacklist *blacklist = NULL;

	int num_clusters = 0;
	char *db_file = NULL;

	log_info (">>> Merge Call step <<<");

	log_info ("Create output dir '%s'", mc->output_dir);
	mkdir_p (mc->output_dir);

	// Merge all databases
	if (mc->output_file != NULL)
		{
			// Use first database
			db_file = xstrdup (mc->output_file);

			log_info ("Connect to database '%s'", db_file);

			// Connect to database
			db = db_connect (db_file);
		}
	else
		{
			// Assemble database output filename
			xasprintf_concat (&db_file, "%s/%s.db",
					mc->output_dir, mc->prefix);

			log_info ("Create and connect to database '%s'", db_file);

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
					log_info ("Index blacklist entries from GTF/GFF3 file '%s'",
							mc->blacklist_region);
					blacklist_index_dump_from_gff (blacklist,
							mc->blacklist_region, mc->filter);
				}
			else
				{
					log_info ("Index blacklist entries from BED file '%s'",
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
			log_info ("Merge all files with '%s'", db_file);
			db_merge (db, array_len (mc->db_files),
					(char **) array_data (mc->db_files));

			// Commit
			db_end_transaction (db);
		}

	// Begin transaction to speed up
	db_begin_transaction (db);

	// Clustering
	log_info ("Run clustering step for '%s'", db_file);
	num_clusters = cluster (cluster_stmt, clustering_stmt,
			mc->epsilon, mc->min_pts, mc->parental_dist, mc->support,
			mc->blacklist_chr, blacklist, mc->padding);

	// Commit
	db_end_transaction (db);

	if (!num_clusters)
		{
			log_warn ("No cluster has been found!");
		}
	else
		{
			// Begin transaction to speed up
			db_begin_transaction (db);

			// Filtering and Annotation
			log_info ("Run retrocopy annotation step for '%s'", db_file);
			retrocopy (retrocopy_stmt, cluster_merging_stmt, mc->near_gene_dist);

			// Commit
			db_end_transaction (db);

			// Begin transaction to speed up
			db_begin_transaction (db);

			// Genotyping
			log_info ("Run genotype annotation step for '%s'", db_file);
			genotype (genotype_stmt, mc->threads, mc->crossing_reads);

			// Commit
			db_end_transaction (db);
		}

	// Cleanup
	xfree (db_file);
	db_finalize (cluster_stmt);
	db_finalize (clustering_stmt);
	db_finalize (blacklist_stmt);
	db_finalize (overlapping_blacklist_stmt);
	db_finalize (retrocopy_stmt);
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
		"Usage: %s merge-call [-h] [-q] [-d] [-l FILE] [-o DIR] [-p STR]\n"
		"       %*c            [-c INT] [-I] [-e INT] [-m INT] [-b STR]\n"
		"       %*c            [-B FILE] [[-T STR] [[-H|S] KEY=VALUE]]\n"
		"       %*c            [-P INT] [-x INT] [-g INT] [-n INT]\n"
		"       %*c            [-t INT] [-C INT] [-i FILE] <FILE> ...\n"
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
		"Input/Output Options:\n"
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
		"\n"
		"SQLite3 Options:\n"
		"   -c, --cache-size           Set SQLite3 cache size in KiB [default:\"%d\"]\n"
		"\n"
		"Clustering Options:\n"
		"   -e, --epsilon              DBSCAN: Maximum distance between two alignments\n"
		"                              inside a cluster [default:\"%d\"]\n"
		"   -m, --min-pts              DBSCAN: Minimum number of points required to form a\n"
		"                              dense region [default:\"%d\"]\n"
		"\n"
		"Filter & Annotation Options:\n"
		"   -b, --blacklist-chr        Avoid clustering from and to this chromosome. This\n"
		"                              option can be passed multiple times [default:\"%s\"]\n"
		"   -B, --blacklist-region     GTF/GFF3/BED blacklisted regions. If the file is in\n"
		"                              GTF/GFF3 format, the user may indicate the 'feature'\n"
		"                              (third column), the 'attribute' (ninth column) and\n"
		"                              its values\n"
		"   -P, --blacklist-padding    Increase the blacklisted regions ranges (left and right)\n"
		"                              by N bases [default:\"%d\"]\n"
		"   -T, --gff-feature          The value of 'feature' (third column) for GTF/GFF3\n"
		"                              file [default:\"%s\"]\n"
		"   -H, --gff-hard-attribute   The 'attribute' (ninth column) for GTF/GFF3\n"
		"                              file. It may be passed in the format key=value\n"
		"                              (e.g. gene_type=pseudogene). Each value will match\n"
		"                              as regex, so 'pseudogene' can capture IG_C_pseudogene,\n"
		"                              IG_V_pseudogene etc. This option can be passed multiple\n"
		"                              times and must be true in all of them\n"
		"   -S, --gff-soft-attribute   Works as 'gff-hard-attribute'. The difference is\n"
		"                              if this option is passed multiple times, it needs\n"
		"                              to be true only once\n"
		"                              [default:\"%s=%s %s=%s\"]\n"
		"   -x, --parental-distance    Minimum distance allowed between a cluster and\n"
		"                              its putative parental gene [default:\"%d\"]\n"
		"   -g, --genotype-support     Minimum number of reads comming from a given source\n"
		"                              (BAM) within a cluster [default:\"%d\"]\n"
		"   -n, --near-gene-distance   Minimum ranked distance between genes in order to\n"
		"                              consider them close [default:\"%d\"]\n"
		"\n"
		"Genotyping Options:\n"
		"   -t, --threads              Number of threads [default:\"%d\"]\n"
		"   -C, --crossing-reads       Minimum number of reads crossing the insertion point\n"
		"                              in order to consider evidence of heterozygosis\n"
		"                              [default:\"%d\"]\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ', pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_CACHE_SIZE, DEFAULT_EPS, DEFAULT_MIN_PTS,
		DEFAULT_BLACKLIST_CHR, DEFAULT_BLACKLIST_PADDING, DEFAULT_GFF_FEATURE, DEFAULT_GFF_ATTRIBUTE1,
		DEFAULT_GFF_ATTRIBUTE_VALUE1, DEFAULT_GFF_ATTRIBUTE2, DEFAULT_GFF_ATTRIBUTE_VALUE2,
		DEFAULT_PARENTAL_DISTANCE, DEFAULT_SUPPORT, DEFAULT_NEAR_GENE_DISTANCE,
		DEFAULT_THREADS, DEFAULT_CROSSING_READS);
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
		.cs               = chr_std_new (),
		.padding          = DEFAULT_BLACKLIST_PADDING,
		.parental_dist    = DEFAULT_PARENTAL_DISTANCE,
		.support          = DEFAULT_SUPPORT,
		.near_gene_dist   = DEFAULT_NEAR_GENE_DISTANCE,
		.crossing_reads   = DEFAULT_CROSSING_READS,
		.threads          = DEFAULT_THREADS
	};
}

static void
merge_call_destroy (MergeCall *mc)
{
	if (mc == NULL)
		return;

	array_free (mc->db_files, 1);
	set_free (mc->blacklist_chr);
	gff_filter_free (mc->filter);
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

	// Test if all alignment files axist
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
		fprintf (stderr, "%s: --blacklist-region '%s': No such file\n", PACKAGE,
				mc->blacklist_region);

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
			fprintf (stderr, "%s: --min_pts must be greater than 2\n", PACKAGE);
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

	if (mc->near_gene_dist < 1)
		{
			fprintf (stderr, "%s: --near-gene-distance must be greater or equal to 1\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (mc->crossing_reads < 1)
		{
			fprintf (stderr, "%s: --crossing-reads must be greater or equal to 1\n", PACKAGE);
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

	// If gff_attribute_values was not set
	if (gff_filter_hard_attribute_size (mc->filter) == 0
			&& gff_filter_soft_attribute_size (mc->filter) == 0)
		{
			gff_filter_insert_soft_attribute (mc->filter, DEFAULT_GFF_ATTRIBUTE1,
					DEFAULT_GFF_ATTRIBUTE_VALUE1);
			gff_filter_insert_soft_attribute (mc->filter, DEFAULT_GFF_ATTRIBUTE2,
					DEFAULT_GFF_ATTRIBUTE_VALUE2);
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
	// was passsed, then set log_level
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
		{"gff-hard-attribute", required_argument, 0, 'H'},
		{"gff-soft-attribute", required_argument, 0, 'S'},
		{"near-gene-distance", required_argument, 0, 'n'},
		{"crossing-reads",     required_argument, 0, 'C'},
		{"threads",            required_argument, 0, 't'},
		{0,                    0,                 0,  0 }
	};

	// Init variables to default values
	MergeCall mc = {};
	merge_call_init (&mc);

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c, i;

	while ((c = getopt_long (argc, argv, "hqdIl:o:p:c:e:m:b:B:P:T:H:S:x:g:n:C:t:i:", opt, &option_index)) >= 0)
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
						mc.near_gene_dist = atoi (optarg);
						break;
					}
				case 'C':
					{
						mc.crossing_reads = atoi (optarg);
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
				case 'H':
					{
						char *scratch = NULL;

						const char *key = strtok_r (optarg,
								"=", &scratch);
						const char *value = strtok_r (NULL,
								"=", &scratch);

						if (key == NULL || value == NULL)
							{
								fprintf (stderr,
										"%s: --gff-hard-attribute KEY=VALUE: "
										"Missing complete pair\n", PACKAGE);
								print_try_help (stderr);
								rc = EXIT_FAILURE; goto Exit;
							}

						gff_filter_insert_hard_attribute (mc.filter,
								key, value);
						break;
					}
				case 'S':
					{
						char *scratch = NULL;

						const char *key = strtok_r (optarg,
								"=", &scratch);
						const char *value = strtok_r (NULL,
								"=", &scratch);

						if (key == NULL || value == NULL)
							{
								fprintf (stderr,
										"%s: --gff-soft-attribute KEY=VALUE: "
										"Missing complete pair\n", PACKAGE);
								print_try_help (stderr);
								rc = EXIT_FAILURE; goto Exit;
							}

						gff_filter_insert_soft_attribute (mc.filter,
								key, value);
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
			run (&mc);
		}

Exit:
	merge_call_destroy (&mc);
	return rc;
}
