#include "config.h"

#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "str.h"
#include "db.h"
#include "log.h"
#include "logger.h"
#include "vcf.h"
#include "make_vcf.h"

#define DEFAULT_PREFIX                "out"
#define DEFAULT_OUTPUT_DIR            "."
#define DEFAULT_LOG_SILENT            0
#define DEFAULT_LOG_LEVEL             LOG_INFO
#define DEFAULT_NEAR_GENE_DIST        10000
#define DEFAULT_ORIENTATION_ERROR     0.05

struct _MakeVcf
{
	// I/O
	const char  *db_file;
	const char  *prefix;
	const char  *output_dir;

	// Log
	Logger      *logger;
	const char  *log_file;
	int          log_level;
	int          silent;

	// Filter
	long         near_gene_dist;
	float        orientation_error;

	// Annotation
	const char  *fasta_file;
};

typedef struct _MakeVcf MakeVcf;

static void
run (MakeVcf *v)
{
	log_trace ("Inside %s", __func__);

	sqlite3 *db = NULL;
	char *vcf_file = NULL;

	// Assemble VCF file
	xasprintf_concat (&vcf_file, "%s/%s.vcf",
			v->output_dir, v->prefix);

	log_info ("Create output dir '%s'", v->output_dir);
	mkdir_p (v->output_dir);

	log_info ("Connect to database %s", v->db_file);
	db = db_connect (v->db_file);

	// Fill options
	VCFOption opt = {
		.near_gene_dist    = v->near_gene_dist,
		.orientation_error = v->orientation_error,
		.fasta_file        = v->fasta_file
	};

	log_info ("Generate VCF file %s", vcf_file);
	vcf (db, vcf_file, &opt);

	log_info ("VCF file '%s' is ready!", vcf_file);
	xfree (vcf_file);
	db_close (db);
}

static void
print_usage (FILE *fp)
{
	int pkg_len = strlen (PACKAGE);
	fprintf (fp,
		"%s\n"
		"\n"
		"Usage: %s make-vcf [-h] [-q] [-d] [-s] [-l FILE] [-o DIR]\n"
		"       %*c          [-p STR] [-n INT] [-e FLOAT]\n"
		"       %*c          [-r FILE] <FILE>\n"
		"\n"
		"Arguments:\n"
		"   One SQLite3 database generated in the 'process-sample'\n"
		"   and 'merge-call' steps\n"
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
		"\n"
		"Filter & Annotation Options:\n"
		"   -n, --near-gene-dist       Minimum distance between genes in order to\n"
		"                              consider them close [default:\"%d\"]\n"
		"   -e, --orientation-error    Maximum error allowed for orientation rho\n"
		"                              [default:\"%.2f\"]\n"
		"   -r, --reference-file       FASTA FILE for the reference genome\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_NEAR_GENE_DIST,
		DEFAULT_ORIENTATION_ERROR);
}

static void
print_try_help (FILE *fp)
{
	fprintf (fp, "Try '%s make-vcf --help' for more information\n",
			PACKAGE);
}

static void
make_vcf_init (MakeVcf *vcf)
{
	*vcf = (MakeVcf) {
		.db_file            = NULL,
		.output_dir         = DEFAULT_OUTPUT_DIR,
		.prefix             = DEFAULT_PREFIX,
		.logger             = NULL,
		.log_file           = NULL,
		.log_level          = DEFAULT_LOG_LEVEL,
		.silent             = DEFAULT_LOG_SILENT,
		.near_gene_dist     = DEFAULT_NEAR_GENE_DIST,
		.orientation_error  = DEFAULT_ORIENTATION_ERROR,
		.fasta_file         = NULL
	};
}

static void
make_vcf_destroy (MakeVcf *vcf)
{
	if (vcf == NULL)
		return;

	logger_free (vcf->logger);

	memset (vcf, 0, sizeof (MakeVcf));
}

static int
make_vcf_validate (MakeVcf *vcf)
{
	int rc = EXIT_SUCCESS;

	/*
	* Validate arguments and mandatory options
	*/

	// DB must exit
	if (!exists (vcf->db_file))
		{
			fprintf (stderr, "%s: SQLite3 database file '%s': No such file\n",
					PACKAGE, vcf->db_file);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*
	* Validate options
	*/

	// If a FASTA file was passed, then it must exist
	if (vcf->fasta_file != NULL && !exists (vcf->fasta_file))
		{
			fprintf (stderr, "%s: --reference-file '%s': No such file\n",
					PACKAGE, vcf->fasta_file);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate near_gene_dist >= 0
	if (vcf->near_gene_dist < 0)
		{
			fprintf (stderr, "%s: --near-gene-dist must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	// Validate orientation_error >= 0.0
	if (vcf->orientation_error < 0.0)
		{
			fprintf (stderr, "%s: --orientation-error must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*
	* Final settings
	*/

	// If it's silent and no log file
	// was passsed, then set log_level
	// to LOG_ERROR - At least print
	// errors
	if (vcf->silent && vcf->log_file == NULL)
		{
			vcf->silent = 0;
			vcf->log_level = LOG_ERROR;
		}

	vcf->logger = logger_new (vcf->log_file, vcf->log_level,
			vcf->silent, 1);

Exit:
	return rc;
}

static void
make_vcf_print (MakeVcf *vcf)
{
	String *msg = NULL;

	msg = string_sized_new (BUFSIZ);

	string_concat_printf (msg, ">> Make VCF step <<\n"
		"\n"
		"#\n"
		"# %s\n"
		"#\n"
		"\n"
		"## Command line parsing with default values\n"
		"\n"
		"# Run %s\n"
		"$ %s make-vcf \\\n"
		"  --output-dir='%s' \\\n"
		"  --prefix='%s' \\\n",
		PACKAGE_STRING, PACKAGE, PACKAGE,
		vcf->output_dir, vcf->prefix);

	if (vcf->log_file != NULL)
		string_concat_printf (msg, "  --log-file='%s' \\\n",
				vcf->log_file);

	if (vcf->log_level <= LOG_DEBUG)
		string_concat_printf (msg, "  --debug \\\n");

	if (vcf->silent)
		string_concat_printf (msg, "  --silent \\\n");

	if (vcf->fasta_file != NULL)
		{
			string_concat_printf (msg,
				"  --reference-file=%s \\\n",
				vcf->fasta_file);
		}

	string_concat_printf (msg,
		"  --near-gene-dist=%li \\\n"
		"  --orientation-error=%.2f\n",
		vcf->near_gene_dist, vcf->orientation_error);

	log_info ("%s", msg->str);
	string_free (msg, 1);
}

int
parse_make_vcf_command_opt (int argc, char **argv)
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
		{"log-file",          required_argument, 0, 'l'},
		{"output-dir",        required_argument, 0, 'o'},
		{"prefix",            required_argument, 0, 'p'},
		{"near-gene-dist",    required_argument, 0, 'n'},
		{"orientation-error", required_argument, 0, 'e'},
		{"reference-file",    required_argument, 0, 'r'},
		{0,                   0,                 0,  0 }
	};

	// Init variables to default values
	MakeVcf vcf = {};
	make_vcf_init (&vcf);

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c = 0;

	while ((c = getopt_long (argc, argv, "hqdl:o:p:n:e:r:", opt, &option_index)) >= 0)
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
						vcf.silent = 1;
						break;
					}
				case 'd':
					{
						vcf.log_level = LOG_DEBUG;
						break;
					}
				case 'l':
					{
						vcf.log_file = optarg;
						break;
					}
				case 'o':
					{
						vcf.output_dir = optarg;
						break;
					}
				case 'p':
					{
						vcf.prefix = optarg;
						break;
					}
				case 'n':
					{
						vcf.near_gene_dist = atoi (optarg);
						break;
					}
				case 'e':
					{
						vcf.orientation_error = atof (optarg);
						break;
					}
				case 'r':
					{
						vcf.fasta_file = optarg;
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

	// Catch argument SQLite3 DB
	vcf.db_file = argv[optind + 1];

	// Validate and init logger
	rc = make_vcf_validate (&vcf);

	// If no error
	if (rc == EXIT_SUCCESS)
		{
			/*
			* RUN FOOLS
			*/
			make_vcf_print (&vcf);
			run (&vcf);
		}

Exit:
	make_vcf_destroy (&vcf);
	return rc;
}
