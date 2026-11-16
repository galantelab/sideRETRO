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
#include <getopt.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "str.h"
#include "db.h"
#include "log.h"
#include "logger.h"
#include "vcf.h"
#include "summary.h"
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
	double       min_qual;
	int          min_mapq;
	int          min_sr;
	int          min_breakpoints;
	int          min_tsd;
	int          max_tsd;
	int          disabled_filters;
	String      *disabled_ids;

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
	char *rnames_file = NULL;
	char *insertion_file = NULL;

	time_t start = time (NULL);
	Summary *s = NULL;
	char b1[SUMMARY_NUM_LEN] = {};
	char b2[SUMMARY_NUM_LEN] = {};

	// Assemble VCF, read names and
	// predicted insertions files
	xasprintf_concat (&vcf_file, "%s/%s.vcf",
			v->output_dir, v->prefix);

	xasprintf_concat (&rnames_file, "%s/%s.rnames.tsv",
			v->output_dir, v->prefix);

	xasprintf_concat (&insertion_file, "%s/%s.insertions.fa",
			v->output_dir, v->prefix);

	log_debug ("Create output dir '%s'", v->output_dir);
	mkdir_p (v->output_dir);

	log_info ("Open the database '%s'", v->db_file);
	db = db_connect (v->db_file);

	// Fill options
	VCFOption opt = {
		.near_gene_dist    = v->near_gene_dist,
		.orientation_error = v->orientation_error,
		.min_qual          = v->min_qual,
		.min_mapq          = v->min_mapq,
		.min_sr            = v->min_sr,
		.min_breakpoints   = v->min_breakpoints,
		.min_tsd           = v->min_tsd,
		.max_tsd           = v->max_tsd,
		.disabled_filters  = v->disabled_filters,
		.fasta_file        = v->fasta_file,
		.rnames_file       = rnames_file,
		.insertion_file    = insertion_file
	};

	vcf (db, vcf_file, &opt);

	s = summary_new ("make-vcf summary", start);

	summary_section (s, "Retrocopies");
	summary_line (s, "VCF records", "%s", summary_count (opt.num_records, b1));
	summary_line (s, "PASS", "%s (%s)", summary_count (opt.num_pass, b1),
			summary_percent (opt.num_pass, opt.num_records, b2));
	summary_line (s, "Predicted sequences", "%s%s",
			summary_count (opt.num_insertions, b1),
			opt.num_insertions == 0 && opt.num_records > 0
				? " (run merge-call with --reference-file)"
				: "");

	summary_section (s, "Output");
	summary_line (s, "VCF", "%s", vcf_file);
	summary_line (s, "Read names", "%s", rnames_file);
	summary_line (s, "Predicted insertions", "%s", insertion_file);

	if (v->fasta_file == NULL)
		summary_text (s, "Without --reference-file, REF is 'N' and TSD is not reported");

	summary_log (s);
	summary_free (s);
	xfree (vcf_file);
	xfree (rnames_file);
	xfree (insertion_file);
	db_close (db);
}

static void
print_usage (FILE *fp)
{
	int pkg_len = strlen (PACKAGE);
	fprintf (fp,
		"%s\n"
		"\n"
		"Usage: %s make-vcf [-h] [-q] [-d] [-l FILE] [-o DIR] [-p STR]\n"
		"       %*c          [-n INT] [-e FLOAT] [-r FILE] [-L FLOAT]\n"
		"       %*c          [-Q INT] [-x STR] <FILE>\n"
		"\n"
		"Write the annotated retrocopies to a VCF and to their companion files\n"
		"\n"
		"Examples:\n"
		"   $ sider vcf in.db\n"
		"   $ sider vcf -r hg38.fa -o result in.db\n"
		"\n"
		"Output:\n"
		"   <prefix>.vcf             The retrocopy insertions and their genotypes\n"
		"   <prefix>.rnames.tsv      The reads supporting each retrocopy, by sample\n"
		"   <prefix>.insertions.fa   The predicted retrocopy sequences\n"
		"\n"
		"Arguments:\n"
		"   SQLite3 database processed by 'merge-call'\n"
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
		"\n"
		"Annotation Options:\n"
		"   -n, --near-gene-dist       Maximum distance between an intergenic insertion\n"
		"                              and a gene to report it as near (NEAR)\n"
		"                              [default:\"%d\"]\n"
		"   -e, --orientation-error    Maximum p-value of the Spearman's rho to report\n"
		"                              the orientation (ORHO) [default:\"%.2f\"]\n"
		"   -r, --reference-file       FASTA file of the reference genome, for the REF\n"
		"                              column, the contigs and the TSD sequence\n"
		"\n"
		"Filter Options:\n"
		"   The FILTER column lists the checks that failed, or PASS: A\n"
		"   retrocopy with its insertion point well defined\n"
		"   -L, --min-qual             LowQual: Minimum QUAL [default:\"%g\"]\n"
		"   -Q, --min-mapq             LowMAPQ: Minimum mean mapping quality of the\n"
		"                              supporting reads [default:\"%d\"]\n"
		"   -S, --min-sr               LowSR: Minimum number of split reads at the\n"
		"                              breakpoints, whatever their side\n"
		"                              [default:\"%d\"]\n"
		"   -B, --min-breakpoints      NoTSD: Minimum number of breakpoints, 1 or 2.\n"
		"                              The target site duplication needs both\n"
		"                              [default:\"%d\"]\n"
		"   -t, --min-tsd              TSDOutlier: Minimum length of the target site\n"
		"                              duplication [default:\"%d\"]\n"
		"   -T, --max-tsd              TSDOutlier: Maximum length of the target site\n"
		"                              duplication [default:\"%d\"]\n"
		"   -x, --no-filter            Disable a filter. It may be passed multiple\n"
		"                              times: LowQual, LowMAPQ, AmbiguousParental,\n"
		"                              LowSR, NoTSD and TSDOutlier\n"
		"\n",
		PACKAGE_STRING, PACKAGE, pkg_len, ' ', pkg_len, ' ',
		DEFAULT_OUTPUT_DIR, DEFAULT_PREFIX, DEFAULT_NEAR_GENE_DIST,
		DEFAULT_ORIENTATION_ERROR, VCF_DEFAULT_MIN_QUAL, VCF_DEFAULT_MIN_MAPQ,
		VCF_DEFAULT_MIN_SR, VCF_DEFAULT_MIN_BKP, VCF_DEFAULT_MIN_TSD,
		VCF_DEFAULT_MAX_TSD);
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
		.min_qual           = VCF_DEFAULT_MIN_QUAL,
		.min_mapq           = VCF_DEFAULT_MIN_MAPQ,
		.min_sr             = VCF_DEFAULT_MIN_SR,
		.min_breakpoints    = VCF_DEFAULT_MIN_BKP,
		.min_tsd            = VCF_DEFAULT_MIN_TSD,
		.max_tsd            = VCF_DEFAULT_MAX_TSD,
		.disabled_filters   = 0,
		.disabled_ids       = string_new (""),
		.fasta_file         = NULL
	};
}

static void
make_vcf_destroy (MakeVcf *vcf)
{
	if (vcf == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	logger_free (vcf->logger);
	string_free (vcf->disabled_ids, 1);

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

	if (vcf->min_qual < 0.0)
		{
			fprintf (stderr, "%s: --min-qual must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (vcf->min_mapq < 0)
		{
			fprintf (stderr, "%s: --min-mapq must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (vcf->min_sr < 0)
		{
			fprintf (stderr, "%s: --min-sr must be greater or equal to 0\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (vcf->min_breakpoints < 1 || vcf->min_breakpoints > 2)
		{
			fprintf (stderr, "%s: --min-breakpoints must be 1 or 2\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (vcf->min_tsd < 1)
		{
			fprintf (stderr, "%s: --min-tsd must be greater or equal to 1\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	if (vcf->max_tsd < vcf->min_tsd)
		{
			fprintf (stderr, "%s: --max-tsd must be greater or equal to --min-tsd\n", PACKAGE);
			rc = EXIT_FAILURE; goto Exit;
		}

	/*
	* Final settings
	*/

	// If it's silent and no log file
	// was passed, then set log_level
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

	if (vcf->disabled_ids->len > 0)
		string_concat_printf (msg, "%s", vcf->disabled_ids->str);

	string_concat_printf (msg,
		"  --near-gene-dist=%li \\\n"
		"  --orientation-error=%.2f \\\n"
		"  --min-qual=%g \\\n"
		"  --min-mapq=%d \\\n"
		"  --min-sr=%d \\\n"
		"  --min-breakpoints=%d \\\n"
		"  --min-tsd=%d \\\n"
		"  --max-tsd=%d\n",
		vcf->near_gene_dist, vcf->orientation_error,
		vcf->min_qual, vcf->min_mapq, vcf->min_sr,
		vcf->min_breakpoints, vcf->min_tsd, vcf->max_tsd);

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
		{"min-qual",          required_argument, 0, 'L'},
		{"min-mapq",          required_argument, 0, 'Q'},
		{"min-sr",            required_argument, 0, 'S'},
		{"min-breakpoints",   required_argument, 0, 'B'},
		{"min-tsd",           required_argument, 0, 't'},
		{"max-tsd",           required_argument, 0, 'T'},
		{"no-filter",         required_argument, 0, 'x'},
		{0,                   0,                 0,  0 }
	};

	// Init variables to default values
	MakeVcf vcf = {};
	make_vcf_init (&vcf);

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c = 0;

	while ((c = getopt_long (argc, argv, "hqdl:o:p:n:e:r:L:Q:S:B:t:T:x:", opt, &option_index)) >= 0)
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
				case 'L':
					{
						vcf.min_qual = atof (optarg);
						break;
					}
				case 'Q':
					{
						vcf.min_mapq = atoi (optarg);
						break;
					}
				case 'S':
					{
						vcf.min_sr = atoi (optarg);
						break;
					}
				case 'B':
					{
						vcf.min_breakpoints = atoi (optarg);
						break;
					}
				case 't':
					{
						vcf.min_tsd = atoi (optarg);
						break;
					}
				case 'T':
					{
						vcf.max_tsd = atoi (optarg);
						break;
					}
				case 'x':
					{
						int flag = vcf_filter_lookup (optarg);

						if (!flag)
							{
								fprintf (stderr, "%s: --no-filter '%s': Unknown filter. "
										"Choose from: LowQual, LowMAPQ, AmbiguousParental, LowSR, "
										"NoTSD, TSDOutlier\n",
										PACKAGE, optarg);
								print_try_help (stderr);
								rc = EXIT_FAILURE; goto Exit;
							}

						vcf.disabled_filters |= flag;
						string_concat_printf (vcf.disabled_ids,
								"  --no-filter=%s \\\n", optarg);
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
