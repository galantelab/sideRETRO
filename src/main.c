#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

#include "htslib/version.h"
#include "sqlite3.h"

#include "process_sample.h"

static void
print_version (FILE *fp)
{
	fprintf (fp,
		"%s\n"
		"\n"
		"Using:\n"
		"  htslib  %s\n"
		"  sqlite3 %s\n",
		PACKAGE_STRING, HTS_VERSION, SQLITE_VERSION);
}

static void
print_usage (FILE *fp)
{
	fprintf (fp,
		"%s\n"
		"\n"
		"Usage: %s [-hv]\n"
		"       %s <command> [options]\n"
		"\n"
		"Structural-Variation Intercorrences Detection\n"
		"Environment for Retrocopies\n"
		"\n"
		"Options:\n"
		"   -h, --help       Show help options\n"
		"   -v, --version    Show current version\n"
		"\n"
		"Commands\n"
		"   process-sample\n"
		"   merge-call\n",
		PACKAGE_STRING, PACKAGE, PACKAGE);
}

static void
print_try_help (FILE *fp)
{
	fprintf (fp, "Try '%s  --help' for more information\n",
			PACKAGE);
}

static int
parse_no_command_opt (int argc, char **argv)
{
	assert (argc > 1 && argv != NULL && *argv != NULL);

	struct option opt[] =
	{
		{"version", no_argument, 0, 'v'},
		{"help",    no_argument, 0, 'h'},
		{0,         0,           0,  0 }
	};

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c;

	while ((c = getopt_long (argc, argv, "vh", opt, &option_index)) >= 0)
		{
			switch (c)
				{
				case 'v':
					{
						print_version (stdout);
						goto Exit;
						break;
					}
				case 'h':
					{
						print_usage (stdout);
						goto Exit;
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

Exit:
	return rc;
}

int
main (int argc, char *argv[])
{
	int rc = EXIT_SUCCESS;

	if (argc == 1)
		print_usage (stdout);
	else if (argv[1][0] == '-')
		rc = parse_no_command_opt (argc, argv);
	else if (!strcmp (argv[1], "process-sample"))
		rc = parse_process_sample_command_opt (argc, argv);
	else if (!strcmp (argv[1], "merge-call"))
		printf ("** Under construction **\n");
	else
		{
			fprintf (stderr, "%s: '%s' is not a valid command\n", PACKAGE, argv[1]);
			print_try_help (stderr);
			rc = EXIT_FAILURE;
		}

	return rc;
}
