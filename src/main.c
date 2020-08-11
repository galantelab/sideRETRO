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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

#include "process_sample.h"
#include "merge_call.h"
#include "make_vcf.h"

static void
print_version (FILE *fp)
{
	fprintf (fp, "%s\n", PACKAGE_STRING);
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
		"A pipeline for detecting\n"
		"Somatic Insertion of DE novo RETROcopies\n"
		"\n"
		"Options:\n"
		"   -h, --help            Show help options\n"
		"   -v, --version         Show current version\n"
		"   -c, --cite            Show citation in BibTeX\n"
		"\n"
		"Commands:\n"
		"   ps,  process-sample   Extract alignments related\n"
		"                         an event of retrocopy\n"
		"   mc,  merge-call       Discover and annotate\n"
		"                         retrocopies\n"
		"   vcf, make-vcf         Generate VCF file with all\n"
		"                         annotate retrocopies\n"
		"\n",
		PACKAGE_STRING, PACKAGE, PACKAGE);
}

static void
print_citation (FILE *fp)
{
	fprintf (fp,
		"@article{10.1093/bioinformatics/btaa689,\n"
		"  author = {Miller, Thiago L A and Orpinelli, Fernanda and Buzzo, "
		"José Leonel L and Galante, Pedro A F},\n"
		"  title = \"{sideRETRO: a pipeline for identifying somatic and "
		"polymorphic insertions of processed pseudogenes or retrocopies}\",\n"
		"  journal = {Bioinformatics},\n"
		"  year = {2020},\n"
		"  month = {07},\n"
		"  issn = {1367-4803},\n"
		"  doi = {10.1093/bioinformatics/btaa689},\n"
		"  url = {https://doi.org/10.1093/bioinformatics/btaa689},\n"
		"  note = {btaa689},\n"
		"}\n");
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
		{"cite",    no_argument, 0, 'c'},
		{"version", no_argument, 0, 'v'},
		{"help",    no_argument, 0, 'h'},
		{0,         0,           0,  0 }
	};

	int rc = EXIT_SUCCESS;
	int option_index = 0;
	int c;

	while ((c = getopt_long (argc, argv, "cvh", opt, &option_index)) >= 0)
		{
			switch (c)
				{
				case 'c':
					{
						print_citation (stdout);
						goto Exit;
						break;
					}
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
	else if (!strcmp (argv[1], "ps") || !strcmp (argv[1], "process-sample"))
		rc = parse_process_sample_command_opt (argc, argv);
	else if (!strcmp (argv[1], "mc") || !strcmp (argv[1], "merge-call"))
		rc = parse_merge_call_command_opt (argc, argv);
	else if (!strcmp (argv[1], "vcf") || !strcmp (argv[1], "make-vcf"))
		rc = parse_make_vcf_command_opt (argc, argv);
	else
		{
			fprintf (stderr, "%s: '%s' is not a valid command\n", PACKAGE, argv[1]);
			print_try_help (stderr);
			rc = EXIT_FAILURE;
		}

	return rc;
}
