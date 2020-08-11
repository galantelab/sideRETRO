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
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/gff.h"

static void
handle_sigabrt (int sig)
{
	if (sig == SIGABRT)
		exit (1);
}

static const char *gff_header =
	"##description: evidence-based annotation of the human genome (GRCh38), version 30 (Ensembl 96)\n"
	"##provider: GENCODE\n"
	"##contact: gencode-help@ebi.ac.uk\n"
	"##format: gff"
	"##date: 2019-03-28";

static const char *gff_body =
	"chr1	.	gene	1000	1100	100	+	.	gene_name=g1;gene_id=ENSG1;transcript_id=ENST1;gene_type=lincRNA;transcript_type=protein_coding;\n"
	"##provider: GENCODE\n"
	"chr1	.	transcript	1	500	.	-	.	gene_name \"ponga\"; gene_desc \"this is only a test\"; type \"joke\";\n"
	"chr1	.	transcript	1000	1100	.	+	.	gene_name=g1;gene_id=ENSG1;transcript_id=ENST1;transcript_type=protein_coding;\n"
	"chr1	.	exon	1000	1100	.	+	1	gene_name=g1;gene_id=ENSG1;transcript_id=ENST1;transcript_type=protein_coding;exon_id=ENSE1;\n"
	"chr1	.	gene	1000	1100	.	+	.	gene_name=g1;gene_id=ENSG1;transcript_id=ENST1;gene_type=lincRNA;ponga=must_fail;\n"
	"chr1	.	gene	1000	1100	.	+	.	gene_name=g1;gene_id=ENSG1;transcript_id=ENST1;transcript_type=protein_coding;\n"
	"chr1	.	gene	1000	1100	.	+	.	gene_name=g1;gene_id=ENSG1;transcript_id=ENST1;gene_type=protein_coding;transcript_type=protein_coding;\n";

static const char *gff_seqname_fatal =
	"##description: evidence-based annotation of the human genome (GRCh38), version 30 (Ensembl 96)\n"
	"		\n";

static const char *gff_source_fatal =
	"chr1		\n";

static const char *gff_feature_fatal =
	"chr1	HAVANA		\n";

static const char *gff_start_fatal =
	"chr1	HAVANA	gene		\n";

static const char *gff_end_fatal =
	"chr1	HAVANA	gene	1		\n";

static const char *gff_score_fatal =
	"chr1	HAVANA	gene	1	2		\n";

static const char *gff_strand_fatal =
	"chr1	HAVANA	gene	1	2	.		\n";

static const char *gff_frame_fatal =
	"chr1	HAVANA	gene	1	2	.	+		\n";

static const char *gff_key_value_fatal =
	"chr1	HAVANA	gene	1	2	.	+	.	gene_name=\n";

static void
create_gff (const char *gff, char *path)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, "%s", gff);

	xfclose (fp);
}

START_TEST (test_gff_header)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_header, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	ck_assert (gff->header != NULL);
	ck_assert (gff_read (gff, entry) == 0);

	gff_close (gff);
	gff_entry_free (entry);

	// For coverage
	gff_close (NULL);
	gff_entry_free (NULL);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_gff_read)
{
	// Our heroes!
	GffFile *gff = NULL;
	GffEntry *entry = NULL;
	GffEntry *dup = NULL;

	// Create gtf file
	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_body, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	ck_assert (gff != NULL && entry != NULL);

	while (gff_read (gff, entry))
		;

	dup = gff_entry_dup (entry);

	// Cleanup
	gff_entry_free (dup);
	gff_entry_free (entry);
	gff_close (gff);
	xunlink (gff_path);
}
END_TEST

START_TEST (test_open_fatal)
{
	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	gff_open_for_reading (gff_path);
}
END_TEST

START_TEST (test_close_fatal)
{
	GffFile *gff;
	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";

	create_gff (gff_header, gff_path);
	gff = gff_open_for_reading (gff_path);

	gff_close (gff);
	gff_close (gff);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_seqname_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_seqname_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_source_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_source_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_feature_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_feature_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_start_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_start_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_end_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_end_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_score_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_score_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_strand_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_strand_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_frame_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_frame_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_key_value_fatal)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_key_value_fatal, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	while (gff_read (gff, entry))
		;

	gff_close (gff);
	gff_entry_free (entry);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_gff_filter)
{
	GffFile *gff = NULL;
	GffEntry *entry = NULL;
	GffFilter *filter = NULL;

	int i = 0;

	char gff_path[] = "/tmp/ponga.gff3.XXXXXX";
	create_gff (gff_body, gff_path);

	gff = gff_open_for_reading (gff_path);
	entry = gff_entry_new ();

	// Only hard attributes
	filter = gff_filter_new ();
	gff_filter_insert_feature (filter, NULL);
	gff_filter_insert_feature (filter, "ponga");
	gff_filter_insert_feature (filter, "gene");
	gff_filter_insert_hard_attribute (filter, "transcript_type", "protein_coding");
	gff_filter_insert_hard_attribute (filter, "gene_id", "ENS");

	while (gff_read_filtered (gff, entry, filter))
		;

	// Rewind
	gff_filter_free (filter);
	filter = gff_filter_new ();
	gff_close (gff);
	gff = gff_open_for_reading (gff_path);

	// Hard and soft attributes
	gff_filter_insert_feature (filter, "gene");
	gff_filter_insert_hard_attribute (filter, "transcript_type", "protein_coding");
	gff_filter_insert_hard_attribute (filter, "gene_id", "ENS");
	gff_filter_insert_soft_attribute (filter, "non_exist", "ponga");
	gff_filter_insert_soft_attribute (filter, "gene_type", "protein_coding");
	gff_filter_insert_soft_attribute (filter, "gene_type", "lincRNA");

	while (gff_read_filtered (gff, entry, filter))
		i++;

	ck_assert_int_eq (i, 2);

	gff_close (gff);
	gff_entry_free (entry);
	gff_filter_free (filter);

	// Coverage
	gff_filter_free (NULL);

	xunlink (gff_path);
}
END_TEST

START_TEST (test_gff_looks_like_gff_file)
{
	const char *filenames[10] = {
		"ponga.gff3",
		"ponga.gtf",
		"ponga.gff",
		"ponga.gff3.gz",
		"ponga.gtf.gz",
		"ponga.gff.gz",
		"ponga.bed",
		"ponga.bed.gz",
		"ponga.gff3.bed.gz",
		"ponga.txt"
	};

	const int true_positives[10] = {
		1, 1, 1, 1, 1, 1, 0, 0, 0, 0
	};

	int i = 0;
	for (; i < 10; i++)
		ck_assert_int_eq (
				gff_looks_like_gff_file (filenames[i]),
				true_positives[i]);
}
END_TEST

Suite *
make_gff_suite (void)
{
	setup_signal (SIGABRT, handle_sigabrt);

	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("GFF");

	/* Core test case */
	tc_core = tcase_create ("Core");

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_test (tc_core, test_gff_header);
	tcase_add_test (tc_core, test_gff_read);
	tcase_add_test (tc_core, test_gff_filter);
	tcase_add_test (tc_core, test_gff_looks_like_gff_file);

	tcase_add_exit_test (tc_abort, test_open_fatal, 1);
	tcase_add_exit_test (tc_abort, test_close_fatal, 1);
	tcase_add_exit_test (tc_abort, test_seqname_fatal, 1);
	tcase_add_exit_test (tc_abort, test_source_fatal, 1);
	tcase_add_exit_test (tc_abort, test_feature_fatal, 1);
	tcase_add_exit_test (tc_abort, test_start_fatal, 1);
	tcase_add_exit_test (tc_abort, test_end_fatal, 1);
	tcase_add_exit_test (tc_abort, test_score_fatal, 1);
	tcase_add_exit_test (tc_abort, test_strand_fatal, 1);
	tcase_add_exit_test (tc_abort, test_frame_fatal, 1);
	tcase_add_exit_test (tc_abort, test_key_value_fatal, 1);

	suite_add_tcase (s, tc_core);
	suite_add_tcase (s, tc_abort);

	return s;
}
