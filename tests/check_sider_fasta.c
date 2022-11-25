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
#include "../src/fasta.h"

static int num_line = 17;
static const char *fasta_body =
	";PONGA test\n"
	"; Comment and more comment\n"
	">chr1 example\n"
	"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
	"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
	"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
	"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
	"; Another comment\n"
	">chr2|example\n"
	"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
	"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
	"; It is possible?\n"
	"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
	"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
	"; Again\n"
	">chrN\n"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";

static const char *fasta_id1 = "chr1";
static const char *fasta_seq1 =
	"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
	"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

static const char *fasta_id2 = "chr2";
static const char *fasta_seq2 =
	"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
	"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

static const char *fasta_id3 = "chrN";
static const char *fasta_seq3 =
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
	"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

static void
create_fasta (const char *fasta, char *path)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, "%s", fasta);

	xfclose (fp);
}

START_TEST (test_fasta_read)
{
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;

	// Create gtf file
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";
	create_fasta (fasta_body, fasta_path);

	fasta = fasta_open_for_reading (fasta_path);
	entry = fasta_entry_new ();

	int i = 0;
	const int loops = 3;

	const char *contigs[] = {
		fasta_id1,
		fasta_id2,
		fasta_id3
	};

	const char *seqs[] = {
		fasta_seq1,
		fasta_seq2,
		fasta_seq3
	};

	while (fasta_read (fasta, entry))
		{
			ck_assert_str_eq (contigs[i], entry->contig->str);
			ck_assert_str_eq (seqs[i], entry->sequence->str);
			i++;
		}

	ck_assert_int_eq (loops, i);
	ck_assert_int_eq (num_line, entry->num_line);

	// Coverage
	fasta_close (NULL);
	fasta_entry_free (NULL);

	fasta_close (fasta);
	fasta_entry_free (entry);
	xunlink (fasta_path);
}
END_TEST

START_TEST (test_empty_file)
{
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";

	create_fasta ("", fasta_path);
	fasta = fasta_open_for_reading (fasta_path);
	entry = fasta_entry_new ();

	ck_assert_int_eq (fasta_read (fasta, entry), 0);

	fasta_entry_free (entry);
	fasta_close (fasta);
	xunlink (fasta_path);
}
END_TEST

START_TEST (test_open_fatal)
{
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";
	fasta_open_for_reading (fasta_path);
}
END_TEST

START_TEST (test_close_fatal)
{
	FastaFile *fasta;
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";

	create_fasta (fasta_body, fasta_path);
	fasta = fasta_open_for_reading (fasta_path);

	fasta_close (fasta);
	fasta_close (fasta);

	xunlink (fasta_path);
}
END_TEST

START_TEST (test_contig_fatal)
{
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";

	create_fasta ("ATCG\n", fasta_path);
	fasta = fasta_open_for_reading (fasta_path);
	entry = fasta_entry_new ();

	while (fasta_read (fasta, entry))
		;

	fasta_entry_free (entry);
	fasta_close (fasta);
	xunlink (fasta_path);
}
END_TEST

START_TEST (test_seq1_fatal)
{
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";

	create_fasta (">PONGA\n", fasta_path);
	fasta = fasta_open_for_reading (fasta_path);
	entry = fasta_entry_new ();

	while (fasta_read (fasta, entry))
		;

	fasta_entry_free (entry);
	fasta_close (fasta);
	xunlink (fasta_path);
}
END_TEST

START_TEST (test_seq2_fatal)
{
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";

	create_fasta (">PONGA1\n>PONGA2\n", fasta_path);
	fasta = fasta_open_for_reading (fasta_path);
	entry = fasta_entry_new ();

	while (fasta_read (fasta, entry))
		;

	fasta_entry_free (entry);
	fasta_close (fasta);
	xunlink (fasta_path);
}
END_TEST

Suite *
make_fasta_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;
	TCase *tc_segfault;

	s = suite_create ("FASTA");

	/* Core test case */
	tc_core = tcase_create ("Core");

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	/* Segfault test case */
	tc_segfault = tcase_create ("Segfault");

	tcase_add_test (tc_core, test_fasta_read);
	tcase_add_test (tc_core, test_empty_file);

	tcase_add_test_raise_signal (tc_abort, test_open_fatal,     SIGABRT);
	tcase_add_test_raise_signal (tc_abort, test_contig_fatal,   SIGABRT);
	tcase_add_test_raise_signal (tc_abort, test_seq1_fatal,     SIGABRT);
	tcase_add_test_raise_signal (tc_abort, test_seq2_fatal,     SIGABRT);

	tcase_add_test_raise_signal (tc_segfault, test_close_fatal, SIGSEGV);

	suite_add_tcase (s, tc_core);
	suite_add_tcase (s, tc_abort);
	suite_add_tcase (s, tc_segfault);

	return s;
}
