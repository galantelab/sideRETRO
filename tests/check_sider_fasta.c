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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"
#include <zlib.h>
#include <sys/stat.h>
#include "../src/hash.h"
#include "../src/list.h"

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

static void
create_truncated_gz (const char *body, char *path)
{
	gzFile gzfp = NULL;
	struct stat st = {};
	int fd = 0;
	int i = 0;

	// A gzip file cut in the middle of its stream
	fd = xmkstemp (path);
	close (fd);

	gzfp = gzopen (path, "wb");
	for (i = 0; i < 2000; i++)
		{
			// Vary the lines, so the stream is long enough
			gzputs (gzfp, body);
			gzprintf (gzfp, "%s%d\n", body[0] == '>' ? "ACGT" : "#", i * 7919);
		}
	gzclose (gzfp);

	ck_assert_int_eq (stat (path, &st), 0);
	ck_assert_int_gt (st.st_size, 1000);
	ck_assert_int_eq (truncate (path, st.st_size / 2), 0);
}

START_TEST (test_truncated_fatal)
{
	FastaFile *fasta = NULL;
	FastaEntry *entry = fasta_entry_new ();
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";

	create_truncated_gz (">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n", fasta_path);
	fasta = fasta_open_for_reading (fasta_path);

	while (fasta_read (fasta, entry))
		;
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

START_TEST (test_fasta_index_regions)
{
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";
	List *regions = list_new ((DestroyNotify) fasta_region_free);
	Hash *contig_len = NULL;
	Hash *idx = NULL;
	char *key = NULL;

	create_fasta (fasta_body, fasta_path);

	// Inside chr1, outside chr1 and an absent contig
	list_append (regions, fasta_region_new (fasta_id1, 1, 4));
	list_append (regions, fasta_region_new (fasta_id1, 1, 1000000));
	list_append (regions, fasta_region_new ("chrPonga", 1, 4));

	idx = fasta_index_regions (fasta_path, regions, &contig_len);

	ck_assert_int_eq (hash_size (idx), 1);

	key = fasta_region_key (fasta_id1, 1, 4);
	ck_assert (hash_lookup (idx, key) != NULL);
	xfree (key);

	ck_assert_int_eq (hash_size (contig_len), 3);
	ck_assert_int_eq (* (long *) hash_lookup (contig_len, fasta_id1),
			strlen (fasta_seq1));

	hash_free (idx);
	hash_free (contig_len);
	list_free (regions);
	xunlink (fasta_path);
}
END_TEST

START_TEST (test_fasta_index_regions_empty_fatal)
{
	char fasta_path[] = "/tmp/ponga.fa.XXXXXX";
	List *regions = list_new (NULL);
	Hash *contig_len = NULL;

	create_fasta ("", fasta_path);
	fasta_index_regions (fasta_path, regions, &contig_len);
}
END_TEST

START_TEST (test_fasta_region_free_null)
{
	fasta_region_free (NULL);
}
END_TEST

static int
check_format_of (int (*check) (const char *), const char *content)
{
	char path[] = "/tmp/ponga_check-XXXXXX";
	FILE *fp = NULL;
	int rc = 0;

	fp = xfdopen (xmkstemp (path), "w");
	xfprintf (fp, "%s", content);
	xfclose (fp);

	rc = check (path);
	xunlink (path);

	return rc;
}

START_TEST (test_fasta_check_format)
{
	// After comments and blank lines, a contig
	ck_assert_int_eq (check_format_of (fasta_check_format,
				">chr1\nACGT\n"), 1);
	ck_assert_int_eq (check_format_of (fasta_check_format,
				";PONGA\n\n>chr1 example\nACGT\n"), 1);

	// Not a FASTA: Text, a contig without name, nothing
	ck_assert_int_eq (check_format_of (fasta_check_format,
				"chr1\t1\t100\n"), 0);
	ck_assert_int_eq (check_format_of (fasta_check_format,
				">\nACGT\n"), 0);
	ck_assert_int_eq (check_format_of (fasta_check_format, ""), 0);
}
END_TEST

Suite *
make_fasta_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("FASTA");

	/* Core test case */
	tc_core = tcase_create ("Core");

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_test (tc_core, test_fasta_check_format);
	tcase_add_test (tc_core, test_fasta_read);
	tcase_add_test (tc_core, test_empty_file);
	tcase_add_test (tc_core, test_fasta_region_free_null);
	tcase_add_test (tc_core, test_fasta_index_regions);

	tcase_add_exit_test (tc_abort, test_open_fatal,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_contig_fatal, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_seq1_fatal,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_seq2_fatal,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_truncated_fatal,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_fasta_index_regions_empty_fatal, EXIT_FAILURE);

	suite_add_tcase (s, tc_core);
	suite_add_tcase (s, tc_abort);

	return s;
}
