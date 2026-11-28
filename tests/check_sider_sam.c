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
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <check.h>
#include "check_sider.h"
#include <sys/resource.h>
#include <signal.h>
#include <unistd.h>

#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/sam.h"

static const char *sam_str =
	"@SQ\tSN:PONGA\tLN:520052\n"
	"@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188	CL:bwa mem -t 1 ponga/ponga.fa ponga.fastq\n"
	"ponga1\t0\tPONGA\t35361\t0\t10M\t*\t0\t0\tATCGATCGAT\t!!!!!!!!!!\tNM:i:0\tMD:Z:52\tAS:i:52\tXS:i:52\n"
	"ponga2\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGAA\t!!!!!!!!!!\tAS:i:0\tXS:i:0\n"
	"ponga3\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGAA\t!!!!!!!!!!\tAS:i:0\tXS:i:0\n"
	"ponga4\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGAA\t!!!!!!!!!!\tAS:i:0\tXS:i:0\n"
	"ponga5\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGAA\t!!!!!!!!!!\tAS:i:0\tXS:i:0\n"
	"ponga6\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGAA\t!!!!!!!!!!\tAS:i:0\tXS:i:0\n"
	"ponga7\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGAA\t!!!!!!!!!!\tAS:i:0\tXS:i:0\n";

static void
create_sam_fd (int fd)
{
	FILE *fp = xfdopen (fd, "w");
	fprintf (fp, sam_str, "");
	xfclose (fp);
}

START_TEST (test_sam_to_bam)
{
	char sam_file1[] = "/tmp/ponga_sam-XXXXXX";
	char sam_file2[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[]  = "/tmp/ponga_bam-XXXXXX";

	int sam_fd1 = xmkstemp (sam_file1);
	create_sam_fd (sam_fd1);

	xmkstemp (sam_file2);
	xmkstemp (bam_file);

	ck_assert_int_eq (sam_to_bam (sam_file1, bam_file), 1);

	if (which ("samtools"))
		{
			char *cmd1 = NULL;
			char *cmd2 = NULL;

			xasprintf (&cmd1, "samtools view -h %s -o %s 2> /dev/null",
					bam_file, sam_file2);
			xasprintf (&cmd2, "diff -q %s %s 2>&1 > /dev/null",
					sam_file1, sam_file2);

			ck_assert_int_eq (system (cmd1), 0);
			ck_assert_int_eq (system (cmd2), 0);

			xfree (cmd1);
			xfree (cmd2);
		}

	xunlink (sam_file1);
	xunlink (sam_file2);
	xunlink (bam_file);
}
END_TEST

START_TEST (test_sam_to_bam_fp)
{
	char sam_file[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[] = "/tmp/ponga_bam-XXXXXX";

	int sam_fd = xmkstemp (sam_file);
	create_sam_fd (sam_fd);

	xmkstemp (bam_file);

	FILE *fp = xfopen (sam_file, "r");
	ck_assert_int_eq (sam_to_bam_fp (fp, bam_file), 1);

	xfclose (fp);

	xunlink (sam_file);
	xunlink (bam_file);
}
END_TEST

START_TEST (test_sam_test_sorted_order)
{
	const char *queryname = "@HD\tVN:1.0\tSO:queryname\n";
	const char *coordinate = "@HD\tVN:1.6\tSO:coordinate\n";
	const char *swapped = "@HD\tSO:queryname\tGO:query\tVN:1.6\n";
	const char *unsorted = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
	const char *no_hd = "@SQ\tSN:chr1\tLN:100\n";

	bam_hdr_t *hdr = NULL;

	hdr = sam_hdr_parse (strlen (queryname), queryname);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "queryname"), 1);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "coordinate"), 0);
	bam_hdr_destroy (hdr);

	hdr = sam_hdr_parse (strlen (coordinate), coordinate);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "coordinate"), 1);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "queryname"), 0);
	bam_hdr_destroy (hdr);

	// The tag SO at any place of @HD
	hdr = sam_hdr_parse (strlen (swapped), swapped);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "queryname"), 1);
	bam_hdr_destroy (hdr);

	// Without SO, or without @HD at all
	hdr = sam_hdr_parse (strlen (unsorted), unsorted);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "queryname"), 0);
	bam_hdr_destroy (hdr);

	hdr = sam_hdr_parse (strlen (no_hd), no_hd);
	ck_assert_int_eq (sam_test_sorted_order (hdr, "queryname"), 0);
	bam_hdr_destroy (hdr);
}
END_TEST

START_TEST (test_sam_strerror)
{
	ck_assert_str_eq (sam_strerror ("ponga.bam", 0), "Failed to open");
	ck_assert_str_eq (sam_strerror ("ponga.bam", ENOENT), strerror (ENOENT));
	ck_assert_str_eq (sam_strerror ("https://ponga.invalid/ponga.bam", EDESTADDRREQ),
			"Could not resolve the host");

	// Only a URL has a host
	ck_assert_str_eq (sam_strerror ("ponga.bam", EDESTADDRREQ),
			strerror (EDESTADDRREQ));
}
END_TEST

static void
write_file (char *path, const char *body)
{
	int fd = xmkstemp (path);
	FILE *fp = xfdopen (fd, "w");
	xfprintf (fp, "%s", body);
	xfclose (fp);
}

START_TEST (test_sam_to_bam_fp_closed_fatal)
{
	// The fd under the stream is closed: It cannot be duplicated
	char sam[] = "/tmp/ponga.sam.XXXXXX";
	FILE *fp = NULL;

	write_file (sam, sam_str);
	fp = xfopen (sam, "r");
	close (fileno (fp));

	sam_to_bam_fp (fp, "/tmp/ponga.bam");
}
END_TEST

START_TEST (test_sam_to_bam_output_fatal)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";

	write_file (sam, sam_str);
	sam_to_bam (sam, "/it-must-not-exist/ponga.bam");
}
END_TEST

START_TEST (test_sam_to_bam_faults)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";
	char bam[] = "/tmp/ponga.bam.XXXXXX";
	char bad[] = "/tmp/ponga.bad.XXXXXX";
	char out[] = "/tmp/ponga.out.XXXXXX";
	int fd = 0;

	write_file (sam, sam_str);

	// Truncated BAM as input: The header cannot be read
	fd = xmkstemp (bam);
	close (fd);
	ck_assert_int_eq (sam_to_bam (sam, bam), 1);
	ck_assert_int_eq (truncate (bam, 20), 0);

	fd = xmkstemp (out);
	close (fd);
	ck_assert_int_eq (sam_to_bam (bam, out), 0);

	// A broken record
	write_file (bad,
		"@SQ\tSN:PONGA\tLN:520052\n"
		"ponga1\tponga\tPONGA\t35361\t0\t10M\t*\t0\t0\tATCGATCGAT\t!!!!!!!!!!\n");
	ck_assert_int_eq (sam_to_bam (bad, out), 0);

	xunlink (sam);
	xunlink (bam);
	xunlink (bad);
	xunlink (out);
}
END_TEST

START_TEST (test_sam_to_bam_full)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";

	// /dev/full as output: The header cannot be written,
	// neither the output stream closed
	write_file (sam, sam_str);
	ck_assert_int_eq (sam_to_bam (sam, "/dev/full"), 0);
	xunlink (sam);
}
END_TEST

START_TEST (test_sam_to_bam_write_fail)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";
	char out[] = "/tmp/ponga.out.XXXXXX";
	struct rlimit limit = {.rlim_cur = 4096, .rlim_max = RLIM_INFINITY};
	FILE *fp = NULL;
	int fd = 0;
	int i = 0;

	// Many records, so the BGZF blocks are flushed while
	// writing them, after the header
	fd = xmkstemp (sam);
	fp = xfdopen (fd, "w");
	xfprintf (fp, "@SQ\tSN:PONGA\tLN:520052\n");
	for (i = 0; i < 20000; i++)
		xfprintf (fp, "ponga%d\t0\tPONGA\t%d\t60\t10M\t*\t0\t0\tATCGATCGAT\t!!!!!!!!!!\n",
				i, i + 1);
	xfclose (fp);

	fd = xmkstemp (out);
	close (fd);

	// The output file cannot grow beyond 4KiB: The header fits,
	// the records do not (EFBIG instead of SIGXFSZ)
	signal (SIGXFSZ, SIG_IGN);
	ck_assert_int_eq (setrlimit (RLIMIT_FSIZE, &limit), 0);

	ck_assert_int_eq (sam_to_bam (sam, out), 0);

	xunlink (sam);
	xunlink (out);
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

START_TEST (test_sam_check_format)
{
	char sam_file[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[] = "/tmp/ponga_bam-XXXXXX";

	// SAM with header, and without it: An alignment
	ck_assert_int_eq (check_format_of (sam_check_format, sam_str), 1);
	ck_assert_int_eq (check_format_of (sam_check_format,
				"ponga1\t0\tPONGA\t35361\t0\t10M\t*\t0\t0\tATCGATCGAT\t!!!!!!!!!!\n"), 1);

	// BAM
	create_sam_fd (xmkstemp (sam_file));
	xmkstemp (bam_file);
	ck_assert_int_eq (sam_to_bam (sam_file, bam_file), 1);
	ck_assert_int_eq (sam_check_format (bam_file), 1);

	// htslib takes any text for SAM: Not FASTA, GTF, nothing
	ck_assert_int_eq (check_format_of (sam_check_format,
				">chr1\nACGT\n"), 0);
	ck_assert_int_eq (check_format_of (sam_check_format,
				"chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_name \"DDX11L1\";\n"), 0);
	ck_assert_int_eq (check_format_of (sam_check_format, ""), 0);

	// It cannot open them: errno tells why
	errno = 0;
	ck_assert_int_eq (sam_check_format ("/tmp/pongatrongaflofa.bam"), -1);
	ck_assert_int_eq (errno, ENOENT);

	errno = 0;
	ck_assert_int_eq (sam_check_format ("ponga://x/ponga.bam"), -1);
	ck_assert_int_eq (errno, EPROTONOSUPPORT);

	xunlink (sam_file);
	xunlink (bam_file);
}
END_TEST

Suite *
make_sam_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_write;
	TCase *tc_abort;

	s = suite_create ("SAM");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_sam_check_format);
	tcase_add_test (tc_core, test_sam_strerror);
	tcase_add_test (tc_core, test_sam_to_bam);
	tcase_add_test (tc_core, test_sam_to_bam_fp);
	tcase_add_test (tc_core, test_sam_test_sorted_order);
	tcase_add_test (tc_core, test_sam_to_bam_faults);
	suite_add_tcase (s, tc_core);

	/*
	* Write failures: htslib (bgzf_close, 1.9 to 1.25) does not release
	* the output stream when its last flush fails - a leak that
	* is not ours. Skip them under valgrind or AddressSanitizer
	* (CK_EXCLUDE_TAGS=htslib-leak)
	*/
	tc_write = tcase_create ("Write");
	tcase_set_tags (tc_write, "htslib-leak");

	tcase_add_test (tc_write, test_sam_to_bam_full);
	tcase_add_test (tc_write, test_sam_to_bam_write_fail);
	suite_add_tcase (s, tc_write);

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_exit_test (tc_abort, test_sam_to_bam_output_fatal,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_sam_to_bam_fp_closed_fatal, EXIT_FAILURE);
	suite_add_tcase (s, tc_abort);

	return s;
}
