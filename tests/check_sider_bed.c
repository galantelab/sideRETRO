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

#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/bed.h"

static const char *bed_header =
	"browser position chr7:127471196-127495720\n"
	"browser hide all\n"
	"track name=\"ColorByStrandDemo\" description=\"Color by strand demonstration\""
	"track name=\"ColorByStrandDemo\" description=\"Color by strand demonstration\"";

static const char *bed_chrom_fatal =
	"browser position chr7:127471196-127495720\n"
	" \n";

static const char *bed_chrom_start_fatal =
	"chr7\n";

static const char *bed_chrom_end_fatal =
	"chr7 127471196\n";

#define BED1 \
	"chr7    127471196  127472363\n" \
	"chr7    127472363  127473530\n" \
	"chr7    127473530  127474697\n" \
	"chr7    127474697  127475864\n" \
	"chr7    127475864  127477031\n" \
	"chr7    127477031  127478198\n" \
	"chr7    127478198  127479365\n" \
	"chr7    127479365  127480532\n" \
	"chr7    127480532  127481699\n"

#define BED2 \
	"chr7    127471196  127472363  Pos1\n" \
	"chr7    127472363  127473530  Pos2\n" \
	"chr7    127473530  127474697  Pos3\n" \
	"chr7    127474697  127475864  Pos4\n" \
	"chr7    127475864  127477031  Neg1\n" \
	"chr7    127477031  127478198  Neg2\n" \
	"chr7    127478198  127479365  Neg3\n" \
	"chr7    127479365  127480532  Pos5\n" \
	"chr7    127480532  127481699  Neg4\n"

#define BED3 \
	"chr7    127471196  127472363  Pos1  0\n" \
	"chr7    127472363  127473530  Pos2  0\n" \
	"chr7    127473530  127474697  Pos3  0\n" \
	"chr7    127474697  127475864  Pos4  0\n" \
	"chr7    127475864  127477031  Neg1  0\n" \
	"chr7    127477031  127478198  Neg2  0\n" \
	"chr7    127478198  127479365  Neg3  0\n" \
	"chr7    127479365  127480532  Pos5  0\n" \
	"chr7    127480532  127481699  Neg4  0\n"

#define BED4 \
	"chr7    127471196  127472363  Pos1  0  +\n" \
	"chr7    127472363  127473530  Pos2  -1  +\n" \
	"chr7    127473530  127474697  Pos3  0  +\n" \
	"chr7    127474697  127475864  Pos4  0  +\n" \
	"chr7    127475864  127477031  Neg1  0  -\n" \
	"chr7    127477031  127478198  Neg2  2000  -\n" \
	"chr7    127478198  127479365  Neg3  0  -\n" \
	"chr7    127479365  127480532  Pos5  0  +\n" \
	"chr7    127480532  127481699  Neg4  0  -\n"

#define BED5 \
	"chr7  127471196  127472363  Pos1  0  +  127471196\n" \
	"chr7  127472363  127473530  Pos2  0  +  127472363\n" \
	"chr7  127473530  127474697  Pos3  0  +  127473530\n" \
	"chr7  127474697  127475864  Pos4  0  +  127474697\n" \
	"chr7  127475864  127477031  Neg1  0  -  127475864\n" \
	"chr7  127477031  127478198  Neg2  0  -  127477031\n" \
	"chr7  127478198  127479365  Neg3  0  -  127478198\n" \
	"chr7  127479365  127480532  Pos5  0  +  127479365\n" \
	"chr7  127480532  127481699  Neg4  0  -  127480532\n"

#define BED6 \
	"chr7  127471196  127472363  Pos1  0  +  127471196  127472363\n" \
	"chr7  127472363  127473530  Pos2  0  +  127472363  127473530\n" \
	"chr7  127473530  127474697  Pos3  0  +  127473530  127474697\n" \
	"chr7  127474697  127475864  Pos4  0  +  127474697  127475864\n" \
	"chr7  127475864  127477031  Neg1  0  -  127475864  127477031\n" \
	"chr7  127477031  127478198  Neg2  0  -  127477031  127478198\n" \
	"chr7  127478198  127479365  Neg3  0  -  127478198  127479365\n" \
	"chr7  127479365  127480532  Pos5  0  +  127479365  127480532\n" \
	"chr7  127480532  127481699  Neg4  0  -  127480532  127481699\n"

#define BED7 \
	"chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0\n" \
	"chr7  127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0\n" \
	"chr7  127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0\n" \
	"chr7  127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0\n" \
	"chr7  127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255,0\n" \
	"chr7  127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255\n" \
	"chr7  127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255\n" \
	"chr7  127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0\n" \
	"chr7  127480532  127481699  Neg4  0  -  127480532  127481699  0,0,256\n"

#define BED8 \
	"track name=pairedReads description=\"Clone Paired Reads\" useScore=1\n" \
	"chr22 1000 5000 cloneA 960 + 1000 5000 0 2\n" \
	"\n" \
	"chr22 2000 6000 cloneB 900 - 2000 6000 0 2\n"

#define BED9 \
	"track name=pairedReads description=\"Clone Paired Reads\" useScore=1\n" \
	"chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488,100\n" \
	"chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399\n"

#define BED10 \
	"track name=pairedReads description=\"Clone Paired Reads\" useScore=1\n" \
	"chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488 0,3512\n" \
	"chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399 0,3601,10\n"

static const char *beds[10] = {BED1, BED2, BED3, BED4, BED5, BED6, BED7, BED8, BED9, BED10};

static void
create_bed (const char *bed, char *path)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, "%s", bed);

	xfclose (fp);
}

START_TEST (test_bed_header)
{
	BedFile *bed = NULL;
	BedEntry *entry = NULL;

	char bed_path[] = "/tmp/ponga.bed.XXXXXX";
	create_bed (bed_header, bed_path);

	bed = bed_open_for_reading (bed_path);
	entry = bed_entry_new ();

	ck_assert (bed->header != NULL);
	ck_assert (bed_read (bed, entry) == 0);

	bed_close (bed);
	bed_entry_free (entry);

	// For coverage
	bed_close (NULL);
	bed_entry_free (NULL);

	xunlink (bed_path);
}
END_TEST

START_TEST (test_beds)
{
	BedFile *bed = NULL;
	BedEntry *entry = NULL;

	char bed_path[] = "/tmp/ponga.bed.XXXXXX";
	create_bed (beds[_i], bed_path);

	bed = bed_open_for_reading (bed_path);
	entry = bed_entry_new ();

	while (bed_read (bed, entry))
		;

	ck_assert_int_gt (entry->num_line, 0);

	bed_close (bed);
	bed_entry_free (entry);

	unlink (bed_path);
}
END_TEST

START_TEST (test_open_fatal)
{
	char bed_path[] = "/tmp/ponga.bed.XXXXXX";
	bed_open_for_reading (bed_path);
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
	BedFile *bed = NULL;
	BedEntry *entry = bed_entry_new ();
	char bed_path[] = "/tmp/ponga.bed.XXXXXX";

	create_truncated_gz ("chr1\t0\t100\tponga\t0\t+\n", bed_path);
	bed = bed_open_for_reading (bed_path);

	while (bed_read (bed, entry))
		;
}
END_TEST

START_TEST (test_chrom_fatal)
{
	BedFile *bed = NULL;
	BedEntry *entry = NULL;

	char bed_path[] = "/tmp/ponga.bed.XXXXXX";
	create_bed (bed_chrom_fatal, bed_path);

	bed = bed_open_for_reading (bed_path);
	entry = bed_entry_new ();

	while (bed_read (bed, entry))
		;

	bed_close (bed);
	bed_entry_free (entry);

	xunlink (bed_path);
}
END_TEST

START_TEST (test_chrom_start_fatal)
{
	BedFile *bed = NULL;
	BedEntry *entry = NULL;

	char bed_path[] = "/tmp/ponga.bed.XXXXXX";
	create_bed (bed_chrom_start_fatal, bed_path);

	bed = bed_open_for_reading (bed_path);
	entry = bed_entry_new ();

	while (bed_read (bed, entry))
		;

	bed_close (bed);
	bed_entry_free (entry);

	xunlink (bed_path);
}
END_TEST

START_TEST (test_chrom_end_fatal)
{
	BedFile *bed = NULL;
	BedEntry *entry = NULL;

	char bed_path[] = "/tmp/ponga.bed.XXXXXX";
	create_bed (bed_chrom_end_fatal, bed_path);

	bed = bed_open_for_reading (bed_path);
	entry = bed_entry_new ();

	while (bed_read (bed, entry))
		;

	bed_close (bed);
	bed_entry_free (entry);

	xunlink (bed_path);
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

START_TEST (test_bed_check_format)
{
	// After comments, 3 fields or more - by tabs or spaces
	ck_assert_int_eq (check_format_of (bed_check_format,
				"#chrom\tstart\tend\nchr1\t0\t100\n"), 1);
	ck_assert_int_eq (check_format_of (bed_check_format,
				"chr1 0 100 gene1 0 +\n"), 1);
	ck_assert_int_eq (check_format_of (bed_check_format,
				"track name=ponga\nchr1\t0\t100\n"), 1);

	// Too few fields, 'chromStart' is not a number, FASTA, nothing
	ck_assert_int_eq (check_format_of (bed_check_format,
				"chr1\t0\n"), 0);
	ck_assert_int_eq (check_format_of (bed_check_format,
				"chr1\tstart\t100\n"), 0);
	ck_assert_int_eq (check_format_of (bed_check_format,
				">chr1\nACGT\n"), 0);
	ck_assert_int_eq (check_format_of (bed_check_format, ""), 0);
}
END_TEST

Suite *
make_bed_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("BED");

	/* Core test case */
	tc_core = tcase_create ("Core");

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_test (tc_core, test_bed_check_format);
	tcase_add_test (tc_core, test_bed_header);
	tcase_add_loop_test (tc_core, test_beds, 0, 10);

	tcase_add_exit_test (tc_abort, test_open_fatal,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_chrom_fatal,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_chrom_start_fatal, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_chrom_end_fatal,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_truncated_fatal,       EXIT_FAILURE);

	suite_add_tcase (s, tc_core);
	suite_add_tcase (s, tc_abort);

	return s;
}
