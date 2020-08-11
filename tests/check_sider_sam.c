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
#include <check.h>
#include "check_sider.h"

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
	char *queryname = "@HD\tVN:1.0\tSO:queryname\n";
	char *coordinate = "@HD\tVN:1.0\tSO:coordinate\n";

	bam_hdr_t hdr;

	hdr.text = queryname;
	hdr.l_text = strlen (queryname);

	ck_assert_int_eq (sam_test_sorted_order (&hdr,
				"queryname"), 1);
	ck_assert_int_eq (sam_test_sorted_order (&hdr,
				"coordinate"), 0);

	hdr.text = coordinate;
	hdr.l_text = strlen (coordinate);

	ck_assert_int_eq (sam_test_sorted_order (&hdr,
				"coordinate"), 1);
	ck_assert_int_eq (sam_test_sorted_order (&hdr,
				"queryname"), 0);
}
END_TEST

Suite *
make_sam_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("SAM");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_sam_to_bam);
	tcase_add_test (tc_core, test_sam_to_bam_fp);
	tcase_add_test (tc_core, test_sam_test_sorted_order);
	suite_add_tcase (s, tc_core);

	return s;
}
