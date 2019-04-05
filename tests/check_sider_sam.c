#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/log.h"
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
	log_set_quiet (1);
	char sam_file1[] = "/tmp/ponga_sam-XXXXXX";
	char sam_file2[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[]  = "/tmp/ponga_bam-XXXXXX";

	int sam_fd1 = xmkstemp (sam_file1);
	int sam_fd2 = xmkstemp (sam_file2);
	int bam_fd  = xmkstemp (bam_file);

	create_sam_fd (sam_fd1);

	ck_assert_int_eq (sam_to_bam (sam_file1, bam_file), 1);

	if (which ("samtools"))
		{
			char *cmd1 = NULL;
			char *cmd2 = NULL;
			int len = 0;

			len = xasprintf (&cmd1, "samtools view -h %s -o %s 2> /dev/null",
					bam_file, sam_file2);
			len = xasprintf (&cmd2, "diff -q %s %s 2>&1 > /dev/null",
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
	log_set_quiet (1);
	char sam_file[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[] = "/tmp/ponga_bam-XXXXXX";

	int sam_fd = xmkstemp (sam_file);
	int bam_fd = xmkstemp (bam_file);

	create_sam_fd (sam_fd);

	FILE *fp = xfopen (sam_file, "r");
	ck_assert_int_eq (sam_to_bam_fp (fp, bam_file), 1);

	xfclose (fp);

	xunlink (sam_file);
	xunlink (bam_file);
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
	suite_add_tcase (s, tc_core);

	return s;
}
