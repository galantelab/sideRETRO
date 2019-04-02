#include "config.h"

#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/log.h"
#include "../src/wrapper.h"
#include "../src/sam.h"

#define HAS_SAMTOOLS (!system ("which samtools 2>&1 > /dev/null"))

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
create_sam (char *filename)
{
	FILE *fp = xfopen (filename, "w");
	fprintf (fp, sam_str);
	xfclose (fp);
}

START_TEST (test_sam_to_bam)
{
	log_set_quiet (1);
	char sam_file1[] = "/tmp/ponga_sam-XXXXXX";
	char sam_file2[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[]  = "/tmp/ponga_bam-XXXXXX";

	mktemp (sam_file1);
	mktemp (sam_file2);
	mktemp (bam_file);

	create_sam (sam_file1);
	ck_assert_int_eq (sam_to_bam (sam_file1, bam_file), 1);

	if (HAS_SAMTOOLS)
		{
			char *cmd1 = NULL;
			char *cmd2 = NULL;

			xasprintf (&cmd1, "samtools view -h %s -o %s 2> /dev/null",
					bam_file, sam_file2);
			xasprintf (&cmd2, "diff -q %s %s 2>&1 > /dev/null",
					sam_file1, sam_file2);

			system (cmd1);
			ck_assert_int_eq (system (cmd2), 0);

			xfree (cmd1);
			xfree (cmd2);
		}

	if (unlink (sam_file1) == -1)
		log_errno_error ("Failed to remove '%s'", sam_file1);

	if (unlink (sam_file2) == -1)
		log_errno_error ("Failed to remove '%s'", sam_file2);

	if (unlink (bam_file) == -1)
		log_errno_error ("Failed to remove '%s'", bam_file);
}
END_TEST

START_TEST (test_sam_to_bam_fp)
{
	log_set_quiet (1);
	char sam_file[] = "/tmp/ponga_sam-XXXXXX";
	char bam_file[] = "/tmp/ponga_bam-XXXXXX";

	mktemp (sam_file);
	mktemp (bam_file);

	create_sam (sam_file);

	FILE *fp = xfopen (sam_file, "r");
	ck_assert_int_eq (sam_to_bam_fp (fp, bam_file), 1);

	xfclose (fp);
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
