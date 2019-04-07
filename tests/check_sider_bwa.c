#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/log.h"
#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/bwa.h"

#define BWA "../bwa/bwa"

static const char *fasta =
	">PONGA\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";

static const char *fastq =
	"@ponga1\n"
	"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	"+\n"
	"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
	"@ponga2\n"
	"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
	"+\n"
	"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
	"@ponga3\n"
	"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
	"+\n"
	"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

static const char *sam[] =
{
	"ponga1\t0\tPONGA\t69\t0\t52M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t"
	"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:0\tMD:Z:52\tAS:i:52\tXS:i:52\n",
	"ponga2\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t"
	"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tAS:i:0\tXS:i:0\n",
	"ponga3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t"
	"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tAS:i:0\tXS:i:0\n"
};

static void
write_file (char *file, const char *content)
{
	int fd;
	FILE *fp = NULL;

	fd = xmkstemp (file);
	fp = xfdopen (fd, "w");

	fprintf (fp, content, "");
	xfclose (fp);
}

static void
make_index (char *file)
{
	int ret;
	int len;
	char *cmd = NULL;

	len = xasprintf (&cmd, "%s index %s 2> /dev/null", BWA, file);
	ret = system (cmd);

	if (ret)
		log_errno_fatal ("bwa index failed");

	xfree (cmd);
}

static int
build_base (char *fasta_file, char *fastq_file)
{
	if (!exists (BWA))
		return 0;

	write_file (fasta_file, fasta);
	write_file (fastq_file, fastq);
	make_index (fasta_file);

	return 1;
}

static void
remove_idx (const char *base)
{
	char *suffix[] = {"amb", "ann", "bwt", "pac", "sa"};
	int i = 0;

	for (; i < sizeof (suffix) / sizeof (char *); i++)
		{
			char *file = NULL;
			int len = 0;

			len = xasprintf (&file, "%s.%s",
					base, suffix[i]);

			xunlink (file);
			xfree (file);
		}
}

START_TEST (test_bwa_mem)
{
	log_set_quiet (1);
	bwa_mem_log_set_level (0);

	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	char fastq_file[] = "/tmp/ponga.fastq.XXXXXX";
	char sam_file[]   = "/tmp/ponga.sam.XXXXXX";
	FILE *fp = NULL;
	char buf[BUFSIZ];
	int i = 0;
	int n = 0;

	if (!build_base (fasta_file, fastq_file))
		return;

	int fd = xmkstemp (sam_file);
	close (fd);

	bwa_mem (fasta_file, fastq_file, NULL, sam_file, 1);
	fp = xfopen (sam_file, "r");

	// skip header
	while (fgets (buf, BUFSIZ, fp) != NULL)
		{
			if (buf[0] != '@')
				break;
		}

	n = sizeof (sam) / sizeof (char *);

	for (; i < n; i++)
		{
			char *s;
			ck_assert_str_eq (buf, sam[i]);
			s = fgets (buf, BUFSIZ, fp);
		}

	xfclose (fp);
	xunlink (fasta_file);
	xunlink (fastq_file);
	xunlink (sam_file);
	remove_idx (fasta_file);
}
END_TEST

Suite *
make_bwa_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("BWA");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_bwa_mem);
	suite_add_tcase (s, tc_core);

	return s;
}
