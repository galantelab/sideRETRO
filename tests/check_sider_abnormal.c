#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/utils.h"
#include "../src/db.h"
#include "../src/abnormal.h"

#define TEST_ABNORMAL_BUFSIZ 64

struct _TestAbnormal
{
	char db_path[TEST_ABNORMAL_BUFSIZ];
	char sam_path[TEST_ABNORMAL_BUFSIZ];
	char gtf_path[TEST_ABNORMAL_BUFSIZ];

	sqlite3 *db;
	sqlite3_stmt *exon_stmt;
	sqlite3_stmt *overlapping_stmt;
	sqlite3_stmt *alignment_stmt;

	ChrStd *cs;

	ExonTree *exon_tree;

	AbnormalArg *arg;
};

typedef struct _TestAbnormal TestAbnormal;

static const char *sam_sorted =
	"@HD\tVN:1.0\tSO:queryname\n"
	"@SQ\tSN:chr1\tLN:248956422\n"
	"@SQ\tSN:chr2\tLN:242193529\n"
	"@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188	CL:bwa mem -t 1 ponga/ponga.fa ponga.fastq\n"
	"E1\t109\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"E1\t157\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\t147\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
	"C2\t97\tchr1\t40\t60\t10M\tchr2\t1\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"C2\t145\tchr2\t1\t60\t10M\tchr1\t40\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"D3\t97\tchr2\t20\t60\t10M\t=\t20000\t19990\tAAAAGGGCCC\t~~~~~~~~~~\n"
	"D3\t145\tchr2\t20000\t60\t10M\t=\t20\t-19990\tCCCCCTTTAG\t~~~~~~~~~~\n"
	"S4\t99\tchr1\t100\t60\t10M\t=\t120\t30\tAAACCCGGGG\t~~~~~~~~~~\n"
	"S4\t147\tchr1\t120\t60\t5M5S\t=\t100\t-30\tGGGCCCCCCC\t~~~~~~~~~~\n"
	"S4\t2195\tchr2\t100\t60\t5H5M\tchr1\t100\t0\tCCCCC\t~~~~~\n";

static const char *sam_unsorted =
	"@SQ\tSN:chr1\tLN:248956422\n"
	"@SQ\tSN:chr2\tLN:242193529\n"
	"@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188	CL:bwa mem -t 1 ponga/ponga.fa ponga.fastq\n"
	"E1\t109\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"E1\t157\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\t147\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
	"C2\t97\tchr1\t40\t60\t10M\tchr2\t1\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"S4\t99\tchr1\t100\t60\t10M\t=\t120\t30\tAAACCCGGGG\t~~~~~~~~~~\n"
	"S4\t147\tchr1\t120\t60\t5M5S\t=\t100\t-30\tGGGCCCCCCC\t~~~~~~~~~~\n"
	"C2\t145\tchr2\t1\t60\t10M\tchr1\t40\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"D3\t97\tchr2\t20\t60\t10M\t=\t20000\t19990\tAAAAGGGCCC\t~~~~~~~~~~\n"
	"D3\t145\tchr2\t20000\t60\t10M\t=\t20\t-19990\tCCCCCTTTAG\t~~~~~~~~~~\n"
	"S4\t2195\tchr2\t100\t60\t5H5M\tchr1\t100\t0\tCCCCC\t~~~~~\n";

static const char *gtf =
	"chr1\t.\texon\t45\t65\t.\t+\t.\t"
	"gene_name \"e1\"; gene_id \"ENG1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; "
	"exon_id \"ENSE1\";\n"
	"chr2\t.\texon\t19990\t20050\t.\t-\t.\t"
	"gene_name \"e2\"; gene_id \"ENG2\"; transcript_id \"t2\"; transcript_type \"protein_coding\"; "
	"exon_id \"ENSE2\";\n";

static sqlite3 *
create_db (char *db_path)
{
	int fd;

	fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
create_file (char *path, const char *txt)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, txt, "");

	xfclose (fp);
}

static void
test_abnormal_init (TestAbnormal *a, const char *sam)
{
	strncpy (a->db_path, "/tmp/ponga.db.XXXXXX",
			TEST_ABNORMAL_BUFSIZ - 1);
	strncpy (a->sam_path, "/tmp/ponga.sam.XXXXXX",
			TEST_ABNORMAL_BUFSIZ - 1);
	strncpy (a->gtf_path, "/tmp/ponga.gtf.XXXXXX",
			TEST_ABNORMAL_BUFSIZ - 1);

	a->db_path[TEST_ABNORMAL_BUFSIZ - 1] = '\0';
	a->sam_path[TEST_ABNORMAL_BUFSIZ - 1] = '\0';
	a->gtf_path[TEST_ABNORMAL_BUFSIZ - 1] = '\0';

	a->db = create_db (a->db_path);
	a->exon_stmt = db_prepare_exon_stmt (a->db);
	a->overlapping_stmt = db_prepare_overlapping_stmt (a->db);
	a->alignment_stmt = db_prepare_alignment_stmt (a->db);

	a->cs = chr_std_new ();

	create_file (a->sam_path, sam);
	create_file (a->gtf_path, gtf);

	a->exon_tree = exon_tree_new (a->exon_stmt,
			a->overlapping_stmt, a->cs);

	exon_tree_index_dump (a->exon_tree, a->gtf_path);

	a->arg = xcalloc (1, sizeof (AbnormalArg));

	*a->arg = (AbnormalArg) {
		.tid = 1,
		.num_threads = 1,
		.sam_file = a->sam_path,
		.exon_tree = a->exon_tree,
		.cs = a->cs,
		.alignment_stmt = a->alignment_stmt,
		.max_distance = 10000,
		.either = 0,
		.exon_frac  = -1,
		.alignment_frac = -1,
	};
}

static void
test_abnormal_destroy (TestAbnormal *a)
{
	if (a == NULL)
		return;

	db_finalize (a->exon_stmt);
	db_finalize (a->overlapping_stmt);
	db_finalize (a->alignment_stmt);

	db_close (a->db);
	chr_std_free (a->cs);

	exon_tree_free (a->exon_tree);

	xfree (a->arg);

	// Remove temp files
	xunlink (a->db_path);
	xunlink (a->sam_path);
	xunlink (a->gtf_path);
}

static sqlite3_stmt *
prepare_alignment_search (sqlite3 *db)
{
	const char sql[] =
		"SELECT qname, chr, type FROM alignment ORDER BY qname ASC, chr ASC";
	return db_prepare (db, sql);
}

START_TEST (test_abnormal_filter_sorted)
{
	// Init AbnormalArg struct and create database
	// and sam files
	TestAbnormal a;
	test_abnormal_init (&a, sam_sorted);

	sqlite3_stmt *search_stmt = NULL;
	const char *qname = NULL;
	const char *chr = NULL;
	int type = 0;
	int i = 0;

	/* TRUE POSITIVE VALUES */
	int alignment_size = 7;

	const char *qnames_with_chr[][2] = {
		{"C2", "chr1"},
		{"C2", "chr2"},
		{"D3", "chr2"},
		{"D3", "chr2"},
		{"S4", "chr1"},
		{"S4", "chr1"},
		{"S4", "chr2"},
	};

	int types[] = {
		ABNORMAL_CHROMOSOME|ABNORMAL_EXONIC,
		ABNORMAL_CHROMOSOME,
		ABNORMAL_DISTANCE,
		ABNORMAL_DISTANCE|ABNORMAL_EXONIC,
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME,
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME,
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME
	};

	// RUN FOOLS
	abnormal_filter (a.arg);

	// Let's get the alignment table values
	search_stmt = prepare_alignment_search (a.db);

	/* TIME TO TEST */
	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			qname = db_column_text (search_stmt, 0);
			ck_assert_str_eq (qname, qnames_with_chr[i][0]);

			chr = db_column_text (search_stmt, 1);
			ck_assert_str_eq (chr, qnames_with_chr[i][1]);

			type = db_column_int (search_stmt, 2);
			ck_assert_int_eq (type, types[i]);
		}

	ck_assert_uint_eq (i, alignment_size);

	// Time to cleanup
	db_finalize (search_stmt);
	test_abnormal_destroy (&a);
}
END_TEST

START_TEST (test_abnormal_filter_unsorted)
{
	// Init AbnormalArg struct and create database
	// and sam files
	TestAbnormal a;
	test_abnormal_init (&a, sam_unsorted);

	sqlite3_stmt *search_stmt = NULL;
	const char *qname = NULL;
	const char *chr = NULL;
	int type = 0;
	int i = 0;

	/* TRUE POSITIVE VALUES */
	int alignment_size = 7;

	const char *qnames_with_chr[][2] = {
		{"C2", "chr1"},
		{"C2", "chr2"},
		{"D3", "chr2"},
		{"D3", "chr2"},
		{"S4", "chr1"},
		{"S4", "chr1"},
		{"S4", "chr2"},
	};

	int types[] = {
		ABNORMAL_CHROMOSOME|ABNORMAL_EXONIC,
		ABNORMAL_CHROMOSOME,
		ABNORMAL_DISTANCE,
		ABNORMAL_DISTANCE|ABNORMAL_EXONIC,
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME,
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME,
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME
	};

	// RUN FOOLS
	abnormal_filter (a.arg);

	// Let's get the alignment table values
	search_stmt = prepare_alignment_search (a.db);

	/* TIME TO TEST */
	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			qname = db_column_text (search_stmt, 0);
			ck_assert_str_eq (qname, qnames_with_chr[i][0]);

			chr = db_column_text (search_stmt, 1);
			ck_assert_str_eq (chr, qnames_with_chr[i][1]);

			type = db_column_int (search_stmt, 2);
			ck_assert_int_eq (type, types[i]);
		}

	ck_assert_uint_eq (i, alignment_size);

	// Time to cleanup
	db_finalize (search_stmt);
	test_abnormal_destroy (&a);
}
END_TEST

Suite *
make_abnormal_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Abnormal");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_abnormal_filter_sorted);
	tcase_add_test (tc_core, test_abnormal_filter_unsorted);
	suite_add_tcase (s, tc_core);

	return s;
}
