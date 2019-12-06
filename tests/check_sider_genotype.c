#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/sam.h"
#include "../src/log.h"
#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/db.h"
#include "../src/genotype.h"

static void
create_file (char *path)
{
	int fd = xmkstemp (path);
	close (fd);
}

static void
populate_sorted_bam (const char *path)
{
	FILE *fp = xfopen (path, "w");

	const char bam[] =
		"@HD\tVN:1.0\tSO:coordinate\n"
		"@SQ\tSN:chr1\tLN:248956422\n"
		"@SQ\tSN:chr2\tLN:248956422\n"
		"@SQ\tSN:chr3\tLN:248956422\n"
		"@SQ\tSN:chr4\tLN:248956422\n"
		"@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188	CL:bwa mem -t 1 ponga/ponga.fa ponga.fastq\n"
		"E1\t109\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"E2\t157\tchr1\t10\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"N1\t99\tchr1\t11\t20\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"N2\t147\tchr1\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
		"N3\t1024\tchr1\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
		"F1\t99\tchr3\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"F2\t99\tchr3\t98\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"M1\t99\tchr3\t99\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"M2\t99\tchr3\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
		"M3\t1024\tchr3\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n";

	xfprintf (fp, "%s", bam);
	xfclose (fp);

	// BAM INDEX work only on BAM, not even SAM!
	fp = xfopen (path, "rb+");
	sam_to_bam_fp (fp, path);
}

static void
populate_unsorted_bam (const char *path)
{
	FILE *fp = xfopen (path, "w");

	const char bam[] =
		"@HD\tVN:1.0\tSO:queryname\n"
		"@SQ\tSN:chr1\tLN:248956422\n"
		"@SQ\tSN:chr2\tLN:248956422\n"
		"@SQ\tSN:chr3\tLN:248956422\n"
		"@SQ\tSN:chr4\tLN:248956422\n"
		"@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188	CL:bwa mem -t 1 ponga/ponga.fa ponga.fastq\n"
		"E1\t99\tchr2\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"E2\t99\tchr2\t98\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"N1\t99\tchr2\t99\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"N2\t99\tchr2\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
		"N3\t1024\tchr2\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
		"F1\t109\tchr4\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"F2\t157\tchr4\t10\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"M1\t99\tchr4\t11\t20\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"M2\t147\tchr4\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
		"M3\t1024\tchr4\t150\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n";

	xfprintf (fp, "%s", bam);

	xfclose (fp);
}

static void
index_bam (const char *path)
{
	if (sam_index_build (path, 0) < 0)
		log_fatal ("Failed to create BAM index for '%s'", path);
}

static sqlite3 *
create_db (char *db_path)
{
	create_file (db_path);
	return db_create (db_path);
}

static void
populate_db (sqlite3 *db, const char *bam1, const char *bam2)
{
	char *sql = NULL;

	// Database dump
	xasprintf (&sql,
		"BEGIN TRANSACTION;\n"
		"INSERT INTO source VALUES (1,1,'%s');\n"
		"INSERT INTO source VALUES (2,1,'%s');\n"
		"INSERT INTO alignment VALUES (1,'q1',97,'chr1',1,100,'100M',100,100,'chr1',1,1,1);\n"
		"INSERT INTO alignment VALUES (2,'q2',97,'chr2',1,100,'100M',100,100,'chr1',1,1,2);\n"
		"INSERT INTO alignment VALUES (3,'q3',97,'chr3',1,100,'100M',100,100,'chr1',1,1,1);\n"
		"INSERT INTO alignment VALUES (4,'q4',97,'chr4',1,100,'100M',100,100,'chr1',1,1,2);\n"
		"INSERT INTO clustering VALUES (1,2,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,2,2,3,100);\n"
		"INSERT INTO clustering VALUES (3,2,3,3,100);\n"
		"INSERT INTO clustering VALUES (4,2,4,3,100);\n"
		"INSERT INTO cluster_merging VALUES (1,1,2);\n"
		"INSERT INTO cluster_merging VALUES (2,2,2);\n"
		"INSERT INTO cluster_merging VALUES (3,3,2);\n"
		"INSERT INTO cluster_merging VALUES (4,4,2);\n"
		"INSERT INTO retrocopy VALUES (1,'chr1',1,200,'PONGA1',1,100,2,1,0.0);\n"
		"INSERT INTO retrocopy VALUES (2,'chr2',1,200,'PONGA2',1,100,2,1,0.0);\n"
		"INSERT INTO retrocopy VALUES (3,'chr3',1,200,'PONGA1',1,100,2,1,0.0);\n"
		"INSERT INTO retrocopy VALUES (4,'chr4',1,200,'PONGA2',1,100,2,1,0.0);\n"
		"COMMIT;", bam1, bam2);

	db_exec (db, sql);
	xfree (sql);
}

START_TEST (test_genotype)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char bam_sorted_file[] = "/tmp/ponga_sorted.bam.XXXXXX";
	char bam_unsorted_file[] = "/tmp/ponga_unsorted.bam.XXXXXX";

	sqlite3 *db = NULL;
	sqlite3_stmt *stmt = NULL;

	db = create_db (db_file);
	create_file (bam_sorted_file);
	create_file (bam_unsorted_file);

	populate_sorted_bam (bam_sorted_file);
	populate_unsorted_bam (bam_unsorted_file);
	populate_db (db, bam_sorted_file, bam_unsorted_file);

	index_bam (bam_sorted_file);
	stmt = db_prepare_genotype_stmt (db);

	genotype (stmt, 2, 0);

	db_finalize (stmt);
	db_close (db);

	xunlink (db_file);
	xunlink (bam_sorted_file);
	xunlink (bam_unsorted_file);
}
END_TEST

Suite *
make_genotype_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Genotype");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_genotype);
	suite_add_tcase (s, tc_core);

	return s;
}
