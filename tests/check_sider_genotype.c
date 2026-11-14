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
		"INSERT INTO alignment VALUES (1,'q1',97,'chr1',1,100,'100M',100,100,'chr1',1,1,1,'*','*');\n"
		"INSERT INTO alignment VALUES (2,'q2',97,'chr2',1,100,'100M',100,100,'chr1',1,1,2,'*','*');\n"
		"INSERT INTO alignment VALUES (3,'q3',97,'chr3',1,100,'100M',100,100,'chr1',1,1,1,'*','*');\n"
		"INSERT INTO alignment VALUES (4,'q4',97,'chr4',1,100,'100M',100,100,'chr1',1,1,2,'*','*');\n"
		"INSERT INTO clustering VALUES (1,2,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,2,2,3,100);\n"
		"INSERT INTO clustering VALUES (3,2,3,3,100);\n"
		"INSERT INTO clustering VALUES (4,2,4,3,100);\n"
		"INSERT INTO cluster_merging VALUES (1,1,2);\n"
		"INSERT INTO cluster_merging VALUES (2,2,2);\n"
		"INSERT INTO cluster_merging VALUES (3,3,2);\n"
		"INSERT INTO cluster_merging VALUES (4,4,2);\n"
		"INSERT INTO alignment VALUES (5,'q5',97,'chr9',1,100,'100M',100,100,'chr1',1,1,1,'*','*');\n"
		"INSERT INTO clustering VALUES (5,2,5,3,100);\n"
		"INSERT INTO cluster_merging VALUES (5,5,2);\n"
		"INSERT INTO retrocopy VALUES (5,'chr9',1,200,'PONGA1',1,100,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (1,'chr1',1,200,'PONGA1',1,100,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (2,'chr2',1,200,'PONGA2',1,100,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (3,'chr3',1,200,'PONGA1',1,100,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (4,'chr4',1,200,'PONGA2',1,100,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
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

	// chr9 is not at the alignment headers, sorted (indexed)
	// nor unsorted (linear search): Its retrocopy still gets a
	// genotype at both sources, with no reference reads
	sqlite3_stmt *search_stmt = db_prepare (db,
			"SELECT COUNT(*), SUM(retrocopy_id = 5),\n"
			"	SUM(retrocopy_id = 5 AND reference_depth = 0)\n"
			"FROM genotype");

	ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
	ck_assert_int_eq (db_column_int (search_stmt, 0), 10);
	ck_assert_int_eq (db_column_int (search_stmt, 1), 2);
	ck_assert_int_eq (db_column_int (search_stmt, 2), 2);
	db_finalize (search_stmt);

	db_finalize (stmt);
	db_close (db);

	xunlink (db_file);
	xunlink (bam_sorted_file);
	xunlink (bam_unsorted_file);
}
END_TEST

static void
run_genotype_on (const char *path)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = create_db (db_file);
	sqlite3_stmt *stmt = NULL;

	populate_db (db, path, path);
	stmt = db_prepare_genotype_stmt (db);

	genotype (stmt, 1, 0);
}

static void
corrupt_second_block (const char *path)
{
	FILE *fp = xfopen (path, "rb+");
	unsigned char head[18];
	long second = 0;

	// BGZF: The block size minus 1 is at bytes 16-17. The
	// header has a block of its own, so the records begin
	// at the second block
	ck_assert_int_eq (fread (head, 1, 18, fp), 18);
	second = (head[16] | (head[17] << 8)) + 1;

	fseek (fp, second + 24, SEEK_SET);
	fputc (0xff, fp);
	fseek (fp, second + 25, SEEK_SET);
	fputc (0x00, fp);
	xfclose (fp);
}

#ifdef HAVE_LD_WRAP
static void
run_genotype_fault (Fault f)
{
	char bam[] = "/tmp/ponga.bam.XXXXXX";

	create_file (bam);
	populate_sorted_bam (bam);
	index_bam (bam);

	fault_inject (f);
	run_genotype_on (bam);
}

START_TEST (test_genotype_bam_init1_fail)
{
	run_genotype_fault (FAULT_BAM_INIT1);
}
END_TEST

START_TEST (test_genotype_query_fail)
{
	run_genotype_fault (FAULT_SAM_ITR_QUERYI);
}
END_TEST

START_TEST (test_genotype_close_fail)
{
	run_genotype_fault (FAULT_HTS_CLOSE);
}
END_TEST
#endif

START_TEST (test_genotype_open_fail)
{
	run_genotype_on ("/it-must-not-exist.bam");
}
END_TEST

START_TEST (test_genotype_header_fail)
{
	char bam[] = "/tmp/ponga.bam.XXXXXX";

	// A BAM truncated after the BGZF header of the first block
	create_file (bam);
	populate_sorted_bam (bam);
	ck_assert_int_eq (truncate (bam, 20), 0);

	run_genotype_on (bam);
}
END_TEST

START_TEST (test_genotype_linear_read_fail)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";
	FILE *fp = NULL;

	// No index: Linear search. The second record is broken
	create_file (sam);
	fp = xfopen (sam, "w");
	xfprintf (fp, "%s",
		"@SQ\tSN:chr1\tLN:248956422\n"
		"R1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
		"R2\tponga\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n");
	xfclose (fp);

	run_genotype_on (sam);
}
END_TEST

START_TEST (test_genotype_indexed_read_fail)
{
	char bam[] = "/tmp/ponga.bam.XXXXXX";

	// Indexed, then its records are corrupted
	create_file (bam);
	populate_sorted_bam (bam);
	index_bam (bam);
	corrupt_second_block (bam);

	run_genotype_on (bam);
}
END_TEST

Suite *
make_genotype_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("Genotype");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_genotype);
	suite_add_tcase (s, tc_core);

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_exit_test (tc_abort, test_genotype_open_fail,         EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_genotype_header_fail,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_genotype_linear_read_fail,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_genotype_indexed_read_fail, EXIT_FAILURE);
#ifdef HAVE_LD_WRAP
	tcase_add_exit_test (tc_abort, test_genotype_bam_init1_fail,    EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_genotype_query_fail,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_genotype_close_fail,        EXIT_FAILURE);
#endif
	suite_add_tcase (s, tc_abort);

	return s;
}
