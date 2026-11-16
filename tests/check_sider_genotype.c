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
#include <math.h>
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
	sqlite3_stmt *intron_stmt = NULL;

	db = create_db (db_file);
	create_file (bam_sorted_file);
	create_file (bam_unsorted_file);

	populate_sorted_bam (bam_sorted_file);
	populate_unsorted_bam (bam_unsorted_file);
	populate_db (db, bam_sorted_file, bam_unsorted_file);

	index_bam (bam_sorted_file);
	stmt = db_prepare_genotype_stmt (db);
	intron_stmt = db_prepare_intron_retention_stmt (db);

	genotype (stmt, intron_stmt, 2, 0);

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
	db_finalize (intron_stmt);
	db_close (db);

	xunlink (db_file);
	xunlink (bam_sorted_file);
	xunlink (bam_unsorted_file);
}
END_TEST

START_TEST (test_genotype_deep)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char sam_file[] = "/tmp/ponga.sam.XXXXXX";
	const int depth = 1200;
	const int depth_mapq0 = 5;

	sqlite3 *db = NULL;
	sqlite3_stmt *stmt = NULL;
	sqlite3_stmt *intron_stmt = NULL;
	FILE *fp = NULL;
	int i = 0;

	db = create_db (db_file);
	create_file (sam_file);

	/*
	* A deep site: (1/2)^depth underflows to zero, and so
	* would log10 of it. The reads with MAPQ 0 (no phred
	* quality filter) must not give log10 (0) either
	*/
	fp = xfopen (sam_file, "w");
	xfprintf (fp,
			"@HD\tVN:1.0\tSO:queryname\n"
			"@SQ\tSN:chr1\tLN:248956422\n");

	for (i = 0; i < depth + depth_mapq0; i++)
		xfprintf (fp, "N%d\t99\tchr1\t95\t%d\t10M\t=\t150\t65\tATCGATCGAT\t~~~~~~~~~~\n",
				i, i < depth ? 60 : 0);

	xfclose (fp);

	populate_db (db, sam_file, sam_file);
	stmt = db_prepare_genotype_stmt (db);
	intron_stmt = db_prepare_intron_retention_stmt (db);

	genotype (stmt, intron_stmt, 1, 0);

	sqlite3_stmt *search_stmt = db_prepare (db,
			"SELECT reference_depth, alternate_depth, ho_ref_likelihood,\n"
			"	he_likelihood, ho_alt_likelihood\n"
			"FROM genotype\n"
			"WHERE retrocopy_id = 1 AND source_id = 1");

	ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
	ck_assert_int_eq (db_column_int (search_stmt, 0), depth + depth_mapq0);
	ck_assert_int_eq (db_column_int (search_stmt, 1), 1);

	for (i = 2; i <= 4; i++)
		ck_assert (isfinite (db_column_double (search_stmt, i)));

	// Homozygous reference is the most likely
	ck_assert (db_column_double (search_stmt, 2) > db_column_double (search_stmt, 3));
	ck_assert (db_column_double (search_stmt, 3) > db_column_double (search_stmt, 4));
	db_finalize (search_stmt);

	db_finalize (stmt);
	db_finalize (intron_stmt);
	db_close (db);

	xunlink (db_file);
	xunlink (sam_file);
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

	genotype (stmt, db_prepare_intron_retention_stmt (db), 1, 0);
}

START_TEST (test_genotype_intron)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char sam_file[] = "/tmp/ponga.sam.XXXXXX";
	char bam_file[] = "/tmp/ponga.bam.XXXXXX";

	sqlite3 *db = NULL;
	sqlite3_stmt *stmt = NULL;
	sqlite3_stmt *intron_stmt = NULL;
	sqlite3_stmt *search_stmt = NULL;
	FILE *fp = NULL;
	char *sql = NULL;
	int source_id = 0;

	/*
	* Retrocopy 1 at chr1:1-200, from GENE at chr2. The union
	* of its exons is 1000-1200, 2000-2100 and 3000-3100, so
	* its introns are 1201-1999 (exons 2 and 3) and 2101-2999
	* (exons 3 and 4). Retrocopy 2 at chr1:300-400 comes from
	* a single exon gene: No intron
	*
	* A1  mate inside the first intron
	* A2  the same, with the mate quality and cigar tags
	* A3  mate reaching the exon 2: Not intronic
	* A4  mate reaching the exon 3 by the length of the read,
	*     but not by its own cigar: Intronic
	* A5  mate with low quality (MQ)
	* A6  low quality alignment
	* A7  supplementary alignment
	* A8  mate inside the second intron
	* A9  mate at another contig
	* A10 duplication
	* A11 mate cigar tag that is not a string: Fall back to
	*     the length of the read
	* B1  retrocopy 2, mate at the intron of the other gene
	*/
	const char sam[] =
		"@HD\tVN:1.0\tSO:coordinate\n"
		"@SQ\tSN:chr1\tLN:248956422\n"
		"@SQ\tSN:chr2\tLN:248956422\n"
		"@SQ\tSN:chr3\tLN:248956422\n"
		"A1\t97\tchr1\t50\t60\t10M\tchr2\t1300\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A2\t97\tchr1\t51\t60\t10M\tchr2\t1500\t0\tATCGATCGAT\t~~~~~~~~~~\tMC:Z:10M\tMQ:i:60\n"
		"A3\t97\tchr1\t52\t60\t10M\tchr2\t1195\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A4\t97\tchr1\t53\t60\t10M\tchr2\t1992\t0\tATCGATCGAT\t~~~~~~~~~~\tMC:Z:2=1X1D1N1M2I4S\n"
		"A5\t97\tchr1\t54\t60\t10M\tchr2\t2500\t0\tATCGATCGAT\t~~~~~~~~~~\tMQ:i:0\n"
		"A6\t97\tchr1\t55\t5\t10M\tchr2\t2500\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A7\t2145\tchr1\t56\t60\t10M\tchr2\t2500\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A8\t97\tchr1\t57\t60\t10M\tchr2\t2500\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A9\t97\tchr1\t58\t60\t10M\tchr3\t1300\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A10\t1121\tchr1\t59\t60\t10M\tchr2\t1300\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"A11\t97\tchr1\t60\t60\t10M\tchr2\t1400\t0\tATCGATCGAT\t~~~~~~~~~~\tMC:i:5\n"
		"B1\t97\tchr1\t350\t60\t10M\tchr2\t1300\t0\tATCGATCGAT\t~~~~~~~~~~\n";

	db = create_db (db_file);
	create_file (sam_file);
	create_file (bam_file);

	// Source 1: SAM, no index (linear search)
	fp = xfopen (sam_file, "w");
	xfprintf (fp, "%s", sam);
	xfclose (fp);

	// Source 2: Indexed BAM
	fp = xfopen (bam_file, "w");
	xfprintf (fp, "%s", sam);
	xfclose (fp);

	fp = xfopen (bam_file, "rb+");
	sam_to_bam_fp (fp, bam_file);
	index_bam (bam_file);

	xasprintf (&sql,
		"BEGIN TRANSACTION;\n"
		"INSERT INTO source VALUES (1,1,'%s');\n"
		"INSERT INTO source VALUES (2,1,'%s');\n"
		"INSERT INTO exon VALUES (1,'GENE','chr2',1000,1100,'+','eg','ee1');\n"
		"INSERT INTO exon VALUES (2,'GENE','chr2',1050,1200,'+','eg','ee2');\n"
		"INSERT INTO exon VALUES (3,'GENE','chr2',2000,2100,'+','eg','ee3');\n"
		"INSERT INTO exon VALUES (4,'GENE','chr2',3000,3100,'+','eg','ee4');\n"
		"INSERT INTO exon VALUES (5,'SINGLE','chr2',5000,5100,'+','es','es1');\n"
		"INSERT INTO alignment VALUES (1,'q1',97,'chr1',1,60,'100M',100,100,'chr2',1,1,1,'*','*');\n"
		"INSERT INTO alignment VALUES (2,'q2',97,'chr1',300,60,'100M',100,100,'chr2',1,1,1,'*','*');\n"
		"INSERT INTO cluster VALUES (1,1,'chr1',1,200,'GENE',31);\n"
		"INSERT INTO cluster VALUES (2,1,'chr1',300,400,'SINGLE',31);\n"
		"INSERT INTO clustering VALUES (1,1,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,1,2,3,100);\n"
		"INSERT INTO cluster_merging VALUES (1,1,1);\n"
		"INSERT INTO cluster_merging VALUES (2,2,1);\n"
		"INSERT INTO retrocopy VALUES (1,'chr1',1,200,'GENE',1,100,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (2,'chr1',300,400,'SINGLE',1,350,2,1,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"COMMIT;", sam_file, bam_file);

	db_exec (db, sql);
	xfree (sql);

	stmt = db_prepare_genotype_stmt (db);
	intron_stmt = db_prepare_intron_retention_stmt (db);

	genotype (stmt, intron_stmt, 2, 20);

	// The same counts with and without the index:
	// A1, A2, A4 and A11 at the first intron, A8 at the second
	search_stmt = db_prepare (db,
			"SELECT source_id, retrocopy_id, gene_name, chr, start, end,\n"
			"	exon_id_left, exon_id_right, reads\n"
			"FROM intron_retention\n"
			"ORDER BY source_id, start");

	for (source_id = 1; source_id <= 2; source_id++)
		{
			ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
			ck_assert_int_eq (db_column_int (search_stmt, 0), source_id);
			ck_assert_int_eq (db_column_int (search_stmt, 1), 1);
			ck_assert_str_eq (db_column_text (search_stmt, 2), "GENE");
			ck_assert_str_eq (db_column_text (search_stmt, 3), "chr2");
			ck_assert_int_eq (db_column_int (search_stmt, 4), 1201);
			ck_assert_int_eq (db_column_int (search_stmt, 5), 1999);
			ck_assert_int_eq (db_column_int (search_stmt, 6), 2);
			ck_assert_int_eq (db_column_int (search_stmt, 7), 3);
			ck_assert_int_eq (db_column_int (search_stmt, 8), 4);

			ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
			ck_assert_int_eq (db_column_int (search_stmt, 0), source_id);
			ck_assert_int_eq (db_column_int (search_stmt, 4), 2101);
			ck_assert_int_eq (db_column_int (search_stmt, 5), 2999);
			ck_assert_int_eq (db_column_int (search_stmt, 6), 3);
			ck_assert_int_eq (db_column_int (search_stmt, 7), 4);
			ck_assert_int_eq (db_column_int (search_stmt, 8), 1);
		}

	ck_assert_int_eq (db_step (search_stmt), SQLITE_DONE);
	db_finalize (search_stmt);

	// Another run cleans the former counts
	genotype (stmt, intron_stmt, 1, 20);

	search_stmt = db_prepare (db, "SELECT COUNT(*) FROM intron_retention");
	ck_assert_int_eq (db_step (search_stmt), SQLITE_ROW);
	ck_assert_int_eq (db_column_int (search_stmt, 0), 4);
	db_finalize (search_stmt);

	db_finalize (stmt);
	db_finalize (intron_stmt);
	db_close (db);

	xunlink (db_file);
	xunlink (sam_file);
	xunlink (bam_file);
}
END_TEST

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
	tcase_add_test (tc_core, test_genotype_deep);
	tcase_add_test (tc_core, test_genotype_intron);
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
