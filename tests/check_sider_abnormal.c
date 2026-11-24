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

#include "../src/wrapper.h"
#include "../src/utils.h"
#include "../src/db.h"
#include "../src/abnormal.h"
#include "../src/sam.h"

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
	"S4\t2195\tchr2\t100\t60\t5H5M\tchr1\t100\t0\tCCCCC\t~~~~~\n"
	"F5\t1123\tchr1\t300\t60\t10M\t=\t320\t30\tAAACCCGGGG\t~~~~~~~~~~\n"
	"F5\t147\tchr1\t320\t60\t5M5S\t=\t300\t-30\tGGGCCCCCCC\t~~~~~~~~~~\n"
	"F5\t2195\tchr2\t300\t60\t5H5M\tchr1\t300\t0\tCCCCC\t~~~~~\n"
	"G6\t97\tchr2\t500\t60\t10M\t=\t30000\t29510\tAAAAGGGCCC\t~~~~~~~~~~\n"
	"G6\t1169\tchr2\t30000\t60\t10M\t=\t500\t-29510\tCCCCCTTTAG\t~~~~~~~~~~\n"
	"X7\t97\tchr1\t200\t60\t10M\tchr2\t200\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"X7\t145\tchr2\t200\t60\t10M\tchr1\t200\t0\tTTTTTGGGGA\t~~~~~~~~~~\n";

static const char *sam_unsorted =
	"@SQ\tSN:chr1\tLN:248956422\n"
	"@SQ\tSN:chr2\tLN:242193529\n"
	"@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188	CL:bwa mem -t 1 ponga/ponga.fa ponga.fastq\n"
	"F5\t1123\tchr1\t300\t60\t10M\t=\t320\t30\tAAACCCGGGG\t~~~~~~~~~~\n"
	"X7\t97\tchr1\t200\t60\t10M\tchr2\t200\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"G6\t97\tchr2\t500\t60\t10M\t=\t30000\t29510\tAAAAGGGCCC\t~~~~~~~~~~\n"
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
	"S4\t2195\tchr2\t100\t60\t5H5M\tchr1\t100\t0\tCCCCC\t~~~~~\n"
	"F5\t147\tchr1\t320\t60\t5M5S\t=\t300\t-30\tGGGCCCCCCC\t~~~~~~~~~~\n"
	"F5\t2195\tchr2\t300\t60\t5H5M\tchr1\t300\t0\tCCCCC\t~~~~~\n"
	"G6\t1169\tchr2\t30000\t60\t10M\t=\t500\t-29510\tCCCCCTTTAG\t~~~~~~~~~~\n"
	"X7\t145\tchr2\t200\t60\t10M\tchr1\t200\t0\tTTTTTGGGGA\t~~~~~~~~~~\n";

/*
 * Exons: e1 (C2), e2 (D3, D7 and D8) and e3 - S4 is exonic by
 * its first read, before its first abnormal alignment (the
 * supplementary). X7 is abnormal, but it does not overlap any
 * exon: It is not stored
 */
static const char *gtf =
	"chr1\t.\texon\t45\t65\t.\t+\t.\t"
	"gene_name \"e1\"; gene_id \"ENG1\"; transcript_id \"t1\"; transcript_type \"protein_coding\"; "
	"exon_id \"ENSE1\";\n"
	"chr2\t.\texon\t19990\t20050\t.\t-\t.\t"
	"gene_name \"e2\"; gene_id \"ENG2\"; transcript_id \"t2\"; transcript_type \"protein_coding\"; "
	"exon_id \"ENSE2\";\n"
	"chr1\t.\texon\t95\t110\t.\t+\t.\t"
	"gene_name \"e3\"; gene_id \"ENG3\"; transcript_id \"t3\"; transcript_type \"protein_coding\"; "
	"exon_id \"ENSE3\";\n";

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

	exon_tree_index_dump (a->exon_tree, a->gtf_path, NULL, NULL, NULL);

	a->arg = xcalloc (1, sizeof (AbnormalArg));

	*a->arg = (AbnormalArg) {
		.tid = 1,
		.inc_step = 1,
		.sam_file = a->sam_path,
		.exon_tree = a->exon_tree,
		.cs = a->cs,
		.alignment_stmt = a->alignment_stmt,
		.max_distance = 10000,
		.either = 0,
		.exon_frac  = -1,
		.alignment_frac = -1,
		.phred_quality = 0,
		.max_base_freq = 1
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
		"SELECT qname, chr, type FROM alignment ORDER BY qname ASC, chr ASC, pos ASC";
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
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME|ABNORMAL_EXONIC,
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
		ABNORMAL_SUPPLEMENTARY|ABNORMAL_CHROMOSOME|ABNORMAL_EXONIC,
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

START_TEST (test_abnormal_filter_seq_qual)
{
	TestAbnormal a;
	sqlite3_stmt *search_stmt = NULL;
	const char *sams[] = {sam_sorted, sam_unsorted};
	int i = 0;
	int j = 0;

	const char sql[] =
		"SELECT seq, qual FROM alignment\n"
		"WHERE qname = 'S4'\n"
		"ORDER BY chr ASC, pos ASC";

	/*
	* The sequence is at the forward strand as it is
	* into the SAM. Only soft clipped alignments keep
	* it: Not '10M', neither '5H5M' (hard clipped bases
	* are not present)
	*/
	const char *seqs[] = {"*", "GGGCCCCCCC", "*"};
	// The base qualities are not stored
	const char *quals[] = {"*", "*", "*"};

	for (i = 0; i < 2; i++)
		{
			test_abnormal_init (&a, sams[i]);
			abnormal_filter (a.arg);

			search_stmt = db_prepare (a.db, sql);

			for (j = 0; db_step (search_stmt) == SQLITE_ROW; j++)
				{
					ck_assert_str_eq (db_column_text (search_stmt, 0), seqs[j]);
					ck_assert_str_eq (db_column_text (search_stmt, 1), quals[j]);
				}

			ck_assert_int_eq (j, 3);

			db_finalize (search_stmt);
			test_abnormal_destroy (&a);
		}
}
END_TEST

static const char *sam_edge =
	"@HD\tVN:1.0\tSO:queryname\n"
	"@SQ\tSN:chr1\tLN:248956422\n"
	"@SQ\tSN:chr2\tLN:242193529\n"
	"D7\t97\tchr2\t20\t60\t290M10S\t=\t20000\t20280\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t*\n"
	"D7\t145\tchr2\t20000\t60\t10M\t=\t20\t-20280\tCCCCCTTTAG\t~~~~~~~~~~\n"
	"D8\t97\tchr2\t30\t60\t5M5S\t=\t20010\t19990\t*\t*\n"
	"D8\t145\tchr2\t20010\t60\t10M\t=\t30\t-19990\tCCCCCTTTAG\t~~~~~~~~~~\n"
	"L9\t97\tchr2\t40\t60\t10M\t=\t40000\t39970\tAAAAAAAAAA\t~~~~~~~~~~\n"
	"L9\t145\tchr2\t40000\t60\t10M\t=\t40\t-39970\tCCCCCTTTAG\t~~~~~~~~~~\n";

static const char *sam_normal =
	"@HD\tVN:1.0\tSO:queryname\n"
	"@SQ\tSN:chr1\tLN:248956422\n"
	"N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\t147\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n";

START_TEST (test_abnormal_filter_edge)
{
	TestAbnormal a;
	AbnormalStat stat = {};
	sqlite3_stmt *stmt = NULL;
	int i = 0;

	/*
	* D7: A 300 bases soft clipped read without quality ('*')
	* D8: A soft clipped read without sequence ('*')
	* L9: A low complexity read (poly-A): Discarded
	*/
	test_abnormal_init (&a, sam_edge);
	a.arg->max_base_freq = 0.75;
	a.arg->hts_threads = 1;
	a.arg->stat = &stat;

	abnormal_filter (a.arg);

	ck_assert_int_eq (stat.alignments, 6);
	ck_assert_int_eq (stat.abnormal, 4);
	ck_assert_int_eq (stat.stored, 4);
	ck_assert_int_eq (stat.sorted, 1);

	stmt = db_prepare (a.db,
			"SELECT qname, flag, length (seq), qual FROM alignment\n"
			"ORDER BY qname ASC, flag ASC");

	/*
	* The base qualities are never stored ('*'), and the
	* sequence only for soft clipped alignments
	*/
	const struct {const char *qname; int flag; int len; const char *qual;} rows[] = {
		{"D7", 97,  300, "*"},
		{"D7", 145, 1,   "*"},
		{"D8", 97,  1,   "*"},
		{"D8", 145, 1,   "*"}
	};

	for (i = 0; db_step (stmt) == SQLITE_ROW; i++)
		{
			ck_assert_int_lt (i, 4);
			ck_assert_str_eq (db_column_text (stmt, 0), rows[i].qname);
			ck_assert_int_eq (db_column_int (stmt, 1), rows[i].flag);
			ck_assert_int_eq (db_column_int (stmt, 2), rows[i].len);
			ck_assert_str_eq (db_column_text (stmt, 3), rows[i].qual);
		}

	ck_assert_int_eq (i, 4);

	db_finalize (stmt);
	test_abnormal_destroy (&a);

	// Only a normal pair: No abnormal alignment
	stat = (AbnormalStat) {};
	test_abnormal_init (&a, sam_normal);
	a.arg->stat = &stat;

	abnormal_filter (a.arg);

	ck_assert_int_eq (stat.alignments, 2);
	ck_assert_int_eq (stat.abnormal, 0);

	test_abnormal_destroy (&a);
}
END_TEST

/*
 * Mapping quality (-Q 20), with the exon e1 at chr1:45-65:
 * S1: A read clipped at chr2 (MAPQ 60), its clip aligned as a
 * supplementary piece with MAPQ 0, and the mate at e1: The piece
 * is left out, and the fragment is kept as a split read
 * P2: The read at chr2 with MAPQ 0: Discarded
 * E3: The mate at e1 with MAPQ 0: Discarded
 * S4: The supplementary piece with MAPQ 0 is the only alignment
 * at an exon: Left out, the fragment is not exonic
 * S5: As S1, but the piece has MAPQ 60: It is a run of A bases
 * (--max-base-freq 0.75), so it is left out as well
 * L6: The read at chr2 is of low complexity: Discarded
 * In the unsorted file, the supplementary pieces come first
 */
static const char *sam_mapq_sorted =
	"@HD\tVN:1.0\tSO:queryname\n"
	"@SQ\tSN:chr1\tLN:248956422\n"
	"@SQ\tSN:chr2\tLN:242193529\n"
	"E3\t97\tchr2\t700\t60\t10M\tchr1\t50\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"E3\t145\tchr1\t50\t0\t10M\tchr2\t700\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"P2\t97\tchr2\t600\t0\t10M\tchr1\t50\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"P2\t145\tchr1\t50\t60\t10M\tchr2\t600\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"L6\t97\tchr2\t400\t60\t10M\tchr1\t50\t0\tAAAAAAAAAT\t~~~~~~~~~~\n"
	"L6\t145\tchr1\t50\t60\t10M\tchr2\t400\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"S1\t97\tchr2\t500\t60\t5M5S\tchr1\t50\t0\tCGATTAAAAA\t~~~~~~~~~~\n"
	"S1\t145\tchr1\t50\t60\t10M\tchr2\t500\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"S1\t2145\tchr2\t9000\t0\t5H5M\tchr1\t50\t0\tAAAAA\t~~~~~\n"
	"S4\t97\tchr2\t800\t60\t5M5S\tchr2\t900\t0\tCGATTAAAAA\t~~~~~~~~~~\n"
	"S4\t145\tchr2\t900\t60\t10M\tchr2\t800\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"S4\t2145\tchr1\t50\t0\t5H5M\tchr2\t900\t0\tAAAAA\t~~~~~\n"
	"S5\t97\tchr2\t1000\t60\t5M5S\tchr1\t50\t0\tCGATTAAAAA\t~~~~~~~~~~\n"
	"S5\t145\tchr1\t50\t60\t10M\tchr2\t1000\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"S5\t2145\tchr2\t9000\t60\t5H5M\tchr1\t50\t0\tAAAAA\t~~~~~\n";

static const char *sam_mapq_unsorted =
	"@SQ\tSN:chr1\tLN:248956422\n"
	"@SQ\tSN:chr2\tLN:242193529\n"
	"S1\t2145\tchr2\t9000\t0\t5H5M\tchr1\t50\t0\tAAAAA\t~~~~~\n"
	"S4\t2145\tchr1\t50\t0\t5H5M\tchr2\t900\t0\tAAAAA\t~~~~~\n"
	"S5\t2145\tchr2\t9000\t60\t5H5M\tchr1\t50\t0\tAAAAA\t~~~~~\n"
	"E3\t145\tchr1\t50\t0\t10M\tchr2\t700\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"L6\t145\tchr1\t50\t60\t10M\tchr2\t400\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"P2\t97\tchr2\t600\t0\t10M\tchr1\t50\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"S1\t145\tchr1\t50\t60\t10M\tchr2\t500\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"S5\t145\tchr1\t50\t60\t10M\tchr2\t1000\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"L6\t97\tchr2\t400\t60\t10M\tchr1\t50\t0\tAAAAAAAAAT\t~~~~~~~~~~\n"
	"S1\t97\tchr2\t500\t60\t5M5S\tchr1\t50\t0\tCGATTAAAAA\t~~~~~~~~~~\n"
	"P2\t145\tchr1\t50\t60\t10M\tchr2\t600\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"E3\t97\tchr2\t700\t60\t10M\tchr1\t50\t0\tAAATTTCCGA\t~~~~~~~~~~\n"
	"S4\t97\tchr2\t800\t60\t5M5S\tchr2\t900\t0\tCGATTAAAAA\t~~~~~~~~~~\n"
	"S4\t145\tchr2\t900\t60\t10M\tchr2\t800\t0\tTTTTTGGGGA\t~~~~~~~~~~\n"
	"S5\t97\tchr2\t1000\t60\t5M5S\tchr1\t50\t0\tCGATTAAAAA\t~~~~~~~~~~\n";

START_TEST (test_abnormal_filter_mapq)
{
	const char *sams[] = {sam_mapq_sorted, sam_mapq_unsorted};
	TestAbnormal a;
	AbnormalStat stat = {};
	sqlite3_stmt *stmt = NULL;
	int i = 0;
	int j = 0;

	const struct {const char *qname; const char *chr; long pos;} rows[] = {
		{"S1", "chr1", 50},
		{"S1", "chr2", 500},
		{"S5", "chr1", 50},
		{"S5", "chr2", 1000}
	};

	for (j = 0; j < 2; j++)
		{
			stat = (AbnormalStat) {};
			test_abnormal_init (&a, sams[j]);
			a.arg->phred_quality = 20;
			a.arg->max_base_freq = 0.75;
			a.arg->hts_threads = 1;
			a.arg->stat = &stat;

			abnormal_filter (a.arg);

			ck_assert_int_eq (stat.stored, 4);

			// Only S1 and S5, without the supplementary pieces,
			// but still marked as supplementary (split reads)
			stmt = db_prepare (a.db,
					"SELECT qname, chr, pos, mapq, type FROM alignment\n"
					"ORDER BY qname ASC, chr ASC");

			for (i = 0; db_step (stmt) == SQLITE_ROW; i++)
				{
					ck_assert_int_lt (i, 4);
					ck_assert_str_eq (db_column_text (stmt, 0), rows[i].qname);
					ck_assert_str_eq (db_column_text (stmt, 1), rows[i].chr);
					ck_assert_int_eq (db_column_int64 (stmt, 2), rows[i].pos);
					ck_assert_int_eq (db_column_int (stmt, 3), 60);
					ck_assert (db_column_int (stmt, 4) & ABNORMAL_SUPPLEMENTARY);
				}

			ck_assert_int_eq (i, 4);

			db_finalize (stmt);
			test_abnormal_destroy (&a);
		}
}
END_TEST

START_TEST (test_abnormal_filter_no_chr)
{
	TestAbnormal a;
	AbnormalStat stat = {};

	/*
	* Mate mapped (no flag 0x8), but without a mate chromosome
	* (RNEXT '*'): Abnormal by chromosome, and the mate
	* chromosome is reported as '*'. The first read is at the
	* exon e1, so the fragment is stored
	*/
	const char *sam_no_chr =
		"@HD\tVN:1.0\tSO:queryname\n"
		"@SQ\tSN:chr1\tLN:248956422\n"
		"S1\t97\tchr1\t50\t60\t10M\t*\t0\t0\tATCGATCGAT\t~~~~~~~~~~\n"
		"S1\t145\tchr1\t9000\t60\t10M\t*\t0\t0\tAAAGGGCCCT\t~~~~~~~~~~\n";

	test_abnormal_init (&a, sam_no_chr);
	a.arg->stat = &stat;

	abnormal_filter (a.arg);

	ck_assert_int_eq (stat.alignments, 2);
	ck_assert_int_eq (stat.abnormal, 2);
	ck_assert_int_eq (stat.stored, 2);

	test_abnormal_destroy (&a);
}
END_TEST

static void
run_abnormal_on (const char *sam_body, const char *path)
{
	TestAbnormal a;

	test_abnormal_init (&a, sam_body);

	if (path != NULL)
		a.arg->sam_file = path;

	abnormal_filter (a.arg);
}

#ifdef HAVE_LD_WRAP
static void
run_abnormal_fault (const char *sam_body, Fault f, int hts_threads)
{
	TestAbnormal a;

	test_abnormal_init (&a, sam_body);
	a.arg->hts_threads = hts_threads;

	fault_inject (f);
	abnormal_filter (a.arg);
}

/* Two fragments: The second one reuses the trash of the first */
static const char *sam_two_frags =
	"@HD\tVN:1.0\tSO:queryname\n"
	"@SQ\tSN:chr1\tLN:248956422\n"
	"N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\t147\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n"
	"N2\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N2\t147\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n";

START_TEST (test_abnormal_threads_fail)
{
	// Not fatal: Go on without decompression threads
	TestAbnormal a;

	test_abnormal_init (&a, sam_normal);
	a.arg->hts_threads = 1;

	fault_inject (FAULT_HTS_SET_THREADS);
	abnormal_filter (a.arg);

	test_abnormal_destroy (&a);
}
END_TEST

START_TEST (test_abnormal_bam_init1_fail)
{
	run_abnormal_fault (sam_normal, FAULT_BAM_INIT1, 0);
}
END_TEST

START_TEST (test_abnormal_sorted_dup_fail)
{
	run_abnormal_fault (sam_sorted, FAULT_BAM_DUP1, 0);
}
END_TEST

START_TEST (test_abnormal_sorted_copy_fail)
{
	run_abnormal_fault (sam_two_frags, FAULT_BAM_COPY1, 0);
}
END_TEST

START_TEST (test_abnormal_unsorted_dup_fail)
{
	run_abnormal_fault (sam_unsorted, FAULT_BAM_DUP1, 0);
}
END_TEST

START_TEST (test_abnormal_sorted_close_fail)
{
	run_abnormal_fault (sam_sorted, FAULT_HTS_CLOSE, 0);
}
END_TEST

START_TEST (test_abnormal_second_read_fail)
{
	TestAbnormal a;
	const char *c = NULL;
	int records = 0;

	// The first reading goes to the end of file (-1):
	// Fail the first read of the second reading
	for (c = sam_unsorted; *c != '\0'; c++)
		if (*c != '@' && (c == sam_unsorted || c[-1] == '\n'))
			records++;

	test_abnormal_init (&a, sam_unsorted);

	fault_inject_after (FAULT_SAM_READ1, records + 1);
	abnormal_filter (a.arg);
}
END_TEST

START_TEST (test_abnormal_rewind_fail)
{
	run_abnormal_fault (sam_unsorted, FAULT_HTS_CLOSE, 0);
}
END_TEST
#endif

START_TEST (test_abnormal_open_fail)
{
	run_abnormal_on (sam_normal, "/it-must-not-exist.bam");
}
END_TEST

START_TEST (test_abnormal_header_fail)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";
	char bam[] = "/tmp/ponga.bam.XXXXXX";
	int fd = 0;

	// A BAM truncated after the BGZF header of the first block
	create_file (sam, sam_normal);
	fd = xmkstemp (bam);
	close (fd);
	ck_assert_int_eq (sam_to_bam (sam, bam), 1);
	ck_assert_int_eq (truncate (bam, 20), 0);

	run_abnormal_on (sam_normal, bam);
}
END_TEST

static const char *sam_broken_sorted =
	"@HD\tVN:1.0\tSO:queryname\n"
	"@SQ\tSN:chr1\tLN:248956422\n"
	"N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\tponga\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n";

static const char *sam_broken_unsorted =
	"@SQ\tSN:chr1\tLN:248956422\n"
	"N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n"
	"N1\tponga\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n";

START_TEST (test_abnormal_sorted_read_fail)
{
	run_abnormal_on (sam_broken_sorted, NULL);
}
END_TEST

START_TEST (test_abnormal_unsorted_read_fail)
{
	run_abnormal_on (sam_broken_unsorted, NULL);
}
END_TEST

Suite *
make_abnormal_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("Abnormal");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_abnormal_filter_sorted);
	tcase_add_test (tc_core, test_abnormal_filter_unsorted);
	tcase_add_test (tc_core, test_abnormal_filter_seq_qual);
	tcase_add_test (tc_core, test_abnormal_filter_edge);
	tcase_add_test (tc_core, test_abnormal_filter_mapq);
	tcase_add_test (tc_core, test_abnormal_filter_no_chr);
#ifdef HAVE_LD_WRAP
	tcase_add_test (tc_core, test_abnormal_threads_fail);
#endif
	suite_add_tcase (s, tc_core);

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_exit_test (tc_abort, test_abnormal_open_fail,          EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_header_fail,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_sorted_read_fail,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_unsorted_read_fail, EXIT_FAILURE);
#ifdef HAVE_LD_WRAP
	tcase_add_exit_test (tc_abort, test_abnormal_bam_init1_fail,     EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_sorted_dup_fail,    EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_sorted_copy_fail,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_unsorted_dup_fail,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_sorted_close_fail,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_rewind_fail,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_abnormal_second_read_fail,   EXIT_FAILURE);
#endif
	suite_add_tcase (s, tc_abort);

	return s;
}
