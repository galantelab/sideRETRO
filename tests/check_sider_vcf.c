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
#include "../src/vcf.h"

static sqlite3 *
create_db (char *db_path)
{
	int fd;

	fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
create_vcf (char *vcf_path)
{
	int fd;

	fd = xmkstemp (vcf_path);
	close (fd);
}

static void
create_fasta (char *fasta_file)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (fasta_file);
	fp = xfdopen (fd, "w");

	xfprintf (fp, "%s",
		">chr10\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		">chr11\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
		">chr12\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		">chr13\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		">chr14\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		">chrY\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
	);

	xfclose (fp);
}

static void
populate_db (sqlite3 *db)
{
	// Database dump
	static const char schema[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO source VALUES (1,1,'PONGA');\n"
		"INSERT INTO exon VALUES (1,'gene1','chr1',1,3000,'+','eg1','ee1');\n"
		"INSERT INTO exon VALUES (13,'gene1','chr1',4000,5000,'+','eg1','ee13');\n"
		"INSERT INTO junction VALUES (1,1,13,5);\n"
		"INSERT INTO evidence VALUES (1,1,1);\n"
		"INSERT INTO evidence VALUES (2,2,2);\n"
		"INSERT INTO evidence VALUES (2,2,4);\n"
		"INSERT INTO evidence VALUES (2,3,8);\n"
		"INSERT INTO insertion VALUES (1,1050,1,7,3000,13,2,2000,'ACGT','GGTTAACC','TTGGCC');\n"
		"INSERT INTO insertion VALUES (8,NULL,NULL,0,150,1,3,NULL,NULL,NULL,NULL);\n"
		"INSERT INTO exon VALUES (2,'gene2_1','chr2',1,3000,'-','eg2','ee2');\n"
		"INSERT INTO exon VALUES (3,'gene2_2','chr2',2000,5000,'-','eg3','ee3');\n"
		"INSERT INTO exon VALUES (4,'gene3_1','chr3',1000,3000,'+','eg4','ee4');\n"
		"INSERT INTO exon VALUES (5,'gene3_2','chr3',5000,8000,'+','eg5','ee5');\n"
		"INSERT INTO exon VALUES (6,'gene4_1','chr4',1000,5000,'+','eg6','ee6');\n"
		"INSERT INTO exon VALUES (7,'gene4_2','chr5',1000,5000,'+','eg7','ee7');\n"
		"INSERT INTO exon VALUES (8,'gene5_1','chr6',1,3000,'-','eg8','ee8');\n"
		"INSERT INTO exon VALUES (9,'gene5_2','chr6',2000,5000,'-','eg9','ee9');\n"
		"INSERT INTO exon VALUES (10,'gene5_3','chr7',10000,13000,'+','eg10','ee10');\n"
		"INSERT INTO exon VALUES (11,'gene5_4','chr7',15000,18000,'+','eg11','ee11');\n"
		"INSERT INTO exon VALUES (12,'gene6','chr10',1,300,'+','eg12','ee12');\n"
		"INSERT INTO alignment VALUES (1,'q1',0x800,'chr10',1,20,'100M10S',100,100,'chr10',1,1,1,'*','*');\n"
		"INSERT INTO alignment VALUES (12,'q1',0x800,'chr10',1,20,'100M10S',100,100,'chr10',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (2,'A00123:45:HV7KJDSXY:1:1101:12345:6789',97,'chr11',1,20,'110M',100,50,'chr11',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (3,'q5',97,'chr11',200,20,'110M',100,50,'chr11',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (4,'q2',0x800,'chr12',250,20,'10H100M',100,100,'chr12',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (5,'q3',0x800,'chr12',200,20,'100M10S',100,50,'chr1',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (6,'q6',97,'chr13',1,20,'100M10S',100,50,'chr13',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (7,'q7',97,'chr13',200,20,'100M10S',100,50,'chr13',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (8,'q8',97,'chr14',1,20,'100M10S',100,50,'chr14',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (9,'q9',97,'chr14',200,20,'100M10S',100,50,'chr14',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (10,'q10',97,'chr14',400,20,'100M10S',100,50,'chr14',1,8,1,'*','*');\n"
		"INSERT INTO alignment VALUES (11,'q11',97,'chr14',500,20,'100M10S',100,50,'chr14',1,8,1,'*','*');\n"
		"INSERT INTO clustering VALUES (1,2,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,2,2,3,100);\n"
		"INSERT INTO clustering VALUES (3,2,3,3,100);\n"
		"INSERT INTO clustering VALUES (4,2,4,3,100);\n"
		"INSERT INTO clustering VALUES (5,2,5,3,100);\n"
		"INSERT INTO clustering VALUES (6,2,6,3,100);\n"
		"INSERT INTO clustering VALUES (7,2,7,3,100);\n"
		"INSERT INTO clustering VALUES (8,2,8,3,100);\n"
		"INSERT INTO clustering VALUES (9,2,9,3,100);\n"
		"INSERT INTO clustering VALUES (10,2,10,3,100);\n"
		"INSERT INTO clustering VALUES (11,2,11,3,100);\n"
		"INSERT INTO cluster VALUES (1,2,'chr10',1,300,'gene1',31);\n"
		"INSERT INTO cluster VALUES (2,2,'chr11',1,300,'gene2_1',31);\n"
		"INSERT INTO cluster VALUES (3,2,'chr11',200,500,'gene2_2',31);\n"
		"INSERT INTO cluster VALUES (4,2,'chr12',1,300,'gene3_1',31);\n"
		"INSERT INTO cluster VALUES (5,2,'chr12',200,500,'gene3_2',31);\n"
		"INSERT INTO cluster VALUES (6,2,'chr13',1,300,'gene4_1',31);\n"
		"INSERT INTO cluster VALUES (7,2,'chr13',200,500,'gene4_2',31);\n"
		"INSERT INTO cluster VALUES (8,2,'chr14',1,300,'gene5_1',31);\n"
		"INSERT INTO cluster VALUES (9,2,'chr14',200,500,'gene5_2',31);\n"
		"INSERT INTO cluster VALUES (10,2,'chr14',400,600,'gene5_3',31);\n"
		"INSERT INTO cluster VALUES (11,2,'chr14',500,700,'gene5_4',31);\n"
		"INSERT INTO cluster_merging VALUES(1,1,2);\n"
		"INSERT INTO cluster_merging VALUES(2,2,2);\n"
		"INSERT INTO cluster_merging VALUES(2,3,2);\n"
		"INSERT INTO cluster_merging VALUES(3,4,2);\n"
		"INSERT INTO cluster_merging VALUES(3,5,2);\n"
		"INSERT INTO cluster_merging VALUES(4,6,2);\n"
		"INSERT INTO cluster_merging VALUES(5,7,2);\n"
		"INSERT INTO cluster_merging VALUES(6,8,2);\n"
		"INSERT INTO cluster_merging VALUES(6,9,2);\n"
		"INSERT INTO cluster_merging VALUES(7,10,2);\n"
		"INSERT INTO cluster_merging VALUES(7,11,2);\n"
		"INSERT INTO cluster_merging VALUES(8,1,2);\n"
		"INSERT INTO retrocopy VALUES(1,'chr10',1,300,'gene1',1,101,2,-1,0,100,2,NULL,0,25,0,'+',6);\n"
		"INSERT INTO retrocopy VALUES(2,'chr11',1,500,'gene2_1/gene2_2',2,350,1,0.0,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES(3,'chr12',1,500,'gene3_1/gene3_2',4,250,2,0.0,0.0,259,3,250,4,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES(4,'chr13',1,300,'gene4_1',8,150,1,0.0,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES(5,'chr13',200,500,'gene4_2',8,350,1,0.0,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES(6,'chr14',1,500,'gene5_1/gene5_2',10,350,1,0.0,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES(7,'chr14',400,700,'gene5_3/gene5_4',12,600,1,0.0,0.0,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES(8,'chrY',1,300,'gene1',1,101,2,1,0,100,1,102,1,30,1,'-',2);\n"
		"INSERT INTO genotype VALUES(1,1,0,0,0.0,0.0,0.0);\n"
		"INSERT INTO genotype VALUES(1,2,0,0,0.0,0.0,0.0);\n"
		"INSERT INTO genotype VALUES(1,3,0,0,0.0,0.0,0.0);\n"
		"INSERT INTO genotype VALUES(1,4,10,10,-33.0,-0.1,-5.0);\n"
		"INSERT INTO genotype VALUES(1,5,10,10,-33.0,-0.1,-5.0);\n"
		"INSERT INTO genotype VALUES(1,6,10,10,-33.0,-0.1,-5.0);\n"
		"INSERT INTO genotype VALUES(1,7,10,10,-33.0,-0.1,-5.0);\n"
		"INSERT INTO genotype VALUES(1,8,10,10,-33.0,-0.1,-5.0);\n"
		"COMMIT;";

	db_exec (db, schema);
}

START_TEST (test_vcf)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	char rnames_file[] = "/tmp/ponga.rnames.XXXXXX";
	char insertion_file[] = "/tmp/ponga.ins.XXXXXX";

	sqlite3 *db = NULL;

	db = create_db (db_file);
	create_vcf (vcf_file);
	create_vcf (rnames_file);
	create_vcf (insertion_file);
	create_fasta (fasta_file);

	// Populate database
	populate_db (db);

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_sr = VCF_DEFAULT_MIN_SR,
		.max_tsd = VCF_DEFAULT_MAX_TSD,
		.fasta_file = fasta_file,
		.rnames_file = rnames_file,
		.insertion_file = insertion_file
	};

	vcf (db, vcf_file, &opt);

	/* TIME TO TEST */
	FILE *fp = NULL;
	char *line = NULL;
	size_t len = 0;
	int i = 0;

	// INFO fields that must be present (or absent)
	// in each retrocopy line, in VCF order. POLARITY
	// comes from ORHO; chr10 and chrY have the poly-A
	// at the other strand (POLYASTRAND): It does not filter
	// (chr10, chr11, chr12, chr13, chr13, chr14, chr14, chrY)
	const char *present[][2] = {
		{";ORHO=-1.000000;POLARITY=-;PG=gene1;PGTYPE=1;PGEXONS=gene1:2;HALLMARKS=POLYA,EEJ;EEJ=1;EEJLIST=ee1-ee13;PG5P=ee1:1050:7;PG3P=ee13:3000:2;SVLEN=2000;LEFT_SVINSSEQ=GGTTAACC;RIGHT_SVINSSEQ=TTGGCC;", ";SR=2;LBKP=100;LSR=2;TSITETR=A:101;POLYA=25;POLYASR=6;POLYASTRAND=+"},
		{";PGEXONS=gene2_1:1,gene2_2:1;",   ";DP="},
		{";SR=7;LBKP=259;LSR=3;RBKP=250;RSR=4;TSDLEN=10;TSITE=DUP;TSD=CCCCCCCCCC", NULL},
		{";DP=",                            NULL},
		{";DP=",                            NULL},
		{";DP=",                            NULL},
		{";DP=",                            NULL},
		{";ORHO=1.000000;POLARITY=+;PG=gene1;PGTYPE=1;PGEXONS=gene1:2;HALLMARKS=POLYA;PG3P=ee1:150:3;DP=", ";SR=2;LBKP=100;LSR=1;RBKP=102;RSR=1;TSDLEN=-1;TSITE=DEL;TSITETR=A:103;POLYA=30;POLYAMIN;POLYASR=2;POLYASTRAND=-"}
	};

	// QUAL and FILTER. Samples without reads give no QUAL. The
	// others: -10 * log10 P(0/0) with GL (-33, -0.1, -5) = 329.0
	const char *qual_filter[] = {
		"\t.\tLowSR;OneBreakpoint\tSVTYPE=INS;",
		"\t.\tAmbiguousParental;LowSR;NoBreakpoint\tSVTYPE=INS;",
		"\t.\tAmbiguousParental\tSVTYPE=INS;",
		"\t329.0\tAmbiguousParental;LowSR;NoBreakpoint\tSVTYPE=INS;",
		"\t329.0\tAmbiguousParental;LowSR;NoBreakpoint\tSVTYPE=INS;",
		"\t329.0\tAmbiguousParental;LowSR;NoBreakpoint\tSVTYPE=INS;",
		"\t329.0\tAmbiguousParental;LowSR;NoBreakpoint\tSVTYPE=INS;",
		"\t329.0\tLowSR\tSVTYPE=INS;"
	};

	const char *absent[] = {
		";RBKP=", ";EEJ", ";POLYA", ";TSD", ";TSD", ";TSD", ";TSD", ";TSD="
	};

	fp = xfopen (vcf_file, "r");

	while (getline (&line, &len, fp) != EOF)
		{
			if (line[0] == '#')
				continue;

			ck_assert_int_lt (i, 8);
			ck_assert_msg (strstr (line, present[i][0]) != NULL,
					"'%s' not found at line %d: %s", present[i][0], i, line);
			ck_assert_msg (present[i][1] == NULL || strstr (line, present[i][1]) != NULL,
					"'%s' not found at line %d: %s", present[i][1], i, line);
			ck_assert_msg (strstr (line, absent[i]) == NULL,
					"'%s' found at line %d: %s", absent[i], i, line);
			ck_assert_msg (strstr (line, qual_filter[i]) != NULL,
					"'%s' not found at line %d: %s", qual_filter[i], i, line);
			i++;
		}

	ck_assert_int_eq (i, 8);
	xfclose (fp);

	// The ID column, in VCF order
	const char *ids[] = {
		"\tRTC000001\t", "\tRTC000002\t", "\tRTC000003\t", "\tRTC000004\t",
		"\tRTC000005\t", "\tRTC000006\t", "\tRTC000007\t", "\tRTC000008\t"
	};

	fp = xfopen (vcf_file, "r");

	for (i = 0; getline (&line, &len, fp) != EOF;)
		{
			if (line[0] == '#')
				continue;

			ck_assert_msg (strstr (line, ids[i]) != NULL,
					"'%s' not found at line %d: %s", ids[i], i, line);
			i++;
		}

	xfclose (fp);

	// Read names by retrocopy: 'sample:read,read'
	// One read per line, sorted by read name for each sample.
	// Read names may have colons (Illumina). The roles come from
	// the evidence of each retrocopy: q1 is split_left only at
	// RTC000001
	const char *rnames[] = {
		"#retrocopy_id\tsample\tread_name\trole\n",
		"RTC000001\tPONGA\tq1\tsplit_left\n",
		"RTC000002\tPONGA\tA00123:45:HV7KJDSXY:1:1101:12345:6789\tsplit_right,polya\n",
		"RTC000002\tPONGA\tq5\tjunction\n",
		"RTC000003\tPONGA\tq2\tdiscordant\n",
		"RTC000003\tPONGA\tq3\tdiscordant\n",
		"RTC000004\tPONGA\tq6\tdiscordant\n",
		"RTC000005\tPONGA\tq7\tdiscordant\n",
		"RTC000006\tPONGA\tq8\tdiscordant\n",
		"RTC000006\tPONGA\tq9\tdiscordant\n",
		"RTC000007\tPONGA\tq10\tdiscordant\n",
		"RTC000007\tPONGA\tq11\tdiscordant\n",
		"RTC000008\tPONGA\tq1\tdiscordant\n"
	};

	const int rnames_size = sizeof (rnames) / sizeof (char *);

	fp = xfopen (rnames_file, "r");

	for (i = 0; getline (&line, &len, fp) != EOF; i++)
		{
			ck_assert_int_lt (i, rnames_size);
			ck_assert_str_eq (line, rnames[i]);
		}

	ck_assert_int_eq (i, rnames_size);
	xfclose (fp);

	// Predicted insertions: Only RTC000001 has one
	const char *fasta[] = {
		">RTC000001 chr10:100 PG=gene1 SVLEN=2000 predicted\n",
		"ACGT\n"
	};

	fp = xfopen (insertion_file, "r");

	for (i = 0; getline (&line, &len, fp) != EOF; i++)
		{
			ck_assert_int_lt (i, 2);
			ck_assert_str_eq (line, fasta[i]);
		}

	ck_assert_int_eq (i, 2);
	xfclose (fp);

	// Thresholds and disabled filters: The disabled ones are
	// neither checked nor declared at the header. One breakpoint
	// is enough (chr10), none is not (chr11); the TSD of 10 bases
	// (chr12) is longer than the maximum of 5, and the deletion
	// of 1 base (chrY) is not checked
	VCFOption opt2 = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file       = fasta_file,
		.min_qual         = 400,
		.min_mapq         = 30,
		.min_sr           = 2,
		.max_tsd          = 5,
		.disabled_filters = VCF_FILTER_AMBIGUOUS_PARENTAL
			| VCF_FILTER_LOW_SR | VCF_FILTER_ONE_BREAKPOINT
	};

	vcf (db, vcf_file, &opt2);

	const char *filters2[] = {
		"\tLowMAPQ\tSVTYPE=INS;",
		"\tLowMAPQ;NoBreakpoint\tSVTYPE=INS;",
		"\tLowMAPQ;TSDOutlier\tSVTYPE=INS;",
		"\tLowQual;LowMAPQ;NoBreakpoint\tSVTYPE=INS;",
		"\tLowQual;LowMAPQ;NoBreakpoint\tSVTYPE=INS;",
		"\tLowQual;LowMAPQ;NoBreakpoint\tSVTYPE=INS;",
		"\tLowQual;LowMAPQ;NoBreakpoint\tSVTYPE=INS;",
		"\tLowQual;LowMAPQ\tSVTYPE=INS;"
	};

	fp = xfopen (vcf_file, "r");

	for (i = 0; getline (&line, &len, fp) != EOF;)
		{
			if (line[0] == '#')
				{
					ck_assert_msg (strstr (line, "ID=AmbiguousParental") == NULL
							&& strstr (line, "ID=LowSR") == NULL
							&& strstr (line, "ID=OneBreakpoint") == NULL,
							"Disabled filter declared: %s", line);
					continue;
				}

			ck_assert_int_lt (i, 8);
			ck_assert_msg (strstr (line, filters2[i]) != NULL,
					"'%s' not found at line %d: %s", filters2[i], i, line);
			i++;
		}

	ck_assert_int_eq (i, 8);
	ck_assert_int_eq (opt2.num_records, 8);
	ck_assert_int_eq (opt2.num_pass, 0);
	xfclose (fp);

	// Filter names
	ck_assert_int_eq (vcf_filter_lookup ("NoBreakpoint"), VCF_FILTER_NO_BREAKPOINT);
	ck_assert_int_eq (vcf_filter_lookup ("OneBreakpoint"), VCF_FILTER_ONE_BREAKPOINT);
	ck_assert_int_eq (vcf_filter_lookup ("NoTSD"), 0);
	ck_assert_int_eq (vcf_filter_lookup ("LowSR"), VCF_FILTER_LOW_SR);
	ck_assert_int_eq (vcf_filter_lookup ("PolarityConflict"), 0);
	ck_assert_int_eq (vcf_filter_lookup ("Ponga"), 0);

	xfree (line);

	db_close (db);
	xunlink (db_file);
	xunlink (fasta_file);
	xunlink (vcf_file);
	xunlink (rnames_file);
	xunlink (insertion_file);
}
END_TEST

static int
grep_vcf (const char *path, const char *needle)
{
	FILE *fp = xfopen (path, "r");
	char *line = NULL;
	size_t len = 0;
	int found = 0;

	while (!found && getline (&line, &len, fp) != EOF)
		found = line[0] != '#' && strstr (line, needle) != NULL;

	xfree (line);
	xfclose (fp);
	return found;
}

/*
 * The record of the retrocopy 'id' has 'needle'
 */
static int
grep_vcf_record (const char *path, const char *id, const char *needle)
{
	FILE *fp = xfopen (path, "r");
	char *line = NULL;
	size_t len = 0;
	int found = 0;

	while (!found && getline (&line, &len, fp) != EOF)
		found = line[0] != '#' && strstr (line, id) != NULL
			&& strstr (line, needle) != NULL;

	xfree (line);
	xfclose (fp);
	return found;
}

START_TEST (test_vcf_edge)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	sqlite3 *db = NULL;
	FILE *fp = NULL;
	int fd = 0;

	/*
	* RTC000001: chr1:1, the first base - POS is 1
	* RTC000002: chr10:500, inside the host gene but not at its
	* exons (INTRONIC), with a TSD of 60 bases (TSDOutlier)
	* RTC000003: chr10:4000, 1000 bases before 'neargene' (NEAR),
	* farther than 2000 bases after the end of 'host' (900)
	* RTC000004: chrZ, absent from the reference: REF is 'N'
	*/
	const char sql[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO source VALUES (1,1,'PONGA');\n"
		"INSERT INTO exon VALUES (1,'host','chr10',100,200,'+','egh','eeh1');\n"
		"INSERT INTO exon VALUES (2,'host','chr10',800,900,'+','egh','eeh2');\n"
		"INSERT INTO exon VALUES (3,'gene1','chr5',1,3000,'+','eg1','ee1');\n"
		"INSERT INTO exon VALUES (4,'neargene','chr10',5000,6000,'-','egn','een');\n"
		"INSERT INTO alignment VALUES (1,'q1',97,'chr1',1,60,'100M',100,100,'chr5',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (2,'q2',97,'chr10',450,60,'100M',100,100,'chr5',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (3,'q3',97,'chr10',3950,60,'100M',100,100,'chr5',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (4,'q4',97,'chrZ',100,60,'100M',100,100,'chr5',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (1,1,'chr1',1,100,'gene1',31);\n"
		"INSERT INTO cluster VALUES (2,1,'chr10',450,550,'gene1',31);\n"
		"INSERT INTO cluster VALUES (3,1,'chr10',3950,4050,'gene1',31);\n"
		"INSERT INTO cluster VALUES (4,1,'chrZ',100,200,'gene1',31);\n"
		"INSERT INTO clustering VALUES (1,1,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,1,2,3,100);\n"
		"INSERT INTO clustering VALUES (3,1,3,3,100);\n"
		"INSERT INTO clustering VALUES (4,1,4,3,100);\n"
		"INSERT INTO cluster_merging VALUES (1,1,1);\n"
		"INSERT INTO cluster_merging VALUES (2,2,1);\n"
		"INSERT INTO cluster_merging VALUES (3,3,1);\n"
		"INSERT INTO cluster_merging VALUES (4,4,1);\n"
		"INSERT INTO retrocopy VALUES (1,'chr1',1,100,'gene1',1,1,2,0,1,NULL,0,1,3,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (2,'chr10',450,550,'gene1',1,500,2,0,1,559,3,500,3,20,0,'+',3);\n"
		"INSERT INTO retrocopy VALUES (3,'chr10',3950,4050,'gene1',1,4000,1,0,1,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (4,'chrZ',100,200,'gene1',1,150,1,0,1,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,1,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO genotype VALUES (1,2,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO genotype VALUES (1,3,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO genotype VALUES (1,4,0,5,-30.0,-1.5,-0.1);\n"
		"COMMIT;";

	fd = xmkstemp (db_file);
	close (fd);
	db = db_create (db_file);
	db_exec (db, sql);

	fd = xmkstemp (fasta_file);
	fp = xfdopen (fd, "w");
	xfprintf (fp, ">chr1\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n>chr10\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
	xfclose (fp);

	fd = xmkstemp (vcf_file);
	close (fd);

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file = fasta_file,
		.near_gene_dist = 2000,
		.orientation_error = 0.05,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ
	};

	vcf (db, vcf_file, &opt);

	ck_assert (grep_vcf (vcf_file, "chr1\t1\tRTC000001\tG\t"));
	ck_assert (grep_vcf (vcf_file, "chr10\t499\tRTC000002\tG\t<INS:ME:RTC>\t"));
	ck_assert (grep_vcf (vcf_file, "\tTSDOutlier\tSVTYPE=INS;"));
	ck_assert (grep_vcf (vcf_file, ";INTRONIC=host;"));
	ck_assert (grep_vcf (vcf_file, ";TSDLEN=60;"));
	ck_assert (grep_vcf (vcf_file, ";TSD=GTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGG;" ));
	ck_assert (grep_vcf (vcf_file, ";NEAR=neargene;"));
	ck_assert (grep_vcf (vcf_file, "chrZ\t149\tRTC000004\tN\t"));

	db_close (db);
	xunlink (db_file);
	xunlink (vcf_file);
	xunlink (fasta_file);
}
END_TEST

START_TEST (test_vcf_parental_chosen)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char insertion_file[] = "/tmp/ponga.ins.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	sqlite3 *db = NULL;
	FILE *fp = NULL;
	char *line = NULL;
	size_t len = 0;
	int found = 0;

	db = create_db (db_file);
	create_vcf (vcf_file);
	create_vcf (insertion_file);
	create_fasta (fasta_file);
	populate_db (db);

	// The first retrocopy (chr10, gene1): As if gene1 had
	// been chosen over the overlapped geneX and geneY
	db_exec (db,
			"UPDATE retrocopy\n"
			"SET parental_gene_name = 'gene1/geneX/geneY', level = 34\n"
			"WHERE id = 1");

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_sr = VCF_DEFAULT_MIN_SR,
		.max_tsd = VCF_DEFAULT_MAX_TSD,
		.fasta_file = fasta_file,
		.insertion_file = insertion_file
	};

	vcf (db, vcf_file, &opt);

	// PG is the chosen gene, the others go to PGALT and the
	// type is still 'overlapped'. It is not ambiguous, and
	// the polarity comes from the strand of the chosen gene
	fp = xfopen (vcf_file, "r");

	while (getline (&line, &len, fp) != EOF)
		{
			if (strncmp (line, "chr10\t", 6))
				continue;

			found++;
			ck_assert_msg (strstr (line,
						";ORHO=-1.000000;POLARITY=-;PG=gene1;PGTYPE=2;"
						"PGALT=geneX,geneY;PGEXONS=") != NULL, "%s", line);
			ck_assert_msg (strstr (line, "\tLowSR;OneBreakpoint\t") != NULL, "%s", line);
		}

	ck_assert_int_eq (found, 1);
	xfclose (fp);

	// The predicted insertion is named after the chosen gene
	fp = xfopen (insertion_file, "r");
	ck_assert_int_ne (getline (&line, &len, fp), EOF);
	ck_assert_str_eq (line, ">RTC000001 chr10:100 PG=gene1 SVLEN=2000 predicted\n");
	xfclose (fp);

	xfree (line);
	db_close (db);
	xunlink (db_file);
	xunlink (vcf_file);
	xunlink (insertion_file);
	xunlink (fasta_file);
}
END_TEST

START_TEST (test_vcf_primer)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	char primer_file[] = "/tmp/ponga.primer.XXXXXX";
	sqlite3 *db = NULL;
	FILE *fp = NULL;
	char buf[BUFSIZ] = {};
	size_t n = 0;
	int fd = 0;

	/*
	* Sequences for primer design, with 10 bases of flank:
	* RTC000001: Both breakpoints (TSD of 8 bases) and both
	* observed sequences. Its parental gene is at '-', and two
	* exons share the same end: The junction is stored with
	* the shorter one, but the flank comes from the longer
	* RTC000002: Only the left breakpoint, at the start of the
	* contig, and no observed sequence. Exons shorter than the
	* flank at its junction
	* RTC000003: Only the right breakpoint, at the end of the
	* contig: The flank is outside it
	* RTC000004: No breakpoints - nothing to write
	* RTC000005: Inconsistent breakpoints - no empty allele
	* RTC000006: Breakpoints with 2 bases between them
	* RTC000007: Contig and exons absent from the reference
	*/
	const char sql[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO source VALUES (1,1,'PONGA');\n"
		"INSERT INTO exon VALUES (1,'geneM','chr2',11,30,'-','egm','eem1');\n"
		"INSERT INTO exon VALUES (2,'geneM','chr2',26,30,'-','egm','eem2');\n"
		"INSERT INTO exon VALUES (3,'geneM','chr2',51,70,'-','egm','eem3');\n"
		"INSERT INTO exon VALUES (4,'geneP','chr2',81,85,'+','egp','eep1');\n"
		"INSERT INTO exon VALUES (5,'geneP','chr2',91,94,'+','egp','eep2');\n"
		"INSERT INTO exon VALUES (6,'geneZ','chrZ',11,30,'+','egz','eez1');\n"
		"INSERT INTO exon VALUES (7,'geneZ','chrZ',51,70,'+','egz','eez2');\n"
		"INSERT INTO junction VALUES (1,2,3,6);\n"
		"INSERT INTO junction VALUES (2,4,5,2);\n"
		"INSERT INTO junction VALUES (7,6,7,3);\n"
		"INSERT INTO insertion VALUES (1,NULL,NULL,0,NULL,NULL,0,NULL,NULL,'ACGTACGT','TTTTGGGG');\n"
		"INSERT INTO insertion VALUES (3,NULL,NULL,0,NULL,NULL,0,NULL,NULL,NULL,'CC');\n"
		"INSERT INTO insertion VALUES (5,NULL,NULL,0,NULL,NULL,0,NULL,NULL,'AA',NULL);\n"
		"INSERT INTO insertion VALUES (7,NULL,NULL,0,NULL,NULL,0,NULL,NULL,'GG','TT');\n"
		"INSERT INTO alignment VALUES (1,'q1',97,'chr1',54,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (1,1,'chr1',54,154,'geneM',31);\n"
		"INSERT INTO clustering VALUES (1,1,1,3,100);\n"
		"INSERT INTO cluster_merging VALUES (1,1,1);\n"
		"INSERT INTO retrocopy VALUES (1,'chr1',54,154,'geneM',1,54,2,0,1,60,3,53,3,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,1,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO alignment VALUES (2,'q2',97,'chr1',6,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (2,1,'chr1',6,106,'geneP',31);\n"
		"INSERT INTO clustering VALUES (2,1,2,3,100);\n"
		"INSERT INTO cluster_merging VALUES (2,2,1);\n"
		"INSERT INTO retrocopy VALUES (2,'chr1',6,106,'geneP',1,6,2,0,1,5,2,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,2,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO alignment VALUES (3,'q3',97,'chr1',195,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (3,1,'chr1',195,295,'geneP',31);\n"
		"INSERT INTO clustering VALUES (3,1,3,3,100);\n"
		"INSERT INTO cluster_merging VALUES (3,3,1);\n"
		"INSERT INTO retrocopy VALUES (3,'chr1',195,295,'geneP',1,195,2,0,1,NULL,0,195,2,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,3,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO alignment VALUES (4,'q4',97,'chr1',120,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (4,1,'chr1',120,220,'geneP',31);\n"
		"INSERT INTO clustering VALUES (4,1,4,3,100);\n"
		"INSERT INTO cluster_merging VALUES (4,4,1);\n"
		"INSERT INTO retrocopy VALUES (4,'chr1',120,220,'geneP',1,120,2,0,1,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,4,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO alignment VALUES (5,'q5',97,'chr1',140,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (5,1,'chr1',140,240,'geneP',31);\n"
		"INSERT INTO clustering VALUES (5,1,5,3,100);\n"
		"INSERT INTO cluster_merging VALUES (5,5,1);\n"
		"INSERT INTO retrocopy VALUES (5,'chr1',140,240,'geneP',1,140,2,0,1,150,3,20,3,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,5,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO alignment VALUES (6,'q6',97,'chr1',103,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (6,1,'chr1',103,203,'geneP',31);\n"
		"INSERT INTO clustering VALUES (6,1,6,3,100);\n"
		"INSERT INTO cluster_merging VALUES (6,6,1);\n"
		"INSERT INTO retrocopy VALUES (6,'chr1',103,203,'geneP',1,103,2,0,1,100,3,103,3,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,6,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO alignment VALUES (7,'q7',97,'chrZ',46,60,'100M',100,100,'chr2',1,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (7,1,'chrZ',46,146,'geneZ',31);\n"
		"INSERT INTO clustering VALUES (7,1,7,3,100);\n"
		"INSERT INTO cluster_merging VALUES (7,7,1);\n"
		"INSERT INTO retrocopy VALUES (7,'chrZ',46,146,'geneZ',1,46,2,0,1,50,3,45,3,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,7,0,5,-30.0,-1.5,-0.1);\n"
		"COMMIT;";

	const char expected[] =
		">RTC000002_empty chr1:1-15 PG=geneP FILTER=LowSR;OneBreakpoint junction=5\n"
		"TTTCCTCATGCAATT\n"
		">RTC000002_eej1 chr2:81-85,91-94 PG=geneP FILTER=LowSR;OneBreakpoint junction=5 strand=+ reads=2\n"
		"CCTGCGGCT\n"
		">RTC000001_empty chr1:43-70 PG=geneM FILTER=PASS junction=10\n"
		"ATAGTAAACCATTTTACGGAGGATACCA\n"
		">RTC000001_left chr1:51-60 PG=geneM FILTER=PASS junction=10\n"
		"CCATTTTACGACGTACGT\n"
		">RTC000001_right chr1:53-62 PG=geneM FILTER=PASS junction=8\n"
		"TTTTGGGGATTTTACGGA\n"
		">RTC000001_eej1 chr2:21-30,51-60 PG=geneM FILTER=PASS junction=10 strand=- reads=6\n"
		"TTCACAAAGCACGAACGGGC\n"
		">RTC000006_empty chr1:91-112 PG=geneP FILTER=PASS junction=10\n"
		"CTAACCTGAGGTAAACCAGGTC\n"
		">RTC000005_left chr1:141-150 PG=geneP FILTER=TSDOutlier junction=10\n"
		"CCTAGCCAAGAA\n";

	fd = xmkstemp (db_file);
	close (fd);
	db = db_create (db_file);
	db_exec (db, sql);

	fd = xmkstemp (fasta_file);
	fp = xfdopen (fd, "w");
	xfprintf (fp, "%s",
		">chr1\n"
		"TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACG\n"
		"GAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC\n"
		"CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA\n"
		"TGACGGATATATATTAAAAA\n"
		">chr2\n"
		"GTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAA\n"
		"GAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGA\n");
	xfclose (fp);

	create_vcf (vcf_file);
	create_vcf (primer_file);

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file = fasta_file,
		.orientation_error = 0.05,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_sr = VCF_DEFAULT_MIN_SR,
		.max_tsd = VCF_DEFAULT_MAX_TSD,
		.primer_file = primer_file,
		.primer_flank = 10
	};

	vcf (db, vcf_file, &opt);

	fp = xfopen (primer_file, "r");
	n = fread (buf, 1, sizeof (buf) - 1, fp);
	buf[n] = '\0';
	xfclose (fp);

	ck_assert_str_eq (buf, expected);
	ck_assert_int_eq (opt.num_primers, 8);

	xunlink (primer_file);

	db_close (db);
	xunlink (db_file);
	xunlink (vcf_file);
	xunlink (fasta_file);
}
END_TEST

START_TEST (test_vcf_intron_retention)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	sqlite3 *db = NULL;
	int fd = 0;

	/*
	* Reads with the mate inside the introns of 'gene1':
	* - intron between ee1 and ee2: 1 read at each sample
	* - intron between ee2 and ee3: 3 reads at the first
	*   sample and 1 at the second one
	* The second retrocopy has no intronic mates
	*/
	const char sql[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO source VALUES (1,1,'PONGA1');\n"
		"INSERT INTO source VALUES (2,1,'PONGA2');\n"
		"INSERT INTO exon VALUES (1,'gene1','chr5',100,200,'+','eg1','ee1');\n"
		"INSERT INTO exon VALUES (2,'gene1','chr5',400,500,'+','eg1','ee2');\n"
		"INSERT INTO exon VALUES (3,'gene1','chr5',800,900,'+','eg1','ee3');\n"
		"INSERT INTO alignment VALUES (1,'q1',97,'chr1',50,60,'100M',100,100,'chr5',100,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (2,'q2',97,'chr2',50,60,'100M',100,100,'chr5',100,2,1,'*','*');\n"
		"INSERT INTO cluster VALUES (1,1,'chr1',50,150,'gene1',31);\n"
		"INSERT INTO cluster VALUES (2,1,'chr2',50,150,'gene1',31);\n"
		"INSERT INTO clustering VALUES (1,1,1,3,100);\n"
		"INSERT INTO clustering VALUES (2,1,2,3,100);\n"
		"INSERT INTO cluster_merging VALUES (1,1,1);\n"
		"INSERT INTO cluster_merging VALUES (2,2,1);\n"
		"INSERT INTO retrocopy VALUES (1,'chr1',50,150,'gene1',1,100,1,0,1,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO retrocopy VALUES (2,'chr2',50,150,'gene1',1,100,1,0,1,NULL,0,NULL,0,NULL,0,NULL,0);\n"
		"INSERT INTO genotype VALUES (1,1,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO genotype VALUES (2,1,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO genotype VALUES (1,2,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO genotype VALUES (2,2,0,5,-30.0,-1.5,-0.1);\n"
		"INSERT INTO intron_retention VALUES (1,1,'gene1','chr5',201,399,1,2,1);\n"
		"INSERT INTO intron_retention VALUES (1,2,'gene1','chr5',201,399,1,2,1);\n"
		"INSERT INTO intron_retention VALUES (1,1,'gene1','chr5',501,799,2,3,3);\n"
		"INSERT INTO intron_retention VALUES (1,2,'gene1','chr5',501,799,2,3,1);\n"
		"COMMIT;";

	fd = xmkstemp (db_file);
	close (fd);
	db = db_create (db_file);
	db_exec (db, sql);

	create_vcf (vcf_file);
	create_fasta (fasta_file);

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file = fasta_file,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_ir = VCF_DEFAULT_MIN_IR
	};

	// An intron needs a sample with enough reads: Two reads
	// from two samples do not report the first intron. IR sums
	// the reads of all samples at the introns reported
	vcf (db, vcf_file, &opt);

	ck_assert (grep_vcf (vcf_file, "chr1\t99\tRTC000001\t"));
	ck_assert (grep_vcf (vcf_file, ";IR=4;IRLIST=ee2-ee3:4;"));
	ck_assert (!grep_vcf (vcf_file, "ee1-ee2:"));

	// A single read is enough: Both introns, by their position
	opt.min_ir = 1;
	vcf (db, vcf_file, &opt);

	ck_assert (grep_vcf (vcf_file, ";IR=6;IRLIST=ee1-ee2:2,ee2-ee3:4;"));

	// No sample has 4 reads at any intron
	opt.min_ir = 4;
	vcf (db, vcf_file, &opt);

	ck_assert (grep_vcf (vcf_file, "chr1\t99\tRTC000001\t"));
	ck_assert (!grep_vcf (vcf_file, "IR="));

	db_close (db);
	xunlink (db_file);
	xunlink (vcf_file);
	xunlink (fasta_file);
}
END_TEST

START_TEST (test_vcf_target_site)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	sqlite3 *db = NULL;

	db = create_db (db_file);
	create_vcf (vcf_file);
	create_fasta (fasta_file);
	populate_db (db);

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file = fasta_file,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_sr = VCF_DEFAULT_MIN_SR,
		.max_tsd = VCF_DEFAULT_MAX_TSD,
		.disabled_filters = VCF_FILTER_LOW_SR
	};

	/*
	* The chrY retrocopy (QUAL 329.0, the only one at 329.0
	* that may pass) with its breakpoints moved: Duplication,
	* blunt insertion and deletion are all valid target sites.
	* Only a duplication too long fails, and one breakpoint
	* or none
	*/
	struct
	{
		const char *update;
		const char *info;
		const char *filter;
		const char *hallmarks;
	} cases[] = {
		{"breakpoint_right = 91, breakpoint_right_sr = 1",
			";TSDLEN=10;TSITE=DUP;", "\t329.0\tPASS\t", ";HALLMARKS=TSD,"},
		{"breakpoint_right = 96, breakpoint_right_sr = 1",
			";TSDLEN=5;TSITE=DUP;", "\t329.0\tPASS\t", ";HALLMARKS=TSD,"},
		{"breakpoint_right = 98, breakpoint_right_sr = 1",
			";TSDLEN=3;TSITE=BLUNT;", "\t329.0\tPASS\t", ";HALLMARKS=POLYA"},
		{"breakpoint_right = 101, breakpoint_right_sr = 1",
			";TSDLEN=0;TSITE=BLUNT;", "\t329.0\tPASS\t"},
		{"breakpoint_right = 106, breakpoint_right_sr = 1",
			";TSDLEN=-5;TSITE=DEL;", "\t329.0\tPASS\t"},
		{"breakpoint_right = 70, breakpoint_right_sr = 1",
			";TSDLEN=31;TSITE=DUP;", "\t329.0\tTSDOutlier\t"},
		{"breakpoint_right = NULL, breakpoint_right_sr = 0",
			";LBKP=100;LSR=1;TSITETR=A:101;POLYA=", "\t329.0\tOneBreakpoint\t"},
		{"breakpoint_left = NULL, breakpoint_left_sr = 0",
			";PG3P=ee1:150:3;DP=", "\t329.0\tNoBreakpoint\t"}
	};

	char sql[256] = {};
	int i = 0;

	for (i = 0; i < (int) (sizeof (cases) / sizeof (cases[0])); i++)
		{
			snprintf (sql, sizeof (sql),
					"UPDATE retrocopy SET %s WHERE id = 8", cases[i].update);
			db_exec (db, sql);

			vcf (db, vcf_file, &opt);

			ck_assert_msg (grep_vcf (vcf_file, cases[i].info),
					"'%s' not found", cases[i].info);
			ck_assert_msg (grep_vcf (vcf_file, cases[i].filter),
					"'%s' not found", cases[i].filter);

			// The TSD is a hallmark only for a duplication
			if (cases[i].hallmarks != NULL)
				ck_assert_msg (grep_vcf_record (vcf_file, "\tRTC000008\t",
							cases[i].hallmarks), "'%s' not found", cases[i].hallmarks);
		}

	// Without both breakpoints, there is no target site
	ck_assert (!grep_vcf (vcf_file, "TSITE=DEL"));

	db_close (db);
	xunlink (db_file);
	xunlink (vcf_file);
	xunlink (fasta_file);
}
END_TEST

START_TEST (test_vcf_tsite_repeat)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	sqlite3 *db = NULL;
	FILE *fp = NULL;
	int fd = 0;

	db = create_db (db_file);
	create_vcf (vcf_file);
	populate_db (db);

	/*
	* chrY (2200 bases) without tandem repeats, but two:
	* - CA x 10 from 95 to 114;
	* - ATTC x 5.25 from 1989 to 2009
	*/
	fd = xmkstemp (fasta_file);
	fp = xfdopen (fd, "w");
	xfprintf (fp, "%s",
		">chrY\n"
		"GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG\n"
		"CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTCACACACACACACACACACATCGGAC\n"
		"TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC\n"
		"GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAA\n"
		"CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT\n"
		"ACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAA\n"
		"GACAGATAGTGCACACGACCGGCGTCGGAGAAACTCTATTTGCCGCCTGACAAGTCAATG\n"
		"CGATCCGTAGGGGCAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAACGATT\n"
		"AACTGATAAATGAGCCCTTTATGACACGGGCATATGACTGGTTTACGATAGTATGTCCAA\n"
		"CGGCGAGCTTTACATTTGCTGTGAGAGGTACAGGGATTAGTGAGAAGCCGTGCGTATCAA\n"
		"TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGC\n"
		"TTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTCCATGAACTTA\n"
		"GCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCA\n"
		"AAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAG\n"
		"CCGGTGACTCCTAATGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCCATAA\n"
		"ATCTGAGCAACCAGCTGAAGCAGGCACGACAGTGCGACATTATATCACTGTGGTAGGTTA\n"
		"GCTTCATCTAATGTCCAACTAGCCGGCCAATTCGCATGATACCTCTCCATCTGACCCAAG\n"
		"ATTGTGCTTGTTCAATTCTTCTTAACGTGATAACAGAATCAAACCTGCCAGGCGGTCGTC\n"
		"GCGGACCTCGGTCGAAGTAGTGGTGCGGATCCAGGGGAACCGTTGACTCAAAAGGAGCTG\n"
		"CCGTCCACCTAACGTGAAGTTCCAAAATCCCAAACCTCTCGAGATATTTATCCAGCAAGG\n"
		"AGTGGCAACGCCCGCTGCTTTAATCGCTACCAAAACGCAAACAAAAGCATACCCAAAAGT\n"
		"ACACGGGTGAGGGAGGTGATATAGTACAGCTACGAAGTATCTGGCGCCTCAATAGGATTA\n"
		"TAGCGGTCTCTCAGGCTGCTTGCCGTCCGGCCCGGCCGCGACACTCCGGTGCAAGCTTAA\n"
		"TTCGTACGTACTTCCCATTGGATCTCGTTTATCGATTAAGCCCGATCTAGGTTCCTAGAG\n"
		"GTTAAATTGGACGTCTTCCCACTCCGTTGCTGCGTGTCTAGGCGGTTTAGCGTAAGCGAA\n"
		"CAGGACCCTGCCTCAGCTCATAAGTCCTTATTCTCTCACGTTGTGTTACGAAAGATTCAC\n"
		"TCGAGGTCGTGTGAGGGTTGGGCTAGCGGCAATTATGAAACTATCACATCACATAAGCGG\n"
		"GCTAGATATAATTTAATCTTAATCCATAAAACACTAGCTCAGCAGTTGAAAAAATGGCTA\n"
		"GGTTCCAGCTTTTGGGGAGACGTCTTTCTGAGGGTCAGCCGTGATTCCGATTCGATTAGA\n"
		"CTGGTCCCCACGGGTCCATGAGTACGAGGAAACTCGGTATCGAGCCTAAAAGTTATAAGG\n"
		"CATCTCGCCCAGGAAAGTAACGACGTATGGGTAGTTCTCCATCACCAGCTATAATGGCTA\n"
		"GCGCACTCTCGTTCCAGGGCGTAGTTACACTGAGCGTGCCATGTCAGCATGCTAGCGTAT\n"
		"CGCCCCCCAATGCCCCGCAATAGGGTAATTCGCCGACGAGTAAGCGTAGATTACACACCC\n"
		"AGGAAACGATTCATTCATTCATTCATTCACTTCATTATAGGTCGTGTAGCGCTAGACAGT\n"
		"CACCTTTAAAGGAAGAATCAGAGGCAAGATCTACGTGGCAGTCTCGTGTTGACGCCTTAG\n"
		"CCGGTGGCGAACAGTATTGACCTGGCCGATGCTAATATTCTGATTTGGGGTTGATTTGCG\n"
		"CTTCAGGCGCTAAAGTGGTTTTGAGTAACATGTCCTTTTG\n");
	xfclose (fp);

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file = fasta_file,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_sr = VCF_DEFAULT_MIN_SR,
		.max_tsd = VCF_DEFAULT_MAX_TSD
	};

	/*
	* The chrY retrocopy with its breakpoints moved. The
	* repeat must touch the site - from one breakpoint to
	* the other, or the single one - or VCF_TSITE_MARGIN
	* bases around it. A NULL 'tr' means no repeat
	*/
	struct
	{
		int left;
		int right;
		const char *tr;
	} cases[] = {
		// Duplication over the repeat
		{100, 91, ";TSITETR=CA:20"},
		// The margin touches the repeat, or not
		{92, 83, ";TSITETR=CA:20"},
		{91, 82, NULL},
		// A single breakpoint
		{110, 0, ";TSITETR=CA:20"},
		{0, 114, ";TSITETR=CA:20"},
		// Far from the repeats
		{300, 291, NULL},
		// A long deletion: A site at each breakpoint
		{300, 2000, ";TSITETR=ATTC:21"},
		// The flank goes beyond the end of the contig
		{2190, 0, NULL},
		// No breakpoint, no target site
		{0, 0, NULL}
	};

	char sql[256] = {};
	int i = 0;

	for (i = 0; i < (int) (sizeof (cases) / sizeof (cases[0])); i++)
		{
			snprintf (sql, sizeof (sql),
					"UPDATE retrocopy SET\n"
					"	breakpoint_left = %d, breakpoint_left_sr = %d,\n"
					"	breakpoint_right = %d, breakpoint_right_sr = %d\n"
					"WHERE id = 8",
					cases[i].left, cases[i].left > 0,
					cases[i].right, cases[i].right > 0);
			db_exec (db, sql);

			vcf (db, vcf_file, &opt);

			ck_assert_msg (grep_vcf_record (vcf_file, "\tRTC000008\t", "chrY\t"),
					"case %d: retrocopy not found", i);

			if (cases[i].tr != NULL)
				ck_assert_msg (grep_vcf_record (vcf_file, "\tRTC000008\t", cases[i].tr),
						"case %d: '%s' not found", i, cases[i].tr);
			else
				ck_assert_msg (!grep_vcf_record (vcf_file, "\tRTC000008\t", ";TSITETR="),
						"case %d: unexpected TSITETR", i);
		}

	db_close (db);
	xunlink (db_file);
	xunlink (vcf_file);
	xunlink (fasta_file);
}
END_TEST

START_TEST (test_vcf_no_genotype_fail)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char vcf_file[] = "/tmp/ponga.vcf.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	sqlite3 *db = NULL;

	db = create_db (db_file);
	create_vcf (vcf_file);
	create_fasta (fasta_file);
	populate_db (db);

	// As a merge-call that stopped before the genotypes
	db_exec (db, "DELETE FROM genotype");

	VCFOption opt = {
		.min_tsd = VCF_DEFAULT_MIN_TSD,
		.fasta_file = fasta_file,
		.min_qual = VCF_DEFAULT_MIN_QUAL,
		.min_mapq = VCF_DEFAULT_MIN_MAPQ,
		.min_sr = VCF_DEFAULT_MIN_SR,
		.max_tsd = VCF_DEFAULT_MAX_TSD
	};

	vcf (db, vcf_file, &opt);
}
END_TEST

Suite *
make_vcf_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("VCF");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_vcf);
	tcase_add_test (tc_core, test_vcf_edge);
	tcase_add_test (tc_core, test_vcf_parental_chosen);
	tcase_add_test (tc_core, test_vcf_primer);
	tcase_add_test (tc_core, test_vcf_intron_retention);
	tcase_add_test (tc_core, test_vcf_target_site);
	tcase_add_test (tc_core, test_vcf_tsite_repeat);
	suite_add_tcase (s, tc_core);

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	tcase_add_exit_test (tc_abort, test_vcf_no_genotype_fail, EXIT_FAILURE);
	suite_add_tcase (s, tc_abort);

	return s;
}
