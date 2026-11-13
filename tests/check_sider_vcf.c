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
		"INSERT INTO alignment VALUES (2,'q4',97,'chr11',1,20,'110M',100,50,'chr11',1,8,1,'*','*');\n"
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
	// in each retrocopy line, in VCF order
	// (chr10, chr11, chr12, chr13, chr13, chr14, chr14, chrY)
	const char *present[][2] = {
		{";ORHO=-1.000000;POLARITY=+;PG=gene1;PGTYPE=1;PGEXONS=gene1:2;EEJ=5;EEJLIST=ee1-ee13;PG5P=ee1:1050:7;PG3P=ee13:3000:2;SVLEN=2000;LEFT_SVINSSEQ=GGTTAACC;RIGHT_SVINSSEQ=TTGGCC;", ";SR=2;LBKP=100;LSR=2;POLYA=25;POLYASR=6"},
		{";PGEXONS=gene2_1:1/gene2_2:1;",   ";DP="},
		{";SR=7;LBKP=259;LSR=3;RBKP=250;RSR=4;TSDLEN=10;TSD=CCCCCCCCCC", NULL},
		{";DP=",                            NULL},
		{";DP=",                            NULL},
		{";DP=",                            NULL},
		{";DP=",                            NULL},
		{";ORHO=1.000000;POLARITY=-;PG=gene1;PGTYPE=1;PGEXONS=gene1:2;PG3P=ee1:150:3;DP=", ";SR=2;LBKP=100;LSR=1;RBKP=102;RSR=1;TSDLEN=-1;POLYA=30;POLYAMIN;POLYASR=2"}
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
	const char *rnames[] = {
		"RTC000001\tPONGA:q1\n",
		"RTC000002\tPONGA:q4,q5\n",
		"RTC000003\tPONGA:q2,q3\n",
		"RTC000004\tPONGA:q6\n",
		"RTC000005\tPONGA:q7\n",
		"RTC000006\tPONGA:q8,q9\n",
		"RTC000007\tPONGA:q10,q11\n",
		"RTC000008\tPONGA:q1\n"
	};

	fp = xfopen (rnames_file, "r");

	for (i = 0; getline (&line, &len, fp) != EOF; i++)
		{
			ck_assert_int_lt (i, 8);
			ck_assert_str_eq (line, rnames[i]);
		}

	ck_assert_int_eq (i, 8);
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

	xfree (line);

	db_close (db);
	xunlink (db_file);
	xunlink (fasta_file);
	xunlink (vcf_file);
	xunlink (rnames_file);
	xunlink (insertion_file);
}
END_TEST

Suite *
make_vcf_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("VCF");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_vcf);
	suite_add_tcase (s, tc_core);

	return s;
}
