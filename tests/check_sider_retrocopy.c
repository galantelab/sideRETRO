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
#include "../src/retrocopy.h"

static sqlite3 *
create_db (char *db_path)
{
	int fd;

	fd = xmkstemp (db_path);
	close (fd);

	return db_create (db_path);
}

static void
create_fasta (char *fasta_file)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (fasta_file);
	fp = xfdopen (fd, "w");

	// chr20 (gene8): The seed of the q22 clip at 1501-1520
	// chr21 (gene9): The reverse complement of the q23 clip
	// (after the poly-T) at the end of the exon, 2983-3000
	xfprintf (fp, "%s",
		">chr20\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"GTACGTACGTGTACGTACGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
		">chr21\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTGCTGCTGCTGCTGCTG\n"
	);

	xfclose (fp);
}

static void
populate_db (sqlite3 *db)
{
	// Database dump
	static const char schema[] =
		"BEGIN TRANSACTION;\n"
		"INSERT INTO exon VALUES (1,'gene1','chr1',1,3000,'+','eg1','ee1');\n"
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
		"INSERT INTO exon VALUES (12,'gene6','chr8',1,3000,'+','eg12','ee12');\n"
		"INSERT INTO exon VALUES (13,'gene7','chr9',1,3000,'+','eg13','ee13');\n"
		"INSERT INTO exon VALUES (14,'gene8','chr20',1,3000,'+','eg14','ee14');\n"
		"INSERT INTO exon VALUES (15,'gene9','chr21',1,3000,'+','eg15','ee15');\n"
		"INSERT INTO exon VALUES (16,'gene10','chr22',1000,1200,'+','eg16','ee16');\n"
		"INSERT INTO exon VALUES (17,'gene10','chr22',2000,2200,'+','eg16','ee17');\n"
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
		"INSERT INTO alignment VALUES (13,'q12',99,'chr15',200,20,'100M10S',110,100,'chr8',1,5,1,'*','*');\n"
		"INSERT INTO alignment VALUES (14,'q13',0x800,'chr15',290,20,'10H100M',100,100,'chr8',1,5,1,'*','*');\n"
		"INSERT INTO alignment VALUES (15,'q14',0x800,'chr15',290,20,'10H100M',100,100,'chr8',1,5,1,'*','*');\n"
		"INSERT INTO alignment VALUES (16,'q15',0x800,'chr15',210,20,'98M2S',100,98,'chr8',1,5,1,'*','*');\n"
		"INSERT INTO alignment VALUES (17,'q16',97,'chr16',300,20,'80M20S',100,80,'chr9',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (18,'q17',97,'chr16',300,20,'80M20S',100,80,'chr9',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (19,'q18',97,'chr16',350,20,'85M15S',100,85,'chr9',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (20,'q19',81,'chr16',420,20,'12S88M',100,88,'chr9',1,2,1,'*','*');\n"
		"INSERT INTO alignment VALUES (21,'q20',81,'chr17',500,20,'30S70M',100,70,'chr20',1,2,1,'GGCCTAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC','*');\n"
		"INSERT INTO alignment VALUES (22,'q21',81,'chr17',500,20,'30S70M',100,70,'chr20',1,2,1,'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC','*');\n"
		"INSERT INTO alignment VALUES (23,'q22',0x800,'chr17',420,20,'70M30S',100,70,'chr20',1,6,1,'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGTACGTACGTGTACGTACGTGTACGTACGT','*');\n"
		"INSERT INTO alignment VALUES (24,'q23',97,'chr18',600,20,'70M30S',100,70,'chr21',1,2,1,'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTCAGCAGCAGCAGCAGCAG','*');\n"
		"INSERT INTO alignment VALUES (25,'q24',97,'chr18',600,20,'70M30S',100,70,'chr21',1,2,1,'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT','*');\n"
		"INSERT INTO alignment VALUES (26,'q25',97,'chr22',1150,20,'50M50S',100,50,'chr19',100,14,1,'*','*');\n"
		"INSERT INTO alignment VALUES (27,'q25',2145,'chr22',2000,20,'50H50M',50,50,'chr19',100,14,1,'*','*');\n"
		"INSERT INTO alignment VALUES (28,'q25',145,'chr19',100,20,'100M',100,100,'chr22',1150,14,1,'*','*');\n"
		"INSERT INTO overlapping VALUES (16,26,1150,50);\n"
		"INSERT INTO overlapping VALUES (17,27,2000,50);\n"
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
		"INSERT INTO clustering VALUES (12,2,13,3,100);\n"
		"INSERT INTO clustering VALUES (12,2,14,3,100);\n"
		"INSERT INTO clustering VALUES (12,2,15,3,100);\n"
		"INSERT INTO clustering VALUES (12,2,16,3,100);\n"
		"INSERT INTO clustering VALUES (13,2,17,3,100);\n"
		"INSERT INTO clustering VALUES (13,2,18,3,100);\n"
		"INSERT INTO clustering VALUES (13,2,19,3,100);\n"
		"INSERT INTO clustering VALUES (13,2,20,3,100);\n"
		"INSERT INTO clustering VALUES (14,2,21,3,100);\n"
		"INSERT INTO clustering VALUES (14,2,22,3,100);\n"
		"INSERT INTO clustering VALUES (14,2,23,3,100);\n"
		"INSERT INTO clustering VALUES (15,2,24,3,100);\n"
		"INSERT INTO clustering VALUES (15,2,25,3,100);\n"
		"INSERT INTO clustering VALUES (16,2,28,3,100);\n"
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
		"INSERT INTO cluster VALUES (12,2,'chr15',200,400,'gene6',31);\n"
		"INSERT INTO cluster VALUES (13,2,'chr16',300,600,'gene7',31);\n"
		"INSERT INTO cluster VALUES (14,2,'chr17',400,600,'gene8',31);\n"
		"INSERT INTO cluster VALUES (15,2,'chr18',600,800,'gene9',31);\n"
		"INSERT INTO cluster VALUES (16,2,'chr19',100,300,'gene10',31);\n"
		"COMMIT;";

	db_exec (db, schema);
}

START_TEST (test_retrocopy)
{
	char db_file[] = "/tmp/ponga.db.XXXXXX";
	char fasta_file[] = "/tmp/ponga.fa.XXXXXX";
	const int near_dist = 3;

	sqlite3 *db = NULL;
	sqlite3_stmt *cluster_merging_stmt = NULL;
	sqlite3_stmt *retrocopy_stmt = NULL;
	sqlite3_stmt *junction_stmt = NULL;
	sqlite3_stmt *insertion_stmt = NULL;

	db = create_db (db_file);
	create_fasta (fasta_file);
	cluster_merging_stmt = db_prepare_cluster_merging_stmt (db);
	retrocopy_stmt = db_prepare_retrocopy_stmt (db);
	junction_stmt = db_prepare_junction_stmt (db);
	insertion_stmt = db_prepare_insertion_stmt (db);

	// Populate database
	populate_db (db);

	retrocopy (retrocopy_stmt,
			cluster_merging_stmt,
			junction_stmt,
			insertion_stmt,
			fasta_file,
			near_dist);

	/* TIME TO TEST */
	sqlite3_stmt *search_stmt = NULL;
	int i = 0;

	const char sql[] =
		"SELECT chr, insertion_point, insertion_point_type,\n"
		"	breakpoint_left, breakpoint_left_sr,\n"
		"	breakpoint_right, breakpoint_right_sr,\n"
		"	polya_len, polya_lower_bound,\n"
		"	COALESCE(polya_strand, '?'), polya_sr\n"
		"FROM retrocopy\n"
		"ORDER BY id ASC";

	// chr10: One clustered split read clipped at the right
	// chr11: No split read - window mean
	// chr12: One split read at each side - clean insertion
	// chr13/chr14: Lone clipped reads without supplementary
	// are noise - window mean
	// chr15: TSD of 10 bases. The 2S clip is ignored
	// chr16: Two clipped reads without supplementary agree
	// at the left. The lone ones at 434 and 420 are ignored
	// chr17: Poly-A at '+'. One read gives the exact length
	// (25) and the other just a lower bound (30)
	// chr18: Poly-A at '-' (run of T after the left flank)
	// chr19: Fragment whose read spans the exon-exon junction
	// of gene10 (exons 16 and 17)
	const char *chrs[] = {
		"chr10", "chr11", "chr12", "chr13", "chr13",
		"chr14", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19"
	};

	long ips[]        = {101, 250, 250, 150, 350, 250, 550, 290, 380, 500, 670, 200};
	int  ip_types[]   = {2,   1,   2,   1,   1,   1,   1,   2,   2,   2,   2,   1};
	long lefts[]      = {100, 0,   249, 0,   0,   0,   0,   299, 379, 489, 669, 0};
	int  left_srs[]   = {1,   0,   1,   0,   0,   0,   0,   1,   2,   1,   2,   0};
	long rights[]     = {0,   0,   250, 0,   0,   0,   0,   290, 0,   500, 0,   0};
	int  right_srs[]  = {0,   0,   1,   0,   0,   0,   0,   2,   0,   2,   0,   0};

	int  polya_lens[]   = {0, 0, 0, 0, 0, 0, 0, 0, 0, 25,  12, 0};
	int  polya_bounds[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,  0};
	const char *polya_strands[] = {"?","?","?","?","?","?","?","?","?","+","-","?"};
	int  polya_srs[]    = {0, 0, 0, 0, 0, 0, 0, 0, 0, 2,   2,  0};

	search_stmt = db_prepare (db, sql);

	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			ck_assert_int_lt (i, 12);
			ck_assert_str_eq (db_column_text (search_stmt, 0), chrs[i]);
			ck_assert_int_eq (db_column_int64 (search_stmt, 1), ips[i]);
			ck_assert_int_eq (db_column_int (search_stmt, 2), ip_types[i]);

			ck_assert_int_eq (db_column_int (search_stmt, 4), left_srs[i]);
			if (left_srs[i] > 0)
				ck_assert_int_eq (db_column_int64 (search_stmt, 3), lefts[i]);
			else
				ck_assert_int_eq (sqlite3_column_type (search_stmt, 3), SQLITE_NULL);

			ck_assert_int_eq (db_column_int (search_stmt, 6), right_srs[i]);
			if (right_srs[i] > 0)
				ck_assert_int_eq (db_column_int64 (search_stmt, 5), rights[i]);
			else
				ck_assert_int_eq (sqlite3_column_type (search_stmt, 5), SQLITE_NULL);

			ck_assert_int_eq (db_column_int (search_stmt, 10), polya_srs[i]);
			ck_assert_int_eq (db_column_int (search_stmt, 8), polya_bounds[i]);
			ck_assert_str_eq (db_column_text (search_stmt, 9), polya_strands[i]);
			if (polya_srs[i] > 0)
				ck_assert_int_eq (db_column_int (search_stmt, 7), polya_lens[i]);
			else
				ck_assert_int_eq (sqlite3_column_type (search_stmt, 7), SQLITE_NULL);
		}

	ck_assert_int_eq (i, 12);
	db_finalize (search_stmt);

	// Only the chr19 retrocopy has an exon-exon junction
	search_stmt = db_prepare (db,
			"SELECT retrocopy_id, exon_id1, exon_id2, reads FROM junction");

	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			ck_assert_int_eq (db_column_int (search_stmt, 0), 12);
			ck_assert_int_eq (db_column_int (search_stmt, 1), 16);
			ck_assert_int_eq (db_column_int (search_stmt, 2), 17);
			ck_assert_int_eq (db_column_int (search_stmt, 3), 1);
		}

	ck_assert_int_eq (i, 1);
	db_finalize (search_stmt);

	// Retrocopy ends into the parental gene:
	// chr17 (rid 10): Retrocopy at '+' of gene8 at '+'. The clip
	// of q22 (left junction) matches chr20:1501, so the 5' end is
	// at 1501. The 3' end is not observed (the clips at the right
	// junction are poly-A), so it is assumed at the exon end
	// (3000): 1500 exon bases + 25 poly-A
	// chr18 (rid 11): Retrocopy at '-' of gene9 at '+'. The clip
	// of q23 after the poly-T matches the reverse complement at
	// chr21:2983-3000, so it is the 3' end at 3000. No 5' end
	search_stmt = db_prepare (db,
			"SELECT retrocopy_id, five_prime_pos, five_prime_exon_id, five_prime_sr,\n"
			"	three_prime_pos, three_prime_exon_id, three_prime_sr, len, seq,\n"
			"	left_seq, right_seq\n"
			"FROM insertion ORDER BY retrocopy_id ASC");

	for (i = 0; db_step (search_stmt) == SQLITE_ROW; i++)
		{
			if (i == 0)
				{
					ck_assert_int_eq (db_column_int (search_stmt, 0), 10);
					ck_assert_int_eq (db_column_int64 (search_stmt, 1), 1501);
					ck_assert_int_eq (db_column_int (search_stmt, 2), 14);
					ck_assert_int_eq (db_column_int (search_stmt, 3), 1);
					ck_assert_int_eq (sqlite3_column_type (search_stmt, 4), SQLITE_NULL);
					ck_assert_int_eq (db_column_int (search_stmt, 6), 0);
					ck_assert_int_eq (db_column_int (search_stmt, 7), 1525);

					// Observed: The clips of q22 (left) and the consensus
					// of q20 and q21 (right): Ties keep the first read
					ck_assert_str_eq (db_column_text (search_stmt, 9),
							"GTACGTACGTGTACGTACGTGTACGTACGT");
					ck_assert_str_eq (db_column_text (search_stmt, 10),
							"GGCCTAAAAAAAAAAAAAAAAAAAAAAAAA");
					ck_assert_str_eq (db_column_text (search_stmt, 8),
							"GTACGTACGTGTACGTACGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAA");
				}
			else
				{
					ck_assert_int_eq (db_column_int (search_stmt, 0), 11);
					ck_assert_int_eq (sqlite3_column_type (search_stmt, 1), SQLITE_NULL);
					ck_assert_int_eq (db_column_int (search_stmt, 3), 0);
					ck_assert_int_eq (db_column_int64 (search_stmt, 4), 3000);
					ck_assert_int_eq (db_column_int (search_stmt, 5), 15);
					ck_assert_int_eq (db_column_int (search_stmt, 6), 1);
					ck_assert_int_eq (sqlite3_column_type (search_stmt, 7), SQLITE_NULL);
					ck_assert_str_eq (db_column_text (search_stmt, 9),
							"TTTTTTTTTTTTCAGCAGCAGCAGCAGCAG");
					ck_assert_int_eq (sqlite3_column_type (search_stmt, 10), SQLITE_NULL);
				}
		}

	ck_assert_int_eq (i, 2);
	db_finalize (search_stmt);

	db_finalize (cluster_merging_stmt);
	db_finalize (retrocopy_stmt);
	db_finalize (junction_stmt);
	db_finalize (insertion_stmt);
	db_close (db);

	xunlink (db_file);
	xunlink (fasta_file);
}
END_TEST

Suite *
make_retrocopy_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Retrocopy");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_retrocopy);
	suite_add_tcase (s, tc_core);

	return s;
}
