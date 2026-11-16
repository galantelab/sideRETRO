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

#ifndef INTRON_H
#define INTRON_H

#include "hash.h"
#include "list.h"
#include "db.h"

/*
 * Intron of a parental gene: A gap into the union of all
 * exons of the gene, so a region that is exon at any
 * transcript is never intron. Coordinates are 1-based
 * inclusive. The exons at its left and right - by the
 * coordinates, whatever the strand of the gene - name it
 */
struct _Intron
{
	char         *gene_name;
	char         *chr;
	long          start;
	long          end;
	int           exon_id_left;
	int           exon_id_right;
};

typedef struct _Intron Intron;

/*
 * Index of the introns of the parental genes of each
 * retrocopy. It is built from the tables exon, cluster
 * and cluster_merging, and it is read only afterwards:
 * The lookup can be shared among threads
 */
struct _IntronTree
{
	Hash         *idx;
	List         *introns;
};

typedef struct _IntronTree IntronTree;

IntronTree * intron_tree_new (sqlite3 *db);

void intron_tree_free (IntronTree *intron_tree);

const Intron * intron_tree_lookup (const IntronTree *intron_tree,
		int retrocopy_id, const char *chr, long low, long high);

void intron_dump (sqlite3_stmt *intron_stmt, const Intron *intron,
		int retrocopy_id, int source_id, int reads);

#define intron_tree_size(t) (list_size ((t)->introns))

#endif /* intron.h */
