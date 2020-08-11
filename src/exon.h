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

#ifndef EXON_H
#define EXON_H

#include "hash.h"
#include "chr.h"
#include "db.h"

struct _ExonTree
{
	Hash         *idx;
	Hash         *cache;
	sqlite3_stmt *exon_stmt;
	sqlite3_stmt *overlapping_stmt;
	ChrStd       *cs;
};

typedef struct _ExonTree ExonTree;

ExonTree * exon_tree_new (sqlite3_stmt *exon_stmt,
		sqlite3_stmt *overlapping_stmt, ChrStd *cs);

void exon_tree_free (ExonTree *exon_tree);

void exon_tree_index_dump (ExonTree *exon_tree, const char *gff_file);

int exon_tree_lookup_dump (ExonTree *exon_tree, const char *chr,
		long low, long high, float exon_overlap_frac,
		float alignment_overlap_frac, int either,
		long alignment_id);

#endif /* exon.h */
