#ifndef EXON_H
#define EXON_H

#include "hash.h"
#include "chr.h"
#include "db.h"

struct _ExonTree
{
	Hash         *idx;
	Hash         *cache;
	sqlite3      *db;
	sqlite3_stmt *exon_stmt;
	sqlite3_stmt *overlapping_stmt;
	ChrStd       *cs;
	long          alignment_id;
};

typedef struct _ExonTree ExonTree;

ExonTree * exon_tree_new (sqlite3 *db, sqlite3_stmt *exon_stmt,
		sqlite3_stmt *overlapping_stmt, ChrStd *cs);

void exon_tree_free (ExonTree *exon_tree);

void exon_tree_index_dump (ExonTree *exon_tree, const char *gff_file);

int exon_tree_lookup_dump (ExonTree *exon_tree, const char *chr,
		long low, long high, float exon_overlap_frac,
		float alignment_overlap_frac, int either,
		long alignment_id);

#endif /* exon.h */
