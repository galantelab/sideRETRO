#ifndef BLACKLIST_H
#define BLACKLIST_H

#include "db.h"
#include "chr.h"
#include "set.h"
#include "hash.h"
#include "gff.h"

struct _Blacklist
{
	Hash         *idx;
	sqlite3_stmt *blacklist_stmt;
	sqlite3_stmt *overlapping_blacklist_stmt;
	ChrStd       *cs;
};

typedef struct _Blacklist Blacklist;

Blacklist * blacklist_new (sqlite3_stmt *blacklist_stmt,
		sqlite3_stmt *overlapping_blacklist_stmt, ChrStd *cs);

void blacklist_free (Blacklist *blacklist);

void blacklist_index_dump_from_gff (Blacklist *blacklist,
		const char *gff_file, const GffFilter *filter);

void blacklist_index_dump_from_bed (Blacklist *blacklist,
		const char *bed_file);

int blacklist_lookup (Blacklist *blacklist, const char *chr,
		long low, long high, long padding, const long cluster_id,
		const long cluster_sid);

#endif /* blacklist.h */
