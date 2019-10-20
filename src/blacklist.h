#ifndef BLACKLIST_H
#define BLACKLIST_H

#include "db.h"
#include "chr.h"
#include "set.h"
#include "hash.h"

struct _Blacklist
{
	Hash         *idx;
	sqlite3_stmt *blacklist_stmt;
	sqlite3_stmt *overlapping_blacklist_stmt;
	ChrStd       *cs;
	long          table_id;
};

typedef struct _Blacklist Blacklist;

Blacklist * blacklist_new (sqlite3_stmt *blacklist_stmt,
		sqlite3_stmt *overlapping_blacklist_stmt, ChrStd *cs);

void blacklist_free (Blacklist *blacklist);

void blacklist_index_dump (Blacklist *blacklist, const char *file,
		const char *feature, const char *attribute,
		const Set *attributes);

int blacklist_lookup (Blacklist *blacklist, const char *chr,
		long low, long high, long padding, const long cluster_id,
		const long cluster_sid);

#endif /* blacklist.h */
