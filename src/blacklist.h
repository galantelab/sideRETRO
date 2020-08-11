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
