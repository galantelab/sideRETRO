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

#ifndef CLUSTER_H
#define CLUSTER_H

#include "db.h"
#include "set.h"
#include "blacklist.h"

enum _ClusterFilter
{
	CLUSTER_FILTER_NONE    = 1,
	CLUSTER_FILTER_CHR     = 2,
	CLUSTER_FILTER_DIST    = 4,
	CLUSTER_FILTER_REGION  = 8,
	CLUSTER_FILTER_SUPPORT = 16,
};

typedef enum _ClusterFilter ClusterFilter;

int cluster (sqlite3_stmt *cluster_stmt, sqlite3_stmt *clustering_stmt,
		const long eps, const int min_pts, const int distance,
		const int support, Set *blacklist_chr, Blacklist *blacklist,
		const int padding);

#endif /* cluster.h */
