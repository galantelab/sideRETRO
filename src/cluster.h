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
		const int support, Set *blacklist_chr,
		Blacklist *blacklist);

#endif /* cluster.h */
