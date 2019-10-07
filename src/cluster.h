#ifndef CLUSTER_H
#define CLUSTER_H

#include "db.h"
#include "set.h"

void cluster (sqlite3_stmt *cluster_stmt, sqlite3_stmt *clustering_stmt,
		const long eps, const int min_pts, const Set *blacklist_chr,
		const int distance);

#endif /* cluster.h */
