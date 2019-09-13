#ifndef CLUSTER_H
#define CLUSTER_H

#include "db.h"
#include "set.h"

void cluster (sqlite3_stmt *clustering_stmt, long eps, int min_pts, Set *blacklist_chr);

#endif /* cluster.h */
