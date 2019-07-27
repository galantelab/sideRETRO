#ifndef CLUSTER_H
#define CLUSTER_H

#include "db.h"

void cluster (sqlite3_stmt *clustering_stmt, long eps, int min_pts);

#endif /* cluster.h */
