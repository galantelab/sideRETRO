#ifndef RECLUSTER_H
#define RECLUSTER_H

#include "db.h"

void recluster (sqlite3_stmt *reclustering_stmt,
		const long eps, const int min_pts,
		const int distance);

#endif /* recluster.h */
