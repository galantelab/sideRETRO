#ifndef DBSCAN_H
#define DBSCAN_H

#include "ibitree.h"
#include "list.h"

enum _Label
{
	UNDEFINED,
	NOISE,
	REACHABLE,
	CORE
};

typedef enum _Label Label;

struct _Point
{
	Label  label;
	int    id;
	int    neighbors;
	long   low;
	long   high;
	void  *data;
};

typedef struct _Point Point;

struct _DBSCAN
{
	List          *points;
	IBiTree       *index;
	DestroyNotify  destroy_data;
};

typedef struct _DBSCAN DBSCAN;

typedef void (*DFunc) (Point *p, void *user_data);

DBSCAN * dbscan_new          (DestroyNotify destroy_data);
void     dbscan_free         (DBSCAN *db);
void     dbscan_insert_point (DBSCAN *db, long low, long high, void *data);
int      dbscan_cluster      (DBSCAN *db, long eps, int min_pts, DFunc func, void *user_data);

#endif /* dbscan.h */
