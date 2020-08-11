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

#include "config.h"

#include "set.h"
#include "ibitree.h"
#include "wrapper.h"
#include "dbscan.h"

void
dbscan_free (DBSCAN *db)
{
	if (db == NULL)
		return;

	Point *p = NULL;
	ListElmt *cur = list_head (db->points);

	for (; cur != NULL; cur = list_next (cur))
		{
			p = list_data (cur);
			if (db->destroy_data != NULL)
				db->destroy_data (p->data);
		}

	list_free (db->points);
	ibitree_free (db->index);
	xfree (db);
}

DBSCAN *
dbscan_new (DestroyNotify destroy_data)
{
	DBSCAN *db = xcalloc (1, sizeof (DBSCAN));

	db->points = list_new ((DestroyNotify) xfree);
	db->index = ibitree_new (NULL);
	db->destroy_data = destroy_data;

	return db;
}

void
dbscan_insert_point (DBSCAN *db,
		long low, long high, void *data)
{
	Point *p = xcalloc (1, sizeof (Point));

	*p = (Point) {
		.low   = low,
		.high  = high,
		.data  = data
	};

	list_append (db->points, p);
	ibitree_insert (db->index, low, high, p);
}

static void
get_point (IBiTreeLookupData *ldata,
		void *user_data)
{
	Point *p = ldata->data;
	List *neighbors = user_data;
	list_append (neighbors, p);
}

static int
range_query (DBSCAN *db, List *neighbors,
		Point *q, long eps)
{
	long center = (q->high + q->low) / 2;
	long low = center - eps;
	return ibitree_lookup (db->index, low > 0 ? low : 1,
			center + eps, -1, -1, 0, get_point, neighbors);
}

static void
process_seed (const Set *seed, DFunc func, void *user_data)
{
	ListElmt *cur = list_head (set_list (seed));
	Point *p = NULL;

	for (; cur != NULL; cur = list_next (cur))
		{
			p = list_data (cur);
			func (p, user_data);
		}
}

static void
set_union_list (Set *set, const List *list)
{
	ListElmt *cur = list_head (list);
	for (; cur != NULL; cur = list_next (cur))
		set_insert (set, list_data (cur));
}

static void
reset_points (List *list)
{
	ListElmt *cur = list_head (list);
	Point *p = NULL;

	for (; cur != NULL; cur = list_next (cur))
		{
			p = list_data (cur);
			p->label = UNDEFINED;
			p->id = 0;
		}
}

int
dbscan_cluster (DBSCAN *db, long eps, int min_pts, DFunc func, void *user_data)
{
	ListElmt *cur = NULL;
	ListElmt *cur_seed = NULL;
	List *neighbors = NULL;
	Set *seed = NULL;
	Point *p = NULL;
	Point *q = NULL;
	int acm = 0;

	int c = 0;
	int n = 0;

	// Set all points to UNDEFINED
	reset_points (db->points);

	cur = list_head (db->points);
	for (; cur != NULL; cur = list_next (cur))
		{
			p = list_data (cur);

			if (p->label != UNDEFINED)
				continue;

			neighbors = list_new (NULL);
			n = range_query (db, neighbors, p, eps);
			p->neighbors = n;

			if (n < min_pts)
				{
					p->label = NOISE;
					list_free (neighbors);
					continue;
				}

			p->label = CORE;
			p->id = ++c;

			seed = set_new (NULL);

			set_union_list (seed, neighbors);
			list_free (neighbors);

			cur_seed = list_head (set_list (seed));

			for (; cur_seed != NULL; cur_seed = list_next (cur_seed))
				{
					q = list_data (cur_seed);
					q->id = c;

					if (q->label == NOISE)
						q->label = REACHABLE;

					if (q->label != UNDEFINED)
						continue;

					q->label = REACHABLE;

					neighbors = list_new (NULL);
					n = range_query (db, neighbors, q, eps);
					q->neighbors = n;

					if (n >= min_pts)
						{
							q->label = CORE;
							set_union_list (seed, neighbors);
						}

					list_free (neighbors);
				}

			process_seed (seed, func, user_data);
			acm++;

			set_free (seed);
		}

	return acm;
}
