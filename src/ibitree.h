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

#ifndef IBITREE_H
#define IBITREE_H

#include "bitree.h"

typedef BiTree IBiTree;

struct _IBiTreeNode
{
	void   *data;
	long    low;
	long    high;
	long    max;
	int     height;
};

typedef struct _IBiTreeNode IBiTreeNode;

IBiTree * ibitree_new    (DestroyNotify destroy_fun);
void      ibitree_free   (IBiTree *tree);
void      ibitree_insert (IBiTree *tree, long low, long high, const void *data);

struct _IBiTreeLookupData
{
	void   *data;
	long    node_low;
	long    node_high;
	long    interval_low;
	long    interval_high;
	long    overlap_pos;
	long    overlap_len;
};

typedef struct _IBiTreeLookupData IBiTreeLookupData;

typedef void (*IBiTreeLookupFunc) (IBiTreeLookupData *ldata, void *user_data);

int ibitree_lookup (IBiTree *tree, long low, long high, float node_overlap_frac,
		float interval_overlap_frac, int either, IBiTreeLookupFunc func, void *user_data);

#define ibitree_size(tree) (bitree_size ((tree)))

#endif /* ibitree.h */
