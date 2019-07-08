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
