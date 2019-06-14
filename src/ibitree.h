#ifndef IBITREE_H
#define IBITREE_H

#include "bitree.h"
#include "types.h"

struct _IBiTreeNode
{
	void   *data;
	long    low;
	long    high;
	long    max;
	int     height;
};

typedef struct _IBiTreeNode IBiTreeNode;

typedef BiTree IBiTree;

IBiTree * ibitree_new    (DestroyNotify destroy_fun);
void      ibitree_free   (IBiTree *tree);
void      ibitree_insert (IBiTree *tree, long low, long high, const void *data);
int       ibitree_lookup (IBiTree *tree, long low, long high, float node_overlap_frac,
		float interval_overlap_frac, int either, Func func, void *user_data);

#define ibitree_size(tree) (bitree_size ((tree)))

#endif /* ibitree.h */
