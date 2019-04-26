#ifndef BITREE_H
#define BITREE_H

#include <stdlib.h>
#include "types.h"

struct _BiTreeNode
{
	void               *data;
	struct _BiTreeNode *left;
	struct _BiTreeNode *right;
};

typedef struct _BiTreeNode BiTreeNode;

struct _BiTree
{
	size_t         size;
	DestroyNotify  destroy_fun;
	BiTreeNode    *root;
};

typedef struct _BiTree BiTree;

BiTree * bitree_new       (DestroyNotify destroy_fun);
void     bitree_free      (BiTree *tree);
void     bitree_rem_left  (BiTree *tree, BiTreeNode *node);
void     bitree_rem_right (BiTree *tree, BiTreeNode *node);
void     bitree_ins_left  (BiTree *tree, BiTreeNode *node, const void *data);
void     bitree_ins_right (BiTree *tree, BiTreeNode *node, const void *data);

enum _BiTreeTraverse
{
	PREORDER = 0,
	INORDER,
	POSTORDER
};

typedef enum _BiTreeTraverse BiTreeTraverse;

void bitree_traverse (BiTreeTraverse traverse, BiTreeNode *node,
		Func func, void *user_data);

#define bitree_size(tree)    ((tree)->size)
#define bitree_root(tree)    ((tree)->root)
#define bitree_is_eob(node)  ((node) == NULL)
#define bitree_is_leaf(node) ((node)->left == NULL && (node)->right == NULL)
#define bitree_data(node)    ((node)->data)
#define bitree_left(node)    ((node)->left)
#define bitree_right(node)   ((node)->right)

#endif /* bitree.h */
