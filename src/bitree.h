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
