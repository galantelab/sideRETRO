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

#include <assert.h>
#include "wrapper.h"
#include "bitree.h"

BiTree *
bitree_new (DestroyNotify destroy_fun)
{
	BiTree *tree = xcalloc (1, sizeof (BiTree));
	tree->destroy_fun = destroy_fun;
	return tree;
}

void
bitree_rem_left (BiTree *tree, BiTreeNode *node)
{
	if (tree == NULL
			|| (tree != NULL && tree->size == 0))
		return;

	BiTreeNode **position = node == NULL
		? &tree->root
		: &node->left;

	if (*position != NULL)
		{
			bitree_rem_left (tree, *position);
			bitree_rem_right (tree, *position);

			if (tree->destroy_fun != NULL)
				tree->destroy_fun ((*position)->data);

			xfree (*position);
			tree->size--;
		}
}

void
bitree_rem_right (BiTree *tree, BiTreeNode *node)
{
	if (tree == NULL
			|| (tree != NULL && tree->size == 0))
		return;

	BiTreeNode **position = node == NULL
		? &tree->root
		: &node->right;

	if (*position != NULL)
		{
			bitree_rem_left (tree, *position);
			bitree_rem_right (tree, *position);

			if (tree->destroy_fun != NULL)
				tree->destroy_fun ((*position)->data);

			xfree (*position);
			tree->size--;
		}
}

void
bitree_free (BiTree *tree)
{
	bitree_rem_left (tree, NULL);
	xfree (tree);
}

void
bitree_ins_left (BiTree *tree, BiTreeNode *node, const void *data)
{
	// Do not allow a NULL element unless the bitree is empty
	// Normally allow insertion only at the end of a branch
	assert ((node == NULL && tree->size == 0)
			|| (node != NULL && node->left == NULL && tree->size > 0));

	BiTreeNode *new_node = NULL;
	BiTreeNode **position = NULL;

	position = node == NULL
		? &tree->root
		: &node->left;

	new_node = xcalloc (1, sizeof (BiTreeNode));
	new_node->data = (void *) data;

	*position = new_node;
	tree->size++;
}

void
bitree_ins_right (BiTree *tree, BiTreeNode *node, const void *data)
{
	// Do not allow a NULL element unless the bitree is empty
	// Normally allow insertion only at the end of a branch
	assert ((node == NULL && tree->size == 0)
			|| (node != NULL && node->right == NULL && tree->size > 0));

	BiTreeNode *new_node = NULL;
	BiTreeNode **position = NULL;

	position = node == NULL
		? &tree->root
		: &node->right;

	new_node = xcalloc (1, sizeof (BiTreeNode));
	new_node->data = (void *) data;

	*position = new_node;
	tree->size++;
}

static void
bitree_preorder (BiTreeNode *node, Func func, void *user_data)
{
	if (!bitree_is_eob (node))
		{
			func (node->data, user_data);
			bitree_preorder (bitree_left (node), func, user_data);
			bitree_preorder (bitree_right (node), func, user_data);
		}
}

static void
bitree_inorder (BiTreeNode *node, Func func, void *user_data)
{
	if (!bitree_is_eob (node))
		{
			bitree_inorder (bitree_left (node), func, user_data);
			func (node->data, user_data);
			bitree_inorder (bitree_right (node), func, user_data);
		}
}

static void
bitree_postorder (BiTreeNode *node, Func func, void *user_data)
{
	if (!bitree_is_eob (node))
		{
			bitree_postorder (bitree_left (node), func, user_data);
			bitree_postorder (bitree_right (node), func, user_data);
			func (node->data, user_data);
		}
}

void
bitree_traverse (BiTreeTraverse traverse, BiTreeNode *node,
		Func func, void *user_data)
{
	assert (node != NULL && func != NULL);
	switch (traverse)
		{
		case PREORDER :
			{
				bitree_preorder (node, func, user_data);
				break;
			}
		case INORDER :
			{
				bitree_inorder (node, func, user_data);
				break;
			}
		case POSTORDER :
			{
				bitree_postorder (node, func, user_data);
				break;
			}
		}
}
