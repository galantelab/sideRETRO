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
bitree_preorder (BiTreeNode *node, BiTFunc func, void *user_data)
{
	if (!bitree_is_eob (node))
		{
			func (node, user_data);
			if (!bitree_is_eob (bitree_left (node)))
				bitree_preorder (bitree_left (node), func, user_data);
			if (!bitree_is_eob (bitree_right (node)))
				bitree_preorder (bitree_right (node), func, user_data);
		}
}

static void
bitree_inorder (BiTreeNode *node, BiTFunc func, void *user_data)
{
	if (!bitree_is_eob (node))
		{
			if (!bitree_is_eob (bitree_left (node)))
				bitree_preorder (bitree_left (node), func, user_data);
			func (node, user_data);
			if (!bitree_is_eob (bitree_right (node)))
				bitree_preorder (bitree_right (node), func, user_data);
		}
}

static void
bitree_postorder (BiTreeNode *node, BiTFunc func, void *user_data)
{
	if (!bitree_is_eob (node))
		{
			if (!bitree_is_eob (bitree_left (node)))
				bitree_preorder (bitree_left (node), func, user_data);
			if (!bitree_is_eob (bitree_right (node)))
				bitree_preorder (bitree_right (node), func, user_data);
			func (node, user_data);
		}
}

void
bitree_traverse (BiTreeTraverse traverse, BiTreeNode *node,
		BiTFunc func, void *user_data)
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
