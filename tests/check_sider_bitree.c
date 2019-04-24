#include "config.h"

#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/list.h"
#include "../src/bitree.h"

static void
bitree_list_append (BiTreeNode *node, void *user_data)
{
	List *list = user_data;
	list_append (list, bitree_data (node));
}

START_TEST (test_bitree)
{
	BiTree *tree = bitree_new (xfree);
	List *list = list_new (NULL);

	bitree_ins_left (tree, NULL, xstrdup ("root"));
	bitree_ins_left (tree, bitree_root (tree), xstrdup ("left"));
	bitree_ins_right (tree, bitree_root (tree), xstrdup ("right"));

	BiTreeNode *left_node = bitree_left (bitree_root (tree));
	bitree_ins_left (tree, left_node, xstrdup ("left_left"));
	bitree_ins_right (tree, left_node, xstrdup ("left_right"));

	BiTreeNode *right_node = bitree_right (bitree_root (tree));
	bitree_ins_left (tree, right_node, xstrdup ("right_left"));
	bitree_ins_right (tree, right_node, xstrdup ("right_right"));

	ck_assert_int_eq (bitree_size (tree), 7);

	bitree_traverse (PREORDER, bitree_root (tree),
			bitree_list_append, list);
	ck_assert_int_eq (list_size (list), 7);

	const char *nodes_preorder[] = {"root", "left", "left_left",
		"left_right", "right", "right_left", "right_right"};

	int i = 0;
	ListElmt *cur = list_head (list);

	for (; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), nodes_preorder[i]);

	list_free (list);
	bitree_free (tree);
}
END_TEST

Suite *
make_bitree_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("BiTree");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_bitree);
	suite_add_tcase (s, tc_core);

	return s;
}
