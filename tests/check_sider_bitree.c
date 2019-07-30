#include "config.h"

#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/list.h"
#include "../src/bitree.h"

static void
catcher (void *data, void *user_data)
{
	List *list = user_data;
	list_append (list, data);
}

START_TEST (test_bitree)
{
	BiTree *tree = bitree_new (xfree);
	List *list1 = list_new (NULL);
	List *list2 = list_new (NULL);
	List *list3 = list_new (NULL);
	ListElmt *cur = NULL;
	int i = 0;

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

	// Test preorder
	bitree_traverse (PREORDER, bitree_root (tree), catcher, list1);
	ck_assert_int_eq (list_size (list1), 7);

	cur = list_head (list1);
	const char *nodes_preorder[] = {"root", "left", "left_left",
		"left_right", "right", "right_left", "right_right"};

	for (i = 0; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), nodes_preorder[i]);

	// Test inorder
	bitree_traverse (INORDER, bitree_root (tree), catcher, list2);
	ck_assert_int_eq (list_size (list2), 7);

	cur = list_head (list2);
	const char *nodes_inorder[] = {"left_left", "left", "left_right",
		"root", "right_left", "right", "right_right"};

	for (i = 0; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), nodes_inorder[i]);

	// Test postorder
	bitree_traverse (POSTORDER, bitree_root (tree), catcher, list3);
	ck_assert_int_eq (list_size (list3), 7);

	cur = list_head (list3);
	const char *nodes_postorder[] = {"left_left", "left_right", "left",
		"right_left", "right_right", "right", "root"};

	for (i = 0; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), nodes_postorder[i]);

	list_free (list1);
	list_free (list2);
	list_free (list3);
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
