#include "config.h"

#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/list.h"
#include "../src/ibitree.h"

static void
catcher (void *data, void *user_data)
{
	List *list = user_data;
	list_append (list, data);
}

START_TEST (test_ibitree)
{
	IBiTree *tree = ibitree_new (xfree);
	List *list = list_new (NULL);
	ListElmt *cur = NULL;
	int acm = 0;
	int i = 0;

	ibitree_insert (tree, 1, 10, xstrdup ("1-10"));
	ibitree_insert (tree, 9, 12, xstrdup ("9-12"));
	ibitree_insert (tree, 5, 10, xstrdup ("5-10"));
	ibitree_insert (tree, 1, 5, xstrdup ("1-5"));
	ibitree_insert (tree, 1, 1, xstrdup ("1-1"));

	ck_assert_int_eq (ibitree_size (tree), 5);

	acm = ibitree_lookup (tree, 6, 20, catcher, list);
	ck_assert_int_eq (acm, 3);
	ck_assert_int_eq (list_size (list), 3);

	cur = list_head (list);
	const char *overlap[] = {"1-10", "5-10", "9-12"};

	for (; cur != NULL; cur = list_next (cur), i++)
		ck_assert_str_eq (list_data (cur), overlap[i]);

	list_free (list);
	ibitree_free (tree);
}
END_TEST

Suite *
make_ibitree_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("IBiTree");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_ibitree);
	suite_add_tcase (s, tc_core);

	return s;
}
