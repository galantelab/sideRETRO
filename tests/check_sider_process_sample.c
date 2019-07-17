#include "config.h"

#include <check.h>
#include "check_sider.h"
#include "../src/sider.h"

START_TEST (test_sider)
{
	ck_assert_int_eq (1, 1);
}
END_TEST


Suite *
make_sider_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("sideRETRO");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_sider);
	suite_add_tcase (s, tc_core);

	return s;
}
