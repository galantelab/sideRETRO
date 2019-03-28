#include "config.h"

#include <check.h>
#include "check_sider.h"

int
main (void)
{
	int number_failed = 0;
	SRunner *sr = NULL;

	sr = srunner_create (make_sider_suite ());
	srunner_add_suite (sr, make_list_suite ());
	srunner_add_suite (sr, make_hash_suite ());
	srunner_add_suite (sr, make_array_suite ());
	srunner_add_suite (sr, make_utils_suite ());
	/*srunner_set_tap (sr, "-");*/

	srunner_run_all (sr, CK_NORMAL);
	number_failed = srunner_ntests_failed (sr);
	srunner_free (sr);

	return number_failed;
}
