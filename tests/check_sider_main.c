#include "config.h"

#include <check.h>
#include "check_sider.h"

int
main (void)
{
	int number_failed = 0;
	SRunner *sr = NULL;

	sr = srunner_create (make_process_sample_suite ());
	srunner_add_suite (sr, make_list_suite ());
	srunner_add_suite (sr, make_hash_suite ());
	srunner_add_suite (sr, make_array_suite ());
	srunner_add_suite (sr, make_utils_suite ());
	srunner_add_suite (sr, make_sam_suite ());
	srunner_add_suite (sr, make_bwa_suite ());
	srunner_add_suite (sr, make_bitree_suite ());
	srunner_add_suite (sr, make_ibitree_suite ());
	srunner_add_suite (sr, make_str_suite ());
	srunner_add_suite (sr, make_db_suite ());
	srunner_add_suite (sr, make_chr_suite ());
	srunner_add_suite (sr, make_exon_suite ());
	srunner_add_suite (sr, make_abnormal_suite ());
	srunner_add_suite (sr, make_gff_suite ());
	srunner_add_suite (sr, make_io_suite ());
	srunner_add_suite (sr, make_dbscan_suite ());
	/*srunner_set_tap (sr, "-");*/

	srunner_run_all (sr, CK_NORMAL);
	number_failed = srunner_ntests_failed (sr);
	srunner_free (sr);

	return number_failed;
}
