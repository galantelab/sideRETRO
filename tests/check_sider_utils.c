#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/utils.h"

START_TEST (test_chomp)
{
	char str[] = "PONGA";
	char str_with_new_line[] = "PONGA\n";

	chomp (str);
	chomp (str_with_new_line);

	ck_assert_str_eq (str, str_with_new_line);
	ck_assert_str_eq (str, "PONGA");
}
END_TEST

START_TEST (test_trim)
{
	char str_left_padding[] = "   PONGA";
	char str_right_padding[] = "PONGA   ";
	char str_left_right_padding[] = "   PONGA   ";

	trim (str_left_padding);
	trim (str_right_padding);
	trim (str_left_right_padding);

	ck_assert_str_eq (str_left_padding, str_right_padding);
	ck_assert_str_eq (str_left_padding, str_left_right_padding);
	ck_assert_str_eq (str_left_padding, "PONGA");
}
END_TEST

START_TEST (test_path_dir)
{
	char path[] = "/home/ponga/ponga.txt";
	char *path_copy = xstrdup (path);

	char *dir = path_dir (path);
	ck_assert_str_eq (path, path_copy);
	ck_assert_str_eq (dir, "/home/ponga");

	free (path_copy);
	free (dir);
}
END_TEST

START_TEST (test_path_file)
{
	char path[] = "/home/ponga/ponga.txt";
	char path_dot[] = "/home/ponga/.ponga.txt";
	char *path_copy = xstrdup (path);
	char *path_dot_copy = xstrdup (path_dot);

	char *file_w_ext = NULL;
	char *file_wo_ext = NULL;
	char *file_dot_w_ext = NULL;
	char *file_dot_wo_ext = NULL;
	char *file_no_change = NULL;

	file_w_ext = path_file (path, 0);
	ck_assert_str_eq (file_w_ext, "ponga.txt");
	ck_assert_str_eq (path, path_copy);

	file_wo_ext = path_file (path, 1);
	ck_assert_str_eq (file_wo_ext, "ponga");

	file_dot_w_ext = path_file (path_dot, 0);
	ck_assert_str_eq (file_dot_w_ext, ".ponga.txt");
	ck_assert_str_eq (path_dot, path_dot_copy);

	file_dot_wo_ext = path_file (path_dot, 1);
	ck_assert_str_eq (file_dot_wo_ext, ".ponga");

	file_no_change = path_file (file_dot_wo_ext, 1);
	ck_assert_str_eq (file_no_change, ".ponga");

	free (path_copy);
	free (path_dot_copy);
	free (file_w_ext);
	free (file_wo_ext);
	free (file_dot_w_ext);
	free (file_dot_wo_ext);
	free (file_no_change);
}
END_TEST

START_TEST (test_which)
{
	ck_assert_int_eq (which ("wc"), 1);
	ck_assert_int_eq (which ("pongadelapenha"), 0);
}
END_TEST

START_TEST (test_exists)
{
	char file[] = "/tmp/pongaXXXXXX";

	int fd = xmkstemp (file);

	ck_assert_int_eq (exists (file), 1);
	ck_assert_int_eq (exists ("pongatrongaflofa"), 0);

	xunlink (file);
}
END_TEST

Suite *
make_utils_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Utils");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_chomp);
	tcase_add_test (tc_core, test_trim);
	tcase_add_test (tc_core, test_path_dir);
	tcase_add_test (tc_core, test_path_file);
	tcase_add_test (tc_core, test_which);
	tcase_add_test (tc_core, test_exists);
	suite_add_tcase (s, tc_core);

	return s;
}
