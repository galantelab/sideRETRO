#include "config.h"

#include <stdlib.h>
#include <check.h>
#include "check_sider.h"

#include "../src/utils.h"
#include "../src/wrapper.h"
#include "../src/gz.c"

static void
handle_sigabrt (int sig)
{
	if (sig == SIGABRT)
		exit (1);
}

static void
create_gz (const char *cnt, char *path)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, "%s", cnt);

	xfclose (fp);
}

static void
create_big_text_gz (char *path)
{
	FILE *fp = NULL;
	int fd, i, j;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	for  (i = 0; i < 10; i++)
		{
			for  (j = 0; j < 10000; j++)
				xfprintf (fp, "ponga");
			xfprintf (fp, "\n");
		}

	xfclose (fp);
}

static void
create_long_line_gz (char *path)
{
	FILE *fp = NULL;
	int fd, i;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	for  (i = 0; i < 10000; i++)
		xfprintf (fp, "ponga");

	xfclose (fp);
}

START_TEST (test_open_fatal)
{
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";
	gz_open_for_reading (gz_path);
}
END_TEST

START_TEST (test_close_fatal)
{
	GzFile *gz;
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";

	create_gz ("PONGA\n", gz_path);
	gz = gz_open_for_reading (gz_path);

	gz_close (gz);
	xunlink (gz_path);

	gz_close (NULL);
	gz_close (gz);
}
END_TEST

START_TEST (test_read_fatal1)
{
	GzFile *gz = NULL;
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";
	char *line = NULL;
	size_t n = 0;

	create_gz ("PONGA\n", gz_path);
	gz = gz_open_for_reading (gz_path);
	gzclose (gz->fp);

	while (gz_getline (gz, &line, &n))
		;

	xfree (line);
	gz_close (gz);
	xunlink (gz_path);
}
END_TEST

START_TEST (test_read_fatal2)
{
	GzFile *gz = NULL;
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";
	char *line = NULL;
	size_t n = 0;

	create_big_text_gz (gz_path);
	gz = gz_open_for_reading (gz_path);

	gz_getline (gz, &line, &n);
	gzclose (gz->fp);

	while (gz_getline (gz, &line, &n))
		;

	xfree (line);
	gz_close (gz);
	xunlink (gz_path);
}
END_TEST

START_TEST (test_read1)
{
	GzFile *gz = NULL;
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";
	char *line = NULL;
	size_t n = 0;

	create_big_text_gz (gz_path);
	gz = gz_open_for_reading (gz_path);

	while (gz_getline (gz, &line, &n))
		;

	xfree (line);
	gz_close (gz);
	xunlink (gz_path);
}
END_TEST

START_TEST (test_read2)
{
	GzFile *gz = NULL;
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";
	char *line = NULL;
	int i = 0;
	int l = 0;
	size_t n = 0;

	const char *gz_cnt_ex =
		"=> 1 ponga\n"
		"=> 2 ponga\n"
		"=> 3 ponga\n"
		"=> 4 ponga\n"
		"=> 5 ponga\n"
		"=> 6 ponga\n"
		"=> 7 ponga\n"
		"=> 8 ponga\n"
		"=> 9 ponga\n"
		"=> 10 ponga\n";

	create_gz (gz_cnt_ex, gz_path);
	gz = gz_open_for_reading (gz_path);

	while (gz_getline (gz, &line, &n))
		{
			sscanf (line, "%*s %d", &l);
			ck_assert_int_eq (l, ++i);
		}

	ck_assert_int_eq (gz_get_num_line (gz), i);
	ck_assert_str_eq (gz_get_filename (gz), gz_path);

	xfree (line);
	gz_close (gz);
	xunlink (gz_path);
}
END_TEST

START_TEST (test_read_long_line)
{
	GzFile *gz = NULL;
	char gz_path[] = "/tmp/ponga.txt.XXXXXX";
	char *line = NULL;
	size_t n = 0;

	create_long_line_gz (gz_path);
	gz = gz_open_for_reading (gz_path);

	while (gz_getline (gz, &line, &n))
		;

	xfree (line);
	gz_close (gz);
	xunlink (gz_path);
}
END_TEST

Suite *
make_gz_suite (void)
{
	setup_signal (SIGABRT, handle_sigabrt);

	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("GZ");

	/* Core test case */
	tc_core = tcase_create ("Core");

	/* Abort test case */
	tc_abort = tcase_create ("Abort");

	suite_add_tcase (s, tc_core);
	tcase_add_test (tc_core, test_read1);
	tcase_add_test (tc_core, test_read2);
	tcase_add_test (tc_core, test_read_long_line);

	suite_add_tcase (s, tc_abort);
	tcase_add_exit_test (tc_abort, test_open_fatal, 1);
	tcase_add_exit_test (tc_abort, test_close_fatal, 1);
	tcase_add_exit_test (tc_abort, test_read_fatal1, 1);
	tcase_add_exit_test (tc_abort, test_read_fatal2, 1);

	return s;
}
