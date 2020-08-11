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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/io.h"

static char *
create_tmpfile (size_t nlines, const char * const *lines)
{
	FILE *fp = NULL;
	int fd;

	char *filename = strdup ("/tmp/tmpfile-XXXXXX");

	fd = mkstemp (filename);
	fp = xfdopen (fd, "w");

	for (size_t i = 0; i < nlines; i++)
		fprintf (fp, "%s\n", lines[i]);

	xfclose (fp);
	return filename;
}

START_TEST (test_read_file_lines)
{
	int i = 0;
	size_t nlines = 5;

	const char *lines[] = {
		"ponga0",
		"ponga1",
		"ponga2",
		"ponga3",
		"ponga4"
	};

	char *filename = create_tmpfile (nlines, lines);

	Array *arr = array_new (xfree);
	read_file_lines (arr, filename);

	ck_assert_int_eq (arr->len, nlines);

	for (i = 0; i < nlines; i++)
		ck_assert_str_eq ((char *) array_get (arr, i), lines[i]);

	array_free (arr, 1);
	xunlink (filename);
	xfree (filename);
}
END_TEST

Suite *
make_io_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("IO");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_read_file_lines);
	suite_add_tcase (s, tc_core);

	return s;
}
