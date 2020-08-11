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

#include <check.h>
#include "check_sider.h"
#include "../src/chr.h"

START_TEST (test_chr_std_lookup)
{
	ChrStd *cs = chr_std_new ();
	const char *chr = NULL;

	const char *chrs[][2] = {
		{"10",     "chr10" },
		{"chrMT",  "chrM"  },
		{"CHr11",  "chr11" },
		{"Chrx",   "chrX"  },
		{"chr21",  "chr21" },
		{"ponga1", "ponga1"},
		{"",       ""      }
	};

	int chrs_size = 7;
	int i = 0;

	for (; i < chrs_size; i++)
		{
			chr = chr_std_lookup (cs, chrs[i][0]);
			ck_assert_str_eq (chr, chrs[i][1]);
		}

	chr_std_free (cs);
}
END_TEST

Suite *
make_chr_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Chr");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_chr_std_lookup);
	suite_add_tcase (s, tc_core);

	return s;
}
