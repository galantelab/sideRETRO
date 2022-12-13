/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2023 Thiago L. A. Miller <tmiller@mochsl.org.br
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

#include "../src/cigar.c"

START_TEST (test_cigar_yylex)
{
	int n_cigar = 4;

	char *cigar_strs[] = {
		"35M",
		"10H100S1024M2I10M20S12H",
		"55M46S",
		"55S46M"
	};

	int n_ops[] = {
		1,
		7,
		2,
		2
	};

	uint32_t ops[][7] = {
		{'M', '0', '0', '0', '0', '0', '0'},
		{'H', 'S', 'M', 'I', 'M', 'S', 'H'},
		{'M', 'S', '0', '0', '0', '0', '0'},
		{'S', 'M', '0', '0', '0', '0', '0'}
	};

	uint32_t lens[][7] = {
		{35,   0,    0, 0,  0,  0,  0},
		{10, 100, 1024, 2, 10, 20, 12},
		{55,  46,    0, 0,  0,  0,  0},
		{55,  46,    0, 0,  0,  0,  0},
	};

	const char *cache = NULL;
	uint32_t val = 0;
	int token = 0;
	int i = 0;
	int j = 0;

	for (i = 0; i < n_cigar; i++)
		{
			j = 0;

			for (token = cigar_yylex (cigar_strs[i], &cache, &val);
					token;
					token = cigar_yylex (NULL, &cache, &val))
				{
					ck_assert_int_ne (token, CIGAR_YY_ERROR);

					if (token == CIGAR_YY_OP)
						ck_assert_uint_eq (val,  ops[i][j++]);
					else
						ck_assert_int_eq (val, lens[i][j]);
				}

			ck_assert_int_eq (j, n_ops[i]);
		}
}

Suite *
make_cigar_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("CIGAR");

	/* Core free test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_cigar_yylex);

	suite_add_tcase (s, tc_core);

	return s;
}
