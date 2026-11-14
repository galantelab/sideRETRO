/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br>
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

#include <string.h>
#include "../src/wrapper.h"
#include "../src/cigar.h"

START_TEST (test_cigar_yylex)
{
	int n_cigar = 9;

	char *cigar_strs[] = {
		"35M",
		"10H100S1024M",
		"55M46S",
		"10M20I",
		"10M5D20M",
		"100M10P",
		"25=10X",
		"101M5N",
		"66M34B"
	};

	int n_ops[] = {
		1,
		3,
		2,
		2,
		3,
		2,
		2,
		2,
		2
	};

	int ops[][3] = {
		{CIGAR_MATCH,     0,                         0},
		{CIGAR_HARD_CLIP, CIGAR_SOFT_CLIP, CIGAR_MATCH},
		{CIGAR_MATCH,     CIGAR_SOFT_CLIP,           0},
		{CIGAR_MATCH,     CIGAR_INS,                 0},
		{CIGAR_MATCH,     CIGAR_DEL,       CIGAR_MATCH},
		{CIGAR_MATCH,     CIGAR_PAD,                 0},
		{CIGAR_EQUAL,     CIGAR_DIFF,                0},
		{CIGAR_MATCH,     CIGAR_REF_SKIP,            0},
		{CIGAR_MATCH,     CIGAR_BACK,                0}
	};

	uint32_t lens[][3] = {
		{35,   0,    0},
		{10, 100, 1024},
		{55,  46,    0},
		{10,  20,    0},
		{10,   5,   20},
		{100, 10,    0},
		{25,  10,    0},
		{101,  5,    0},
		{66,  34,    0}
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
					ck_assert_int_ne (token, CIGAR_ERROR);

					if (token == CIGAR_LENGTH)
						ck_assert_uint_eq (val, lens[i][j]);
					else
						ck_assert_int_eq (token,  ops[i][j++]);
				}

			ck_assert_int_eq (j, n_ops[i]);
		}
}

START_TEST (test_cigar_yyparser)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;

	cigar = cigar_yyparser ("10H5S80M3I10M2D5M1S", &n_cigar);

	ck_assert_ptr_ne (cigar, NULL);
	ck_assert_uint_eq (n_cigar, 8);

	ck_assert_int_eq (cigar_op (cigar[0]), CIGAR_HARD_CLIP);
	ck_assert_uint_eq (cigar_oplen (cigar[0]), 10);
	ck_assert_int_eq (cigar_op (cigar[1]), CIGAR_SOFT_CLIP);
	ck_assert_uint_eq (cigar_oplen (cigar[1]), 5);
	ck_assert_int_eq (cigar_op (cigar[2]), CIGAR_MATCH);
	ck_assert_uint_eq (cigar_oplen (cigar[2]), 80);
	ck_assert_int_eq (cigar_op (cigar[5]), CIGAR_DEL);
	ck_assert_int_eq (cigar_op (cigar[7]), CIGAR_SOFT_CLIP);
	ck_assert_uint_eq (cigar_oplen (cigar[7]), 1);

	// Same encoding of BAM: 'M' == 0
	ck_assert_uint_eq (cigar[2], (80 << 4) | 0);

	xfree (cigar);
}
END_TEST

START_TEST (test_cigar_yyparser_unavailable)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 1;

	cigar = cigar_yyparser ("*", &n_cigar);

	ck_assert_ptr_ne (cigar, NULL);
	ck_assert_uint_eq (n_cigar, 0);

	xfree (cigar);
}
END_TEST

START_TEST (test_cigar_yyparser_fail)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	int i = 0;

	const char *cigar_strs[] = {
		"",
		"M",
		"10M5",
		"10M5Z",
		"10M5H10M",
		"5S10M5S10M",
		"10M5S5H5M",
		"268435456M",
		"99999999999M"
	};

	for (i = 0; i < sizeof (cigar_strs) / sizeof (char *); i++)
		{
			cigar = cigar_yyparser (cigar_strs[i], &n_cigar);
			ck_assert_msg (cigar == NULL, "cigar '%s' should fail", cigar_strs[i]);
		}

	// The max length is ok
	cigar = cigar_yyparser ("268435455M", &n_cigar);
	ck_assert_ptr_ne (cigar, NULL);
	ck_assert_uint_eq (cigar_oplen (cigar[0]), 268435455);
	xfree (cigar);
}
END_TEST

START_TEST (test_cigar_iter)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	CigarIter iter = {};
	CigarOp op = {};
	int i = 0;

	// 10H 5S 10M 2I 5M 3D 4M 100N 6M 3S 7H
	cigar = cigar_yyparser ("10H5S10M2I5M3D4M100N6M3S7H", &n_cigar);
	ck_assert_ptr_ne (cigar, NULL);

	CigarTokenType types[] = {
		CIGAR_HARD_CLIP, CIGAR_SOFT_CLIP, CIGAR_MATCH,
		CIGAR_INS, CIGAR_MATCH, CIGAR_DEL, CIGAR_MATCH,
		CIGAR_REF_SKIP, CIGAR_MATCH, CIGAR_SOFT_CLIP,
		CIGAR_HARD_CLIP
	};

	int  qpos[] = {0,    0,    5,   15,   17,   22,   22,   26,   26,   32,   35};
	int  qlen[] = {0,    5,   10,    2,    5,    0,    4,    0,    6,    3,    0};
	long rpos[] = {1000, 1000, 1000, 1010, 1010, 1015, 1018, 1022, 1122, 1128, 1128};
	long rlen[] = {0,    0,   10,    0,    5,    3,    4,  100,    6,    0,    0};

	cigar_iter_init (&iter, cigar, n_cigar, 1000);

	while (cigar_iter_next (&iter, &op))
		{
			ck_assert_int_lt (i, n_cigar);
			ck_assert_int_eq (op.type, types[i]);
			ck_assert_int_eq (op.qpos, qpos[i]);
			ck_assert_int_eq (op.qlen, qlen[i]);
			ck_assert_int_eq (op.rpos, rpos[i]);
			ck_assert_int_eq (op.rlen, rlen[i]);
			i++;
		}

	ck_assert_int_eq (i, n_cigar);

	// Empty cigar
	cigar_iter_init (&iter, NULL, 0, 1000);
	ck_assert_int_eq (cigar_iter_next (&iter, &op), 0);

	xfree (cigar);
}
END_TEST

START_TEST (test_cigar_op_seq)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	CigarIter iter = {};
	CigarOp op = {};
	char *slice = NULL;
	int i = 0;

	const char *seq = "AAAAACCCCCGGTTTTTA";

	const char *slices[] = {
		NULL,
		"AAAAA",
		"CCCCC",
		"GG",
		NULL,
		"TTTTT",
		"A"
	};

	cigar = cigar_yyparser ("3H5S5M2I3D5M1S", &n_cigar);
	ck_assert_ptr_ne (cigar, NULL);

	cigar_iter_init (&iter, cigar, n_cigar, 1);

	while (cigar_iter_next (&iter, &op))
		{
			slice = cigar_op_seq (&op, seq);

			if (slices[i] == NULL)
				ck_assert_ptr_eq (slice, NULL);
			else
				ck_assert_str_eq (slice, slices[i]);

			xfree (slice);
			i++;
		}

	// Sequence shorter than the cigar
	cigar_iter_init (&iter, cigar, n_cigar, 1);
	cigar_iter_next (&iter, &op);
	cigar_iter_next (&iter, &op);
	ck_assert_ptr_eq (cigar_op_seq (&op, "AAA"), NULL);

	xfree (cigar);
}
END_TEST

START_TEST (test_cigar_soft_clip)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	uint32_t left = 0;
	uint32_t right = 0;
	int i = 0;

	const char *cigar_strs[] = {
		"100M",
		"20S80M",
		"80M20S",
		"10H20S60M10S5H",
		"15S",
		"5S40M5I40M10S",
		"40M2D40M",
		"10S80M5H"
	};

	// The clip at the end stops at the first operation
	// that is not a clip (I, D, M)
	uint32_t lefts[]  = {0, 20,  0, 20, 15, 5, 0, 10};
	uint32_t rights[] = {0,  0, 20, 10,  0, 10, 0,  0};

	for (i = 0; i < sizeof (cigar_strs) / sizeof (char *); i++)
		{
			cigar = cigar_yyparser (cigar_strs[i], &n_cigar);
			ck_assert_ptr_ne (cigar, NULL);

			cigar_soft_clip (cigar, n_cigar, &left, &right);

			ck_assert_uint_eq (left, lefts[i]);
			ck_assert_uint_eq (right, rights[i]);

			xfree (cigar);
		}
}
END_TEST

START_TEST (test_cigar_clip)
{
	uint32_t *cigar = NULL;
	uint32_t n_cigar = 0;
	uint32_t left = 0;
	uint32_t right = 0;
	int i = 0;

	const char *cigar_strs[] = {
		"100M",
		"20S80M",
		"80M20H",
		"10H20S60M10S5H",
		"15H",
		"5S40M5I40M10S",
		"40M2D40M",
		"10S80M5H"
	};

	// The clip at the end stops at the first operation
	// that is not a clip (I, D, M)
	uint32_t lefts[]  = {0, 20,  0, 30, 15, 5, 0, 10};
	uint32_t rights[] = {0,  0, 20, 15,  0, 10, 0,  5};

	for (i = 0; i < sizeof (cigar_strs) / sizeof (char *); i++)
		{
			cigar = cigar_yyparser (cigar_strs[i], &n_cigar);
			ck_assert_ptr_ne (cigar, NULL);

			cigar_clip (cigar, n_cigar, &left, &right);

			ck_assert_uint_eq (left, lefts[i]);
			ck_assert_uint_eq (right, rights[i]);

			xfree (cigar);
		}
}
END_TEST

Suite *
make_cigar_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("CIGAR");

	/* Core free test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_cigar_yylex);
	tcase_add_test (tc_core, test_cigar_yyparser);
	tcase_add_test (tc_core, test_cigar_yyparser_unavailable);
	tcase_add_test (tc_core, test_cigar_yyparser_fail);
	tcase_add_test (tc_core, test_cigar_iter);
	tcase_add_test (tc_core, test_cigar_op_seq);
	tcase_add_test (tc_core, test_cigar_soft_clip);
	tcase_add_test (tc_core, test_cigar_clip);

	suite_add_tcase (s, tc_core);

	return s;
}
