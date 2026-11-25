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
#include <string.h>
#include "check_sider.h"

#include "../src/tandem.h"

static int
find (const char *seq, long from, long to, TandemRepeat *tr)
{
	return tandem_repeat_find (seq, strlen (seq), from, to, tr);
}

START_TEST (test_tandem_dinucleotide)
{
	TandemRepeat tr = {};

	// AC x 8 (16 bases) at offset 4
	const char *seq = "GTTGACACACACACACACACGTTGAGT";

	ck_assert_int_eq (find (seq, 10, 10, &tr), 1);
	ck_assert_int_eq (tr.start, 4);
	ck_assert_int_eq (tr.span, 16);
	ck_assert_int_eq (tr.unit_len, 2);
	ck_assert_str_eq (tr.unit, "AC");

	// It must overlap the range: Before and after it
	ck_assert_int_eq (find (seq, 0, 3, &tr), 0);
	ck_assert_int_eq (find (seq, 20, 26, &tr), 0);

	// The range may cover it only by its edges
	ck_assert_int_eq (find (seq, 0, 4, &tr), 1);
	ck_assert_int_eq (find (seq, 19, 26, &tr), 1);
}
END_TEST

START_TEST (test_tandem_homopolymer)
{
	TandemRepeat tr = {};

	// A run of 12 'A': The unit is 'A', not 'AA' or 'AAA'
	const char *seq = "GCTGAAAAAAAAAAAACGTC";

	ck_assert_int_eq (find (seq, 0, 19, &tr), 1);
	ck_assert_int_eq (tr.start, 4);
	ck_assert_int_eq (tr.span, 12);
	ck_assert_int_eq (tr.unit_len, 1);
	ck_assert_str_eq (tr.unit, "A");
}
END_TEST

START_TEST (test_tandem_min_copies)
{
	TandemRepeat tr = {};

	// TTTTC with 2.6 copies (13 bases): Long enough, but
	// too few copies. Its runs of 'T' are short too
	ck_assert_int_eq (find ("GACGTTTTCTTTTCTTTGCAG", 0, 20, &tr), 0);

	// GCC with 3.67 copies (11 bases): A repeat
	ck_assert_int_eq (find ("ATGGCCGCCGCCGCATGT", 0, 17, &tr), 1);
	ck_assert_int_eq (tr.span, 11);
	ck_assert_str_eq (tr.unit, "GCC");
}
END_TEST

START_TEST (test_tandem_min_span)
{
	TandemRepeat tr = {};

	// TA x 4.5 (9 bases): Enough copies, too short
	ck_assert_int_eq (find ("GCGTATATATATGCG", 0, 14, &tr), 0);

	// TA x 5 (10 bases): A repeat
	ck_assert_int_eq (find ("GCGTATATATATAGCG", 0, 15, &tr), 1);
	ck_assert_int_eq (tr.span, 10);
	ck_assert_str_eq (tr.unit, "TA");
}
END_TEST

START_TEST (test_tandem_longest)
{
	TandemRepeat tr = {};

	// TTCA x 5 (20 bases) and a run of 10 'G': The longest
	// wins, whatever its unit
	const char *seq = "GGGGGGGGGGTTCATTCATTCATTCATTCAC";

	ck_assert_int_eq (find (seq, 0, 30, &tr), 1);
	ck_assert_int_eq (tr.start, 10);
	ck_assert_int_eq (tr.span, 20);
	ck_assert_str_eq (tr.unit, "TTCA");

	// Only the run of 'G' overlaps
	ck_assert_int_eq (find (seq, 0, 5, &tr), 1);
	ck_assert_int_eq (tr.span, 10);
	ck_assert_str_eq (tr.unit, "G");

	// A hexanucleotide: AAGGGT x 3 (18 bases)
	ck_assert_int_eq (find ("CAAGGGTAAGGGTAAGGGTC", 0, 19, &tr), 1);
	ck_assert_int_eq (tr.span, 18);
	ck_assert_str_eq (tr.unit, "AAGGGT");
}
END_TEST

START_TEST (test_tandem_edges)
{
	TandemRepeat tr = {};

	// At the start and at the end of the sequence
	ck_assert_int_eq (find ("CACACACACACAGTTG", 0, 15, &tr), 1);
	ck_assert_int_eq (tr.start, 0);
	ck_assert_int_eq (tr.span, 12);

	ck_assert_int_eq (find ("GTTGCACACACACACA", 0, 15, &tr), 1);
	ck_assert_int_eq (tr.start, 4);
	ck_assert_int_eq (tr.span, 12);

	// A sequence shorter than any repeat, and an empty one
	ck_assert_int_eq (find ("ACGT", 0, 3, &tr), 0);
	ck_assert_int_eq (find ("", 0, 0, &tr), 0);
}
END_TEST

START_TEST (test_tandem_unknown)
{
	TandemRepeat tr = {};

	// Unknown bases are not a repeat, and break one
	ck_assert_int_eq (find ("ACNNNNNNNNNNNNACGT", 0, 17, &tr), 0);
	ck_assert_int_eq (find ("GCACACANCACACAGT", 0, 15, &tr), 0);

	// Neither without repeats at all
	ck_assert_int_eq (find ("ACGTTGCAAGTCCTAGGATC", 0, 19, &tr), 0);
}
END_TEST

Suite *
make_tandem_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("Tandem");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_tandem_dinucleotide);
	tcase_add_test (tc_core, test_tandem_homopolymer);
	tcase_add_test (tc_core, test_tandem_min_copies);
	tcase_add_test (tc_core, test_tandem_min_span);
	tcase_add_test (tc_core, test_tandem_longest);
	tcase_add_test (tc_core, test_tandem_edges);
	tcase_add_test (tc_core, test_tandem_unknown);
	suite_add_tcase (s, tc_core);

	return s;
}
