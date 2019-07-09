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
