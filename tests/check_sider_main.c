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

#include <stdlib.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/log.h"

#define LOG_DEBUG_KEY "LOG_DEBUG"

int
main (void)
{
	const char *log_debug = NULL;
	int number_failed = 0;
	SRunner *sr = NULL;

	log_set_quiet (1);

	log_debug = secure_getenv (LOG_DEBUG_KEY);
	if (log_debug != NULL && atoi (log_debug))
		log_set_quiet (0);

	sr = srunner_create (make_process_sample_suite ());
	srunner_add_suite (sr, make_list_suite ());
	srunner_add_suite (sr, make_hash_suite ());
	srunner_add_suite (sr, make_array_suite ());
	srunner_add_suite (sr, make_utils_suite ());
	srunner_add_suite (sr, make_sam_suite ());
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
	srunner_add_suite (sr, make_cluster_suite ());
	srunner_add_suite (sr, make_wrapper_suite ());
	srunner_add_suite (sr, make_db_merge_suite ());
	srunner_add_suite (sr, make_set_suite ());
	srunner_add_suite (sr, make_correlation_suite ());
	srunner_add_suite (sr, make_bed_suite ());
	srunner_add_suite (sr, make_blacklist_suite ());
	srunner_add_suite (sr, make_retrocopy_suite ());
	srunner_add_suite (sr, make_dedup_suite ());
	srunner_add_suite (sr, make_genotype_suite ());
	srunner_add_suite (sr, make_fasta_suite ());
	srunner_add_suite (sr, make_vcf_suite ());
	srunner_add_suite (sr, make_gz_suite ());
	srunner_add_suite (sr, make_graph_suite ());
	srunner_add_suite (sr, make_debrujin_suite ());
	srunner_add_suite (sr, make_floyd_warshall_suite ());
	srunner_add_suite (sr, make_hungarian_suite ());
	srunner_add_suite (sr, make_graph_unipath_suite ());
	srunner_add_suite (sr, make_graph_enumerate_suite ());
	/*srunner_set_tap (sr, "-");*/

	srunner_run_all (sr, CK_NORMAL);
	number_failed = srunner_ntests_failed (sr);
	srunner_free (sr);

	return number_failed;
}
