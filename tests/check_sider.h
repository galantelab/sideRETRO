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

#ifndef CHECK_SIDER_H
#define CHECK_SIDER_H

Suite * make_list_suite           (void);
Suite * make_hash_suite           (void);
Suite * make_array_suite          (void);
Suite * make_utils_suite          (void);
Suite * make_sam_suite            (void);
Suite * make_bitree_suite         (void);
Suite * make_ibitree_suite        (void);
Suite * make_str_suite            (void);
Suite * make_db_suite             (void);
Suite * make_chr_suite            (void);
Suite * make_exon_suite           (void);
Suite * make_abnormal_suite       (void);
Suite * make_gff_suite            (void);
Suite * make_io_suite             (void);
Suite * make_process_sample_suite (void);
Suite * make_dbscan_suite         (void);
Suite * make_cluster_suite        (void);
Suite * make_wrapper_suite        (void);
Suite * make_db_merge_suite       (void);
Suite * make_set_suite            (void);
Suite * make_correlation_suite    (void);
Suite * make_bed_suite            (void);
Suite * make_blacklist_suite      (void);
Suite * make_retrocopy_suite      (void);
Suite * make_dedup_suite          (void);
Suite * make_genotype_suite       (void);
Suite * make_fasta_suite          (void);
Suite * make_vcf_suite            (void);
Suite * make_gz_suite             (void);
Suite * make_graph_suite          (void);
Suite * make_debrujin_suite       (void);
Suite * make_floyd_warshall_suite (void);
Suite * make_hungarian_suite      (void);
Suite * make_graph_unipath_suite  (void);

#endif /* check_sider.h */
