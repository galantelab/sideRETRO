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

#ifndef VCF_H
#define VCF_H

#include "db.h"

struct _VCFOption
{
	long        near_gene_dist;
	float       orientation_error;
	const char *fasta_file;

	// Read names and predicted insertions files
	const char *rnames_file;
	const char *insertion_file;

	// Filters: Thresholds and the disabled ones
	// (VCFFilter flags)
	double      min_qual;
	int         min_mapq;
	int         disabled_filters;

	// Output: Number of retrocopies written to
	// the VCF and to the insertions file
	int         num_records;
	int         num_insertions;
	int         num_pass;
};

typedef struct _VCFOption VCFOption;

/*
 * FILTER column: Each flag is a check a canonical
 * retrocopy should pass. A check that cannot be
 * evaluated (e.g. TSD without both breakpoints)
 * does not fail
 */
enum _VCFFilter
{
	VCF_FILTER_LOW_QUAL           = 1,
	VCF_FILTER_LOW_MAPQ           = 2,
	VCF_FILTER_AMBIGUOUS_PARENTAL = 4,
	VCF_FILTER_NO_TSD             = 8,
	VCF_FILTER_TSD_OUTLIER        = 16,
	VCF_FILTER_NO_POLYA           = 32,
	VCF_FILTER_NO_JUNCTION        = 64,
	VCF_FILTER_POLARITY_CONFLICT  = 128
};

typedef enum _VCFFilter VCFFilter;

#define VCF_DEFAULT_MIN_QUAL      20.0
#define VCF_DEFAULT_MIN_MAPQ      20
#define VCF_TSD_MAX_LEN           50
#define VCF_TSD_MIN_LEN           -20

int vcf_filter_lookup (const char *id);

void vcf (sqlite3 *db, const char *output_file, VCFOption *opt);

#endif /* vcf.h */
