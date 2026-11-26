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

	// Sequences for primer design: It needs the
	// reference genome. Bases of each flank
	const char *primer_file;
	int         primer_flank;

	// Filters: Thresholds and the disabled ones
	// (VCFFilter flags)
	double      min_qual;
	int         min_mapq;
	int         min_sr;
	int         min_tsd;
	int         max_tsd;
	int         min_pgmates;
	int         disabled_filters;

	// Primary chromosomes (NonPrimary): An extended
	// regular expression. NULL for the default
	const char *primary_chr;

	// Annotation: Minimum number of reads of a
	// sample to report an intron as retained
	int         min_ir;

	// Output: Number of retrocopies written to
	// the VCF and to the insertions file
	int         num_records;
	int         num_insertions;
	int         num_primers;
	int         num_pass;
};

typedef struct _VCFOption VCFOption;

/*
 * FILTER column: Each flag is a check that failed.
 * A PASS retrocopy has enough split reads and both
 * breakpoints, which describe the target site: A
 * duplication, a blunt insertion or a deletion. Only
 * the duplication length is checked - a long one comes
 * from breakpoints of unrelated reads. And it needs a
 * hallmark of retrotransposition: An exon-exon junction,
 * the poly-A tail or a duplication outside a tandem
 * repeat. It must be at a primary chromosome - not at
 * an alt, decoy, unplaced or HLA contig - and a copy
 * longer than a read: Some mates inside the parental
 * exons
 */
enum _VCFFilter
{
	VCF_FILTER_LOW_QUAL           = 1,
	VCF_FILTER_LOW_MAPQ           = 2,
	VCF_FILTER_AMBIGUOUS_PARENTAL = 4,
	VCF_FILTER_LOW_SR             = 8,
	VCF_FILTER_NO_BREAKPOINT      = 16,
	VCF_FILTER_ONE_BREAKPOINT     = 32,
	VCF_FILTER_TSD_OUTLIER        = 64,
	VCF_FILTER_NO_HALLMARK        = 128,
	VCF_FILTER_NON_PRIMARY        = 256,
	VCF_FILTER_NO_PARENTAL_MATE   = 512
};

typedef enum _VCFFilter VCFFilter;

#define VCF_DEFAULT_MIN_QUAL      20.0
#define VCF_DEFAULT_MIN_MAPQ      20
#define VCF_DEFAULT_MIN_SR        3
#define VCF_DEFAULT_MIN_TSD       5
#define VCF_DEFAULT_MAX_TSD       25
#define VCF_DEFAULT_MIN_PGMATES   1
#define VCF_DEFAULT_PRIMER_FLANK  200
#define VCF_DEFAULT_MIN_IR        2

// Chromosomes 1-22, X and Y, with or without 'chr'
#define VCF_DEFAULT_PRIMARY_CHR   "^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$"

int vcf_filter_lookup (const char *id);

void vcf (sqlite3 *db, const char *output_file, VCFOption *opt);

#endif /* vcf.h */
