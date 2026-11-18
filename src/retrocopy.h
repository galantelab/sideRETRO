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

#ifndef RETROCOPY_H
#define RETROCOPY_H

#include "db.h"

/*
 * Types of retrocopy insertion context
 * - RETROCOPY_PASS => When there is no
 *   overlapping among other clusters;
 * - RETROCOPY_OVERLAPPED_PARENTALS
 *   => Clusters overlap and their parental
 *   genes also overlaps each other (or
 *   they are the same gene);
 * - RETROCOPY_NEAR_PARENTALS => Clusters
 *   overlap and their parental genes are
 *   next to each other;
 * - RETROCOPY_HOTSPOT => Clusters overlap
 *   and their parental genes are not near
 *   or are at different chromosomes: They
 *   are not merged, so there is a retrocopy
 *   for each parental gene at the same
 *   segment, whatever the individuals
 *   carrying them;
 * - RETROCOPY_AMBIGUOUS => Reserved: It is
 *   never set. It was meant for the clusters
 *   of a hotspot which cannot be told apart
 *   by the individuals carrying them;
 * - RETROCOPY_PARENTAL_CHOSEN => Overlapped
 *   parental genes, but the evidence points
 *   to one of them: It is the first name of
 *   the parental genes, followed by the
 *   others ('chosen/other'). Set together
 *   with RETROCOPY_OVERLAPPED_PARENTALS;
 * - RETROCOPY_TSITE_MERGED => Two retrocopies of
 *   the same parental gene, each one with a single
 *   breakpoint, facing each other: The flanks of
 *   a deletion at the target site longer than the
 *   cluster window. They were merged into one, and
 *   the flag is set along with the other levels;
 * P.S: Clusters are merged in
 *      RETROCOPY_OVERLAPPED_PARENTALS
 *      and RETROCOPY_NEAR_PARENTALS
 */
enum _RetrocopyLevel
{
	RETROCOPY_PASS                 = 1,
	RETROCOPY_OVERLAPPED_PARENTALS = 2,
	RETROCOPY_NEAR_PARENTALS       = 4,
	RETROCOPY_HOTSPOT              = 8,
	RETROCOPY_AMBIGUOUS            = 16,
	RETROCOPY_PARENTAL_CHOSEN      = 32,
	RETROCOPY_TSITE_MERGED         = 64
};

typedef enum _RetrocopyLevel RetrocopyLevel;


/*
 * Types of insertion points calculation
 * - RETROCOPY_INSERTION_POINT_WINDOW_MEAN =>
 *   insertion point is calculated on the
 *   window range mean;
 * - RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE =>
 *   insertion point is calculated on the supplementary
 *   reads mode;
 */
enum _RetrocopyInsertionPoint
{
	RETROCOPY_INSERTION_POINT_WINDOW_MEAN        = 1,
	RETROCOPY_INSERTION_POINT_SUPPLEMENTARY_MODE = 2
};

typedef enum _RetrocopyInsertionPoint RetrocopyInsertionPoint;

/*
 * Minimum clipped length of an alignment to be
 * used as a breakpoint candidate. Avoid noise
 * from the ends of the reads
 */
#define RETROCOPY_MIN_CLIP 10

/*
 * Minimum run of A (or T) at the clipped sequence
 * of a read crossing the breakpoint to be called
 * as the poly-A tail
 */
#define RETROCOPY_MIN_POLYA 6

/*
 * Maximum number of bases between the flank and the
 * poly-A tail: The breakpoint may be shifted by the
 * aligner, or untemplated bases added at the junction
 */
#define RETROCOPY_MAX_POLYA_OFFSET 3

/*
 * Choice among overlapped parental genes: Minimum
 * number of fragments explained by one gene alone,
 * and maximum p-value of the orientation to compare
 * its strand with the poly-A strand
 */
#define RETROCOPY_MIN_EXCLUSIVE 3
#define RETROCOPY_ORIENTATION_ALPHA 0.05

/*
 * Maximum distance between the end of a piece of a read
 * spanning an exon-exon junction and the exon boundary:
 * Inside the exon (bases left unaligned by mismatches near
 * the junction) and outside it, into the intron (sequence
 * shared by the exon end and the intron before the other
 * exon, aligned at both pieces)
 */
#define RETROCOPY_JUNCTION_MAX_INNER 5
#define RETROCOPY_JUNCTION_MAX_OUTER 10

/*
 * Deletion at the target site longer than the cluster
 * window: Default maximum length, and the number of
 * fragments of a sample at one flank, with none at the
 * other, that tells two different events
 */
#define RETROCOPY_DEFAULT_MAX_TSD_DEL 50000
#define RETROCOPY_TSITE_MIN_CONFLICT 5

/*
 * Seed used to search the clipped sequence at the
 * junction into the parental exons: Its maximum
 * length and the minimum length allowed
 */
/*
 * Role of a read in the evidence of a retrocopy. A read
 * may have more than one role. A read with none is a
 * discordant read
 */
enum _RetrocopyEvidence
{
	RETROCOPY_EVIDENCE_SPLIT_LEFT  = 1,
	RETROCOPY_EVIDENCE_SPLIT_RIGHT = 2,
	RETROCOPY_EVIDENCE_POLYA       = 4,
	RETROCOPY_EVIDENCE_JUNCTION    = 8
};

typedef enum _RetrocopyEvidence RetrocopyEvidence;

#define RETROCOPY_SEED_LEN 20
#define RETROCOPY_MIN_SEED 12

void retrocopy (sqlite3_stmt *retrocopy_stmt,
		sqlite3_stmt *cluster_merging_stmt,
		sqlite3_stmt *junction_stmt,
		sqlite3_stmt *insertion_stmt,
		sqlite3_stmt *evidence_stmt,
		const char *fasta_file,
		int near_gene_dist,
		long max_tsd_del);

#endif /* retrocopy.h */
