#ifndef RETROCOPY_H
#define RETROCOPY_H

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
 *   or are at different chromosomes; Each
 *   cluster needs to be present at least
 *   at one individual at a time (mutually
 *   exclusive);
 * - RETROCOPY_AMBIGUOUS => Clusters overlap
 *   and their parental genes are not near
 *   or are at different chromosomes - just
 *   like in RETROCOPY_HOTSPOT - but there is
 *   no evidence of presence of each cluster
 *   at least at one individual;
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
	RETROCOPY_AMBIGUOUS            = 16
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

void retrocopy (sqlite3_stmt *retrocopy_stmt,
		sqlite3_stmt *cluster_merging_stmt);

#endif /* retrocopy.h */
