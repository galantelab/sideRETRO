#ifndef RETROCOPY_H
#define RETROCOPY_H

/*
 * Types of retrocopy insertion context
 * - RETROCOPY_PASS => When there is no
 *   overlapping among other clusters;
 * - RETROOPY_OVERLAPPED_PARENTALS
 *   => Clusters overlap and their parental
 *   genes also overlaps each other;
 * - RETROOPY_NEAR_PARENTALS => Clusters
 *   overlap and their parental genes are
 *   next to each other;
 * - RETROOPY_HOTSPOT => Clusters overlap
 *   and their parental genes are not near
 *   or are at different chromosomes; Each
 *   cluster needs to be present at least
 *   at one individual at a time (mutually
 *   exclusive);
 * - RETROOPY_AMBIGUOUS => Clusters overlap
 *   and their parental genes are not near
 *   or are at different chromosomes - just
 *   like in RETROOPY_HOTSPOT - but there is
 *   no evidence of presence of each cluster
 *   at least at one individual;
 * P.S: Clusters are not merged in
 *      RETROCOPY_PASS and RETROOPY_HOTSPOT;
 *      Retrocopy can accumulate flags for
 *      RETROOPY_HOTSPOT and RETROOPY_AMBIGUOUS;
 */
enum _RetrocopyType
{
	RETROCOPY_PASS                = 1,
	RETROOPY_OVERLAPPED_PARENTALS = 2,
	RETROOPY_NEAR_PARENTALS       = 4,
	RETROOPY_HOTSPOT              = 8,
	RETROOPY_AMBIGUOUS            = 16
};

typedef enum _RetrocopyType RetrocopyType;

#endif /* retrocopy.h */
