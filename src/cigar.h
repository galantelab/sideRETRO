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

#ifndef CIGAR_H
#define CIGAR_H

#include <stdint.h>

/*
 * The parsed cigar follows the same encoding of
 * the BAM format (and htslib): The lower 4 bits
 * hold the operation index into CIGAR_STR and the
 * upper 28 bits hold the operation length
 */

#define CIGAR_STR        "MIDNSHP=XB"
#define CIGAR_STR_SAFE   CIGAR_STR "??????"
#define CIGAR_SHIFT      4
#define CIGAR_MASK       0xf
#define CIGAR_MAX_LEN    (UINT32_MAX >> CIGAR_SHIFT)

/*
 * 2 bits for each operation in CIGAR_STR order:
 * bit 1 consumes query and bit 2 consumes reference
 */
#define CIGAR_TYPE       0x3C1A7

#define cigar_idx(t)     ((t) - CIGAR_MATCH)
#define cigar_chr(t)     (CIGAR_STR_SAFE[cigar_idx(t)])

#define cigar_opidx(c)   ((c) & CIGAR_MASK)
#define cigar_op(c)      (cigar_opidx(c) + CIGAR_MATCH)
#define cigar_oplen(c)   ((c) >> CIGAR_SHIFT)
#define cigar_opchr(c)   (CIGAR_STR_SAFE[cigar_opidx(c)])

#define cigar_type(t)          ((CIGAR_TYPE >> (cigar_idx(t) << 1)) & 3)
#define cigar_consumes_query(t) (cigar_type(t) & 1)
#define cigar_consumes_ref(t)   (cigar_type(t) & 2)

enum _CigarTokenType
{
	CIGAR_MATCH      = 258,
	CIGAR_INS        = 259,
	CIGAR_DEL        = 260,
	CIGAR_REF_SKIP   = 261,
	CIGAR_SOFT_CLIP  = 262,
	CIGAR_HARD_CLIP  = 263,
	CIGAR_PAD        = 264,
	CIGAR_EQUAL      = 265,
	CIGAR_DIFF       = 266,
	CIGAR_BACK       = 267,
	CIGAR_LENGTH     = 268,
	CIGAR_ERROR      = 269
};

typedef enum _CigarTokenType CigarTokenType;

/*
 * Cigar navigator: Walk through the cigar operations
 * keeping track of the coordinates into the query
 * (read sequence, SEQ) and into the reference
 */

struct _CigarOp
{
	CigarTokenType  type;
	uint32_t        len;

	// 0-based offset into the read sequence (SEQ)
	// where the operation starts, and how many bases
	// of the read it consumes
	int             qpos;
	int             qlen;

	// Reference position where the operation starts
	// (same base as the 'pos' passed to the iterator),
	// and how many bases of the reference it consumes
	long            rpos;
	long            rlen;
};

typedef struct _CigarOp CigarOp;

struct _CigarIter
{
	/* PRIVATE  */
	const uint32_t *cigar;
	uint32_t        n_cigar;
	uint32_t        idx;
	int             qpos;
	long            rpos;
};

typedef struct _CigarIter CigarIter;

int        cigar_yylex      (const char *cigar_str, const char **cache, uint32_t *yyval);
uint32_t * cigar_yyparser   (const char *cigar_str, uint32_t *n_cigar);

void       cigar_iter_init  (CigarIter *iter, const uint32_t *cigar, uint32_t n_cigar, long pos);
int        cigar_iter_next  (CigarIter *iter, CigarOp *op);

char     * cigar_op_seq     (const CigarOp *op, const char *seq);

void       cigar_soft_clip  (const uint32_t *cigar, uint32_t n_cigar, uint32_t *left, uint32_t *right);
void       cigar_clip       (const uint32_t *cigar, uint32_t n_cigar, uint32_t *left, uint32_t *right);

#endif /* cigar.h */
