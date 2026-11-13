/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2023 Thiago L. A. Miller <tmiller@mochsl.org.br
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

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "log.h"
#include "wrapper.h"
#include "cigar.h"

int
cigar_yylex (const char *cigar_str, const char **cache,
		uint32_t *yyval)
{
	assert (cache != NULL);
	assert (yyval != NULL);

	const char *c = NULL;
	int token = 0;

	// Initiate cache
	if (cigar_str != NULL)
		*cache = cigar_str;

	// To facilitate access pointer
	c = *cache;

	// EOF
	if (*c == '\0')
		return 0;

	if (isdigit (*c))
		{
			uint32_t i = *c - '0';
			int overflow = 0;

			while (isdigit (*(++c)))
				{
					// Keep consuming the digits, but avoid
					// wrapping around uint32_t
					if (i > (UINT32_MAX - 9) / 10)
						overflow = 1;
					else
						i = (10 * i) + *c - '0';
				}

			// Unget
			c--;

			if (overflow)
				{
					log_error ("CIGAR operation length overflow");
					token = CIGAR_ERROR;
				}
			else
				{
					*yyval = i;
					token = CIGAR_LENGTH;
				}
		}
	else
		{
			switch (*c)
				{
				case 'M': {token = CIGAR_MATCH;     break;}
				case 'I': {token = CIGAR_INS;       break;}
				case 'D': {token = CIGAR_DEL;       break;}
				case 'N': {token = CIGAR_REF_SKIP;  break;}
				case 'S': {token = CIGAR_SOFT_CLIP; break;}
				case 'H': {token = CIGAR_HARD_CLIP; break;}
				case 'P': {token = CIGAR_PAD;       break;}
				case '=': {token = CIGAR_EQUAL;     break;}
				case 'X': {token = CIGAR_DIFF;      break;}
				case 'B': {token = CIGAR_BACK;      break;}
				default:
					{
						log_error ("Invalid CIGAR operation: '%c'", *c);
						token = CIGAR_ERROR;
					}
				}
		}

	*cache = c + 1;
	return token;
}

uint32_t *
cigar_yyparser (const char *cigar_str, uint32_t *n_cigar)
{
	assert (cigar_str != NULL);
	assert (n_cigar != NULL);

	// For cigar_yylex
	const char *cache = NULL;
	uint32_t yyval = 0;
	int token = 0;
	int prev_token = 0;

	// Operations
	uint32_t op_pair = 0;
	uint32_t op_len = 0;

	// Object
	uint32_t n_ops = 0;
	uint32_t alloc = 8;
	uint32_t *cigar = xcalloc (alloc, sizeof (uint32_t));

	// Initial state
	int init_op = 0;
	int soft_end = 0;
	int hard_end = 0;

	// SAM: '*' means the cigar is unavailable
	if (cigar_str[0] == '*' && cigar_str[1] == '\0')
		{
			*n_cigar = 0;
			return cigar;
		}

	if (cigar_str[0] == '\0')
		{
			log_error ("Empty cigar");
			goto Error;
		}

	for (token = cigar_yylex (cigar_str, &cache, &yyval);
			token;
			token = cigar_yylex (NULL, &cache, &yyval))
		{
			if (token == CIGAR_ERROR)
				{
					log_error ("Error parsing cigar '%s'", cigar_str);
					goto Error;
				}
			else if (token == CIGAR_LENGTH)
				{
					if (yyval > CIGAR_MAX_LEN)
						{
							log_error ("Operation length %u is greater than %u at cigar '%s'",
								yyval, CIGAR_MAX_LEN, cigar_str);
							goto Error;
						}

					op_len = yyval;
					init_op = 1;
				}
			else
				{
					if (!init_op)
						{
							log_error ("Missing length to operation '%c' at cigar '%s'",
								cigar_chr (token), cigar_str);
							goto Error;
						}

					if (hard_end)
						{
							log_error ("Premature ending at cigar '%s'",
								cigar_str);
							goto Error;
						}

					if (token == CIGAR_HARD_CLIP && n_ops > 0)
						{
							hard_end = 1;
						}
					else if (token == CIGAR_SOFT_CLIP
							&& (n_ops > 0 && prev_token != CIGAR_HARD_CLIP))
						{
							soft_end = 1;
						}
					else if (soft_end)
						{
							log_error ("Premature ending at cigar '%s'",
								cigar_str);
							goto Error;
						}

					// Set operation: 4 bits = char and 28 bits = len
					op_pair = (op_len << CIGAR_SHIFT) | cigar_idx (token);
					cigar[n_ops++] = op_pair;

					if (n_ops >= alloc)
						cigar = xrealloc (cigar, sizeof (uint32_t) * (alloc += 8));

					prev_token = token;
					init_op = 0;
				}
		}

	if (init_op)
		{
			log_error ("Truncated cigar '%s'", cigar_str);
			goto Error;
		}

	*n_cigar = n_ops;
	return xrealloc (cigar, sizeof (uint32_t) * n_ops);

Error:
	xfree (cigar);
	return NULL;
}

void
cigar_iter_init (CigarIter *iter, const uint32_t *cigar,
		uint32_t n_cigar, long pos)
{
	assert (iter != NULL);
	assert (cigar != NULL || n_cigar == 0);

	*iter = (CigarIter) {
		.cigar   = cigar,
		.n_cigar = n_cigar,
		.idx     = 0,
		.qpos    = 0,
		.rpos    = pos
	};
}

int
cigar_iter_next (CigarIter *iter, CigarOp *op)
{
	assert (iter != NULL);
	assert (op != NULL);

	uint32_t c = 0;

	// End of cigar
	if (iter->idx >= iter->n_cigar)
		return 0;

	c = iter->cigar[iter->idx++];

	op->type = cigar_op (c);
	op->len = cigar_oplen (c);

	// Where the operation starts
	op->qpos = iter->qpos;
	op->rpos = iter->rpos;

	// How many bases it consumes. Hard clip and
	// padding consume neither query nor reference
	op->qlen = cigar_consumes_query (op->type) ? op->len : 0;
	op->rlen = cigar_consumes_ref (op->type) ? op->len : 0;

	// Move forward
	iter->qpos += op->qlen;
	iter->rpos += op->rlen;

	return 1;
}

char *
cigar_op_seq (const CigarOp *op, const char *seq)
{
	assert (op != NULL);
	assert (seq != NULL);

	char *slice = NULL;

	// Operation does not consume the query: D, N, H, P, B
	if (op->qlen == 0)
		return NULL;

	if ((size_t) (op->qpos + op->qlen) > strlen (seq))
		{
			log_error ("Operation '%c' at %d (len %d) is beyond the sequence length %zu",
					cigar_chr (op->type), op->qpos, op->qlen, strlen (seq));
			return NULL;
		}

	slice = xcalloc (op->qlen + 1, sizeof (char));
	memcpy (slice, seq + op->qpos, op->qlen);

	return slice;
}

static void
__cigar_clip (const uint32_t *cigar, uint32_t n_cigar,
		int hard, uint32_t *left, uint32_t *right)
{
	int op = 0;
	int i = 0;
	int j = 0;

	*left = 0;
	*right = 0;

	// Clip at the beginning. Skip hard clip
	// if only soft clip is required
	for (i = 0; i < (int) n_cigar; i++)
		{
			op = cigar_op (cigar[i]);

			if (op == CIGAR_HARD_CLIP && !hard)
				continue;

			if (op != CIGAR_SOFT_CLIP && op != CIGAR_HARD_CLIP)
				break;

			*left += cigar_oplen (cigar[i]);
		}

	// Clip at the end - do not cross
	// the beginning operations
	for (j = (int) n_cigar - 1; j > i; j--)
		{
			op = cigar_op (cigar[j]);

			if (op == CIGAR_HARD_CLIP && !hard)
				continue;

			if (op != CIGAR_SOFT_CLIP && op != CIGAR_HARD_CLIP)
				break;

			*right += cigar_oplen (cigar[j]);
		}
}

void
cigar_soft_clip (const uint32_t *cigar, uint32_t n_cigar,
		uint32_t *left, uint32_t *right)
{
	assert (cigar != NULL || n_cigar == 0);
	assert (left != NULL);
	assert (right != NULL);

	__cigar_clip (cigar, n_cigar, 0, left, right);
}

void
cigar_clip (const uint32_t *cigar, uint32_t n_cigar,
		uint32_t *left, uint32_t *right)
{
	assert (cigar != NULL || n_cigar == 0);
	assert (left != NULL);
	assert (right != NULL);

	__cigar_clip (cigar, n_cigar, 1, left, right);
}
