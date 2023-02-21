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

#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <assert.h>
#include "log.h"
#include "wrapper.h"

#include "cigar.h"

#define CIGAR_STR        "MIDNSHP=XB"
#define CIGAR_STR_SAFE   CIGAR_STR "??????"
#define CIGAR_SHIFT      4
#define CIGAR_MASK       0xf
#define CIGAR_TYPE       0x3C1A7

#define cigar_idx(t)     ((t) - CIGAR_MATCH)
#define cigar_chr(t)     (CIGAR_STR_SAFE[cigar_idx(t)])

#define cigar_opidx(c)   ((c) & CIGAR_MASK)
#define cigar_op(c)      (cigar_opidx(c) + CIGAR_MATCH)
#define cigar_oplen(c)   ((c) >> CIGAR_SHIFT)
#define cigar_opchr(c)   (CIGAR_STR_SAFE[cigar_opidx(c)])


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

			while (isdigit (*(++c)))
				i = (10 * i) + *c - '0';

			// Unget
			c--;

			*yyval = i;
			token = CIGAR_LENGTH;
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
