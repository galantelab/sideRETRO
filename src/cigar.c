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

enum _CigarYYTokenType
{
	CIGAR_YY_OP    = 258,
	CIGAR_YY_OPLEN = 259,
	CIGAR_YY_ERROR = 260
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
			token = CIGAR_YY_OPLEN;

			goto Return;
		}

	switch (*c)
		{
			case 'M':
			case 'I':
			case 'D':
			case 'N':
			case 'S':
			case 'H':
			case 'P':
			case '=':
			case 'X':
			case 'B':
				{
					*yyval = *c;
					token = CIGAR_YY_OP;

					goto Return;
				}
			default:
				{
					log_error ("Invalid CIGAR operation: '%c'", *c);
					token = CIGAR_YY_ERROR;
				}
		}

Return:
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

	// Operations
	uint32_t *cigar_o = NULL;
	uint32_t op = 0;
	int32_t n_ops = 0;
	char op_char = 0;
	char prev_op_char = 0;
	uint32_t op_len = 0;

	// Initial state
	int state = CIGAR_YY_OPLEN;
	int soft_end = 0;
	int hard_end = 0;

	for (token = cigar_yylex (cigar_str, &cache, &yyval);
			token;
			token = cigar_yylex (NULL, &cache, &yyval))
		{
			if (token == CIGAR_YY_OP)
				{
					op_char = (char) yyval;

					if (state == CIGAR_YY_OPLEN)
						{
							log_error ("Missing length to operation '%c' at cigar '%s'",
									op_char, cigar_str);
							goto Error;
						}

					if (hard_end)
						{
							log_error ("Premature ending at cigar '%s'",
									cigar_str);
							goto Error;
						}

					switch (op_char)
						{
						case 'H':
							{
								if (n_ops > 0)
									hard_end = 1;

								break;
							}
						case 'S':
							{
								if (n_ops > 0 && prev_op_char != 'H')
									soft_end = 1;

								break;
							}
						default:
							{
								if (soft_end)
									{
										log_error ("Premature ending at cigar '%s'",
												cigar_str);
										goto Error;
									}
							}
						}

					n_ops++;
					prev_op_char = op_char;
					state = CIGAR_YY_OPLEN;
				}
			else if (token == CIGAR_YY_OPLEN)
				{
					op_len = yyval;
					state = CIGAR_YY_OP;
				}
			else
			{
				log_error ("Error parsing cigar '%s'", cigar_str);
				goto Error;
			}
		}

	if (state == CIGAR_YY_OP)
		{
			log_error ("Truncated cigar '%s'", cigar_str);
			goto Error;
		}

	return cigar_o;

Error:
	xfree (cigar_o);
	return NULL;
}
