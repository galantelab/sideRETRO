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

#include "config.h"

#include <string.h>
#include <assert.h>
#include "tandem.h"

/*
 * The run of the unit of 'k' bases starting at 'i': It
 * goes on while each base repeats the one 'k' bases
 * before. An unknown base (N) breaks it. Return the end
 * of the run (one past the last base)
 */
static long
tandem_run_end (const char *seq, long len, long i, int k)
{
	long j = 0;

	for (j = i; j < i + k; j++)
		if (j >= len || seq[j] == 'N')
			return j;

	while (j < len && seq[j] != 'N' && seq[j] == seq[j - k])
		j++;

	return j;
}

/*
 * Find the longest perfect tandem repeat at 'seq' (upper
 * case, 'len' bases) that overlaps the bases from 'from'
 * to 'to' (offsets, inclusive). Between two of the same
 * length, the one of the shortest unit wins: 'AAAAAAAAAA'
 * is a run of 'A', not of 'AA'. Return 1 and fill 'tr' if
 * found, or 0
 */
int
tandem_repeat_find (const char *seq, long len, long from, long to,
		TandemRepeat *tr)
{
	assert (seq != NULL && tr != NULL);

	long i = 0;
	long j = 0;
	long span = 0;
	int k = 0;
	int found = 0;

	memset (tr, 0, sizeof (TandemRepeat));

	for (k = 1; k <= TANDEM_MAX_UNIT; k++)
		{
			i = 0;

			while (i < len)
				{
					j = tandem_run_end (seq, len, i, k);
					span = j - i;

					if (span >= TANDEM_MIN_COPIES * k
							&& span >= TANDEM_MIN_SPAN
							&& i <= to && (j - 1) >= from
							&& span > tr->span)
						{
							tr->start = i;
							tr->span = span;
							tr->unit_len = k;
							memcpy (tr->unit, seq + i, k);
							tr->unit[k] = '\0';
							found = 1;
						}

					/*
					* A run longer than its unit: The runs starting
					* inside it, but at its last unit, end at the
					* same base and are shorter
					*/
					i = span > k ? j - k + 1 : i + 1;
				}
		}

	return found;
}
