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

#ifndef TANDEM_H
#define TANDEM_H

/*
 * A perfect tandem repeat: A unit of 1 to TANDEM_MAX_UNIT
 * bases repeated at least TANDEM_MIN_COPIES times, side by
 * side, along at least TANDEM_MIN_SPAN bases. The last copy
 * may be partial, as 'CACACAC'
 */
#define TANDEM_MAX_UNIT   6
#define TANDEM_MIN_COPIES 3
#define TANDEM_MIN_SPAN   10

struct _TandemRepeat
{
	// Offset of the first base at the sequence
	long start;
	// Number of bases covered
	long span;
	// The unit, as it is at the start
	int  unit_len;
	char unit[TANDEM_MAX_UNIT + 1];
};

typedef struct _TandemRepeat TandemRepeat;

int tandem_repeat_find (const char *seq, long len, long from, long to,
		TandemRepeat *tr);

#endif /* tandem.h */
