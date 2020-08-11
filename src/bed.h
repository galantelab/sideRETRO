/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2020 Thiago L. A. Miller <tmiller@mochsl.org.br
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

#ifndef BED_H
#define BED_H

#include <stdlib.h>
#include "gz.h"

struct _BedFile
{
	GzFile       *gz;
	const char   *header;
	char         *buf;
	size_t        buf_size;
	unsigned int  eof:1;
};

typedef struct _BedFile BedFile;

struct _BedEntry
{
	size_t        chrom_size;
	size_t        name_size;
	char         *chrom;
	long          chrom_start;
	long          chrom_end;
	char         *name;
	int           score;
	char          strand;
	long          thick_start;
	long          thick_end;
	int           rgb[3];
	size_t        num_blocks;
	int           block_count;
	int          *block_sizes;
	int          *block_starts;
	size_t        num_field;
	size_t        num_line;
};

typedef struct _BedEntry BedEntry;

BedFile *  bed_open_for_reading (const char *path);
void       bed_close            (BedFile *bed);

BedEntry * bed_entry_new        (void);
void       bed_entry_free       (BedEntry *entry);

int        bed_read             (BedFile *bed, BedEntry *entry);

#endif /* bed.h */
