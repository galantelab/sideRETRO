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

#ifndef FASTA_H
#define FASTA_H

#include <stdlib.h>
#include "gz.h"
#include "str.h"

struct _FastaFile
{
	GzFile       *gz;
	char         *buf;
	size_t        buf_size;
	unsigned int  eof:1;
};

typedef struct _FastaFile FastaFile;

struct _FastaEntry
{
	String       *contig;
	String       *sequence;
	size_t        num_line;
};

typedef struct _FastaEntry FastaEntry;

FastaFile *  fasta_open_for_reading (const char *path);
void         fasta_close            (FastaFile *fasta);
FastaEntry * fasta_entry_new        (void);
void         fasta_entry_free       (FastaEntry *entry);
int          fasta_read             (FastaFile *fasta, FastaEntry *entry);

#endif /* fasta.h */
