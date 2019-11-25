#ifndef FASTA_H
#define FASTA_H

#include <stdlib.h>
#include <zlib.h>
#include "str.h"

struct _FastaFile
{
	gzFile        fp;
	size_t        buf_alloc;
	char         *buf;
	const char   *filename;
	size_t        num_line;
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
