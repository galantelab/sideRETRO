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
