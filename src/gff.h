#ifndef GFF_H
#define GFF_H

#include <stdlib.h>
#include <zlib.h>

struct _GffFile
{
	gzFile        fp;
	const char   *filename;
	const char   *header;
	char         *buf;
	size_t        buf_size;
	size_t        num_line;
	unsigned int  eof:1;
};

typedef struct _GffFile GffFile;

struct _GffAttribute
{
	size_t        key_size;
	size_t        value_size;
	char         *key;
	char         *value;
};

typedef struct _GffAttribute GffAttribute;

struct _GffEntry
{
	size_t        seqname_size;
	size_t        source_size;
	size_t        feature_size;
	char         *seqname;
	char         *source;
	char         *feature;
	size_t        start;
	size_t        end;
	float         score;
	char          strand;
	short int     frame;
	size_t        attributes_size;
	size_t        num_attributes;
	GffAttribute *attributes;
	size_t        num_line;
};

typedef struct _GffEntry GffEntry;

GffFile    * gff_open           (const char *path, const char *mode);
void         gff_close          (GffFile *gff);
int          gff_read           (GffFile *gff, GffEntry *entry);
GffEntry   * gff_entry_new      (void);
GffEntry   * gff_entry_copy     (const GffEntry *entry);
void         gff_entry_free     (GffEntry *entry);
const char * gff_attribute_find (GffEntry *entry, const char *key);

#define gff_attribute_get(entry,i) (&(entry)->attributes[(i)])

#endif /* gff.h */
