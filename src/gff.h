#ifndef GFF_H
#define GFF_H

#include <stdlib.h>
#include <zlib.h>
#include "hash.h"

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

struct _GffFilter
{
	const char *feature;
	Hash       *hard_attributes;
	Hash       *soft_attributes;
	Hash       *values_regex;
};

typedef struct _GffFilter GffFilter;

GffFile    * gff_open_for_reading               (const char *path);
void         gff_close                          (GffFile *gff);
int          gff_read                           (GffFile *gff, GffEntry *entry);
GffEntry   * gff_entry_new                      (void);
GffEntry   * gff_entry_dup                      (const GffEntry *from);
void         gff_entry_copy                     (GffEntry *to, const GffEntry *from);
void         gff_entry_free                     (GffEntry *entry);
const char * gff_attribute_find                 (const GffEntry *entry, const char *key);
GffFilter  * gff_filter_new                     (const char *feature);
void         gff_filter_free                    (GffFilter *filter);
void         gff_filter_insert_hard_attribute   (GffFilter *filter, const char *key, const char *value);
void         gff_filter_insert_soft_attribute   (GffFilter *filter, const char *key, const char *value);
int          gff_read_filtered                  (GffFile *gff, GffEntry *entry, const GffFilter *filter);

#define gff_attribute_get(entry,i) (&(entry)->attributes[(i)])

#endif /* gff.h */
