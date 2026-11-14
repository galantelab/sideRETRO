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

#ifndef GFF_H
#define GFF_H

#include <stdlib.h>
#include "hash.h"
#include "list.h"
#include "gz.h"

struct _GffFile
{
	GzFile       *gz;
	const char   *header;
	char         *buf;
	size_t        buf_size;
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

/*
 * Each attribute group comes from one 'key=value' option and
 * holds alternative key/value pairs: 'k1=v1|v2|k2=v3' means
 * k1=v1 OR k1=v2 OR k2=v3. An entry passes the filter if it
 * matches the feature (if any) and at least one pair of EACH
 * group. Values match exactly
 */
struct _GffFilter
{
	const char *feature;
	List       *attributes;
};

typedef struct _GffFilter GffFilter;

int          gff_looks_like_gff_file            (const char *filename);
GffFile    * gff_open_for_reading               (const char *path);
void         gff_close                          (GffFile *gff);
int          gff_read                           (GffFile *gff, GffEntry *entry);
GffEntry   * gff_entry_new                      (void);
GffEntry   * gff_entry_dup                      (const GffEntry *from);
void         gff_entry_copy                     (GffEntry *to, const GffEntry *from);
void         gff_entry_free                     (GffEntry *entry);
const char * gff_attribute_find                 (const GffEntry *entry, const char *key);
GffFilter  * gff_filter_new                     (void);
void         gff_filter_insert_feature          (GffFilter *filter, const char *feature);
void         gff_filter_free                    (GffFilter *filter);
int          gff_filter_insert_attribute        (GffFilter *filter, const char *key_values);
int          gff_filter_match                   (const GffFilter *filter, const GffEntry *entry);
int          gff_read_filtered                  (GffFile *gff, GffEntry *entry, const GffFilter *filter);

#define gff_attribute_get(entry,i) (&(entry)->attributes[(i)])
#define gff_filter_attribute_size(filter) (list_size ((filter)->attributes))

#endif /* gff.h */
