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

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "gff.h"

#define GFF_ATTRSIZ 8

int
gff_looks_like_gff_file (const char *filename)
{
	assert (filename != NULL);

	regex_t reg;
	int rc = 0;

	const char regex_str[] = "\\.(gff|gff3|gtf)(\\.gz)*$";

	if (regcomp (&reg, regex_str, REG_EXTENDED|REG_NOSUB) != 0)
		log_fatal ("Error regcomp for '%s'", regex_str);

	rc = regexec (&reg, filename, 0, (regmatch_t *) NULL, 0);

	regfree (&reg);
	return rc == 0;
}

GffEntry *
gff_entry_new (void)
{
	GffEntry *entry = xcalloc (1, sizeof (GffEntry));
	return entry;
}

void
gff_entry_copy (GffEntry *to, const GffEntry *from)
{
	int i = 0;

	// Release the previous values, if any
	xfree (to->seqname);
	xfree (to->source);
	xfree (to->feature);

	to->seqname = xstrdup (from->seqname);
	to->source = xstrdup (from->source);
	to->feature = xstrdup (from->feature);

	to->seqname_size = strlen (from->seqname) + 1;
	to->source_size = strlen (from->source) + 1;
	to->feature_size = strlen (from->feature) + 1;

	to->start = from->start;
	to->end = from->end;
	to->score = from->score;
	to->strand = from->strand;
	to->frame = from->frame;

	for (i = 0; i < to->attributes_size; i++)
		{
			xfree (to->attributes[i].key);
			xfree (to->attributes[i].value);
		}

	to->attributes = xrealloc (to->attributes,
			from->num_attributes * sizeof (GffAttribute));

	for (i = 0; i < from->num_attributes; i++)
		{
			to->attributes[i].key =
				xstrdup (from->attributes[i].key);
			to->attributes[i].value =
				xstrdup (from->attributes[i].value);
			to->attributes[i].key_size =
				strlen (from->attributes[i].key) + 1;
			to->attributes[i].value_size =
				strlen (from->attributes[i].value) + 1;
		}

	to->num_attributes = from->num_attributes;
	to->attributes_size = from->num_attributes;
}

GffEntry *
gff_entry_dup (const GffEntry *from)
{
	GffEntry *to = gff_entry_new ();
	gff_entry_copy (to, from);
	return to;
}

void
gff_entry_free (GffEntry *entry)
{
	if (entry == NULL)
		return;

	int i = 0;

	for (; i < entry->attributes_size; i++)
		{
			xfree (entry->attributes[i].key);
			xfree (entry->attributes[i].value);
		}

	xfree (entry->attributes);
	xfree (entry->seqname);
	xfree (entry->source);
	xfree (entry->feature);

	xfree (entry);
}

static void
gff_get_header (GffFile *gff)
{
	char *header = NULL;
	int rc = 0;

	// Catch header '##' and ignore blank lines
	while ((rc = gz_getline (gff->gz, &gff->buf, &gff->buf_size))
			&& (gff->buf[0] == '#' || gff->buf[0] == '\n'))
		{
			if (gff->buf[0] == '#')
				header = xstrdup_concat (header, gff->buf);
		}

	// Reached end of file
	if (!rc)
		gff->eof = 1;

	gff->header = chomp (header);
}

GffFile *
gff_open_for_reading (const char *path)
{
	assert (path != NULL);

	GffFile *gff = NULL;

	gff = xcalloc (1, sizeof (GffFile));
	gff->gz = gz_open_for_reading (path);

	gff_get_header (gff);

	return gff;
}

void
gff_close (GffFile *gff)
{
	if (gff == NULL)
		return;

	gz_close (gff->gz);

	xfree ((void *) gff->header);
	xfree (gff->buf);
	xfree (gff);
}

const char *
gff_attribute_find (const GffEntry *entry, const char *key)
{
	assert (entry != NULL && key != NULL);

	const char *value = NULL;
	const GffAttribute *cur = NULL;
	int i = 0;

	for (; i < entry->num_attributes; i++)
		{
			cur = &entry->attributes[i];
			if (!strcmp (cur->key, key))
				{
					value = cur->value;
					break;
				}
		}

	return value;
}

int
gff_read (GffFile *gff, GffEntry *entry)
{
	assert (gff != NULL && entry != NULL);

	int rc, i;
	char *token, *subtoken;
	char *attr_key, *attr_value;
	char *saveptr1, *saveptr2, *saveptr3;

	token = subtoken = attr_key = attr_value = NULL;
	saveptr1 = saveptr2 = saveptr3 = NULL;

	// end of file
	if (gff->eof)
		return 0;

	gff->buf = chomp (gff->buf);
	entry->num_line = gz_get_num_line (gff->gz);

	token = strtok_r (gff->buf, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'seqname' (field 1) at line %zu", entry->num_line);

	entry->seqname_size = entry_set (&entry->seqname,
			entry->seqname_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'source' (field 2) at line %zu", entry->num_line);

	entry->source_size = entry_set (&entry->source,
			entry->source_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'feature' (field 3) at line %zu", entry->num_line);

	entry->feature_size = entry_set (&entry->feature,
			entry->feature_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'start' (field 4) at line %zu", entry->num_line);

	entry->start = atol (token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'end' (field 5) at line %zu", entry->num_line);

	entry->end = atol (token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'score' (field 6) at line %zu", entry->num_line);

	entry->score = strcmp (token, ".")
		? atof (token)
		: -1;

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'strand' (field 7) at line %zu", entry->num_line);

	entry->strand = *token;

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("Missing 'frame' (field 8) at line %zu", entry->num_line);

	entry->frame = strcmp (token, ".")
		? atoi (token)
		: -1 ;

	token = strtok_r (NULL, "", &saveptr1);
	entry->num_attributes = i = 0;

	if (token != NULL)
		{
			for (subtoken = strtok_r (token, ";", &saveptr2); subtoken != NULL;
					subtoken = strtok_r (NULL, ";", &saveptr2))
				{
					attr_key = strtok_r (subtoken, " =", &saveptr3);
					attr_value = strtok_r (NULL, "", &saveptr3);

					if (attr_value == NULL)
						log_fatal ("Missing value for attribute %d (%s) at line %zu",
								i, attr_key, entry->num_line);

					if ((i + 1) > entry->attributes_size)
						entry->attributes_size = buf_expand ((void **) &entry->attributes,
								sizeof (GffAttribute), entry->attributes_size, GFF_ATTRSIZ);

					attr_key = trim (attr_key);

					entry->attributes[i].key_size = entry_set (&entry->attributes[i].key,
							entry->attributes[i].key_size, attr_key);

					attr_value = trimc (trim (attr_value), '"');

					entry->attributes[i].value_size = entry_set (&entry->attributes[i].value,
							entry->attributes[i].value_size, attr_value);

					i++;
				}
		}

	entry->num_attributes = i;

	/*
	* get next entry
	* ignore comments and blank lines
	*/
	while ((rc = gz_getline (gff->gz, &gff->buf, &gff->buf_size))
			&& (gff->buf[0] == '#' || gff->buf[0] == '\n'))
		;

	// Reached end of file
	if (!rc)
		gff->eof = 1;

	return 1;
}

struct _GffPair
{
	char *key;
	char *value;
};

typedef struct _GffPair GffPair;

static void
gff_pair_free (GffPair *pair)
{
	if (pair == NULL)
		return;

	xfree (pair->key);
	xfree (pair->value);
	xfree (pair);
}

static void
gff_group_free (List *group)
{
	list_free (group);
}

GffFilter *
gff_filter_new (void)
{
	GffFilter *filter = xcalloc (1, sizeof (GffFilter));

	filter->attributes = list_new ((DestroyNotify) gff_group_free);

	return filter;
}

void
gff_filter_free (GffFilter *filter)
{
	if (filter == NULL)
		return;

	xfree ((void *) filter->feature);
	list_free (filter->attributes);

	xfree (filter);
}

void
gff_filter_insert_feature (GffFilter *filter, const char *feature)
{
	xfree ((void *) filter->feature);
	filter->feature = feature != NULL
		? xstrdup (feature)
		: NULL;
}

int
gff_filter_insert_attribute (GffFilter *filter, const char *key_values)
{
	assert (filter != NULL && key_values != NULL);

	List *group = NULL;
	GffPair *pair = NULL;
	char *copy = NULL;
	char *scratch = NULL;
	char *token = NULL;
	char *eq = NULL;
	const char *key = NULL;
	int rc = 1;

	/*
	* Parse 'k1=v1|v2|k2=v3': A token with '=' sets the key,
	* a token without it is another value for the last key
	*/
	group = list_new ((DestroyNotify) gff_pair_free);
	copy = xstrdup (key_values);

	for (token = strtok_r (copy, "|", &scratch); token != NULL;
			token = strtok_r (NULL, "|", &scratch))
		{
			eq = strchr (token, '=');

			if (eq != NULL)
				{
					*eq = '\0';
					key = token;
					token = eq + 1;
				}

			if (key == NULL || *key == '\0' || *token == '\0')
				{
					rc = 0;
					break;
				}

			pair = xcalloc (1, sizeof (GffPair));
			pair->key = xstrdup (key);
			pair->value = xstrdup (token);

			list_append (group, pair);
		}

	if (rc && list_size (group) > 0)
		list_append (filter->attributes, group);
	else
		{
			rc = 0;
			list_free (group);
		}

	xfree (copy);
	return rc;
}

static inline int
gff_value_match (const char *value, const char *expected)
{
	size_t len = strlen (expected);
	const char *cur = value;
	const char *next = NULL;

	/*
	* Exact match. GFF3 may pack multiple values separated
	* by comma (e.g. 'tag=basic,CCDS'), so test each one
	*/
	while (cur != NULL)
		{
			next = strchr (cur, ',');

			if ((next == NULL && !strcmp (cur, expected))
					|| (next != NULL && (size_t) (next - cur) == len
						&& !strncmp (cur, expected, len)))
				return 1;

			cur = next != NULL ? next + 1 : NULL;
		}

	return 0;
}

static inline int
gff_group_match (const List *group, const GffEntry *entry)
{
	const ListElmt *cur = NULL;
	const GffPair *pair = NULL;
	const GffAttribute *attr = NULL;
	size_t i = 0;

	// Any pair of the group. A key may appear more than
	// once into the attributes (e.g. 'tag' in GTF)
	for (cur = list_head (group); cur != NULL; cur = list_next (cur))
		{
			pair = list_data (cur);

			for (i = 0; i < entry->num_attributes; i++)
				{
					attr = gff_attribute_get (entry, i);

					if (!strcmp (attr->key, pair->key)
							&& gff_value_match (attr->value, pair->value))
						return 1;
				}
		}

	return 0;
}

int
gff_filter_match (const GffFilter *filter, const GffEntry *entry)
{
	assert (filter != NULL && entry != NULL);

	const ListElmt *cur = NULL;

	// Return 0 if there is a feature to be filtered
	// and it does not match
	if (filter->feature != NULL
			&& strcmp (filter->feature, entry->feature))
		return 0;

	// All groups must match. Without groups, any entry
	// passes and it works as a regular 'gff_read'
	for (cur = list_head (filter->attributes); cur != NULL; cur = list_next (cur))
		if (!gff_group_match (list_data (cur), entry))
			return 0;

	return 1;
}

int
gff_read_filtered (GffFile *gff, GffEntry *entry, const GffFilter *filter)
{
	assert (filter != NULL);

	while (gff_read (gff, entry))
		if (gff_filter_match (filter, entry))
			return 1;

	return 0;
}
