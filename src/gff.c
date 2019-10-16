#include "config.h"

#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "set.h"
#include "gff.h"

#define GFF_BUFSIZ  128
#define GFF_ATTRSIZ 8

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

	to->num_attributes = from->num_attributes;
	to->attributes_size = from->num_attributes;
	to->attributes = xcalloc (from->num_attributes,
			sizeof (GffAttribute));

	for (; i < from->num_attributes; i++)
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

static inline size_t
nearest_pow (size_t num)
{
	size_t n = 1;

	while (n < num && n > 0)
		n <<= 1;

	return n ? n : num;
}

static size_t
gff_buf_expand (void **buf, size_t size,
		size_t old_nmemb, size_t length)
{
	size_t final_nmemb = nearest_pow (old_nmemb + length);
	*buf = xrealloc (*buf, size * final_nmemb);
	memset (*buf + old_nmemb * size, 0,
			size * (final_nmemb - old_nmemb));
	return final_nmemb;
}

static size_t
gff_entry_set (char **buf, size_t buf_size, const char *entry)
{
	size_t entry_size = strlen (entry);

	if (entry_size >= buf_size)
		buf_size = gff_buf_expand ((void **) buf, sizeof (char),
				buf_size, entry_size - buf_size + 1);

	*buf = strncpy (*buf, entry, buf_size);
	return buf_size;
}

static int
gff_getline (GffFile *gff)
{
	const char *err_msg = NULL;
	char *line = NULL;
	size_t len = 0;
	int rc = 0;

	line = gzgets (gff->fp, gff->buf, gff->buf_size);
	if (line == NULL)
		{
			err_msg = gzerror (gff->fp, &rc);
			if (rc != Z_OK)
				log_fatal ("Could not read entry from '%s': %s",
						gff->filename, err_msg);
			else
				return 0;
		}

	len = strlen (line);

	while (line[len - 1] != '\n')
		{
			size_t old_size = gff->buf_size;
			gff->buf_size = gff_buf_expand ((void **) &gff->buf,
					sizeof (char), gff->buf_size, GFF_BUFSIZ);

			line = gzgets (gff->fp, &gff->buf[old_size - 1],
					gff->buf_size - old_size + 1);
			if (line == NULL)
				{
					err_msg = gzerror (gff->fp, &rc);
					if (rc != Z_OK)
						log_fatal ("Could not read entry from '%s': %s",
								gff->filename, err_msg);
					else
						return 0;
				}

			len = strlen (line);
		}

	return 1;
}

static void
gff_get_header (GffFile *gff)
{
	char *header = NULL;
	int rc = 0;

	// Catch header '##' and ignore blank lines
	while ((rc = gff_getline (gff)) &&
			(gff->buf[0] == '#' || gff->buf[0] == '\n'))
		{
			if (gff->buf[0] == '#')
				header = xstrdup_concat (header, gff->buf);
			gff->num_line++;
		}

	// Reached end of file
	if (!rc)
		gff->eof = 1;
	else
		gff->num_line++;

	gff->header = chomp (header);
}

GffFile *
gff_open_for_reading (const char *path)
{
	assert (path != NULL);

	GffFile *gff = NULL;
	gzFile fp = NULL;

	gff = xcalloc (1, sizeof (GffFile));
	fp = gzopen (path, "rb");

	if (fp == NULL)
		log_errno_fatal ("Could not open '%s' for reading", path);

	gff->fp = fp;
	gff->filename = xstrdup (path);

	gff->buf = xcalloc (GFF_BUFSIZ, sizeof (char));
	gff->buf_size = GFF_BUFSIZ;

	gff_get_header (gff);
	return gff;
}

void
gff_close (GffFile *gff)
{
	if (gff == NULL)
		return;

	int rc;

	rc = gzclose (gff->fp);
	if (rc != Z_OK)
		log_fatal ("Could not close file '%s': %s", gff->filename,
				gzerror (gff->fp, &rc));

	xfree ((void *) gff->filename);
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
	char *saveptr1, *saveptr2;

	token = subtoken = saveptr1 = saveptr2 = NULL;

	// end of file
	if (gff->eof)
		return 0;

	gff->buf = chomp (gff->buf);
	entry->num_line = gff->num_line;

	token = strtok_r (gff->buf, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'seqname' (field 1) at line %zu", entry->num_line);

	entry->seqname_size = gff_entry_set (&entry->seqname,
			entry->seqname_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'source' (field 2) at line %zu", entry->num_line);

	entry->source_size = gff_entry_set (&entry->source,
			entry->source_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'feature' (field 3) at line %zu", entry->num_line);

	entry->feature_size = gff_entry_set (&entry->feature,
			entry->feature_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'start' (field 4) at line %zu", entry->num_line);

	entry->start = atol (token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'end' (field 5) at line %zu", entry->num_line);

	entry->end = atol (token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'score' (field 6) at line %zu", entry->num_line);

	entry->score = strcmp (token, ".")
		? atof (token)
		: -1;

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'strand' (field 7) at line %zu", entry->num_line);

	entry->strand = *token;

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'frame' (field 8) at line %zu", entry->num_line);

	entry->frame = strcmp (token, ".")
		? atoi (token)
		: -1 ;

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->num_attributes = i = 0;

	if (token != NULL)
		{
			subtoken = strtok_r (token, ";= ", &saveptr2);
			while (subtoken != NULL)
				{
					subtoken = trimc (trim (subtoken), '"');

					if ((i + 1) > entry->attributes_size)
						entry->attributes_size = gff_buf_expand ((void **) &entry->attributes,
								sizeof (GffAttribute), entry->attributes_size, GFF_ATTRSIZ);

					entry->attributes[i].key_size = gff_entry_set (&entry->attributes[i].key,
							entry->attributes[i].key_size, subtoken);

					subtoken = strtok_r (NULL, ";= ", &saveptr2);
					subtoken = trimc (trim (subtoken), '"');
					if (subtoken == NULL)
						log_fatal ("missing value for attribute %d (%s) at line %zu",
								i, entry->attributes[i].key, entry->num_line);

					entry->attributes[i].value_size = gff_entry_set (&entry->attributes[i].value,
							entry->attributes[i].value_size, subtoken);

					subtoken = strtok_r (NULL, ";= ", &saveptr2);
					i++;
				}
		}

	entry->num_attributes = i;

	/*
	* get next entry
	* ignore comments and blank lines
	*/
	while ((rc = gff_getline (gff)) &&
			(gff->buf[0] == '#' || gff->buf[0] == '\n'))
		gff->num_line++;

	// Reached end of file
	if (!rc)
		gff->eof = 1;
	else
		gff->num_line++;

	return 1;
}

static void
regex_free (regex_t *reg)
{
	if (reg == NULL)
		return;

	regfree (reg);
	xfree (reg);
}

GffFilter *
gff_filter_new (const char *feature)
{
	assert (feature != NULL);

	GffFilter *filter = xcalloc (1, sizeof (GffFilter));

	filter->feature = xstrdup (feature);
	filter->hard_attributes = hash_new (xfree,
			(DestroyNotify) set_free);
	filter->soft_attributes = hash_new (xfree,
			(DestroyNotify) set_free);
	filter->values_regex = hash_new (NULL,
			(DestroyNotify) regex_free);

	return filter;
}

void
gff_filter_free (GffFilter *filter)
{
	if (filter == NULL)
		return;

	xfree ((void *) filter->feature);
	hash_free (filter->hard_attributes);
	hash_free (filter->soft_attributes);
	hash_free (filter->values_regex);

	xfree (filter);
}

static inline void
gff_filter_insert_value_regex (GffFilter *filter,
		const char *value)
{
	regex_t *reg = xcalloc (1, sizeof (regex_t));

	if (regcomp (reg, value, REG_EXTENDED|REG_NOSUB) != 0)
		log_fatal ("Error regcomp for '%s'", value);

	hash_insert (filter->values_regex, value, reg);
}

static inline int
gff_filter_match_value_regex (const GffFilter *filter,
		const char *value_key, const char *value)
{
	regex_t *reg = hash_lookup (filter->values_regex, value_key);
	return reg != NULL
		&& regexec (reg, value, 0, (regmatch_t *) NULL, 0) == 0;
}

void
gff_filter_insert_hard_attribute (GffFilter *filter,
		const char *key, const char *value)
{
	assert (filter != NULL && key != NULL && value != NULL);

	Set *values = hash_lookup (filter->hard_attributes, key);

	if (values == NULL)
		{
			values = set_new_full (str_hash, str_equal, xfree);
			hash_insert (filter->hard_attributes, xstrdup (key), values);
		}

	if (set_insert (values, xstrdup (value)))
		gff_filter_insert_value_regex (filter, value);
}

void
gff_filter_insert_soft_attribute (GffFilter *filter,
		const char *key, const char *value)
{
	assert (filter != NULL && key != NULL && value != NULL);

	Set *values = hash_lookup (filter->soft_attributes, key);

	if (values == NULL)
		{
			values = set_new_full (str_hash, str_equal, xfree);
			hash_insert (filter->soft_attributes, xstrdup (key), values);
		}

	if (set_insert (values, xstrdup (value)))
		gff_filter_insert_value_regex (filter, value);
}

static inline int
gff_filter_lookup_hard_attributes (const GffFilter *filter, const GffEntry *entry)
{
	HashIter iter;
	ListElmt *cur = NULL;

	const Set *values = NULL;
	const char *key = NULL;
	const char *value = NULL;

	hash_iter_init (&iter, filter->hard_attributes);
	while (hash_iter_next (&iter, (void **) &key, (void **) &values))
		{
			value = gff_attribute_find (entry, key);
			if (value == NULL)
				return 0;

			cur = list_head (set_list (values));
			for (; cur != NULL; cur = list_next (cur))
				if (!gff_filter_match_value_regex (filter, list_data (cur), value))
					return 0;
		}

	return 1;
}

static inline int
gff_filter_lookup_soft_attributes (const GffFilter *filter, const GffEntry *entry)
{
	HashIter iter;
	ListElmt *cur = NULL;

	const Set *values = NULL;
	const char *key = NULL;
	const char *value = NULL;

	hash_iter_init (&iter, filter->soft_attributes);
	while (hash_iter_next (&iter, (void **) &key, (void **) &values))
		{
			value = gff_attribute_find (entry, key);
			if (value == NULL)
				continue;

			cur = list_head (set_list (values));
			for (; cur != NULL; cur = list_next (cur))
				if (gff_filter_match_value_regex (filter, list_data (cur), value))
					return 1;
		}

	return 0;
}

static inline int
gff_filter_lookup (const GffFilter *filter, const GffEntry *entry)
{
	// Return 0 if there is a feature to be filtered
	// and it does not match
	if (filter->feature != NULL
			&& strcmp (filter->feature, entry->feature))
		return 0;

	// Filter hard attributes returns 1 if all tests succeed
	// or if there is no hard attributes
	if (!gff_filter_lookup_hard_attributes (filter, entry))
		return 0;

	// Soft filtering returns 1 if at least one test succeed
	// and returns 0 if there is no soft attributes. So, test
	// now for soft_attributes values and return 1 if none
	// was passed. It also returns 1 if no filter was set, and
	// this way, works as a regular 'gff_read'
	if (hash_size (filter->soft_attributes) == 0)
		return 1;

	// Finally, return the soft attributes filtering if all
	// above succeed or was not set at all
	return gff_filter_lookup_soft_attributes (filter, entry);
}

int
gff_read_filtered (GffFile *gff, GffEntry *entry, const GffFilter *filter)
{
	assert (filter != NULL);

	while (gff_read (gff, entry))
		if (gff_filter_lookup (filter, entry))
			return 1;

	return 0;
}
