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
		log_fatal ("missing 'seqname' (field 1) at line %zu", entry->num_line);

	entry->seqname_size = entry_set (&entry->seqname,
			entry->seqname_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'source' (field 2) at line %zu", entry->num_line);

	entry->source_size = entry_set (&entry->source,
			entry->source_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'feature' (field 3) at line %zu", entry->num_line);

	entry->feature_size = entry_set (&entry->feature,
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
						log_fatal ("missing value for attribute %d (%s) at line %zu",
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

static void
regex_free (regex_t *reg)
{
	if (reg == NULL)
		return;

	regfree (reg);
	xfree (reg);
}

GffFilter *
gff_filter_new (void)
{
	GffFilter *filter = xcalloc (1, sizeof (GffFilter));

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

void
gff_filter_insert_feature (GffFilter *filter, const char *feature)
{
	xfree ((void *) filter->feature);
	filter->feature = feature != NULL
		? xstrdup (feature)
		: NULL;
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
