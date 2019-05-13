#include "config.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "gff.h"

#define GFF_BUFSIZ  128
#define GFF_ATTRSIZ 8

GffFile *
gff_open (const char *path, const char *mode)
{
	assert (path != NULL && mode != NULL
			&& *mode != 'w' && *mode != 'a');

	GffFile *gff = NULL;
	gzFile fp = NULL;

	gff = xcalloc (1, sizeof (GffFile));
	fp = gzopen (path, mode);

	if (fp == NULL)
		{
			if (*mode == 'w' || *mode == 'a')
				log_errno_fatal ("Could not open '%s' for writing", path);
			else
				log_errno_fatal ("Could not open '%s' for reading", path);
		}

	gff->fp = fp;
	gff->filename = xstrdup (path);
	gff->buf = xcalloc (GFF_BUFSIZ, sizeof (char));
	gff->buf_size = GFF_BUFSIZ;

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
	xfree (gff->buf);
	xfree (gff);
}

GffEntry *
gff_entry_new (void)
{
	GffEntry *entry = xcalloc (1, sizeof (GffEntry));
	return entry;
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

int
gff_read (GffFile *gff, GffEntry *entry)
{
	assert (gff != NULL && entry != NULL);

	int rc, i;
	char *token, *subtoken;
	char *saveptr1, *saveptr2;

	// ignore comments and header
	while ((rc = gff_getline (gff)) && gff->buf[0] == '#')
		;

	// end of file
	if (!rc)
		return 0;

	gff->buf = chomp (gff->buf);

	token = strtok_r (gff->buf, "\t", &saveptr1);
	entry->seqname_size = gff_entry_set (&entry->seqname,
			entry->seqname_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->source_size = gff_entry_set (&entry->source,
			entry->source_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->feature_size = gff_entry_set (&entry->feature,
			entry->feature_size, token);

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->start = atol (token);

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->end = atol (token);

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->score = strcmp (token, ".")
		? atof (token)
		: -1;

	token = strtok_r (NULL, "\t", &saveptr1);
	entry->strand = *token;

	token = strtok_r (NULL, "\t", &saveptr1);
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

					entry->attributes[i].value_size = gff_entry_set (&entry->attributes[i].value,
							entry->attributes[i].value_size, subtoken);

					subtoken = strtok_r (NULL, ";= ", &saveptr2);
					i++;
				}
		}

	entry->num_attributes = i;
	gff->num_line++;

	return 1;
}
