#include "config.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "bed.h"

#define BED_BUFSIZ 128

BedEntry *
bed_entry_new (void)
{
	BedEntry *entry = xcalloc (1, sizeof (BedEntry));
	return entry;
}

void
bed_entry_free (BedEntry *entry)
{
	if (entry == NULL)
		return;

	xfree (entry->chrom);
	xfree (entry->name);
	xfree (entry->block_sizes);
	xfree (entry->block_starts);

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
bed_buf_expand (void **buf, size_t size,
		size_t old_nmemb, size_t length)
{
	size_t final_nmemb = nearest_pow (old_nmemb + length);
	*buf = xrealloc (*buf, size * final_nmemb);
	memset (*buf + old_nmemb * size, 0,
			size * (final_nmemb - old_nmemb));
	return final_nmemb;
}

static size_t
bed_entry_set (char **buf, size_t buf_size, const char *entry)
{
	size_t entry_size = strlen (entry);

	if (entry_size >= buf_size)
		buf_size = bed_buf_expand ((void **) buf, sizeof (char),
				buf_size, entry_size - buf_size + 1);

	*buf = strncpy (*buf, entry, buf_size);
	return buf_size;
}

static int
bed_getline (BedFile *bed)
{
	const char *err_msg = NULL;
	char *line = NULL;
	size_t len = 0;
	int rc = 0;

	line = gzgets (bed->fp, bed->buf, bed->buf_size);
	if (line == NULL)
		{
			err_msg = gzerror (bed->fp, &rc);
			if (rc != Z_OK)
				log_fatal ("Could not read entry from '%s': %s",
						bed->filename, err_msg);
			else
				return 0;
		}

	len = strlen (line);

	while (line[len - 1] != '\n')
		{
			size_t old_size = bed->buf_size;
			bed->buf_size = bed_buf_expand ((void **) &bed->buf,
					sizeof (char), bed->buf_size, BED_BUFSIZ);

			line = gzgets (bed->fp, &bed->buf[old_size - 1],
					bed->buf_size - old_size + 1);
			if (line == NULL)
				{
					err_msg = gzerror (bed->fp, &rc);
					if (rc != Z_OK)
						log_fatal ("Could not read entry from '%s': %s",
								bed->filename, err_msg);
					else
						return 0;
				}

			len = strlen (line);
		}

	return 1;
}

static inline int
bed_is_header (BedFile *bed)
{
	return strstr (bed->buf, "browser") == bed->buf
		|| strstr (bed->buf, "track") == bed->buf;
}

static void
bed_get_header (BedFile *bed)
{
	char *header = NULL;
	int rc = 0;

	// Catch header 'browser', 'track' and ignore blank lines
	while ((rc = bed_getline (bed)) && (bed->buf[0] == '\n'
				|| bed_is_header (bed)))
		{
			if (bed->buf[0] != '\n')
				header = xstrdup_concat (header, bed->buf);
			bed->num_line++;
		}

	// Reached end of file
	if (!rc)
		bed->eof = 1;
	else
		bed->num_line++;

	bed->header = chomp (header);
}

int
bed_read (BedFile *bed, BedEntry *entry)
{
	assert (bed != NULL && entry != NULL);

	int rc, i, num_field;
	char *token, *subtoken;
	char *saveptr1, *saveptr2;

	token = subtoken = saveptr1 = saveptr2 = NULL;
	rc = i = num_field = 0;

	// end of file
	if (bed->eof)
		return 0;

	bed->buf = chomp (bed->buf);
	entry->num_line = bed->num_line;

	token = strtok_r (bed->buf, "\t ", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'chrom' (field 1) at line %zu", entry->num_line);

	entry->chrom_size = bed_entry_set (&entry->chrom,
			entry->chrom_size, token);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'chromStart' (field 2) at line %zu", entry->num_line);

	entry->chrom_start = atol (token);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'chromEnd' (field 3) at line %zu", entry->num_line);

	entry->chrom_end = atol (token);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	entry->name_size = bed_entry_set (&entry->name,
			entry->name_size, token);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	entry->score = atoi (token);
	if (entry->score < 0 || entry->score > 1000)
		log_warn ("field 5 'score' out of bound at line %zu", entry->num_line);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	entry->strand = *token;

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	entry->thick_start = atol (token);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	entry->thick_end = atol (token);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	subtoken = strtok_r (token, ",", &saveptr2);
	for (i = 0; subtoken != NULL; i++)
		{
			if (i < 3)
				{
					entry->rgb[i] = atoi (subtoken);
					if (entry->rgb[i] < 0 || entry->rgb[i] > 255)
						log_warn ("field 9 'rgb' [%d] out of bound at line %zu",
								i, entry->num_line);
				}

			subtoken = strtok_r (NULL, ",", &saveptr2);
		}

	if (i > 3)
		log_warn ("field 9 'rgb' has more than 3 values");

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	entry->block_count = atoi (token);
	
	if (entry->block_count > entry->num_blocks)
		{
			size_t size = 0;

			size = bed_buf_expand ((void **) &entry->block_sizes,
					sizeof (int), entry->num_blocks, entry->block_count);

			size = bed_buf_expand ((void **) &entry->block_starts,
					sizeof (int), entry->num_blocks, entry->block_count);

			entry->num_blocks = size;
		}

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	subtoken = strtok_r (token, ",", &saveptr2);
	for (i = 0; subtoken != NULL; i++)
		{
			if (i < entry->block_count)
				entry->block_sizes[i] = atol (subtoken);

			subtoken = strtok_r (NULL, ",", &saveptr2);
		}

	if (i > entry->block_count)
		log_warn ("field 11 'blockSizes' has more values than 'blockCount' at line %zu",
				entry->num_line);

	num_field++;

	token = strtok_r (NULL, "\t ", &saveptr1);
	if (token == NULL)
		goto EOE;

	subtoken = strtok_r (token, ",", &saveptr2);
	for (i = 0; subtoken != NULL; i++)
		{
			if (i < entry->block_count)
				entry->block_starts[i] = atol (subtoken);

			subtoken = strtok_r (NULL, ",", &saveptr2);
		}

	if (i > entry->block_count)
		log_warn ("field 12 'blockStarts' has more values than 'blockCount' at line %zu",
				entry->num_line);

	num_field++;

EOE:
	// Update the number of fields
	entry->num_field = num_field;

	/*
	* get next entry
	* ignore comments and blank lines
	*/
	while ((rc = bed_getline (bed)) &&
			(bed->buf[0] == '\n' || bed_is_header (bed)))
		bed->num_line++;

	// Reached end of file
	if (!rc)
		bed->eof = 1;
	else
		bed->num_line++;

	return 1;
}

BedFile *
bed_open_for_reading (const char *path)
{
	assert (path != NULL);

	BedFile *bed = NULL;
	gzFile fp = NULL;

	bed = xcalloc (1, sizeof (BedFile));
	fp = gzopen (path, "rb");

	if (fp == NULL)
		log_errno_fatal ("Could not open '%s' for reading", path);

	bed->fp = fp;
	bed->filename = xstrdup (path);

	bed->buf = xcalloc (BED_BUFSIZ, sizeof (char));
	bed->buf_size = BED_BUFSIZ;

	bed_get_header (bed);
	return bed;
}

void
bed_close (BedFile *bed)
{
	if (bed == NULL)
		return;

	int rc;

	rc = gzclose (bed->fp);
	if (rc != Z_OK)
		log_fatal ("Could not close file '%s': %s", bed->filename,
				gzerror (bed->fp, &rc));

	xfree ((void *) bed->filename);
	xfree ((void *) bed->header);
	xfree (bed->buf);
	xfree (bed);
}
