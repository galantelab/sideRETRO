#include "config.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "bed.h"

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
	while ((rc = gz_getline (bed->gz, &bed->buf, &bed->buf_size))
			&& (bed->buf[0] == '\n' || bed_is_header (bed)))
		{
			if (bed->buf[0] != '\n')
				header = xstrdup_concat (header, bed->buf);
		}

	// Reached end of file
	if (!rc)
		bed->eof = 1;

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
	entry->num_line = gz_get_num_line (bed->gz);

	token = strtok_r (bed->buf, "\t ", &saveptr1);
	if (token == NULL)
		log_fatal ("missing 'chrom' (field 1) at line %zu", entry->num_line);

	entry->chrom_size = entry_set (&entry->chrom,
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

	entry->name_size = entry_set (&entry->name,
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

			size = buf_expand ((void **) &entry->block_sizes,
					sizeof (int), entry->num_blocks, entry->block_count);

			size = buf_expand ((void **) &entry->block_starts,
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
	while ((rc = gz_getline (bed->gz, &bed->buf, &bed->buf_size))
			&& (bed->buf[0] == '\n' || bed_is_header (bed)))
		;

	// Reached end of file
	if (!rc)
		bed->eof = 1;

	return 1;
}

BedFile *
bed_open_for_reading (const char *path)
{
	assert (path != NULL);

	BedFile *bed = NULL;

	bed = xcalloc (1, sizeof (BedFile));
	bed->gz = gz_open_for_reading (path);

	bed_get_header (bed);

	return bed;
}

void
bed_close (BedFile *bed)
{
	if (bed == NULL)
		return;

	gz_close (bed->gz);

	xfree ((void *) bed->header);
	xfree (bed->buf);

	xfree (bed);
}
