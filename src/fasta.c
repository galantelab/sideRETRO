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

#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "log.h"
#include "fasta.h"

#define GZ_BUFSIZ             1 << 20   //   1MiB
#define FASTA_SEQ_BUFSIZ    100 << 20   // 100MiB
#define FASTA_CONTIG_BUFSIZ 128

/*
 * The first line of a FASTA file - after comments and
 * blank lines - is the header of a contig: '>contig'
 */
int
fasta_check_format (const char *path)
{
	assert (path != NULL);

	char *line = NULL;
	int rc = 0;

	line = gz_first_line (path, ";");
	rc = line != NULL && line[0] == '>' && line[1] != '\0';

	xfree (line);
	return rc;
}

FastaEntry *
fasta_entry_new (void)
{
	FastaEntry *entry = xcalloc (1, sizeof (FastaEntry));

	entry->sequence = string_sized_new (FASTA_SEQ_BUFSIZ);
	entry->contig = string_sized_new (FASTA_CONTIG_BUFSIZ);

	return entry;
}

void
fasta_entry_free (FastaEntry *entry)
{
	if (entry == NULL)
		return;

	string_free (entry->sequence, 1);
	string_free (entry->contig, 1);

	xfree (entry);
}

FastaFile *
fasta_open_for_reading (const char *path)
{
	assert (path != NULL);

	FastaFile *fasta = NULL;

	fasta = xcalloc (1, sizeof (FastaFile));
	fasta->gz = gz_open_for_reading (path);

	gzbuffer (gz_get_fp (fasta->gz), GZ_BUFSIZ);

	return fasta;
}

void
fasta_close (FastaFile *fasta)
{
	if (fasta == NULL)
		return;

	gz_close (fasta->gz);
	xfree (fasta->buf);
	xfree (fasta);
}

static inline int
fasta_get_contig (FastaFile *fasta, FastaEntry *entry)
{
	int rc = 0;
	char *token = NULL;
	char *saveptr = NULL;

	if (fasta->buf == NULL || fasta->buf[0] != '>')
		{
			// Skip comments and blank lines
			while ((rc = gz_getline (fasta->gz, &fasta->buf, &fasta->buf_size))
					&& (fasta->buf[0] == ';' || fasta->buf[0] == '\n'))
				;

			entry->num_line = gz_get_num_line (fasta->gz);

			// Empty file?
			if (!rc)
				return 0;
		}

	fasta->buf = chomp (fasta->buf);
	token = strtok_r (fasta->buf, " |", &saveptr);

	if (token == NULL || fasta->buf[0] != '>')
		log_fatal ("Missing 'contig' at line %zu of '%s'",
				entry->num_line, gz_get_filename (fasta->gz));

	string_set (entry->contig, token + 1);

	return 1;
}

static inline int
fasta_get_seq (FastaFile *fasta, FastaEntry *entry)
{
	size_t num_line_old = gz_get_num_line (fasta->gz);

	// Reset
	string_reset (entry->sequence);

	while (gz_getline (fasta->gz, &fasta->buf, &fasta->buf_size))
		{
			fasta->buf = chomp (fasta->buf);

			// Skip comments and blank lines
			if (fasta->buf[0] == ';' || fasta->buf[0] == '\n')
				{
					// Update num_line_old in order to avoid
					// >ID1
					// ;
					// >ID2
					num_line_old++;
					continue;
				}

			if (fasta->buf[0] == '>')
				{
					if ((num_line_old + 1) == gz_get_num_line (fasta->gz))
						log_fatal ("Missing sequence for contig %s at line %zu of '%s'",
								entry->contig->str, num_line_old + 1, gz_get_filename (fasta->gz));
					break;
				}

			string_concat (entry->sequence, fasta->buf);
		}

	entry->num_line = gz_get_num_line (fasta->gz);

	if (num_line_old == entry->num_line)
		log_fatal ("Missing sequence for contig %s at line %zu of '%s'",
				entry->contig->str, entry->num_line, gz_get_filename (fasta->gz));

	return fasta->buf[0] == '>' ? 1 : 0;
}

int
fasta_read (FastaFile *fasta, FastaEntry *entry)
{
	assert (fasta != NULL && entry != NULL);

	// End of file or empty
	if (fasta->eof || !fasta_get_contig (fasta, entry))
		{
			fasta->eof = 1;
			return 0;
		}

	if (!fasta_get_seq (fasta, entry))
		fasta->eof = 1;

	return 1;
}

FastaRegion *
fasta_region_new (const char *chr, long start, long end)
{
	assert (chr != NULL && start > 0 && end >= start);

	FastaRegion *region = xcalloc (1, sizeof (FastaRegion));

	*region = (FastaRegion) {
		.chr   = xstrdup (chr),
		.start = start,
		.end   = end
	};

	return region;
}

void
fasta_region_free (FastaRegion *region)
{
	if (region == NULL)
		return;

	xfree (region->chr);
	xfree (region);
}

char *
fasta_region_key (const char *chr, long start, long end)
{
	char *key = NULL;
	xasprintf (&key, "%s:%li-%li", chr, start, end);
	return key;
}

static void
fasta_region_list_free (List *list)
{
	list_free (list);
}

Hash *
fasta_index_regions (const char *path, const List *regions, Hash **contig_len)
{
	log_trace ("inside '%s'", __func__);
	assert (path != NULL && regions != NULL);

	Hash *idx = NULL;
	Hash *by_chr = NULL;
	List *list = NULL;
	ListElmt *cur = NULL;
	FastaRegion *region = NULL;
	FastaFile *fasta = NULL;
	FastaEntry *entry = NULL;
	char *seq = NULL;
	long *len_alloc = NULL;
	long len = 0;

	/*
	* Read the genome only once and keep just the
	* requested regions: 'chr:start-end' => sequence.
	* Group the regions by contig before
	*/
	by_chr = hash_new_full (str_hash, str_equal, NULL,
			(DestroyNotify) fasta_region_list_free);

	for (cur = list_head (regions); cur != NULL; cur = list_next (cur))
		{
			region = list_data (cur);
			list = hash_lookup (by_chr, region->chr);

			if (list == NULL)
				{
					list = list_new (NULL);
					hash_insert (by_chr, region->chr, list);
				}

			list_append (list, region);
		}

	idx = hash_new_full (str_hash, str_equal, xfree, xfree);

	if (contig_len != NULL)
		*contig_len = hash_new_full (str_hash, str_equal, xfree, xfree);

	fasta = fasta_open_for_reading (path);
	entry = fasta_entry_new ();

	while (fasta_read (fasta, entry))
		{
			if (contig_len != NULL)
				{
					len_alloc = xcalloc (1, sizeof (long));
					*len_alloc = entry->sequence->len;
					hash_insert (*contig_len, xstrdup (entry->contig->str),
							len_alloc);
				}

			list = hash_lookup (by_chr, entry->contig->str);
			if (list == NULL)
				continue;

			for (cur = list_head (list); cur != NULL; cur = list_next (cur))
				{
					region = list_data (cur);

					if (region->end > (long) entry->sequence->len)
						{
							log_warn ("Region %s:%li-%li is outside contig range %zu",
									region->chr, region->start, region->end,
									entry->sequence->len);
							continue;
						}

					len = region->end - region->start + 1;
					seq = xcalloc (len + 1, sizeof (char));
					memcpy (seq, entry->sequence->str + region->start - 1, len);
					upper (seq);

					hash_insert (idx, fasta_region_key (region->chr,
								region->start, region->end), seq);
				}
		}

	if (contig_len != NULL && !hash_size (*contig_len))
		log_fatal ("FASTA file '%s' has no entries", path);

	fasta_entry_free (entry);
	fasta_close (fasta);
	hash_free (by_chr);

	return idx;
}
