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

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "log.h"
#include "fasta.h"

#define GZ_BUFSIZ             1 << 20   //   1MiB
#define FASTA_SEQ_BUFSIZ    100 << 20   // 100MiB
#define FASTA_CONTIG_BUFSIZ 128

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
		log_fatal ("missing 'contig' at line %zu",
				entry->num_line);

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
						log_fatal ("missing sequence for contig %s at line %zu",
								entry->contig->str, num_line_old + 1);
					break;
				}

			string_concat (entry->sequence, fasta->buf);
		}

	entry->num_line = gz_get_num_line (fasta->gz);

	if (num_line_old == entry->num_line)
		log_fatal ("missing sequence for contig %s at line %zu",
				entry->contig->str, entry->num_line);

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
