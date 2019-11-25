#include "config.h"

#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "log.h"
#include "fasta.h"

#define GZ_BUFSIZ             1 << 20   //   1MiB
#define FASTA_BUFSIZ          1 << 10   //   1KiB
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
	gzFile fp = NULL;

	fasta = xcalloc (1, sizeof (FastaFile));
	fp = gzopen (path, "rb");

	if (fp == NULL)
		log_errno_fatal ("Could not open '%s' for reading", path);

	gzbuffer (fp, GZ_BUFSIZ);

	fasta->fp = fp;
	fasta->filename = xstrdup (path);

	fasta->buf = xcalloc (FASTA_BUFSIZ, sizeof (char));
	fasta->buf_alloc = FASTA_BUFSIZ;

	return fasta;
}

void
fasta_close (FastaFile *fasta)
{
	if (fasta == NULL)
		return;

	int rc;

	rc = gzclose (fasta->fp);
	if (rc != Z_OK)
		log_fatal ("Could not close file '%s': %s", fasta->filename,
				gzerror (fasta->fp, &rc));

	xfree ((void *) fasta->filename);
	xfree (fasta->buf);

	xfree (fasta);
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
fasta_buf_expand (void **buf, size_t size,
		size_t old_nmemb, size_t length)
{
	size_t final_nmemb = nearest_pow (old_nmemb + length);
	*buf = xrealloc (*buf, size * final_nmemb);
	memset (*buf + old_nmemb * size, 0,
			size * (final_nmemb - old_nmemb));
	return final_nmemb;
}

static inline int
fasta_getline (FastaFile *fasta)
{
	const char *err_msg = NULL;
	char *line = NULL;
	size_t len = 0;
	int rc = 0;

	line = gzgets (fasta->fp, fasta->buf, fasta->buf_alloc);
	if (line == NULL)
		{
			err_msg = gzerror (fasta->fp, &rc);
			if (rc != Z_OK)
				log_fatal ("Could not read entry from '%s': %s",
						fasta->filename, err_msg);
			else
				return 0;
		}

	len = strlen (line);

	while (line[len - 1] != '\n')
		{
			size_t old_alloc = fasta->buf_alloc;
			fasta->buf_alloc = fasta_buf_expand ((void **) &fasta->buf,
					sizeof (char), fasta->buf_alloc, FASTA_BUFSIZ);

			line = gzgets (fasta->fp, &fasta->buf[old_alloc - 1],
					fasta->buf_alloc - old_alloc + 1);
			if (line == NULL)
				{
					err_msg = gzerror (fasta->fp, &rc);
					if (rc != Z_OK)
						log_fatal ("Could not read entry from '%s': %s",
								fasta->filename, err_msg);
					else
						return 0;
				}

			len = strlen (line);
		}

	return 1;
}

static inline int
fasta_get_contig (FastaFile *fasta, FastaEntry *entry)
{
	int rc = 0;
	char *token = NULL;
	char *saveptr = NULL;

	if (fasta->buf[0] != '>')
		{
			// Skip comments
			while ((rc = fasta_getline (fasta)))
				{
					fasta->num_line++;
					if (fasta->buf[0] != ';')
						break;
				}

			// Empty file?
			if (!rc)
				return 0;
		}

	fasta->buf = chomp (fasta->buf);
	token = strtok_r (fasta->buf, " |", &saveptr);

	if (token == NULL || fasta->buf[0] != '>')
		log_fatal ("missing 'contig' at line %zu",
				fasta->num_line);

	string_set (entry->contig, token + 1);
	entry->num_line = fasta->num_line;

	return 1;
}

static inline int
fasta_get_seq (FastaFile *fasta, FastaEntry *entry)
{
	size_t num_line_old = fasta->num_line;

	// Reset
	string_reset (entry->sequence);

	while (fasta_getline (fasta))
		{
			fasta->buf = chomp (fasta->buf);
			fasta->num_line++;

			// Skip comments
			if (fasta->buf[0] == ';')
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
					if ((num_line_old + 1) == fasta->num_line)
						log_fatal ("missing sequence for contig %s at line %zu",
								entry->contig->str, fasta->num_line);
					break;
				}

			string_concat (entry->sequence, fasta->buf);
		}

	if (num_line_old == fasta->num_line)
		log_fatal ("missing sequence for contig %s at line %zu",
				entry->contig->str, fasta->num_line);

	entry->num_line = fasta->num_line;
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
