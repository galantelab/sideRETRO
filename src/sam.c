#include "config.h"

#include "htslib/hts.h"
#include "htslib/hfile.h"
#include <fcntl.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "log.h"
#include "sam.h"

int
sam_to_bam_fp (FILE *fp, const char *output_file)
{
	log_trace ("inside '%s'", __func__);
	assert (fp != NULL && output_file != NULL);

	hFILE *hfile = NULL;
	samFile *in = NULL;
	samFile *out = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *align = NULL;
	int success = 1;
	int rc;
	int fd_copy;

	// Copy the fd in order to cleanup
	// the memory associated with both
	// filehandles
	fd_copy = fcntl (fileno (fp), F_DUPFD, 0);
	if (fd_copy < 0)
		log_errno_fatal ("Failed to duplicate fd");

	// Create a hFILE* from the copied fd
	hfile = hdopen (fd_copy, "rb");
	if (hfile == NULL)
		log_errno_fatal ("Failed to create hFILE");

	// Finally, create a samFile* from
	// a opened hFILE*
	in = hts_hopen (hfile, "", "rb");
	if (in == NULL)
		log_errno_fatal ("Failed to create samFile");

	// Output BAM file
	out = sam_open (output_file, "wb");
	if (out == NULL)
		log_errno_fatal ("Failed to open '%s' for writing",
				output_file);

	// Get the sam header
	hdr = sam_hdr_read (in);
	if (hdr == NULL)
		{
			log_error ("Failed to read sam header");
			success = 0; goto Exit;
		}

	// Write the header from input
	if (sam_hdr_write (out, hdr) < 0)
		{
			log_error ("Failed to write sam header");
			success = 0; goto Exit;
		}

	align = bam_init1 ();
	if (align == NULL)
		log_errno_fatal ("Failed to create bam_init1");

	while (1)
		{
			// Returns -1 or < -1 in case of
			// error
			rc = sam_read1 (in, hdr, align);
			if (rc < 0)
				break;

			// Returns -1 in case of error
			if (sam_write1 (out, hdr, align) < 0)
				{
					log_error ("Failed to write sam alignment");
					success = 0; goto Exit;
				}
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		{
			log_error ("Failed to read sam alignment");
			success = 0; goto Exit;
		}

Exit:

	if (hdr != NULL)
		bam_hdr_destroy (hdr);

	if (align != NULL)
		bam_destroy1 (align);

	if (in != NULL)
		{
			if (sam_close (in) < 0)
				{
					log_errno_error ("Failed to close input stream");
					success = 0;
				}
		}

	if (out != NULL)
		{
			if (sam_close (out) < 0)
				{
					log_errno_error ("Failed to close output stream");
					success = 0;
				}
		}

	return success;
}

int
sam_to_bam (const char *input_file, const char *output_file)
{
	log_trace ("inside '%s'", __func__);
	assert (input_file != NULL);

	int success = 1;
	FILE *fp = NULL;

	fp = xfopen (input_file, "r");
	success = sam_to_bam_fp (fp, output_file);

	xfclose (fp);
	return success;
}

int
sam_test_sorted_order (const bam_hdr_t *hdr, const char *value)
{
	assert (hdr != NULL && value != NULL);

	char *sorted_by = NULL;
	float version = 0;
	int success = 0;

	sorted_by = xcalloc (hdr->l_text, sizeof (char));

	if (sscanf (hdr->text, "@HD\tVN:%f\tSO:%s",
				&version, sorted_by) == 2)
		{
			chomp (sorted_by);
			success = !strcmp (sorted_by, value);
		}

	xfree (sorted_by);
	return success;
}
