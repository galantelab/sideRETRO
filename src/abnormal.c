#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <assert.h>
#include "list.h"
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "abnormal.h"

static volatile long abnormal_acm = 0;

static inline void
push (const bam1_t *align, List *in, List *trash)
{
	if (trash->size)
		{
			ListElmt *head = list_head (trash);

			if (bam_copy1 (list_data (head), align) == NULL)
				log_fatal ("Failed to copy sam alignment");

			list_remove_link (trash, head);
			list_append_link (in, head);
		}
	else
		{
			bam1_t *align_copy = bam_dup1 (align);

			if (align_copy == NULL)
				log_fatal ("Failed to duplicate sam alignment");

			list_append (in, align_copy);
		}
}

static inline void
clean (List *in, List *trash)
{
	ListElmt *cur = list_head (in);

	while (cur != NULL)
		{
			ListElmt *next = list_next (cur);
			list_remove_link (in, cur);
			list_append_link (trash, cur);
			cur = next;
		}
}

static inline int
inside_fragment (const bam1_t *align, const List *in)
{
	return !strcmp (bam_get_qname (align),
			bam_get_qname ((bam1_t *) list_data (list_head (in))));
}

static int
test_sorted_queryname (const bam_hdr_t *hdr)
{
	assert (hdr != NULL);

	char *sorted_by = NULL;
	float version = 0;
	int success = 0;

	sorted_by = xcalloc (hdr->l_text, sizeof (char));

	if (sscanf (hdr->text, "@HD\tVN:%f\tSO:%s", &version, sorted_by) == 2)
		{
			chomp (sorted_by);
			success = !strcmp (sorted_by, "queryname");
		}

	xfree (sorted_by);
	return success;
}

static void
dump_if_abnormal (samFile *out, bam_hdr_t *hdr, List *stack)
{
	ListElmt *cur = NULL;
	bam1_t *align = NULL;
	int is_abnormal = 0;

	for (cur = list_head (stack); cur != NULL; cur = list_next (cur))
		{
			align = list_data (cur);

			/*
			* abnormal alignment must be:
			* - paired-end
			* - mapped
			* - mate mapped
			*/
			if (!(align->core.flag & 1)
					|| (align->core.flag & 0x4)
					|| (align->core.flag & 0x8))
				{
					is_abnormal = 0;
					break;
				}

			/*
			* - supplementary
			* - or reads at different chromosomes
			* - or distance between the pairs is bigger
			*   than ABNORMAL_DISTANCE_CUTOFF
			*/
			if ((align->core.flag & 0x800)
						|| (align->core.tid != align->core.mtid)
						|| (abs (align->core.pos - align->core.mpos) > ABNORMAL_DISTANCE_CUTOFF))
				{
					is_abnormal = 1;
				}
		}

	if (is_abnormal)
		{
			abnormal_acm++;
			for (cur = list_head (stack); cur != NULL; cur = list_next (cur))
				{
					align = list_data (cur);
					if (sam_write1 (out, hdr, align) < 0)
						log_fatal ("Failed to write sam abnormal alignment");
				}
		}
}

void
abnormal_filter (const char *input_file, const char *output_file)
{
	log_trace ("Inside %s", __func__);
	assert (input_file != NULL && output_file != NULL);

	samFile *in = NULL;
	samFile *out = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *align = NULL;
	List *stack = NULL;
	List *cache = NULL;
	int rc = 0;

	in = sam_open (input_file, "rb");
	if (in == NULL)
		log_errno_fatal ("Failed to open '%s' for reading", input_file);

	hdr = sam_hdr_read (in);
	if (hdr == NULL)
		log_fatal ("Failed to read sam header");

	// The SAM/BAM file must be sorted by queryname
	assert (test_sorted_queryname (hdr));

	// Filtered sam file
	out = sam_open (output_file, "wb");
	if (out == NULL)
		log_errno_fatal ("Failed to open '%s' for writing", output_file);

	// Write the header from input
	if (sam_hdr_write (out, hdr) < 0)
		log_fatal ("Failed to write sam header");

	align = bam_init1 ();
	if (align == NULL)
		log_errno_fatal ("Failed to create bam1_t");

	// Keep all reads from the same fragment
	// into the list
	stack = list_new ((DestroyNotify) bam_destroy1);
	cache = list_new ((DestroyNotify) bam_destroy1);

	log_info ("Searching for abnormal alignments into '%s'",
			input_file);

	while (1)
		{
			// Returns -1 or < -1 in case of
			// error
			rc = sam_read1 (in, hdr, align);
			if (rc < 0)
				break;

			if (list_size (stack) && !inside_fragment (align, stack))
				{
					dump_if_abnormal (out, hdr, stack);
					clean (stack, cache);
				}

			push (align, stack, cache);
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_fatal ("Failed to read sam alignment");

	// The stack is filled after the parsing,
	// so it is late in relation to the file loop.
	// Therefore, test the last bunch of
	// alignments
	dump_if_abnormal (out, hdr, stack);

	// Just print the amount of abnormal alignments
	log_info ("Found %lu abnormal alignments", abnormal_acm);

	if (sam_close (in) < 0)
		log_errno_fatal ("Failed to close input stream");

	if (sam_close (out) < 0)
		log_errno_fatal ("Failed to close output stream");

	bam_hdr_destroy (hdr);
	bam_destroy1 (align);

	list_free (stack);
	list_free (cache);
}
