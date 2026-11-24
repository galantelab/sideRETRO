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

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sam.h"
#include "list.h"
#include "hash.h"
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "str.h"
#include "abnormal.h"

struct _AbnormalFilter
{
	int            tid;
	int            inc_step;
	const char    *sam_file;
	ExonTree      *exon_tree;
	ChrStd        *cs;
	sqlite3_stmt  *alignment_stmt;
	int            phred_quality;
	int            queryname_sorted;
	int            max_distance;
	float          max_base_freq;
	int            either;
	float          exon_frac;
	float          alignment_frac;
	int            hts_threads;
	AbnormalStat  *stat;
	samFile       *in;
	bam_hdr_t     *hdr;
	bam1_t        *align;
	String        *cigar;
	char          *seq;
	int            seq_alloc;
	long           alignment_id;
	long           alignment_acm;
	long           abnormal_acm;
	long           exonic_acm;
	long           stored_acm;
};

typedef struct _AbnormalFilter AbnormalFilter;

static void
abnormal_filter_open (AbnormalFilter *argf)
{
	// Open SAM/BAM/CRAM file
	argf->in = sam_open (argf->sam_file, "rb");
	if (argf->in == NULL)
		log_errno_fatal ("Failed to open '%s' for reading",
				argf->sam_file);

	// Use the spare threads to decompress BAM/CRAM
	if (argf->hts_threads > 0
			&& hts_set_threads (argf->in, argf->hts_threads) < 0)
		log_error ("Failed to set %d decompression threads for '%s'",
				argf->hts_threads, argf->sam_file);

	// Get the header
	argf->hdr = sam_hdr_read (argf->in);
	if (argf->hdr == NULL)
		log_fatal ("Failed to read sam header from '%s'",
				argf->sam_file);
}

static void
abnormal_filter_init (AbnormalFilter *argf)
{
	// Open SAM/BAM/CRAM file and get the header
	abnormal_filter_open (argf);

	// Alloc alignment fields struct
	argf->align = bam_init1 ();
	if (argf->align == NULL)
		log_errno_fatal ("Failed to create bam1_t for '%s'",
				argf->sam_file);

	// If the user did not say the file is sorted,
	// then test if the header ponts to that
	if (!argf->queryname_sorted
			&& sam_test_sorted_order (argf->hdr, "queryname"))
		argf->queryname_sorted = 1;

	// Alloc String cigar to avoid
	// free and realloc often
	argf->cigar = string_sized_new (128);

	// Alloc the sequence buffer for the
	// same reason. They grow if needed
	argf->seq_alloc = 256;
	argf->seq = xcalloc (argf->seq_alloc, sizeof (char));

	// Init alignment_id to its thread id
	// Whenever it is needed to update its value,
	// sum the number of threads - in order to
	// avoid database insertion chocking and
	// constraints
	argf->alignment_id = argf->tid;

	// Zero the accumulators
	argf->alignment_acm = 0;
	argf->abnormal_acm = 0;
	argf->exonic_acm = 0;
	argf->stored_acm = 0;
}

static void
abnormal_filter_destroy (AbnormalFilter *argf)
{
	if (argf == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	if (sam_close (argf->in) < 0)
		log_errno_fatal ("Failed to close input stream for '%s'",
				argf->sam_file);

	bam_hdr_destroy (argf->hdr);
	bam_destroy1 (argf->align);

	string_free (argf->cigar, 1);

	xfree (argf->seq);
}

static int
is_base_overly_freq (const bam1_t *align, float freq)
{
	uint8_t *seq = NULL;
	int bases[5] = {};
	int i = 0;

	seq = bam_get_seq (align);

	for (i = 0; i < align->core.l_qseq; i++)
		bases[seq_nt16_int[bam_seqi (seq, i)]]++;

	for (i = 0; i < 5; i++)
		if (((float) bases[i] / align->core.l_qseq) > freq)
			return 1;

	return 0;
}

int
abnormal_quality_control (const bam1_t *align, int phred_quality,
		float max_base_freq)
{
	/*
	* abnormal alignment must be:
	* - paired-end
	* - mapped
	* - mate mapped
	* - not a duplication
	* - phred-quality
	* - base frequency
	*/
	return (align->core.flag & 0x1)
		&& !(align->core.flag & 0x4)
		&& !(align->core.flag & 0x8)
		&& !(align->core.flag & 0x400)
		&& (align->core.qual >= phred_quality)
		&& !is_base_overly_freq (align, max_base_freq);
}

static inline int
abnormal_classifier (const bam1_t *align, int max_distance,
		int phred_quality, float max_base_freq,
		AbnormalType *type)
{
	if (!abnormal_quality_control (align, phred_quality, max_base_freq))
		return 0;

	// Default value
	*type = ABNORMAL_NONE;

	/*
	* - reads at different chromosomes
	* - or distance between the pairs is bigger
	*   than max_distance
	*/
	if (align->core.tid != align->core.mtid)
		{
			*type |= ABNORMAL_CHROMOSOME;
		}
	else if (abs (align->core.pos - align->core.mpos)
			> max_distance)
		{
			*type |= ABNORMAL_DISTANCE;
		}

	/*
	* - supplementary
	*/
	if (align->core.flag & 0x800)
		{
			*type |= ABNORMAL_SUPPLEMENTARY;
		}

	return 1;
}

/*
 * The base qualities are not stored: No step reads them
 * and they took a quarter of the database. The column is
 * kept with the SAM '*' (not stored), so the databases of
 * former versions can still be merged
 */
#define ABNORMAL_NO_QUAL "*"

static void
fill_seq (AbnormalFilter *argf, const bam1_t *align)
{
	uint8_t *seq = NULL;
	int len = 0;
	int i = 0;

	len = align->core.l_qseq;

	// Grow buffers if needed: [len] == '\0'
	if (len + 1 > argf->seq_alloc)
		{
			argf->seq_alloc = len + 1;
			argf->seq = xrealloc (argf->seq, argf->seq_alloc);
		}

	// SAM: '*' means the sequence is not stored
	if (len == 0)
		{
			strcpy (argf->seq, "*");
			return;
		}

	seq = bam_get_seq (align);

	// The sequence is always at the forward strand
	// of the reference - as in the SAM format
	for (i = 0; i < len; i++)
		argf->seq[i] = seq_nt16_str[bam_seqi (seq, i)];

	argf->seq[len] = '\0';
}

/*
 * The same exon lookup of dump_alignment, without
 * dumping: Does the alignment overlap any exon?
 */
static int
overlaps_exon (AbnormalFilter *argf, const bam1_t *align)
{
	const char *chr = NULL;
	const char *chr_std = NULL;
	int rlen = 0;

	chr = align->core.tid > -1 // LCOV_EXCL_LINE: htslib unmaps any read without chromosome
		? argf->hdr->target_name[align->core.tid]
		: "*";

	chr_std = chr_std_lookup (argf->cs, chr);
	rlen = bam_cigar2rlen (align->core.n_cigar, bam_get_cigar (align));

	return exon_tree_lookup (argf->exon_tree, chr_std,
			align->core.pos + 1, align->core.pos + (rlen < 1 ? 1 : rlen),
			argf->exon_frac, argf->alignment_frac, argf->either) > 0;
}

/*
 * A supplementary alignment below the minimum mapping quality
 * or of low complexity. The clip of a read crossing the poly-A
 * junction is often aligned as a short supplementary piece at a
 * random place: A run of A bases, with MAPQ 0. The fragment would
 * lose all its alignments with it - the read at the insertion
 * site included. So it is left out: It does not fail its
 * fragment, but it is not stored either
 */
static inline int
abnormal_left_out (const AbnormalFilter *argf, const bam1_t *align)
{
	return (align->core.flag & 0x800)
		&& (align->core.qual < argf->phred_quality
				|| is_base_overly_freq (align, argf->max_base_freq));
}

/*
 * Classify an alignment of a fragment: A left out alignment
 * passes the other checks, and its type still marks the
 * fragment - as a split read, for instance
 */
static inline int
abnormal_fragment_classifier (AbnormalFilter *argf, const bam1_t *align,
		AbnormalType *type)
{
	const int left_out = abnormal_left_out (argf, align);

	return abnormal_classifier (align, argf->max_distance,
			left_out ? 0 : argf->phred_quality,
			left_out ? 1.0 : argf->max_base_freq, type);
}

static void
dump_alignment (AbnormalFilter *argf, const bam1_t *align,
		AbnormalType type)
{
	uint32_t *cigar = NULL;
	int n_cigar = 0;
	const char *chr = NULL;
	const char *chr_std = NULL;
	const char *chr_next = NULL;
	const char *chr_std_next = NULL;
	const char *qname = NULL;
	int qlen = 0;
	int rlen = 0;
	int len = 0;
	int acm = 0;
	int soft_clip = 0;

	cigar = bam_get_cigar (align);
	argf->cigar = string_clear (argf->cigar);

	for (n_cigar = 0; n_cigar < align->core.n_cigar; n_cigar++)
		{
			argf->cigar = string_concat_printf (argf->cigar, "%d%c",
					bam_cigar_oplen (cigar[n_cigar]),
					bam_cigar_opchr (cigar[n_cigar]));

			if (bam_cigar_op (cigar[n_cigar]) == BAM_CSOFT_CLIP)
				soft_clip = 1;
		}

	qlen = bam_cigar2qlen (align->core.n_cigar, cigar);
	rlen = bam_cigar2rlen (align->core.n_cigar, cigar);
	len = rlen < 1 ? 1 : rlen;

	qname = bam_get_qname (align);

	/*
	* Get the sequence as string - only for soft clipped
	* alignments: The clipped bases are the only ones read
	* (breakpoints, poly-A tail and the retrocopy ends).
	* Otherwise, SAM '*' (not stored)
	*/
	if (soft_clip)
		fill_seq (argf, align);
	else
		strcpy (argf->seq, "*");

	chr = align->core.tid > -1 // LCOV_EXCL_LINE: htslib unmaps any read without chromosome
		? argf->hdr->target_name[align->core.tid]
		: "*";

	chr_next = align->core.mtid > -1
		? argf->hdr->target_name[align->core.mtid]
		: "*";

	// Standardize chromosomes
	chr_std = chr_std_lookup (argf->cs, chr);
	chr_std_next = chr_std_lookup (argf->cs, chr_next);

	// Dump overlapping exon with alignment
	acm = exon_tree_lookup_dump (argf->exon_tree, chr_std,
			align->core.pos + 1, align->core.pos + len,
			argf->exon_frac, argf->alignment_frac,
			argf->either, argf->alignment_id);

	if (acm > 0)
		{
			type |= ABNORMAL_EXONIC;
			argf->exonic_acm++;
			log_debug ("Alignment [%li] %s %s:%li overlaps %d exons",
					argf->alignment_id, qname, chr_std,
					(long int) align->core.pos + 1, acm);
		}

	log_debug ("Dump abnormal alignment [%li] %s %d %s:%li type %d",
			argf->alignment_id, qname, align->core.flag, chr_std,
			(long int) align->core.pos + 1, type);

	db_insert_alignment (argf->alignment_stmt,
			argf->alignment_id, qname, align->core.flag,
			chr_std, align->core.pos + 1, align->core.qual,
			argf->cigar->str, qlen, rlen, chr_std_next,
			align->core.mpos + 1, type, argf->tid,
			argf->seq, ABNORMAL_NO_QUAL);

	argf->stored_acm++;

	// sum the number of files - in order to
	// avoid database insertion chocking and
	// constraints
	argf->alignment_id += argf->inc_step;
}

static inline void
dump_stack_if_abnormal (AbnormalFilter *argf, const List *stack)
{
	const ListElmt *cur = NULL;
	const bam1_t *align = NULL;
	AbnormalType rtype = ABNORMAL_NONE;
	AbnormalType type = ABNORMAL_NONE;

	for (cur = list_head (stack); cur != NULL;
			cur = list_next (cur))
		{
			align = list_data (cur);

			if (!abnormal_fragment_classifier (argf, align, &rtype))
				return;

			type |= rtype;
		}

	if (type == ABNORMAL_NONE)
		return;

	argf->abnormal_acm += stack->size;

	/*
	* A fragment without any alignment overlapping an exon
	* is never used: The clustering needs a mate at the
	* parental gene. Do not store it
	*/
	for (cur = list_head (stack); cur != NULL;
			cur = list_next (cur))
		if (!abnormal_left_out (argf, list_data (cur))
				&& overlaps_exon (argf, list_data (cur)))
			break;

	if (cur == NULL)
		return;

	for (cur = list_head (stack); cur != NULL;
			cur = list_next (cur))
		{
			align = list_data (cur);
			if (!abnormal_left_out (argf, align))
				dump_alignment (argf, align, type);
		}
}

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

static void
parse_sorted_sam (AbnormalFilter *argf)
{
	// Keep all reads from the same fragment
	// into the list
	List *stack = list_new ((DestroyNotify) bam_destroy1);
	List *cache = list_new ((DestroyNotify) bam_destroy1);

	int rc = 0;

	while ((rc = sam_read1 (argf->in, argf->hdr, argf->align)) >= 0)
		{
			argf->alignment_acm++;

			if (list_size (stack)
					&& !inside_fragment (argf->align, stack))
				{
					dump_stack_if_abnormal (argf, stack);
					clean (stack, cache);
				}

			push (argf->align, stack, cache);
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_errno_fatal ("Failed to read sam alignment from '%s'",
				argf->sam_file);

	// The stack is filled after the parsing,
	// so it is late in relation to the file loop.
	// Therefore, test the last bunch of
	// alignments
	dump_stack_if_abnormal (argf, stack);

	// Clean
	list_free (stack);
	list_free (cache);
}

static void
sam_rewind (AbnormalFilter *argf)
{
	if (sam_close (argf->in) < 0)
		log_errno_fatal ("Failed to rewind '%s'",
				argf->sam_file);

	bam_hdr_destroy (argf->hdr);

	// Open again, get the header and
	// set the decompression threads
	abnormal_filter_open (argf);
}

/*
 * For the unsorted mode, each abnormal fragment keeps
 * the type of the abnormality, if some alignment
 * failed the classifier and the index (into the file)
 * of the first abnormal alignment
 */
struct _AbnormalFragment
{
	AbnormalType   type;
	int            failed;
	int            exonic;
	long           first;
	List          *buffer;
};

typedef struct _AbnormalFragment AbnormalFragment;

static AbnormalFragment *
abnormal_fragment_new (AbnormalType type, long first)
{
	AbnormalFragment *frag = xcalloc (1, sizeof (AbnormalFragment));

	frag->type = type;
	frag->first = first;

	return frag;
}

static void
abnormal_fragment_free (AbnormalFragment *frag)
{
	if (frag == NULL)
		return; // LCOV_EXCL_LINE: static: never NULL

	list_free (frag->buffer);
	xfree (frag);
}

static void
abnormal_fragment_fail (AbnormalFragment *frag)
{
	frag->failed = 1;

	// Release the alignments - if any
	list_free (frag->buffer);
	frag->buffer = NULL;
}

static void
parse_unsorted_sam (AbnormalFilter *argf)
{
	int rc = 0;
	int pass = 0;
	long idx = 0;
	const char *name = NULL;
	const ListElmt *cur = NULL;
	bam1_t *align_copy = NULL;
	AbnormalType type = 0;
	AbnormalFragment *frag = NULL;
	Hash *abnormal_frags = NULL;

	/*
	* An abnormal fragment must have at least one abnormal
	* alignment and all its alignments must pass the classifier.
	* When the file is not sorted by queryname, the alignments
	* of the same fragment may be anywhere into the file,
	* so we need to read it twice:
	*
	* First reading: index all fragments with an abnormal
	* alignment, saving the position (idx) of the first one.
	* Any alignment that fails the classifier after that
	* position marks the fragment as failed.
	*
	* Second reading: only the alignments before the first
	* abnormal alignment are unknown. Keep them into a buffer
	* until reach the first abnormal alignment. At this point,
	* all failures are known, so dump the buffer and the
	* remaining alignments - or discard them all.
	*
	* A fragment without any alignment overlapping an exon is
	* never used - the clustering needs a mate at the parental
	* gene - so it is not dumped. Its alignments from the first
	* abnormal one on are tested at the first reading, and the
	* ones before it at the second reading, into the buffer.
	*/
	abnormal_frags = hash_new (xfree,
			(DestroyNotify) abnormal_fragment_free);

	log_debug ("Index all abnormal fragments from '%s'", argf->sam_file);

	/*First reading*/
	for (idx = 0; (rc = sam_read1 (argf->in, argf->hdr, argf->align)) >= 0; idx++)
		{
			argf->alignment_acm++;

			pass = abnormal_fragment_classifier (argf, argf->align, &type);

			frag = hash_lookup (abnormal_frags,
					bam_get_qname (argf->align));

			if (frag != NULL)
				{
					// Alignment after the first abnormal one
					if (!pass)
						frag->failed = 1;
					else
						frag->type |= type;
				}
			else if (pass && type != ABNORMAL_NONE)
				{
					name = xstrdup (bam_get_qname (argf->align));
					frag = abnormal_fragment_new (type, idx);
					hash_insert (abnormal_frags, name, frag);
				}

			if (frag != NULL && !frag->failed && !frag->exonic
					&& !abnormal_left_out (argf, argf->align))
				frag->exonic = overlaps_exon (argf, argf->align);
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_errno_fatal ("Failed to read sam alignment from '%s'",
				argf->sam_file);

	// Read file again to filter and dump
	sam_rewind (argf);

	log_debug ("Filter and catch all indexed abnormal fragments from '%s'",
			argf->sam_file);

	/*Second reading*/
	for (idx = 0; (rc = sam_read1 (argf->in, argf->hdr, argf->align)) >= 0; idx++)
		{
			frag = hash_lookup (abnormal_frags,
					bam_get_qname (argf->align));

			if (frag == NULL || frag->failed)
				continue;

			if (idx < frag->first)
				{
					// Unknown alignment: test and keep it. It is not
					// abnormal - otherwise it would be the first - so it
					// is never a left out (supplementary) alignment
					pass = abnormal_fragment_classifier (argf, argf->align, &type);

					if (!pass)
						{
							abnormal_fragment_fail (frag);
							continue;
						}

					align_copy = bam_dup1 (argf->align);
					if (align_copy == NULL)
						log_fatal ("Failed to duplicate sam alignment");

					if (frag->buffer == NULL)
						frag->buffer = list_new ((DestroyNotify) bam_destroy1);

					list_append (frag->buffer, align_copy);

					if (!frag->exonic)
						frag->exonic = overlaps_exon (argf, align_copy);

					continue;
				}

			// Reached the first abnormal alignment: all
			// failures are known - dump the buffer
			if (idx == frag->first && frag->buffer != NULL)
				{
					for (cur = list_head (frag->buffer); cur != NULL;
							cur = list_next (cur))
						{
							if (frag->exonic)
								dump_alignment (argf, list_data (cur), frag->type);
							argf->abnormal_acm++;
						}

					list_free (frag->buffer);
					frag->buffer = NULL;
				}

			// From the first abnormal alignment on, all
			// alignments passed the classifier - otherwise
			// the fragment would be marked as failed. The left
			// out ones are not stored
			if (frag->exonic && !abnormal_left_out (argf, argf->align))
				dump_alignment (argf, argf->align, frag->type);
			argf->abnormal_acm++;
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_errno_fatal ("Failed to read sam alignment from '%s'",
				argf->sam_file);

	// Clean
	hash_free (abnormal_frags);
}

void
abnormal_filter (AbnormalArg *arg)
{
	assert (arg != NULL && arg->sam_file != NULL
			&& arg->alignment_stmt != NULL && arg->exon_tree
			&& arg->cs && arg->tid >= 0 && arg->inc_step > 0
			&& arg->phred_quality >= 0 && arg->max_base_freq > 0);

	AbnormalFilter argf = {};
	memcpy (&argf, arg, sizeof (AbnormalArg));

	// Opon SAM/BAM/CRAM
	// Allocate resources
	abnormal_filter_init (&argf);

	// Let's make this work!
	if (argf.queryname_sorted)
		{
			log_info ("Search for abnormal alignments in '%s' "
					"(grouped by read name: one reading)", argf.sam_file);
			parse_sorted_sam (&argf);
		}
	else
		{
			log_info ("Search for abnormal alignments in '%s' "
					"(not grouped by read name: two readings)", argf.sam_file);
			parse_unsorted_sam (&argf);
		}

	// Just print the amount of abnormal alignments
	if (argf.abnormal_acm > 0)
		log_info ("Found %li abnormal alignments out of %li in '%s': "
			"%li (%.2f%%) overlap exons. Stored %li: The fragments with "
			"an alignment overlapping an exon",
			argf.abnormal_acm, argf.alignment_acm, argf.sam_file, argf.exonic_acm,
			(float) (argf.exonic_acm * 100) / argf.abnormal_acm,
			argf.stored_acm);
	else
		log_warn ("No abnormal alignment found out of %li in '%s'",
				argf.alignment_acm, argf.sam_file);

	// Counts for the final summary
	if (argf.stat != NULL)
		*argf.stat = (AbnormalStat) {
			.alignments = argf.alignment_acm,
			.abnormal   = argf.abnormal_acm,
			.exonic     = argf.exonic_acm,
			.stored     = argf.stored_acm,
			.sorted     = argf.queryname_sorted
		};

	// Cleanup
	abnormal_filter_destroy (&argf);
}
