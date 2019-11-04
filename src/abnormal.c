#include "config.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
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
	samFile       *in;
	bam_hdr_t     *hdr;
	bam1_t        *align;
	String        *cigar;
	long           alignment_id;
	long           alignment_acm;
	long           abnormal_acm;
	long           exonic_acm;
};

typedef struct _AbnormalFilter AbnormalFilter;

static void
abnormal_filter_init (AbnormalFilter *argf)
{
	// Open SAM/BAM file
	argf->in = sam_open (argf->sam_file, "rb");
	if (argf->in == NULL)
		log_errno_fatal ("Failed to open '%s' for reading",
				argf->sam_file);

	// Get the header
	argf->hdr = sam_hdr_read (argf->in);
	if (argf->hdr == NULL)
		log_fatal ("Failed to read sam header from '%s'",
				argf->sam_file);

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
}

static void
abnormal_filter_destroy (AbnormalFilter *argf)
{
	if (argf == NULL)
		return;

	if (sam_close (argf->in) < 0)
		log_errno_fatal ("Failed to close input stream for '%s'",
				argf->sam_file);

	bam_hdr_destroy (argf->hdr);
	bam_destroy1 (argf->align);

	string_free (argf->cigar, 1);
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

static inline int
abnormal_classifier (const bam1_t *align, int max_distance,
		int phred_quality, float max_base_freq,
		AbnormalType *type)
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
	if (!(align->core.flag & 0x1)
			|| (align->core.flag & 0x4)
			|| (align->core.flag & 0x8)
			|| (align->core.flag & 0x400)
			|| (align->core.qual < phred_quality)
			|| is_base_overly_freq (align, max_base_freq))
		{
			return 0;
		}

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

	cigar = bam_get_cigar (align);
	argf->cigar = string_clear (argf->cigar);

	for (n_cigar = 0; n_cigar < align->core.n_cigar; n_cigar++)
		argf->cigar = string_concat_printf (argf->cigar, "%d%c",
				bam_cigar_oplen (cigar[n_cigar]),
				bam_cigar_opchr (cigar[n_cigar]));

	qlen = bam_cigar2qlen (align->core.n_cigar, cigar);
	rlen = bam_cigar2rlen (align->core.n_cigar, cigar);
	len = rlen < 1 ? 1 : rlen;

	qname = bam_get_qname (align);

	chr = align->core.tid > -1
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
			log_debug ("Alignment [%li] %s %s:%d overlaps %d exons",
					argf->alignment_id, qname, chr_std,
					align->core.pos + 1, acm);
		}

	log_debug ("Dump abnormal alignment [%li] %s %d %s:%d type %d",
			argf->alignment_id, qname, align->core.flag, chr_std,
			align->core.pos + 1, type);

	db_insert_alignment (argf->alignment_stmt,
			argf->alignment_id, qname, align->core.flag,
			chr_std, align->core.pos + 1, align->core.qual,
			argf->cigar->str, qlen, rlen, chr_std_next,
			align->core.mpos + 1, type, argf->tid);

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

			if (!abnormal_classifier (align, argf->max_distance,
						argf->phred_quality, argf->max_base_freq,
						&rtype))
				return;

			type |= rtype;
		}

	if (type != ABNORMAL_NONE)
		{
			for (cur = list_head (stack); cur != NULL;
					cur = list_next (cur))
				{
					align = list_data (cur);
					dump_alignment (argf, align, type);
					argf->abnormal_acm++;
				}
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
sam_rewind (samFile *in, bam_hdr_t **hdr)
{
	int rc = 0;

	rc = in->is_bgzf
		? bgzf_seek (in->fp.bgzf, 0L, SEEK_SET)
		: hseek (in->fp.hfile, 0L, SEEK_SET);

	if (rc < 0)
		log_errno_fatal ("Failed to rewind SAM/BAM file");

	bam_hdr_destroy (*hdr);
	*hdr = sam_hdr_read (in);

	if (*hdr == NULL)
		log_fatal ("Failed to read sam header in rewinding");
}

static void
purge_blacklist (const char *name, Hash *ids)
{
	if (hash_remove (ids, name))
		log_debug ("Remove blacklisted alignment '%s'",
				name);
}

static void
parse_unsorted_sam (AbnormalFilter *argf)
{
	int rc = 0;
	const char *name = NULL;
	AbnormalType type = 0;
	AbnormalType *type_copy = NULL;
	Hash *abnormal_ids = NULL;
	List *blacklist_ids = NULL;

	// All abnormal alignments are keeped
	// into a hash if the BAM/SAM is not sorted
	// by queryname
	abnormal_ids = hash_new (xfree, xfree);

	// Keep all alignments, whose at least
	// one read not pass the constraints
	blacklist_ids = list_new (xfree);

	log_debug ("Index all fragment ids from '%s'", argf->sam_file);

	while ((rc = sam_read1 (argf->in, argf->hdr, argf->align)) >= 0)
		{
			argf->alignment_acm++;

			if (!abnormal_classifier (argf->align, argf->max_distance,
						argf->phred_quality, argf->max_base_freq, &type))
				{
					name = xstrdup (bam_get_qname (argf->align));
					list_append (blacklist_ids, name);
				}
			else if (type != ABNORMAL_NONE)
				{
					type_copy = hash_lookup (abnormal_ids,
							bam_get_qname (argf->align));

					if (type_copy == NULL)
						{
							name = xstrdup (bam_get_qname (argf->align));
							type_copy = xcalloc (1, sizeof (AbnormalType));
							hash_insert (abnormal_ids, name, type_copy);
						}

					*type_copy |= type;
				}
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_errno_fatal ("Failed to read sam alignment from '%s'",
				argf->sam_file);

	// Read file twice in order to catch all
	// supplementary alignments
	sam_rewind (argf->in, &argf->hdr);

	log_debug ("Remove blacklisted fragments from '%s'",
			argf->sam_file);

	// Remove blacklisted alignments
	list_foreach (blacklist_ids, (Func) purge_blacklist,
			abnormal_ids);

	log_debug ("Catch all indexed abnormal fragments from '%s'",
			argf->sam_file);

	// Get all reads from indexed fragments
	while ((rc = sam_read1 (argf->in, argf->hdr, argf->align)) >= 0)
		{
			type_copy = hash_lookup (abnormal_ids,
					bam_get_qname (argf->align));

			if (type_copy != NULL)
				{
					dump_alignment (argf, argf->align, *type_copy);
					argf->abnormal_acm++;
				}
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_errno_fatal ("Failed to read sam alignment from '%s'",
				argf->sam_file);

	// Clean
	hash_free (abnormal_ids);
	list_free (blacklist_ids);
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

	// Opon SAM/BAM
	// Allocate resources
	abnormal_filter_init (&argf);

	log_info ("Searching for abnormal alignments into '%s'",
			argf.sam_file);

	// Let's make this work!
	if (argf.queryname_sorted)
		{
			log_info ("Parsing 'sorted file' mode");
			parse_sorted_sam (&argf);
		}
	else
		{
			log_info ("Parsing 'unsorted file' mode");
			parse_unsorted_sam (&argf);
		}

	log_info ("Processed %li alignments for '%s'",
			argf.alignment_acm, argf.sam_file);

	// Just print the amount of abnormal alignments
	if (argf.abnormal_acm > 0)
		log_info ("Found %li abnormal alignments for '%s': "
			"%li abnormal alignments fall inside some exonic region (%.2f%%)",
			argf.abnormal_acm, argf.sam_file, argf.exonic_acm,
			(float) (argf.exonic_acm * 100) / argf.abnormal_acm);
	else
		log_info ("File '%s' has no abnormal alignments", argf.sam_file);

	// Cleanup
	abnormal_filter_destroy (&argf);
}
