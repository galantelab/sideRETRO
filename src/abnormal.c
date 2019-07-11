#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sam.h"
#include "list.h"
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "str.h"
#include "abnormal.h"

struct _AbnormalFilter
{
	int            tid;
	int            num_threads;
	const char    *sam_file;
	ExonTree      *exon_tree;
	ChrStd        *cs;
	sqlite3       *db;
	sqlite3_stmt  *alignment_stmt;
	int            either;
	float          node_overlap_frac;
	float          interval_overlap_frac;
	samFile       *in;
	bam_hdr_t     *hdr;
	bam1_t        *align;
	List          *stack;
	List          *cache;
	String        *cigar;
	long           alignment_id;
	long           fragment_acm;
	long           abnormal_acm;
	long           exonic_acm;
};

typedef struct _AbnormalFilter AbnormalFilter;

static void
abnormal_filter_init (AbnormalFilter *argf)
{
	argf->in = sam_open (argf->sam_file, "rb");
	if (argf->in == NULL)
		log_errno_fatal ("Failed to open '%s' for reading",
				argf->sam_file);

	argf->hdr = sam_hdr_read (argf->in);
	if (argf->hdr == NULL)
		log_fatal ("Failed to read sam header from '%s'",
				argf->sam_file);

	// The SAM/BAM file must be sorted by queryname
	assert (sam_test_sorted_order (argf->hdr, "queryname"));

	argf->align = bam_init1 ();
	if (argf->align == NULL)
		log_errno_fatal ("Failed to create bam1_t for '%s'",
				argf->sam_file);

	argf->cigar = string_sized_new (128);

	// Keep all reads from the same fragment
	// into the list
	argf->stack = list_new ((DestroyNotify) bam_destroy1);
	argf->cache = list_new ((DestroyNotify) bam_destroy1);

	// Init alignment_id to its thread id
	// Whenever it is needed to update its value,
	// sum the number of threads - in order to
	// avoid database insertion chocking and
	// constraints
	argf->alignment_id = argf->tid;

	// Zero the accumulators
	argf->fragment_acm = 0;
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

	list_free (argf->stack);
	list_free (argf->cache);

	char *str = NULL;
	str = string_free (argf->cigar, 1);
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
dump_alignment (AbnormalFilter *argf, int type)
{
	ListElmt *cur = NULL;
	bam1_t *align = NULL;
	uint32_t *cigar = NULL;
	int n_cigar = 0;
	const char *chr = NULL;
	const char *chr_std = NULL;
	const char *chr_next = NULL;
	const char *chr_std_next = NULL;
	const char *qname = NULL;
	int qlen = 0;
	int rlen = 0;
	int acm = 0;
	int align_type = 0;

	for (cur = list_head (argf->stack); cur != NULL;
			cur = list_next (cur))
		{
			align = list_data (cur);

			cigar = bam_get_cigar (align);
			argf->cigar = string_clear (argf->cigar);

			for (n_cigar = 0; n_cigar < align->core.n_cigar; n_cigar++)
				argf->cigar = string_concat_printf (argf->cigar, "%d%c",
						bam_cigar_oplen (cigar[n_cigar]),
						bam_cigar_opchr (cigar[n_cigar]));

			qlen = bam_cigar2qlen (align->core.n_cigar, cigar);
			rlen = bam_cigar2rlen (align->core.n_cigar, cigar);

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

			// Reset align_type to type
			align_type = type;

			// Dump overlapping exon with alignment
			acm = exon_tree_lookup_dump (argf->exon_tree, chr_std,
					align->core.pos + 1, align->core.pos + rlen,
					argf->node_overlap_frac, argf->interval_overlap_frac,
					argf->either, argf->alignment_id);

			if (acm > 0)
				{
					align_type |= ABNORMAL_EXONIC;
					argf->exonic_acm++;
					log_debug ("Alignment %s %s:%d overlaps %d exons",
							qname, chr_std, align->core.pos + 1, acm);
				}

			log_debug ("Dump abnormal alignment %s %d %s:%d type %d",
					qname, align->core.flag, chr_std, align->core.pos + 1,
					type);

			db_insert_alignment (argf->db, argf->alignment_stmt,
					argf->alignment_id, qname, align->core.flag,
					chr_std, align->core.pos + 1, align->core.qual,
					argf->cigar->str, qlen, rlen, chr_std_next,
					align->core.mpos + 1, align_type,
					argf->tid);

			// sum the number of threads - in order to
			// avoid database insertion chocking and
			// constraints
			argf->alignment_id += argf->num_threads;
		}
}

static void
dump_if_abnormal (AbnormalFilter *argf)
{
	ListElmt *cur = NULL;
	bam1_t *align = NULL;
	int type = ABNORMAL_NONE;

	argf->fragment_acm++;

	for (cur = list_head (argf->stack); cur != NULL;
			cur = list_next (cur))
		{
			align = list_data (cur);

			/*
			* abnormal alignment must be:
			* - paired-end
			* - mapped
			* - mate mapped
			*/
			if (!(align->core.flag & 0x1)
					|| (align->core.flag & 0x4)
					|| (align->core.flag & 0x8))
				{
					return;
				}

			/*
			* - supplementary
			*/
			if (align->core.flag & 0x800)
				{
					if (!(type & ABNORMAL_SUPPLEMENTARY))
						type |= ABNORMAL_SUPPLEMENTARY;
				}

			/*
			* - reads at different chromosomes
			* - or distance between the pairs is bigger
			*   than ABNORMAL_DISTANCE_CUTOFF
			*/
			if (align->core.tid != align->core.mtid)
				{
					if (!(type & ABNORMAL_CHROMOSOME))
						type |= ABNORMAL_CHROMOSOME;
				}
			else if (abs (align->core.pos - align->core.mpos)
					> ABNORMAL_DISTANCE_CUTOFF)
				{
					if (!(type & ABNORMAL_DISTANCE))
						type |= ABNORMAL_DISTANCE;
				}
		}

	if (type != ABNORMAL_NONE)
		{
			dump_alignment (argf, type);
			argf->abnormal_acm++;
		}
}

void
abnormal_filter (AbnormalArg *arg)
{
	assert (arg != NULL && arg->sam_file != NULL && arg->db != NULL
			&& arg->alignment_stmt != NULL && arg->exon_tree && arg->cs
			&& arg->tid >= 0 && arg->num_threads > 0);

	int rc = 0;
	AbnormalFilter argf = {};
	memcpy (&argf, arg, sizeof (AbnormalArg));

	// Opon SAM/BAM
	// Allocate resources
	abnormal_filter_init (&argf);

	log_info ("Searching for abnormal alignments into '%s'",
			argf.sam_file);

	while (1)
		{
			// Returns -1 or < -1 in case of
			// error
			rc = sam_read1 (argf.in, argf.hdr, argf.align);
			if (rc < 0)
				break;

			if (list_size (argf.stack)
					&& !inside_fragment (argf.align, argf.stack))
				{
					dump_if_abnormal (&argf);
					clean (argf.stack, argf.cache);
				}

			push (argf.align, argf.stack, argf.cache);
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_fatal ("Failed to read sam alignment from '%s'",
				argf.sam_file);

	// The stack is filled after the parsing,
	// so it is late in relation to the file loop.
	// Therefore, test the last bunch of
	// alignments
	dump_if_abnormal (&argf);

	// Just print the amount of abnormal alignments
	if (argf.abnormal_acm > 0)
		log_info ("Found %li abnormal alignments for '%s': "
			"%li abnormal alignments falls inside some exonic region (%.2f%)",
			argf.abnormal_acm, argf.sam_file, argf.exonic_acm,
			(float) (argf.exonic_acm * 100) / argf.abnormal_acm);
	else
		log_info ("File '%s' has no abnormal alignments", argf.sam_file);

	// Cleanup
	abnormal_filter_destroy (&argf);
}
