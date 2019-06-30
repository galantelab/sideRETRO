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
#include "str.h"
#include "hash.h"
#include "gff.h"
#include "ibitree.h"
#include "db.h"
#include "abnormal.h"

struct _AbnormalArg
{
	samFile       *in;
	bam_hdr_t     *hdr;
	bam1_t        *align;
	sqlite3       *db;
	sqlite3_stmt  *exon_stmt;
	sqlite3_stmt  *alignment_stmt;
	sqlite3_stmt  *overlapping_stmt;
	Hash          *tree_idx;
	List          *stack;
	List          *cache;
	String        *cigar;
	int            alignment_id;
	int            fragment_acm;
	int            abnormal_acm;
	int            exonic_acm;
};

typedef struct _AbnormalArg AbnormalArg;

static void
index_dump_gff_file (AbnormalArg *arg, const char *gff_file)
{
	log_trace ("Inside %s", __func__);

	GffFile *gff = gff_open (gff_file, "rb");
	GffEntry *entry = gff_entry_new ();

	Hash *exons_found = hash_new (HASH_MEDIUM_SIZE,
			xfree, NULL);

	IBiTree *tree = NULL;

	int rc = 0;
	int table_id = 0;
	int *alloc_id = NULL;
	int hash_true = 1;

	const char *gene_name = NULL;
	const char *gene_id = NULL;
	const char *exon_id = NULL;
	const char *transcript_type = NULL;
	char strand[1];

	while (gff_read (gff, entry))
		{
			transcript_type = gff_attribute_find (entry, "transcript_type");

			if (!strcmp (entry->feature, "exon")
					&& transcript_type != NULL
					&& !strcmp (transcript_type, "protein_coding"))
				{
					gene_name = gff_attribute_find (entry, "gene_name");
					gene_id = gff_attribute_find (entry, "gene_id");
					exon_id = gff_attribute_find (entry, "exon_id");

					if (gene_name == NULL || gene_id == NULL || exon_id == NULL)
						{
							log_warn ("Missing gene_name|gene_id|exon_id at line %zu",
								entry->num_line);
							continue;
						}

					if (hash_contains (exons_found, exon_id))
						continue;
					else
						hash_insert (exons_found,
								xstrdup (exon_id), &hash_true);

					log_debug ("Index exon from gene '%s' at %s:%zu-%zu", gene_name,
							entry->seqname, entry->start, entry->end);

					strand[0] = entry->strand;
					alloc_id = xcalloc (1, sizeof (int));
					* (int *) alloc_id = ++table_id;

					tree = hash_lookup (arg->tree_idx, entry->seqname);

					if (tree == NULL)
						{
							tree = ibitree_new (xfree);
							rc = hash_insert (arg->tree_idx,
									xstrdup (entry->seqname), tree);
						}

					ibitree_insert (tree, entry->start, entry->end, alloc_id);

					db_insert_exon (arg->db, arg->exon_stmt, table_id,
							gene_name, entry->seqname, entry->start, entry->end,
							strand, gene_id, exon_id);
				}
		}

	hash_free (exons_found);
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

static AbnormalArg *
abnormal_arg_new (const char *sam_file, const char *db_path)
{
	log_trace ("Inside %s", __func__);

	AbnormalArg *arg = xcalloc (1, sizeof (AbnormalArg));

	arg->in = sam_open (sam_file, "rb");
	if (arg->in == NULL)
		log_errno_fatal ("Failed to open '%s' for reading",
				sam_file);

	arg->hdr = sam_hdr_read (arg->in);
	if (arg->hdr == NULL)
		log_fatal ("Failed to read sam header");

	// The SAM/BAM file must be sorted by queryname
	assert (test_sorted_queryname (arg->hdr));

	arg->align = bam_init1 ();
	if (arg->align == NULL)
		log_errno_fatal ("Failed to create bam1_t");

	// Connect to database
	arg->db = db_create (db_path);
	arg->exon_stmt = db_prepare_exon_stmt (arg->db);
	arg->alignment_stmt = db_prepare_alignment_stmt (arg->db);
	arg->overlapping_stmt = db_prepare_overlapping_stmt (arg->db);

	// Increase the cache size to 1GiB
	db_cache_size (arg->db, 100000);

	// Begin transaction to speed up
	db_begin_transaction (arg->db);

	// Keep all reads from the same fragment
	// into the list
	arg->stack = list_new ((DestroyNotify) bam_destroy1);
	arg->cache = list_new ((DestroyNotify) bam_destroy1);

	arg->cigar = string_sized_new (128);

	arg->tree_idx = hash_new (HASH_SMALL_SIZE, xfree,
			(DestroyNotify) ibitree_free);

	return arg;
}

static void
abnormal_arg_free (AbnormalArg *arg)
{
	assert (arg != NULL);

	if (sam_close (arg->in) < 0)
		log_errno_fatal ("Failed to close input stream");

	bam_hdr_destroy (arg->hdr);
	bam_destroy1 (arg->align);

	list_free (arg->stack);
	list_free (arg->cache);

	char *str = NULL;
	str = string_free (arg->cigar, 1);

	hash_free (arg->tree_idx);

	db_end_transaction (arg->db);
	db_finalize (arg->db, arg->exon_stmt);
	db_finalize (arg->db, arg->alignment_stmt);
	db_finalize (arg->db, arg->overlapping_stmt);
	db_close (arg->db);

	xfree (arg);
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
dump_if_overlaps_exon (void *data, void *user_data)
{
	const int *exon_id = data;
	AbnormalArg *arg = user_data;

	log_debug ("Dump overlapping: exon_id %d, alignment_id %d",
			*exon_id, arg->alignment_id);

	db_insert_overlapping (arg->db, arg->overlapping_stmt,
			*exon_id, arg->alignment_id);

	arg->exonic_acm++;
}

static void
dump_alignment (AbnormalArg *arg, int type)
{
	ListElmt *cur = NULL;
	IBiTree *tree = NULL;
	bam1_t *align = NULL;
	uint32_t *cigar = NULL;
	int n_cigar = 0;
	const char *chr = NULL;
	const char *chr_next = NULL;
	const char *qname = NULL;
	int qlen = 0;
	int rlen = 0;
	int acm = 0;

	for (cur = list_head (arg->stack); cur != NULL;
			cur = list_next (cur))
		{
			align = list_data (cur);

			cigar = bam_get_cigar (align);
			arg->cigar = string_clear (arg->cigar);

			for (n_cigar = 0; n_cigar < align->core.n_cigar; n_cigar++)
				arg->cigar = string_concat_printf (arg->cigar, "%d%c",
						bam_cigar_oplen (cigar[n_cigar]),
						bam_cigar_opchr (cigar[n_cigar]));

			qlen = bam_cigar2qlen (align->core.n_cigar, cigar);
			rlen = bam_cigar2rlen (align->core.n_cigar, cigar);

			qname = bam_get_qname (align);

			chr = align->core.tid > -1
				? arg->hdr->target_name[align->core.tid]
				: "*";

			chr_next = align->core.mtid > -1
				? arg->hdr->target_name[align->core.mtid]
				: "*";

			log_debug ("Dump abnormal alignment %s %d %s:%ld",
					qname, align->core.flag, chr, align->core.pos);

			db_insert_alignment (arg->db, arg->alignment_stmt,
					++(arg->alignment_id), qname, align->core.flag,
					chr, align->core.pos + 1, align->core.qual,
					arg->cigar->str, qlen, rlen, chr_next,
					align->core.mpos + 1, type);

			tree = hash_lookup (arg->tree_idx, chr);

			if (tree != NULL)
				{
					acm = ibitree_lookup (tree, align->core.pos + 1,
							align->core.pos + rlen + 1, 0.5, 0.5, 0,
							dump_if_overlaps_exon, arg);

					log_debug ("Alignment %s %s:%ld type %d overlaps %d exons",
							qname, chr, align->core.pos + 1, type, acm);
				}
		}
}

static void
dump_if_abnormal (AbnormalArg *arg)
{
	ListElmt *cur = NULL;
	bam1_t *align = NULL;
	int type = ABNORMAL_NONE;

	arg->fragment_acm++;

	for (cur = list_head (arg->stack); cur != NULL;
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
		}

	for (cur = list_head (arg->stack); cur != NULL
			&& type == ABNORMAL_NONE; cur = list_next (cur))
		{
			align = list_data (cur);

			/*
			* - supplementary
			*/
			if (align->core.flag & 0x800)
				type += ABNORMAL_SUPPLEMENTARY;

			/*
			* - reads at different chromosomes
			* - or distance between the pairs is bigger
			*   than ABNORMAL_DISTANCE_CUTOFF
			*/
			if (align->core.tid != align->core.mtid)
				type += ABNORMAL_CHROMOSOME;
			else if (abs (align->core.pos - align->core.mpos)
					> ABNORMAL_DISTANCE_CUTOFF)
				type += ABNORMAL_DISTANCE;
		}

	if (type != ABNORMAL_NONE)
		{
			log_debug ("Dump abnormal fragment of type %d", type);
			dump_alignment (arg, type);
			arg->abnormal_acm++;
		}
}

void
abnormal_filter (const char *sam_file,
		const char *gff_file, const char *db_path)
{
	log_trace ("Inside %s", __func__);
	assert (sam_file != NULL && gff_file != NULL
			&& db_path != NULL);

	int rc = 0;
	AbnormalArg *arg = NULL;

	// Opon SAM/BAM
	// Connect to database
	// Allocate resources
	arg = abnormal_arg_new (sam_file, db_path);

	// Index protein coding genes into the database
	// and its exons into an intervalar tree by
	// chromosome
	log_info ("Index annotation file '%s'", gff_file);
	index_dump_gff_file (arg, gff_file);

	log_info ("Searching for abnormal alignments into '%s'",
			sam_file);

	while (1)
		{
			// Returns -1 or < -1 in case of
			// error
			rc = sam_read1 (arg->in, arg->hdr, arg->align);
			if (rc < 0)
				break;

			if (list_size (arg->stack)
					&& !inside_fragment (arg->align, arg->stack))
				{
					dump_if_abnormal (arg);
					clean (arg->stack, arg->cache);
				}

			push (arg->align, arg->stack, arg->cache);
		}

	// Catch if it ocurred an error
	// in reading from input
	if (rc < -1)
		log_fatal ("Failed to read sam alignment");

	// The stack is filled after the parsing,
	// so it is late in relation to the file loop.
	// Therefore, test the last bunch of
	// alignments
	dump_if_abnormal (arg);

	// Just print the amount of abnormal alignments
	log_info ("Found %d abnormal alignments\n"
			"%d abnormal alignments falls inside some exonic region"
			" (%.2f%)", arg->abnormal_acm, arg->exonic_acm,
			(float) (arg->exonic_acm * 100) / arg->abnormal_acm);

	// Cleanup
	abnormal_arg_free (arg);
}
