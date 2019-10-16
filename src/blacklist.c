#include "config.h"

#include <assert.h>
#include <regex.h>
#include "wrapper.h"
#include "log.h"
#include "gff.h"
#include "bed.h"
#include "ibitree.h"
#include "str.h"
#include "blacklist.h"

struct _BlacklistData
{
	Blacklist *blacklist;
	long       cluster_id;
};

typedef struct _BlacklistData BlacklistData;

Blacklist *
blacklist_new (sqlite3_stmt *blacklist_stmt,
		sqlite3_stmt *overlapping_blacklist_stmt,
		ChrStd *cs)
{
	assert (blacklist_stmt != NULL
			&& overlapping_blacklist_stmt != NULL
			&& cs != NULL);

	Blacklist *blacklist = xcalloc (1, sizeof (Blacklist));

	blacklist->blacklist_stmt = blacklist_stmt;
	blacklist->overlapping_blacklist_stmt =
		overlapping_blacklist_stmt;

	blacklist->cs = cs;

	blacklist->idx = hash_new (xfree,
			(DestroyNotify) ibitree_free);

	return blacklist;
}

void
blacklist_free (Blacklist *blacklist)
{
	if (blacklist == NULL)
		return;

	hash_free (blacklist->idx);
	xfree (blacklist);
}

static inline void
init_attribute_regex (regex_t *reg, const Set *attributes)
{
	String *regex_s = NULL;
	ListElmt *cur = NULL;

	regex_s = string_new (NULL);
	cur = list_head (set_list (attributes));

	string_concat_printf (regex_s, "(%s", (char *) list_data (cur));

	for (cur = list_next (cur); cur != NULL; cur = list_next (cur))
		string_concat_printf (regex_s, "|%s", (char *) list_data (cur));

	string_concat (regex_s, ")");

	log_debug ("GFF/GTF attributes regex: %s", regex_s->str);

	if (regcomp (reg, regex_s->str, REG_EXTENDED|REG_NOSUB) != 0)
		log_fatal ("Error regcomp for '%s'", regex_s->str);

	string_free (regex_s, 1);
}

static inline int
match_attribute (regex_t *reg, const char *str)
{
	return regexec (reg, str, 0, (regmatch_t *) NULL, 0) == 0;
}

static void
blacklist_index_dump_from_gff (Blacklist *blacklist, const char *gff_file,
		const char *feature, const char *attribute, const Set *attributes)
{
	assert (blacklist != NULL && gff_file != NULL && feature != NULL
			&& attribute != NULL && attributes != NULL
			&& set_size (attributes) > 0);

	// Compile regex
	regex_t reg;
	init_attribute_regex (&reg, attributes);

	// Open for reading gff file
	GffFile *gff = gff_open_for_reading (gff_file);
	GffEntry *entry = gff_entry_new ();

	IBiTree *tree = NULL;

	long table_id = blacklist->table_id;
	long *alloc_id = NULL;

	const char *chr_std = NULL;
	const char *gene_name = NULL;
	const char *attribute_value = NULL;

	while (gff_read (gff, entry))
		{
			attribute_value = gff_attribute_find (entry, attribute);

			if (!strcmp (entry->feature, feature)
					&& attribute_value != NULL
					&& match_attribute (&reg, attribute_value))
				{
					gene_name = gff_attribute_find (entry, "gene_name");

					if (gene_name == NULL)
						gene_name = "blacklist";

					chr_std = chr_std_lookup (blacklist->cs, entry->seqname);

					log_debug ("Index blacklist '%s' at %s:%zu-%zu", gene_name,
							chr_std, entry->start, entry->end);

					alloc_id = xcalloc (1, sizeof (long));
					* (long *) alloc_id = ++table_id;

					tree = hash_lookup (blacklist->idx, chr_std);

					if (tree == NULL)
						{
							tree = ibitree_new (xfree);
							hash_insert (blacklist->idx,
									xstrdup (chr_std), tree);
						}

					ibitree_insert (tree, entry->start, entry->end,
							alloc_id);

					db_insert_blacklist (blacklist->blacklist_stmt,
							table_id, gene_name, chr_std, entry->start,
							entry->end);
				}
		}

	blacklist->table_id = table_id;

	gff_entry_free (entry);
	gff_close (gff);
}

static void
blacklist_index_dump_from_bed (Blacklist *blacklist, const char *bed_file)
{
	assert (blacklist != NULL && bed_file != NULL);

	// Open for reading bed file
	BedFile *bed = bed_open_for_reading (bed_file);
	BedEntry *entry = bed_entry_new ();

	IBiTree *tree = NULL;

	long table_id = blacklist->table_id;
	long *alloc_id = NULL;

	const char *chr_std = NULL;
	const char *name = NULL;

	while (bed_read (bed, entry))
		{
			name = entry->num_field > 3 && entry->name != NULL
				? entry->name
				: "blacklist";

			chr_std = chr_std_lookup (blacklist->cs, entry->chrom);

			log_debug ("Index blacklist '%s' at %s:%zu-%zu", name,
					chr_std, entry->chrom_start, entry->chrom_end);

			alloc_id = xcalloc (1, sizeof (long));
			* (long *) alloc_id = ++table_id;

			tree = hash_lookup (blacklist->idx, chr_std);

			if (tree == NULL)
				{
					tree = ibitree_new (xfree);
					hash_insert (blacklist->idx,
							xstrdup (chr_std), tree);
				}

			ibitree_insert (tree, entry->chrom_start,
					entry->chrom_end, alloc_id);

			db_insert_blacklist (blacklist->blacklist_stmt,
					table_id, name, chr_std, entry->chrom_start,
					entry->chrom_end);
		}

	blacklist->table_id = table_id;

	bed_entry_free (entry);
	bed_close (bed);
}

static void
dump_if_overlaps_blacklist (IBiTreeLookupData *ldata,
		void *user_data)
{
	const long *blacklist_id = ldata->data;
	BlacklistData *data = user_data;

	log_debug ("Dump overlapping blacklist [%li] %li-%li with cluster [%li] %li-%li at %li-%li",
			*blacklist_id, ldata->node_low, ldata->node_high, data->cluster_id,
			ldata->interval_low, ldata->interval_high, ldata->overlap_pos,
			ldata->overlap_pos + ldata->overlap_len - 1);

	db_insert_overlapping_blacklist (data->blacklist->overlapping_blacklist_stmt, *blacklist_id,
			data->cluster_id, ldata->overlap_pos, ldata->overlap_len);
}

int
blacklist_lookup (Blacklist *blacklist, const char *chr,
		long low, long high, long padding,
		const long cluster_id)
{
	assert (blacklist != NULL && chr != NULL
			&& padding >= 0);

	int acm = 0;
	IBiTree *tree = NULL;

	tree = hash_lookup (blacklist->idx, chr);

	if (tree != NULL)
		{
			low -= padding;
			BlacklistData data = {blacklist, cluster_id};
			acm = ibitree_lookup (tree, low > 0 ? low : 0, high + padding,
					-1, -1, 0, dump_if_overlaps_blacklist, &data);
		}

	return acm;
}

static inline int
is_gff_file (const char *file)
{
	regex_t reg;
	const char regex_str[] = "\\.(gff|gff3|gtf)(\\.gz)*$";

	if (regcomp (&reg, regex_str, REG_EXTENDED|REG_NOSUB) != 0)
		log_fatal ("Error regcomp for '%s'", regex_str);

	return regexec (&reg, file, 0, (regmatch_t *) NULL, 0) == 0;
}

void
blacklist_index_dump (Blacklist *blacklist, const char *file,
		const char *feature, const char *attribute,
		const Set *attributes)
{
	assert (blacklist != NULL && file != NULL);

	if (is_gff_file (file))
		{
			log_info ("Index blacklist entries from GFF/GTF file '%s'", file);
			blacklist_index_dump_from_gff (blacklist, file, feature, attribute, attributes);
		}
	else
		{
			log_info ("Index blacklist entries from BED file '%s'", file);
			blacklist_index_dump_from_bed (blacklist, file);
		}
}
