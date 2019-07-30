#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <check.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/gff.h"

static const char *gtf_header =
	"##description: evidence-based annotation of the human genome (GRCh38), version 30 (Ensembl 96)\n"
	"##provider: GENCODE\n"
	"##contact: gencode-help@ebi.ac.uk\n"
	"##format: gtf\n"
	"##date: 2019-03-28";

static const char *gtf_body =
	"chr1\t.\tgene\t1\t1000\t.\t+\t.\t"
	"gene_name \"g1\"; gene_id \"ENSG1\"; transcript_id \"ENST1\"; transcript_type \"protein_coding\";\n"
	"chr1\t.\ttranscript\t1\t1000\t.\t+\t.\t"
	"gene_name \"g1\"; gene_id \"ENSG1\"; transcript_id \"ENST1\"; transcript_type \"protein_coding\";\n"
	"chr1\t.\texon\t1\t100\t.\t+\t.\t"
	"gene_name \"g1\"; gene_id \"ENSG1\"; transcript_id \"ENST1\"; transcript_type \"protein_coding\"; "
	"exon_id \"ENSE1\"\n";

static const int gtf_num_lines = 8;

static void
create_gtf (char *path)
{
	FILE *fp = NULL;
	int fd;

	fd = xmkstemp (path);
	fp = xfdopen (fd, "w");

	xfprintf (fp, gtf_header, "");
	xfprintf (fp, "\n");
	xfprintf (fp, gtf_body, "");

	xfclose (fp);
}

START_TEST (test_gff_read)
{
	// Our heroes!
	GffFile *fp = NULL;
	GffEntry *e = NULL;
	GffEntry *dup = NULL;

	// Create gtf file
	char gtf_path[] = "/tmp/ponga.gtf.XXXXXX";
	create_gtf (gtf_path);

	int i = 0;

	/* TRUE POSITIVE */
	int gtf_size = 3;
	const char *chrs[] = {"chr1", "chr1", "chr1"};
	const char *types[] = {"gene", "transcript", "exon"};
	const char *gene_ids[] = {"ENSG1", "ENSG1", "ENSG1"};
	const char *transcript_ids[] = {"ENST1", "ENST1", "ENST1"};
	const char *exon_ids[] = {NULL, NULL, "ENSE1"};
	const int starts[] = {1, 1, 1};
	const int ends[] = {1000, 1000, 100};

	// Open gtf_path and allocate GffEntry
	fp = gff_open (gtf_path, "rb");
	e = gff_entry_new ();

	/* TEST GTF HEADER */
	ck_assert (fp != NULL && e != NULL);
	ck_assert_str_eq (gtf_header, fp->header);

	/* TIME TO TEST THE GTF BODY */
	while (gff_read (fp, e))
		{
			ck_assert_str_eq (chrs[i], e->seqname);
			ck_assert_str_eq (types[i], e->feature);

			ck_assert_int_eq (starts[i], e->start);
			ck_assert_int_eq (ends[i], e->end);

			ck_assert_str_eq (gff_attribute_find (e, "gene_id"),
					gene_ids[i]);
			ck_assert_str_eq (gff_attribute_find (e, "transcript_id"),
					transcript_ids[i]);

			if (exon_ids[i] != NULL)
				ck_assert_str_eq (gff_attribute_find (e, "exon_id"),
						exon_ids[i]);

			i++;
		}

	ck_assert_int_eq (e->num_line, gtf_num_lines);
	ck_assert_int_eq (i, gtf_size);

	dup = gff_entry_dup (e);

	ck_assert_str_eq (dup->seqname, e->seqname);
	ck_assert_str_eq (dup->feature, e->feature);

	ck_assert_int_eq (dup->start, e->start);
	ck_assert_int_eq (dup->end, e->end);

	ck_assert_str_eq (gff_attribute_find (dup, "gene_id"),
			gff_attribute_find (e, "gene_id"));
	ck_assert_str_eq (gff_attribute_find (dup, "transcript_id"),
			gff_attribute_find (e, "transcript_id"));

	// Cleanup
	gff_entry_free (dup);
	gff_entry_free (e);
	gff_close (fp);
	xunlink (gtf_path);
}
END_TEST

Suite *
make_gff_suite (void)
{
	Suite *s;
	TCase *tc_core;

	s = suite_create ("GFF");

	/* Core test case */
	tc_core = tcase_create ("Core");

	tcase_add_test (tc_core, test_gff_read);
	suite_add_tcase (s, tc_core);

	return s;
}
