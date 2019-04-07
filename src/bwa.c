#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include <bwa/bntseq.h>
#include <bwa/bwa.h>
#include <bwa/bwamem.h>
#include <bwa/kvec.h>
#include <bwa/kseq.h>
#include <bwa/utils.h>
KSEQ_DECLARE(gzFile)

#include "log.h"
#include "wrapper.h"
#include "bwa.h"

// It's necessary, because there is no <kopen.h>
void * kopen  (const char *fn, int *_fd);
int    kclose (void *a);

static void
bwa_fprint_sam_hdr (FILE *fp, const bntseq_t *bns,
		const char *rg_line)
{
	int i;
	for (i = 0; i < bns->n_seqs; ++i)
		xfprintf (fp, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name,
				bns->anns[i].len);

	if (rg_line)
		xfprintf (fp, "%s\n", rg_line);

	xfprintf (fp, "@PG\tID:%s\tPN:%s\tVN:%s\n", "bwa",
			PACKAGE_NAME, PACKAGE_VERSION);

	xfflush (fp);
}

void
bwa_mem_log_set_level (int level)
{
	bwa_verbose = level;
}

// Modified based on main_mem in fastmap.c
void
bwa_mem (const char *db, const char *read, const char *mate,
		const char *out, int n_threads)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && read != NULL && out != NULL);

	mem_opt_t *opt;
	int fd, fd2, i, n, copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	bwaidx_t *idx;
	char *rg_line = 0;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;
	mem_pestat_t *pes0 = 0;
	FILE *fpo;

	opt = mem_opt_init ();
	opt->n_threads = n_threads > 0 ? n_threads : 1;

	bwa_fill_scmat (opt->a, opt->b, opt->mat);
	if ((idx = bwa_idx_load (db, BWA_IDX_ALL)) == 0)
		log_fatal ("Failed to load index from '%s'", db);

	ko = kopen (read, &fd);
	if (ko == 0)
		log_errno_fatal ("kopen failed");

	fp = gzdopen (fd, "r");
	ks = kseq_init (fp);

	if (mate != NULL)
		{
		if (opt->flag & MEM_F_PE)
			{
				if (bwa_verbose >= 2)
					fprintf (stderr,
							"[W::%s] when '-p' is in use, the second query file will be ignored.\n",
							__func__);
			}
		else
			{
				ko2 = kopen (mate, &fd2);
				if (ko2 == 0)
					log_errno_fatal ("kopen failed");

				fp2 = gzdopen (fd2, "r");
				ks2 = kseq_init (fp2);
				opt->flag |= MEM_F_PE;
			}
		}

	fpo = xopen (out, "w");
	bwa_fprint_sam_hdr (fpo, idx->bns, rg_line);

	while ((seqs = bseq_read (opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0)
		{
			if ((opt->flag & MEM_F_PE) && (n&1) == 1)
				{
					if (bwa_verbose >= 2)
						fprintf (stderr,
								"[W::%s] odd number of reads in the PE mode; last read dropped\n",
								__func__);
					n = n >> 1 << 1;
				}

			if (!copy_comment)
				{
					for (i = 0; i < n; ++i)
						{
							xfree (seqs[i].comment);
							seqs[i].comment = 0;
						}
				}

			mem_process_seqs (opt, idx->bwt, idx->bns, idx->pac,
			n_processed, n, seqs, pes0);
			n_processed += n;

			for (i = 0; i < n; ++i)
				{
					xfputs (seqs[i].sam, fpo);
					xfree (seqs[i].name);
					xfree (seqs[i].comment);
					xfree (seqs[i].seq);
					xfree (seqs[i].qual);
					xfree (seqs[i].sam);
				}

			xfree (seqs);
		}

	xfree (opt);
	bwa_idx_destroy (idx);
	kseq_destroy (ks);
	err_fclose (fpo);
	err_gzclose (fp);
	kclose (ko);

	if (ks2)
		{
			kseq_destroy (ks2);
			err_gzclose (fp2);
			kclose (ko2);
		}
}
