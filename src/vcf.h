#ifndef VCF_H
#define VCF_H

#include "db.h"

struct _VCFOption
{
	long        near_gene_distance;
	float       alpha_error;
	const char *fasta_file;
};

typedef struct _VCFOption VCFOption;

void vcf (sqlite3 *db, const char *output_file, VCFOption *opt);

#endif /* vcf.h */
