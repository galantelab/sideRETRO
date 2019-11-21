#ifndef VCF_H
#define VCF_H

#include "db.h"

void vcf (sqlite3 *db, const char *fasta_file, const char *output_file);

#endif /* vcf.h */
