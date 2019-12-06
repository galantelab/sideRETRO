#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "db.h"

void genotype (sqlite3_stmt *genotype_stmt, int threads, int phred_quality);

#endif /* genotype.h */
