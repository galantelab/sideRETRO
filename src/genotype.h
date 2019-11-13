#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "db.h"

enum _Zygosity
{
	HOMOZYGOUS   = 0,
	HETEROZYGOUS
};

typedef enum _Zygosity Zygosity;

void genotype (sqlite3_stmt *genotype_stmt, int threads, int crossing_reads);

#endif /* genotype.h */
