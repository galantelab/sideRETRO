#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "db.h"

enum _Zygosity
{
	HOMOZYGOUS   = 0,
	HETEROZYGOUS
};

typedef enum _Zygosity Zygosity;

void genotype (sqlite3 *db, int threads);

#endif /* genotype.h */
