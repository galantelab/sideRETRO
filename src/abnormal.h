#ifndef ABNORMAL_H
#define ABNORMAL_H

#include "db.h"
#include "chr.h"
#include "exon.h"

#define ABNORMAL_DISTANCE_CUTOFF 10000

#define ABNORMAL_NONE            0
#define ABNORMAL_DISTANCE        1
#define ABNORMAL_CHROMOSOME      2
#define ABNORMAL_SUPPLEMENTARY   4
#define ABNORMAL_EXONIC          8

struct _AbnormalArg
{
	int            tid;
	int            num_threads;
	const char    *sam_file;
	ExonTree      *exon_tree;
	ChrStd        *cs;
	sqlite3       *db;
	sqlite3_stmt  *alignment_stmt;
	int            either;
	float          exon_frac;
	float          alignment_frac;
};

typedef struct _AbnormalArg AbnormalArg;

void abnormal_filter (AbnormalArg *arg);

#endif /* abnormal.h */
