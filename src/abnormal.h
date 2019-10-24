#ifndef ABNORMAL_H
#define ABNORMAL_H

#include "db.h"
#include "chr.h"
#include "exon.h"

enum _AbnormalType
{
	ABNORMAL_NONE          = 0,
	ABNORMAL_DISTANCE      = 1,
	ABNORMAL_CHROMOSOME    = 2,
	ABNORMAL_SUPPLEMENTARY = 4,
	ABNORMAL_EXONIC        = 8
};

typedef enum _AbnormalType AbnormalType;

struct _AbnormalArg
{
	int            tid;
	int            inc_step;
	const char    *sam_file;
	ExonTree      *exon_tree;
	ChrStd        *cs;
	sqlite3_stmt  *alignment_stmt;
	int            phred_quality;
	int            queryname_sorted;
	int            max_distance;
	float          max_base_freq;
	int            either;
	float          exon_frac;
	float          alignment_frac;
};

typedef struct _AbnormalArg AbnormalArg;

void abnormal_filter (AbnormalArg *arg);

#endif /* abnormal.h */
