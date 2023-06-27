/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2020 Thiago L. A. Miller <tmiller@mochsl.org.br
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
	const char    *output_dir;
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
