/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br>
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

#ifndef SUMMARY_H
#define SUMMARY_H

#include <time.h>
#include "str.h"

/*
 * Summary printed at the end of each command: A framed
 * block with a title, sections and 'key  value' lines
 */

#define SUMMARY_WIDTH    68
#define SUMMARY_KEY_LEN  24
#define SUMMARY_NUM_LEN  32

struct _Summary
{
	String *msg;
	time_t  start;
};

typedef struct _Summary Summary;

Summary    * summary_new          (const char *title, time_t start);
void         summary_free         (Summary *s);

void         summary_section      (Summary *s, const char *name);
void         summary_line         (Summary *s, const char *key, const char *fmt, ...)
	__attribute__((format (printf, 3, 4)));
void         summary_text         (Summary *s, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

void         summary_log          (Summary *s);

const char * summary_count        (long n, char buf[SUMMARY_NUM_LEN]);
const char * summary_percent      (long part, long total, char buf[SUMMARY_NUM_LEN]);
const char * summary_elapsed      (time_t start, char buf[SUMMARY_NUM_LEN]);

#endif /* summary.h */
