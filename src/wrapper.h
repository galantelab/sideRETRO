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

#ifndef WRAPPER_H
#define WRAPPER_H

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

void * xmalloc    (size_t size);
void * xcalloc    (size_t nmemb, size_t size);
void * xrealloc   (void *ptr, size_t size);
char * xstrdup    (const char *str);
void   xfree      (void *ptr);

FILE * xfopen     (const char *path, const char *mode);
FILE * xfdopen    (int fd, const char *mode);
void   xfclose    (FILE *fp);
void   xfflush    (FILE *fp);

FILE * xpopen     (const char *cmd, const char *mode);
int    xpclose    (FILE *pp);

void   xunlink    (const char *file);
int    xmkstemp   (char *template);
void   xmkdir     (const char *pathname, int mode);
void   xsigaction (int sig, const struct sigaction *restrict act,
		struct sigaction *restrict oact);

void   xfputs     (const char *str, FILE *fp);

int    xvasprintf (char **strp, const char *fmt, va_list ap)
	__attribute__((format (printf, 2, 0)));
int    xasprintf  (char **strp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));
int    xsnprintf  (char *str, size_t size, const char *fmt, ...)
	__attribute__((format (printf, 3, 4)));
int    xvsnprintf (char *str, size_t size, const char *fmt, va_list ap)
	__attribute__((format (printf, 3, 0)));
int    xfprintf   (FILE *fp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

#endif /* wrapper.h */
