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

#ifndef UTILS_H
#define UTILS_H

int fequal         (const double a, const double b);
int equalstring    (const void *a, const void *b);
int casequalstring (const void *a, const void *b);
int cmpstringp     (const void *p1, const void *p2);
int casecmpstringp (const void *p1, const void *p2);

char * chomp       (char *str);
char * trim        (char *str);
char * trimc       (char *str, int c);

char * upper       (char *s);
char * lower       (char *s);

char * path_dir    (const char *path);
char * path_file   (const char *path, int rm_ext);

int    which       (const char *cmd);
int    exists      (const char *file);
void   mkdir_p     (const char *path);

char * xstrdup_concat   (char *dest, const char *src);
int    xasprintf_concat (char **strp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

void setup_signal (int sig, void (*handler)(int));

size_t buf_expand (void **buf, size_t size,
		size_t old_nmemb, size_t length);
size_t entry_set  (char **buf, size_t buf_size,
		const char *entry);

#endif /* utils.h */
