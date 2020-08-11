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

#ifndef STR_H
#define STR_H

#include <stdlib.h>

struct _String
{
	char   *str;
	size_t  len;
	size_t  alloc;
};

typedef struct _String String;

String * string_new           (const char *str);
String * string_sized_new     (size_t size);
char   * string_free          (String *s, int free_segment);

String * string_set           (String *s, const char *str);
String * string_clear         (String *s);
String * string_concat        (String *s, const char *str);

#define string_reset(s) ((s)->str[0] = '\0', (s)->len = 0)

String * string_printf        (String *s, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));
String * string_concat_printf (String *s, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

#endif /* string.h */
