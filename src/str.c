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

#include "config.h"

#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "str.h"

static inline void
string_maybe_expand (String *s, size_t len)
{
	if (len > s->alloc)
		{
			len = nearest_pow (len);
			s->str = xrealloc (s->str, sizeof (char) * len);
			memset (s->str + s->alloc, 0,
					sizeof (char) * (len - s->alloc));
			s->alloc = len;
		}
}

String *
string_set (String *s, const char *str)
{
	assert (s != NULL);

	if (str != NULL)
		{
			size_t str_len = strlen (str);
			string_maybe_expand (s, str_len + 1);

			memcpy (s->str, str, sizeof (char) * (str_len + 1));
			s->len = str_len;
		}

	return s;
}

String *
string_clear (String *s)
{
	assert (s != NULL);

	memset (s->str, 0, s->alloc * sizeof (char));
	s->len = 0;

	return s;
}

String *
string_concat (String *s, const char *str)
{
	assert (s != NULL);

	if (str != NULL)
		{
			size_t str_len = strlen (str);
			string_maybe_expand (s, s->len + str_len + 1);

			memcpy (s->str + s->len, str,
					sizeof (char) * (str_len + 1));
			s->len += str_len;
		}

	return s;
}

String *
string_printf (String *s, const char *fmt, ...)
{
	assert (s != NULL);

	size_t str_len = 0;
	char one_char[1];
	va_list ap, copy_ap;

	va_start (ap, fmt);
	va_copy (copy_ap, ap);

	str_len = xvsnprintf (one_char, 1, fmt, ap);
	string_maybe_expand (s, str_len + 1);

	str_len = xvsnprintf (s->str, str_len + 1, fmt, copy_ap);
	s->len = str_len;

	va_end (ap);
	va_end (copy_ap);

	return s;
}

String *
string_concat_printf (String *s, const char *fmt, ...)
{
	assert (s != NULL);

	size_t str_len = 0;
	char one_char[1];
	va_list ap, copy_ap;

	va_start (ap, fmt);
	va_copy (copy_ap, ap);

	str_len = xvsnprintf (one_char, 1, fmt, ap);
	string_maybe_expand (s, s->len + str_len + 1);

	str_len = xvsnprintf (s->str + s->len,
			str_len + 1, fmt, copy_ap);
	s->len += str_len;

	va_end (ap);
	va_end (copy_ap);

	return s;
}

String *
string_new (const char *str)
{
	String *s = xcalloc (1, sizeof (String));

	s->str = xcalloc (1, sizeof (char));
	s->alloc = 1;

	return string_set (s, str);
}

String *
string_sized_new (size_t size)
{
	size_t alloc = size > 0 ? size : 1;
	String *s = xcalloc (1, sizeof (String));

	s->str = xcalloc (alloc, sizeof (char));
	s->alloc = alloc;

	return s;
}

char *
string_free (String *s, int free_segment)
{
	if (s == NULL)
		return NULL;

	char *str = NULL;

	if (free_segment)
		xfree (s->str);
	else
		str = s->str;

	xfree (s);
	return str;
}
