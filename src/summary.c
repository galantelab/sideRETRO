/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br
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

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "summary.h"

static void
summary_rule (Summary *s, const char *title)
{
	int len = 0;
	int left = 0;
	int right = 0;
	int i = 0;

	if (title == NULL)
		{
			for (i = 0; i < SUMMARY_WIDTH; i++)
				string_concat (s->msg, "=");

			string_concat (s->msg, "\n");
			return;
		}

	// '=== title ==='
	len = strlen (title) + 2;
	left = (SUMMARY_WIDTH - len) / 2;
	right = SUMMARY_WIDTH - len - left;

	for (i = 0; i < left; i++)
		string_concat (s->msg, "=");

	string_concat_printf (s->msg, " %s ", title);

	for (i = 0; i < right; i++)
		string_concat (s->msg, "=");

	string_concat (s->msg, "\n");
}

Summary *
summary_new (const char *title, time_t start)
{
	assert (title != NULL);

	Summary *s = xcalloc (1, sizeof (Summary));

	s->msg = string_sized_new (BUFSIZ);
	s->start = start;

	// The log line prefix comes before the first line
	string_concat (s->msg, "\n");
	summary_rule (s, title);

	return s;
}

void
summary_free (Summary *s)
{
	if (s == NULL)
		return;

	string_free (s->msg, 1);
	xfree (s);
}

void
summary_section (Summary *s, const char *name)
{
	assert (s != NULL && name != NULL);
	string_concat_printf (s->msg, "\n %s\n", name);
}

void
summary_line (Summary *s, const char *key, const char *fmt, ...)
{
	assert (s != NULL && key != NULL && fmt != NULL);

	va_list ap;
	char *value = NULL;

	va_start (ap, fmt);
	xvasprintf (&value, fmt, ap);
	va_end (ap);

	string_concat_printf (s->msg, "   %-*s %s\n",
			SUMMARY_KEY_LEN, key, value);

	xfree (value);
}

void
summary_text (Summary *s, const char *fmt, ...)
{
	assert (s != NULL && fmt != NULL);

	va_list ap;
	char *text = NULL;

	va_start (ap, fmt);
	xvasprintf (&text, fmt, ap);
	va_end (ap);

	string_concat_printf (s->msg, "   %s\n", text);

	xfree (text);
}

void
summary_log (Summary *s)
{
	assert (s != NULL);

	char buf[SUMMARY_NUM_LEN] = {};

	string_concat_printf (s->msg, "\n   %-*s %s\n",
			SUMMARY_KEY_LEN, "Elapsed time",
			summary_elapsed (s->start, buf));

	summary_rule (s, NULL);

	log_info ("%s", s->msg->str);
}

const char *
summary_count (long n, char buf[SUMMARY_NUM_LEN])
{
	char digits[SUMMARY_NUM_LEN] = {};
	int len = 0;
	int i = 0;
	int j = 0;

	// 1234567 => 1,234,567
	len = snprintf (digits, SUMMARY_NUM_LEN, "%ld", n < 0 ? -n : n);

	if (n < 0)
		buf[j++] = '-';

	for (i = 0; i < len; i++)
		{
			if (i > 0 && (len - i) % 3 == 0)
				buf[j++] = ',';

			buf[j++] = digits[i];
		}

	buf[j] = '\0';
	return buf;
}

const char *
summary_percent (long part, long total, char buf[SUMMARY_NUM_LEN])
{
	if (total <= 0)
		snprintf (buf, SUMMARY_NUM_LEN, "-");
	else
		snprintf (buf, SUMMARY_NUM_LEN, "%.1f%%",
				(double) part * 100 / total);

	return buf;
}

const char *
summary_elapsed (time_t start, char buf[SUMMARY_NUM_LEN])
{
	long secs = (long) difftime (time (NULL), start);

	if (secs < 0)
		secs = 0;

	if (secs >= 3600)
		snprintf (buf, SUMMARY_NUM_LEN, "%ldh %02ldm %02lds",
				secs / 3600, (secs % 3600) / 60, secs % 60);
	else if (secs >= 60)
		snprintf (buf, SUMMARY_NUM_LEN, "%ldm %02lds",
				secs / 60, secs % 60);
	else
		snprintf (buf, SUMMARY_NUM_LEN, "%lds", secs);

	return buf;
}
