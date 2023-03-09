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

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <libgen.h>
#include <string.h>
#include <signal.h>
#include <float.h>
#include <math.h>
#include "wrapper.h"
#include "utils.h"

int
fequal (const double a, const double b)
{
	return fabs (a - b) < DBL_EPSILON;
}

int
equalstring (const void *a, const void *b)
{
	return !strcmp (a, b);
}

int
casequalstring (const void *a, const void *b)
{
	return !strcasecmp (a, b);
}

int
cmpstringp (const void *p1, const void *p2)
{
	return strcmp (* (char * const *) p1, * (char * const *) p2);
}

int
casecmpstringp (const void *p1, const void *p2)
{
	return strcasecmp (* (char * const *) p1, * (char * const *) p2);
}

char *
chomp (char *str)
{
	if (str == NULL)
		return NULL;

	size_t len = strlen (str);

	if (len && str[len - 1] == '\n')
		str[len - 1] = '\0';

	return str;
}

char *
trimc (char *str, int c)
{
	if (str == NULL)
		return NULL;

	size_t start, end;
	size_t len = strlen (str);

	// Empty string
	if (len == 0)
		return str;

	// Leading chars 'c'
	for (start = 0; start < len && str[start] == c; start++)
		;

	// All 'c' characters
	if (start == len)
		{
			str[0] = '\0';
			return str;
		}

	// Trailing chars 'c'
	for (end = len - 1; end >= start && str[end] == c; end--)
		;

	memmove (str, str + start, sizeof (char) * (end - start + 1));
	str[end - start + 1] = '\0';

	return str;
}

char *
trim (char *str)
{
	return trimc (str, ' ');
}

char *
upper (char *s)
{
	int i = 0;
	while (s[i])
		{
			if (s[i] >= 'a' && s[i] <= 'z')
				s[i] -= 32;
			i++;
		}
	return s;
}

char *
lower (char *s)
{
	int i = 0;
	while (s[i])
		{
			if (s[i] >= 'A' && s[i] <= 'Z')
				s[i] += 32;
			i++;
		}
	return s;
}

char *
path_dir (const char *path)
{
	if (path == NULL)
		return NULL;

	char *dir = NULL;
	char *path_copy = NULL;

	path_copy = xstrdup (path);
	dir = xstrdup (dirname (path_copy));

	xfree (path_copy);
	return dir;
}

char *
path_file (const char *path, int rm_ext)
{
	if (path == NULL)
		return NULL;

	char *file = NULL;
	char *path_copy = NULL;

	path_copy = xstrdup (path);
	file = xstrdup (basename (path_copy));

	if (rm_ext)
		{
			char *dot = strrchr (file, '.');
			if (dot != NULL && dot != file)
				*dot = '\0';
		}

	xfree (path_copy);
	return file;
}

char *
xstrdup_concat (char *dest, const char *src)
{
	size_t len_dest = 0;
	size_t len_src = 0;

	len_src = strlen (src);
	if (dest != NULL)
		len_dest = strlen (dest);

	dest = xrealloc (dest, sizeof (char) * (len_dest + len_src + 1));

	memcpy (dest + len_dest, src, sizeof (char) * len_src);
	dest[len_dest + len_src] = '\0';

	return dest;
}

int
xasprintf_concat (char **strp, const char *fmt, ...)
{
	if (strp == NULL)
		return 0;

	char *tmp_str1 = NULL;
	char *tmp_str2 = NULL;
	int ret = 0;
	va_list ap;

	va_start (ap, fmt);
	tmp_str1 = *strp;

	ret = xvasprintf (strp, fmt, ap);

	if (tmp_str1 != NULL)
		{
			tmp_str2 = *strp;
			ret = xasprintf (strp, "%s%s", tmp_str1, tmp_str2);
		}

	xfree (tmp_str1);
	xfree (tmp_str2);
	va_end (ap);

	return ret;
}

int
which (const char *cmd)
{
	char *paths = NULL;
	char *path = NULL;
	char *scratch = NULL;
	char paths_copy[BUFSIZ];
	char exec[BUFSIZ];
	int found = 0;

	paths = secure_getenv ("PATH");
	if (paths == NULL)
		return found;

	strncpy (paths_copy, paths, BUFSIZ);
	paths_copy[BUFSIZ - 1] = '\0';

	path = strtok_r (paths_copy, ":", &scratch);

	while (path != NULL)
		{
			struct stat sb;
			snprintf (exec, BUFSIZ, "%s/%s", path, cmd);

			if (stat (exec, &sb) == 0 && sb.st_mode & S_IXUSR)
				{
					found = 1;
					break;
				}

			path = strtok_r (NULL, ":", &scratch);
		}

	return found;
}

int
exists (const char *file)
{
	struct stat sb;
	return stat (file, &sb) == 0 && sb.st_mode & S_IFREG;
}

void
mkdir_p (const char *path)
{
	/* Adapted from http://stackoverflow.com/a/2336245/119527 */
	char path_copy[BUFSIZ];
	char *p;

	strncpy (path_copy, path, BUFSIZ);
	path_copy[BUFSIZ - 1] = '\0';

	/* Iterate the string */
	for (p = path_copy + 1; *p; p++)
		{
			if (*p == '/')
				{
					/* Temporarily truncate */
					*p = '\0';

					xmkdir (path_copy, S_IRWXU);

					*p = '/';
				}
		}

	xmkdir (path_copy, S_IRWXU);
}

void
setup_signal (int sig, void (*handler)(int))
{
	struct sigaction action;

	action.sa_handler = handler;
	sigemptyset (&action.sa_mask);
	action.sa_flags = 0;

	xsigaction (sig, &action, NULL);
}

size_t
buf_expand (void **buf, size_t size,
		size_t old_nmemb, size_t length)
{
	size_t final_nmemb = nearest_pow (old_nmemb + length);
	*buf = xrealloc (*buf, size * final_nmemb);
	memset (*buf + old_nmemb * size, 0,
			size * (final_nmemb - old_nmemb));
	return final_nmemb;
}

size_t
entry_set (char **buf, size_t buf_size, const char *entry)
{
	size_t entry_size = strlen (entry);

	if (entry_size >= buf_size)
		buf_size = buf_expand ((void **) buf, sizeof (char),
				buf_size, entry_size - buf_size + 1);

	memcpy (*buf, entry, sizeof (char) * entry_size);
	(*buf)[entry_size] = '\0';

	return buf_size;
}
