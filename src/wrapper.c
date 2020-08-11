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

#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "log.h"
#include "wrapper.h"

char *
xstrdup (const char *str)
{
	char *ret = strdup (str);

	if (ret == NULL)
		log_errno_fatal ("strdup failed");

	return ret;
}

void *
xmalloc (size_t size)
{
	void *ret = malloc (size);

	if (ret == NULL && !size)
		ret = malloc (1);

	if (ret == NULL)
		log_errno_fatal ("malloc failed");

	return ret;
}

int
xvasprintf (char **strp, const char *fmt, va_list ap)
{
	int ret = vasprintf (strp, fmt, ap);

	if (ret < 0)
		log_errno_fatal ("vasprintf failed");

	return ret;
}

int
xasprintf (char **strp, const char *fmt, ...)
{
	int ret = 0;
	va_list ap;

	va_start (ap, fmt);
	ret = vasprintf (strp, fmt, ap);

	if (ret < 0)
		log_errno_fatal ("asprintf failed");

	va_end (ap);
	return ret;
}

int
xsnprintf (char *str, size_t size, const char *fmt, ...)
{
	int ret = 0;
	va_list ap;

	va_start (ap, fmt);
	ret = vsnprintf (str, size, fmt, ap);

	if (ret < 0)
		log_errno_fatal ("snprintf failed");

	va_end (ap);
	return ret;
}

int
xvsnprintf (char *str, size_t size, const char *fmt, va_list ap)
{
	int ret = 0;

	ret = vsnprintf (str, size, fmt, ap);

	if (ret < 0)
		log_errno_fatal ("vsnprintf failed");

	return ret;
}

void *
xcalloc (size_t nmemb, size_t size)
{
	void *ret = calloc (nmemb, size);

	if (ret == NULL && (!nmemb || !size))
		ret = calloc (1, 1);

	if (ret == NULL)
		log_errno_fatal ("calloc failed");

	return ret;
}

void *
xrealloc (void *ptr, size_t size)
{
	void *ret = realloc (ptr, size);

	if (ret == NULL && !size)
		ret = realloc (ret, 1);

	if (ret == NULL)
		log_errno_fatal ("realloc failed");

	return ret;
}

void
xfree (void *ptr)
{
	if (ptr == NULL)
		return;
	free (ptr);
}

FILE *
xfopen (const char *path, const char *mode)
{
	FILE *fp = fopen (path, mode);
	if (fp != NULL)
		return fp;

	if (*mode && mode[1] == '+')
		log_errno_fatal ("Could not open '%s' for reading and writing", path);
	else if (*mode == 'w' || *mode == 'a')
		log_errno_fatal ("Could not open '%s' for writing", path);
	else
		log_errno_fatal ("Could not open '%s' for reading", path);
}

FILE *
xfdopen (int fd, const char *mode)
{
	FILE *fp = fdopen (fd, mode);
	if (fp != NULL)
		return fp;

	if (*mode && mode[1] == '+')
		log_errno_fatal ("Could not open fd '%d' for reading and writing", fd);
	else if (*mode == 'w' || *mode == 'a')
		log_errno_fatal ("Could not open fd '%d' for writing", fd);
	else
		log_errno_fatal ("Could not open fd '%d' for reading", fd);
}

void
xfclose (FILE *fp)
{
	if (fp == NULL)
		return;

	if (fclose (fp) == EOF)
		log_errno_fatal ("Could not close file stream");
}

FILE *
xpopen (const char *cmd, const char *mode)
{
	FILE *pp = popen (cmd, mode);
	if (pp != NULL)
		return pp;

	if (*mode && mode[0] == 'w')
		log_errno_fatal ("Could not open pipe for writing to '%s'", cmd);
	else
		log_errno_fatal ("Could not open pipe for reading from '%s'", cmd);
}

int
xpclose (FILE *pp)
{
	if (pp == NULL)
		return EOF;

	int stat = pclose (pp);

	if (stat == EOF)
		log_errno_fatal ("Could not close pipe");

	return WEXITSTATUS (stat);
}

void
xunlink (const char *file)
{
	if (unlink (file) == -1)
		log_errno_fatal ("Could not remove file '%s'", file);
}

int
xmkstemp (char *template)
{
	int fd = mkstemp (template);
	if (fd == -1)
		log_errno_fatal ("Could not create temporary file '%s'", template);
	return fd;
}

int
xfprintf (FILE *fp, const char *fmt, ...)
{
	int ret = 0;
	va_list ap;

	va_start (ap, fmt);
	ret = vfprintf (fp, fmt, ap);

	if (ret < 0)
		log_errno_fatal ("fprintf failed");

	return ret;
}

void
xfputs (const char *str, FILE *fp)
{
	if (fputs (str, fp) == EOF)
		log_errno_fatal ("fputs failed");
}

void
xfflush (FILE *fp)
{
	if (fflush (fp) == EOF)
		log_errno_fatal ("fflush failed");
}

void
xmkdir (const char *path, int mode)
{
	errno = 0;
	if (mkdir (path, mode) == -1)
		if (errno != EEXIST)
			log_errno_fatal ("mkdir failed to create dir '%s'", path);
}

void
xsigaction (int sig, const struct sigaction *restrict act,
		struct sigaction *restrict oact)
{
	if (sigaction (sig, act, oact) == -1)
		log_errno_fatal ("sigaction failed");
}
