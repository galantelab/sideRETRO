#include "config.h"

#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include "log.h"
#include "wrapper.h"

#define die(fmt,...) \
	do { \
		log_fatal (fmt ": %s", ##__VA_ARGS__, strerror (errno)); \
		abort (); \
	} while (0)

char *
xstrdup (const char *str)
{
	char *ret = strdup (str);

	if (ret == NULL)
		die ("strdup failed");

	return ret;
}

void *
xmalloc (size_t size)
{
	void *ret = malloc (size);

	if (ret == NULL && !size)
		ret = malloc (1);

	if (ret == NULL)
		die ("malloc failed");

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
		die ("asprintf failed");

	va_end (ap);
	return ret;
}

void *
xcalloc (size_t nmemb, size_t size)
{
	void *ret = calloc (nmemb, size);

	if (ret == NULL && (!nmemb || !size))
		ret = calloc (1, 1);

	if (ret == NULL)
		die ("calloc failed");

	return ret;
}

void *
xrealloc (void *ptr, size_t size)
{
	void *ret = realloc (ptr, size);

	if (ret == NULL && !size)
		ret = realloc (ret, 1);

	if (ret == NULL)
		die ("realloc failed");

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
	while (1)
		{
			FILE *fp = fopen (path, mode);
			if (fp)
				return fp;

			if (errno == EINTR)
				continue;

			if (*mode && mode[1] == '+')
				die ("Could not open '%s' for reading and writing", path);
			else if (*mode == 'w' || *mode == 'a')
				die ("Could not open '%s' for writing", path);
			else
				die ("Could not open '%s' for reading", path);
		}
}

void
xfclose (FILE *fp)
{
	if (fp == NULL)
		return;

	if (fclose (fp) == EOF)
		die ("Could not close file stream");
}
