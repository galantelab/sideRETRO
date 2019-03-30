#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <ctype.h>
#include <string.h>
#include "wrapper.h"
#include "utils.h"

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
trim (char *str)
{
	if (str == NULL)
		return NULL;

	size_t start, end;
	size_t len = strlen (str);

	// Empty string
	if (len == 0)
		return str;

	for (start = 0; start < len && isspace (str[start]); start++)
		;

	// All spaces
	if (start == len)
		{
			str[0] = '\0';
			return str;
		}

	for (end = len - 1; end >= start && isspace (str[end]); end--)
		;

	memmove (str, str + start, sizeof (char) * (end - start + 1));
	str[end - start + 1] = '\0';

	return str;
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
			char *dot = NULL;

			if (*file == '.')
				{
					dot = strrchr (file, '.');
					if (dot != NULL && dot != file)
						*dot = '\0';
				}
			else
				{
					dot = strchr (file, '.');
					if (dot != NULL)
						*dot = '\0';
				}
		}

	xfree (path_copy);
	return file;
}
