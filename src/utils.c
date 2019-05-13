#include "config.h"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <libgen.h>
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
trimc (char *str, int c)
{
	if (str == NULL)
		return NULL;

	size_t start, end;
	size_t len = strlen (str);

	// Empty string
	if (len == 0)
		return str;

	for (start = 0; start < len && str[start] == c; start++)
		;

	// All 'c' characters
	if (start == len)
		{
			str[0] = '\0';
			return str;
		}

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
