#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
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
