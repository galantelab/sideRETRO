#include "config.h"

#include <string.h>
#include "wrapper.h"

char *
xstrdup (const char *str)
{
	char *ret = strdup (str);

	if (ret == NULL)
		{
		}

	return ret;
}

void *
xmalloc (size_t size)
{
	void *ret = malloc (size);

	if (ret == NULL && !size)
		ret = malloc (1);

	if (ret == NULL)
		{
		}

	return ret;
}

void *
xcalloc (size_t nmemb, size_t size)
{
	void *ret = calloc (nmemb, size);

	if (ret == NULL && (!nmemb || !size))
		ret = calloc (1, 1);

	if (ret == NULL)
		{
		}

	return ret;
}

void *
xrealloc (void *ptr, size_t size)
{
	void *ret = realloc (ptr, size);

	if (ret == NULL && !size)
		ret = realloc (ret, 1);

	if (ret == NULL)
		{
		}

	return ret;
}
