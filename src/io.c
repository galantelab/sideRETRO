#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wrapper.h"
#include "utils.h"
#include "log.h"
#include "io.h"

void
read_file_lines (Array *arr, const char *filename)
{
	assert (arr != NULL && filename != NULL);

	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	size_t nread;

	errno = 0;

	fp = xfopen (filename, "r");

	while (!feof (fp))
		{
			nread = getline (&line, &len, fp);

			// End of file or error
			if (nread == EOF)
				{
					if (errno)
						log_errno_fatal ("Failed to read '%s'",
								filename);
					break;
				}

			chomp (line);
			trim (line);

			// Line is empty
			if (*line == '\0')
				continue;

			char *line_copy = xstrdup (line);
			array_add (arr, line_copy);
		}

	xfree (line);
	xfclose (fp);
}
