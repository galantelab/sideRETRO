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
