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

#ifndef GZ_H
#define GZ_H

#include <stdlib.h>
#include <zlib.h>

struct _GzFile
{
	gzFile        fp;
	const char   *filename;
	char         *buf;
	size_t        buf_size;
	size_t        num_line;
};

typedef struct _GzFile GzFile;

GzFile     * gz_open_for_reading (const char *path);
void         gz_close            (GzFile *gz);
int          gz_getline          (GzFile *gz, char **lineptr, size_t *n);

#define      gz_get_fp(gz)       ((gz)->fp)
#define      gz_get_filename(gz) ((gz)->filename)
#define      gz_get_num_line(gz) ((gz)->num_line)

#endif /* gz.h */
