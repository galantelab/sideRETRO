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
