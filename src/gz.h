#ifndef GZ_H
#define GZ_H

#include <stdlib.h>

typedef struct _GzFile GzFile;

GzFile     * gz_open_for_reading (const char *path);
void         gz_close            (GzFile *gz);
const char * gz_get_filename     (const GzFile *gz);
size_t       gz_get_num_line     (const GzFile *gz);
int          gz_getline          (GzFile *gz, char **lineptr, size_t *n);

#endif /* gz.h */
