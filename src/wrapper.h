#ifndef WRAPPER_H
#define WRAPPER_H

#include <stdio.h>
#include <stdlib.h>

void * xmalloc    (size_t size);
void * xcalloc    (size_t nmemb, size_t size);
void * xrealloc   (void *ptr, size_t size);
char * xstrdup    (const char *str);
void   xfree      (void *ptr);

FILE * xfopen     (const char *path, const char *mode);
void   xfclose    (FILE *fp);
FILE * xpopen     (const char *cmd, const char *mode);
int    xpclose    (FILE *pp);

int    xvasprintf (char **strp, const char *fmt, va_list ap)
	__attribute__((format (printf, 2, 0)));
int    xasprintf  (char **strp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

#endif /* wrapper.h */
