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
FILE * xfdopen    (int fd, const char *mode);
void   xfclose    (FILE *fp);
void   xfflush    (FILE *fp);

FILE * xpopen     (const char *cmd, const char *mode);
int    xpclose    (FILE *pp);

void   xunlink    (const char *file);
int    xmkstemp   (char *template);

void   xfputs     (const char *str, FILE *fp);

int    xvasprintf (char **strp, const char *fmt, va_list ap)
	__attribute__((format (printf, 2, 0)));
int    xasprintf  (char **strp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));
int    xfprintf (FILE *fp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

#endif /* wrapper.h */
