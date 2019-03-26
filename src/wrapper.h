#ifndef WRAPPER_H
#define WRAPPER_H

#include <stdlib.h>

void * xmalloc   (size_t size);
void * xcalloc   (size_t nmemb, size_t size);
void * xrealloc  (void *ptr, size_t size);
char * xstrdup   (const char *str);
int    xasprintf (char **strp, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

#endif /* wrapper.h */
