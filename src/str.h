#ifndef STR_H
#define STR_H

#include <stdlib.h>

struct _String
{
	char   *str;
	size_t  len;
	size_t  alloc;
};

typedef struct _String String;

String * string_new           (const char *str);
String * string_sized_new     (size_t size);
char   * string_free          (String *s, int free_segment);

String * string_set           (String *s, const char *str);
String * string_clear         (String *s);
String * string_concat        (String *s, const char *str);

String * string_printf        (String *s, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));
String * string_concat_printf (String *s, const char *fmt, ...)
	__attribute__((format (printf, 2, 3)));

#endif /* string.h */
