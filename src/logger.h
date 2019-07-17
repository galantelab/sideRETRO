#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <pthread.h>

struct _Logger
{
	FILE            *fp;
	pthread_mutex_t *lock;
	int              level;
	int              silent;
	int              color;
};

typedef struct _Logger Logger;

Logger * logger_new  (const char *file, int level, int silent, int color);
void     logger_free (Logger *l);

#endif /* logger.h */
