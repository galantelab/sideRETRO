#include "config.h"

#include <assert.h>
#include "wrapper.h"
#include "log.h"
#include "logger.h"

static void
logger_lock_print (void *udata, int lock)
{
	pthread_mutex_t *l = udata;
	if (lock)
		pthread_mutex_lock (l);
	else
		pthread_mutex_unlock (l);
}

static void
logger_set (Logger *l)
{
	assert (l != NULL);

	log_set_fp (l->fp);
	log_set_udata (l->lock);

	log_set_level (l->level);
	log_set_quiet (l->silent);
	log_set_color (l->color);

	log_set_lock (logger_lock_print);
}

Logger *
logger_new  (const char *file, int level, int silent, int color)
{
	FILE *fp = NULL;
	Logger *l = xcalloc (1, sizeof (Logger));
	pthread_mutex_t *lock = xcalloc (1, sizeof (pthread_mutex_t));

	if (pthread_mutex_init (lock, NULL) != 0)
		log_errno_fatal ("Failed to create pthread mutex");

	if (file != NULL)
		fp = xfopen (file, "w");

	*l = (Logger) {
		.fp     = fp,
		.lock   = lock,
		.level  = level,
		.silent = silent,
		.color  = color
	};

	logger_set (l);
	return l;
}

void
logger_free (Logger *l)
{
	if (l == NULL)
		return;

	xfclose (l->fp);
	pthread_mutex_destroy (l->lock);
	xfree (l->lock);
	xfree (l);
}
