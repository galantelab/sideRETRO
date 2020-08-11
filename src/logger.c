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
