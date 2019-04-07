#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#include "sam.h"
#include "bwa.h"
#include "wrapper.h"
#include "log.h"
#include "utils.h"
#include "align.h"

void
align_with_bwa_mem (const char *db, const char *read,
		const char *mate, const char *output, int n_threads)
{
	log_trace ("Inside %s", __func__);
	assert (db != NULL && read != NULL && output != NULL);

	pid_t child_id;
	int fildes[2];

	log_debug ("Create pipe between 'bwa_mem' and 'sam_to_bam'");
	if (pipe (fildes) == -1)
		log_errno_fatal ("pipe failed");

	log_debug ("Fork process: Run bwa_mem in child job");
	child_id = fork ();
	if (child_id == -1)
		log_errno_fatal ("fork failed");

	if (child_id == 0)
		{
			// Child - align with bwa mem
			log_debug ("Inside child: Prepare to run 'bwa_mem'");

			FILE *fp = NULL;
			char *dir = NULL;
			char *log_path = NULL;

			n_threads = n_threads < 1 ? 1 : n_threads;

			dir = path_dir (output);
			xasprintf_concat (&log_path, "%s/bwa.log", dir);

			fp = xfopen (log_path, "w");

			close (fildes[0]);

			if (mate == NULL)
				log_info ("Run bwa mem -t %d %s %s",
						n_threads, db, read);
			else
				log_info ("Run bwa mem -t %d %s %s %s",
						n_threads, db, read, mate);

			log_info ("For details, see: %s", log_path);

			if (dup2 (fildes[1], 1) == -1)
				log_errno_fatal ("dup2 failed");

			if (dup2 (fileno (fp), 2) == -1)
				log_errno_fatal ("dup2 failed");

			log_debug ("Run bwa wrapper 'bwa_mem'");
			bwa_mem (db, read, mate, "-", n_threads);

			close (fildes[1]);
			xfclose (fp);
			xfree (dir);
			xfree (log_path);

			log_info ("bwa_mem finished. Everything all right");
			exit (EXIT_SUCCESS);
		}
	else
		{
			// Parent - reads from pipe
			log_debug ("Inside parent: Read from child and write to bam format");

			FILE *fp = NULL;
			int ret = 0;
			clock_t start, end;
			double cpu_time = 0;
			pid_t ret_id;

			close (fildes[1]);
			fp = xfdopen (fildes[0], "r");

			log_debug ("Read from bwa_mem and write to bam format");
			start = clock ();

			ret = sam_to_bam_fp (fp, output);
			if (!ret)
				log_fatal ("bam writing failed");

			ret_id = waitpid (child_id, &ret, 0);
			if (ret_id == -1 && errno != EINTR)
				log_errno_fatal ("waitpid failed");

			end = clock ();
			cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
			log_info ("bwa_mem | sam_to_bam: elapsed time %fs", cpu_time);

			ret = !WIFEXITED (ret);
			log_info ("bwa returned status '%d'", ret);

			if (ret)
				log_fatal ("bwa failed!");

			xfclose (fp);
		}
}
