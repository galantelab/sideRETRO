/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br
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

/*
 * Fault injection: The error paths of sideRETRO are tested with
 * real faults whenever possible (/dev/full, invalid paths, broken
 * files, SQLite limits). Out of memory and pipe failures cannot be
 * produced for real, so the check executable is linked with
 * '-Wl,--wrap=...' (HAVE_LD_WRAP): Each wrapped function fails once
 * when a test arms it, and calls the real function otherwise
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <wchar.h>
#include <unistd.h>
#include <assert.h>
#include <check.h>
#include <sqlite3.h>
#include "check_sider.h"

#include "../src/wrapper.h"
#include "../src/array.h"
#include "../src/io.h"
#include "../src/db.h"

#ifdef HAVE_LD_WRAP

#include <regex.h>
#include <pthread.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/hfile.h"
#include "htslib/thread_pool.h"
#include "../src/log.h"
#include "../src/logger.h"
#include "../src/gff.h"
#include "../src/gz.h"
#include "../src/sam.h"

/*
 * Armed failures: Countdown of calls until the failure.
 * Zero means disarmed
 */
static int armed[FAULT_NUM] = {};

void
fault_inject_after (Fault f, int calls)
{
	assert (f >= 0 && f < FAULT_NUM && calls >= 0);
	armed[f] = calls + 1;
}

void
fault_inject (Fault f)
{
	fault_inject_after (f, 0);
}

static inline int
fault_fire (Fault f, int err)
{
	if (armed[f] == 0 || --armed[f] > 0)
		return 0;

	errno = err;
	return 1;
}

void      * __real_malloc           (size_t size);
void      * __real_calloc           (size_t nmemb, size_t size);
void      * __real_realloc          (void *ptr, size_t size);
char      * __real_strdup           (const char *s);
FILE      * __real_popen            (const char *cmd, const char *mode);
int         __real_pclose           (FILE *pp);
int         __real_regcomp          (regex_t *preg, const char *regex, int cflags);
int         __real_pthread_mutex_init (pthread_mutex_t *m, const pthread_mutexattr_t *a);
int         __real_gzclose          (gzFile file);
bam1_t    * __real_bam_init1        (void);
bam1_t    * __real_bam_dup1         (const bam1_t *bsrc);
bam1_t    * __real_bam_copy1        (bam1_t *bdst, const bam1_t *bsrc);
int         __real_hts_close        (htsFile *fp);
hts_itr_t * __real_sam_itr_queryi   (const hts_idx_t *idx, int tid, int beg, int end);
int         __real_hts_set_threads  (htsFile *fp, int n);
int         __real_sam_hdr_write    (samFile *fp, const bam_hdr_t *h);
int         __real_sam_read1        (samFile *fp, bam_hdr_t *h, bam1_t *b);
hFILE     * __real_hdopen           (int fd, const char *mode);
htsFile   * __real_hts_hopen        (hFILE *fp, const char *fn, const char *mode);

void *
__wrap_malloc (size_t size)
{
	return fault_fire (FAULT_MALLOC, ENOMEM)
		? NULL
		: __real_malloc (size);
}

void *
__wrap_calloc (size_t nmemb, size_t size)
{
	return fault_fire (FAULT_CALLOC, ENOMEM)
		? NULL
		: __real_calloc (nmemb, size);
}

void *
__wrap_realloc (void *ptr, size_t size)
{
	return fault_fire (FAULT_REALLOC, ENOMEM)
		? NULL
		: __real_realloc (ptr, size);
}

char *
__wrap_strdup (const char *s)
{
	return fault_fire (FAULT_STRDUP, ENOMEM)
		? NULL
		: __real_strdup (s);
}

FILE *
__wrap_popen (const char *cmd, const char *mode)
{
	return fault_fire (FAULT_POPEN, EMFILE)
		? NULL
		: __real_popen (cmd, mode);
}

int
__wrap_regcomp (regex_t *preg, const char *regex, int cflags)
{
	return fault_fire (FAULT_REGCOMP, ENOMEM)
		? REG_ESPACE
		: __real_regcomp (preg, regex, cflags);
}

int
__wrap_pthread_mutex_init (pthread_mutex_t *m, const pthread_mutexattr_t *a)
{
	return fault_fire (FAULT_PTHREAD_MUTEX_INIT, ENOMEM)
		? ENOMEM
		: __real_pthread_mutex_init (m, a);
}

bam1_t *
__wrap_bam_init1 (void)
{
	return fault_fire (FAULT_BAM_INIT1, ENOMEM)
		? NULL
		: __real_bam_init1 ();
}

bam1_t *
__wrap_bam_dup1 (const bam1_t *bsrc)
{
	return fault_fire (FAULT_BAM_DUP1, ENOMEM)
		? NULL
		: __real_bam_dup1 (bsrc);
}

bam1_t *
__wrap_bam_copy1 (bam1_t *bdst, const bam1_t *bsrc)
{
	return fault_fire (FAULT_BAM_COPY1, ENOMEM)
		? NULL
		: __real_bam_copy1 (bdst, bsrc);
}

hts_itr_t *
__wrap_sam_itr_queryi (const hts_idx_t *idx, int tid, int beg, int end)
{
	return fault_fire (FAULT_SAM_ITR_QUERYI, ENOMEM)
		? NULL
		: __real_sam_itr_queryi (idx, tid, beg, end);
}

int
__wrap_hts_set_threads (htsFile *fp, int n)
{
	return fault_fire (FAULT_HTS_SET_THREADS, EAGAIN)
		? -1
		: __real_hts_set_threads (fp, n);
}

int
__wrap_sam_hdr_write (samFile *fp, const bam_hdr_t *h)
{
	return fault_fire (FAULT_SAM_HDR_WRITE, EIO)
		? -1
		: __real_sam_hdr_write (fp, h);
}

int
__wrap_sam_read1 (samFile *fp, bam_hdr_t *h, bam1_t *b)
{
	// Less than -1: Read error (-1 is the end of file)
	return fault_fire (FAULT_SAM_READ1, EIO)
		? -4
		: __real_sam_read1 (fp, h, b);
}

hFILE *
__wrap_hdopen (int fd, const char *mode)
{
	if (fault_fire (FAULT_HDOPEN, ENOMEM))
		{
			close (fd);
			return NULL;
		}

	return __real_hdopen (fd, mode);
}

htsFile *
__wrap_hts_hopen (hFILE *fp, const char *fn, const char *mode)
{
	if (fault_fire (FAULT_HTS_HOPEN, ENOMEM))
		{
			hclose_abruptly (fp);
			return NULL;
		}

	return __real_hts_hopen (fp, fn, mode);
}

/*
 * The closing functions close for real, so nothing
 * leaks, and then report the failure
 */

int
__wrap_pclose (FILE *pp)
{
	int rc = __real_pclose (pp);
	return fault_fire (FAULT_PCLOSE, ECHILD) ? -1 : rc;
}

int
__wrap_gzclose (gzFile file)
{
	int rc = __real_gzclose (file);
	return fault_fire (FAULT_GZCLOSE, EIO) ? Z_ERRNO : rc;
}

int
__wrap_hts_close (htsFile *fp)
{
	int rc = __real_hts_close (fp);
	return fault_fire (FAULT_HTS_CLOSE, EIO) ? -1 : rc;
}

/* Out of memory: sideRETRO aborts with a message */

START_TEST (test_xmalloc_oom)
{
	fault_inject (FAULT_MALLOC);
	xmalloc (16);
}
END_TEST

START_TEST (test_xcalloc_oom)
{
	fault_inject (FAULT_CALLOC);
	xcalloc (1, 16);
}
END_TEST

START_TEST (test_xrealloc_oom)
{
	void *p = xmalloc (16);
	fault_inject (FAULT_REALLOC);
	xrealloc (p, 32);
}
END_TEST

START_TEST (test_xstrdup_oom)
{
	fault_inject (FAULT_STRDUP);
	xstrdup ("ponga");
}
END_TEST

START_TEST (test_xalloc_zero)
{
	/*
	* A zero size allocation may return NULL: Then
	* allocate one byte, so the caller never gets NULL
	*/
	void *p = NULL;

	fault_inject (FAULT_MALLOC);
	p = xmalloc (0);
	ck_assert_ptr_ne (p, NULL);
	xfree (p);

	fault_inject (FAULT_CALLOC);
	p = xcalloc (0, 16);
	ck_assert_ptr_ne (p, NULL);
	xfree (p);

	// realloc (p, 0) is never called: It is undefined
	p = xmalloc (16);
	p = xrealloc (p, 0);
	ck_assert_ptr_ne (p, NULL);
	xfree (p);
}
END_TEST

/* Pipes */

START_TEST (test_xpopen_r_fail)
{
	fault_inject (FAULT_POPEN);
	xpopen ("true", "r");
}
END_TEST

START_TEST (test_xpopen_w_fail)
{
	fault_inject (FAULT_POPEN);
	xpopen ("true", "w");
}
END_TEST

START_TEST (test_xpclose_fail)
{
	FILE *pp = xpopen ("true", "r");
	fault_inject (FAULT_PCLOSE);
	xpclose (pp);
}
END_TEST

/* System and libraries */

START_TEST (test_regcomp_fail)
{
	fault_inject (FAULT_REGCOMP);
	gff_looks_like_gff_file ("ponga.gff3");
}
END_TEST

START_TEST (test_mutex_init_fail)
{
	fault_inject (FAULT_PTHREAD_MUTEX_INIT);
	logger_new (NULL, LOG_INFO, 1, 0);
}
END_TEST

START_TEST (test_gzclose_fail)
{
	char path[] = "/tmp/ponga.txt.XXXXXX";
	int fd = xmkstemp (path);
	GzFile *gz = NULL;

	close (fd);
	gz = gz_open_for_reading (path);
	xunlink (path);

	fault_inject (FAULT_GZCLOSE);
	gz_close (gz);
}
END_TEST

/* sam.c: A SAM to convert */

static const char *fault_sam =
	"@SQ\tSN:chr1\tLN:248956422\n"
	"R1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n";

static void
fault_sam_to_bam (Fault f)
{
	char sam[] = "/tmp/ponga.sam.XXXXXX";
	char bam[] = "/tmp/ponga.bam.XXXXXX";
	FILE *fp = NULL;
	int fd = 0;

	fd = xmkstemp (sam);
	fp = xfdopen (fd, "w");
	xfprintf (fp, "%s", fault_sam);
	xfclose (fp);

	fd = xmkstemp (bam);
	close (fd);

	fault_inject (f);
	ck_assert_int_eq (sam_to_bam (sam, bam), 0);

	xunlink (sam);
	xunlink (bam);
}

START_TEST (test_sam_hdopen_fail)
{
	fault_sam_to_bam (FAULT_HDOPEN);
}
END_TEST

START_TEST (test_sam_hts_hopen_fail)
{
	fault_sam_to_bam (FAULT_HTS_HOPEN);
}
END_TEST

START_TEST (test_sam_bam_init1_fail)
{
	fault_sam_to_bam (FAULT_BAM_INIT1);
}
END_TEST

START_TEST (test_sam_hdr_write_fail)
{
	// Not fatal: sam_to_bam returns false
	fault_sam_to_bam (FAULT_SAM_HDR_WRITE);
}
END_TEST

START_TEST (test_sam_close_fail)
{
	// The input stream is closed first
	fault_sam_to_bam (FAULT_HTS_CLOSE);
}
END_TEST

#endif /* HAVE_LD_WRAP */

/* Real faults */

START_TEST (test_xpclose_null)
{
	ck_assert_int_eq (xpclose (NULL), EOF);
}
END_TEST

START_TEST (test_xunlink_fail)
{
	xunlink ("/it-must-not-exist/ponga");
}
END_TEST

START_TEST (test_xmkstemp_fail)
{
	char template[] = "/it-must-not-exist/ponga.XXXXXX";
	xmkstemp (template);
}
END_TEST

/*
 * A wide string with a lone surrogate cannot be converted
 * to a multibyte string: The printf family fails (EILSEQ)
 */
static const wchar_t bad_wide[] = {0xd800, 0};

START_TEST (test_xsnprintf_fail)
{
	char buf[16];
	xsnprintf (buf, sizeof (buf), "%ls", bad_wide);
}
END_TEST

static void
call_xvsnprintf (const char *fmt, ...)
{
	char buf[16];
	va_list ap;

	va_start (ap, fmt);
	xvsnprintf (buf, sizeof (buf), fmt, ap);
	va_end (ap);
}

START_TEST (test_xvsnprintf_fail)
{
	call_xvsnprintf ("%ls", bad_wide);
}
END_TEST

START_TEST (test_xasprintf_fail)
{
	char *s = NULL;
	xasprintf (&s, "%ls", bad_wide);
}
END_TEST

static void
call_xvasprintf (const char *fmt, ...)
{
	char *s = NULL;
	va_list ap;

	va_start (ap, fmt);
	xvasprintf (&s, fmt, ap);
	va_end (ap);
}

START_TEST (test_xvasprintf_fail)
{
	call_xvasprintf ("%ls", bad_wide);
}
END_TEST

/* /dev/full: Every write fails with ENOSPC */

static FILE *
open_full (int buffered)
{
	FILE *fp = xfopen ("/dev/full", "w");

	if (!buffered)
		setvbuf (fp, NULL, _IONBF, 0);

	return fp;
}

START_TEST (test_xfprintf_full)
{
	xfprintf (open_full (0), "%s", "ponga");
}
END_TEST

START_TEST (test_xfputs_full)
{
	xfputs ("ponga", open_full (0));
}
END_TEST

START_TEST (test_xfflush_full)
{
	FILE *fp = open_full (1);
	fprintf (fp, "%s", "ponga");
	xfflush (fp);
}
END_TEST

START_TEST (test_xfclose_full)
{
	FILE *fp = open_full (1);
	fprintf (fp, "%s", "ponga");
	xfclose (fp);
}
END_TEST

/* Reading a directory: getline fails with EISDIR */

START_TEST (test_read_file_lines_dir)
{
	Array *arr = array_new (xfree);
	read_file_lines (arr, "/");
}
END_TEST

/* SQLite: Real failures, each with its hint */

static sqlite3 *
fault_db (char *path)
{
	int fd = xmkstemp (path);
	close (fd);
	return db_create (path);
}

START_TEST (test_db_close_busy)
{
	// A statement not finalized keeps the database busy
	char path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = fault_db (path);

	db_prepare (db, "SELECT * FROM schema");
	db_close (db);
}
END_TEST

START_TEST (test_db_finalize_error)
{
	// finalize returns the error of the last step
	char path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = fault_db (path);
	sqlite3_stmt *stmt = NULL;

	db_exec (db, "CREATE TABLE ponga (id INTEGER PRIMARY KEY)");
	db_exec (db, "INSERT INTO ponga VALUES (1)");

	stmt = db_prepare (db, "INSERT INTO ponga VALUES (1)");
	ck_assert_int_ne (sqlite3_step (stmt), SQLITE_DONE);
	db_finalize (stmt);
}
END_TEST

START_TEST (test_db_bind_null_range)
{
	char path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = fault_db (path);
	sqlite3_stmt *stmt = db_prepare (db, "SELECT ?1");

	db_bind_null (stmt, 0);
}
END_TEST

START_TEST (test_db_full)
{
	// No more pages: SQLITE_FULL
	char path[] = "/tmp/ponga.db.XXXXXX";
	sqlite3 *db = fault_db (path);
	int i = 0;

	db_exec (db, "PRAGMA max_page_count = 20");
	db_exec (db, "CREATE TABLE ponga (s TEXT)");

	for (i = 0; i < 1000; i++)
		db_exec (db, "INSERT INTO ponga VALUES (hex (randomblob (1000)))");
}
END_TEST

START_TEST (test_db_notadb)
{
	// Not a database: SQLITE_NOTADB
	char path[] = "/tmp/ponga.db.XXXXXX";
	int fd = xmkstemp (path);
	FILE *fp = xfdopen (fd, "w");

	xfprintf (fp, "%s", "I am not a database, but I am long enough to be read as one. "
			"I am not a database, but I am long enough to be read as one.\n");
	xfclose (fp);

	db_exec (db_open (path, SQLITE_OPEN_READWRITE), "SELECT * FROM schema");
}
END_TEST

/*
 * I/O error: A VFS that wraps the default one and fails
 * all writes. It is registered as the default VFS only at
 * the forked test process
 */

static sqlite3_io_methods failing_methods;
static const sqlite3_io_methods *real_methods = NULL;
static sqlite3_vfs failing_vfs;
static sqlite3_vfs *real_vfs = NULL;

static int
failing_write (sqlite3_file *f, const void *buf, int amt, sqlite3_int64 off)
{
	return SQLITE_IOERR_WRITE;
}

static int
failing_open (sqlite3_vfs *vfs, const char *name, sqlite3_file *f,
		int flags, int *out_flags)
{
	int rc = real_vfs->xOpen (real_vfs, name, f, flags, out_flags);

	if (rc == SQLITE_OK && f->pMethods != NULL)
		{
			real_methods = f->pMethods;
			failing_methods = *real_methods;
			failing_methods.xWrite = failing_write;
			f->pMethods = &failing_methods;
		}

	return rc;
}

START_TEST (test_db_ioerr)
{
	char path[] = "/tmp/ponga.db.XXXXXX";
	int fd = xmkstemp (path);
	sqlite3 *db = NULL;

	close (fd);

	real_vfs = sqlite3_vfs_find (NULL);
	failing_vfs = *real_vfs;
	failing_vfs.zName = "failing";
	failing_vfs.xOpen = failing_open;
	sqlite3_vfs_register (&failing_vfs, 1);

	db = db_open (path, SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE);
	db_exec (db, "CREATE TABLE ponga (id INTEGER)");
}
END_TEST

Suite *
make_fault_suite (void)
{
	Suite *s;
	TCase *tc_core;
	TCase *tc_abort;

	s = suite_create ("Fault");

	tc_core = tcase_create ("Core");
	tc_abort = tcase_create ("Abort");

	tcase_add_test (tc_core, test_xpclose_null);

#ifdef HAVE_LD_WRAP
	tcase_add_test (tc_core, test_xalloc_zero);

	tcase_add_exit_test (tc_abort, test_xmalloc_oom,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xcalloc_oom,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xrealloc_oom,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xstrdup_oom,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xpopen_r_fail, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xpopen_w_fail, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xpclose_fail,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_regcomp_fail,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_mutex_init_fail,    EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_gzclose_fail,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_sam_hdopen_fail,    EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_sam_hts_hopen_fail, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_sam_bam_init1_fail, EXIT_FAILURE);

	tcase_add_test (tc_core, test_sam_hdr_write_fail);
	tcase_add_test (tc_core, test_sam_close_fail);
#endif

	tcase_add_exit_test (tc_abort, test_xunlink_fail,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xmkstemp_fail,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xsnprintf_fail,      EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xvsnprintf_fail,     EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xasprintf_fail,      EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xvasprintf_fail,     EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xfprintf_full,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xfputs_full,         EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xfflush_full,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_xfclose_full,        EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_read_file_lines_dir, EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_db_close_busy,       EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_db_finalize_error,   EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_db_bind_null_range,  EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_db_full,             EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_db_notadb,           EXIT_FAILURE);
	tcase_add_exit_test (tc_abort, test_db_ioerr,            EXIT_FAILURE);

	suite_add_tcase (s, tc_core);
	suite_add_tcase (s, tc_abort);

	return s;
}
