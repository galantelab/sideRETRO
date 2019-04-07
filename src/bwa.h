#ifndef BWA_H
#define BWA_H

void bwa_mem_log_set_level (int level);

void bwa_mem (const char *db, const char *read, const char *mate,
		const char *out, int n_threads);

#endif /* bwa.h */
