#ifndef SAM_SORT_H
#define SAM_SORT_H

#include <stdlib.h>

int sam_sort (int is_by_qname, const char *fn,
		const char *prefix, size_t max_mem, int n_threads);

#endif /* sam_sort.h */
