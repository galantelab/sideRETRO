#ifndef SAM_H
#define SAM_H

#include <htslib/sam.h>

int sam_to_bam_fp (FILE *fp, const char *output_file);
int sam_to_bam    (const char *input_file, const char *output_file);
int sam_test_sorted_order (const bam_hdr_t *hdr, const char *value);

#endif /* sam.h */
