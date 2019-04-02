#ifndef SAM_H
#define SAM_H

int sam_to_bam_fp (FILE *fp, const char *output_file);
int sam_to_bam    (const char *input_file, const char *output_file);

#endif /* sam.h */
