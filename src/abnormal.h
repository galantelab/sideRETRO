#ifndef ABNORMAL_H
#define ABNORMAL_H

#define ABNORMAL_DISTANCE_CUTOFF 10000

#define ABNORMAL_NONE            0
#define ABNORMAL_DISTANCE        1
#define ABNORMAL_CHROMOSOME      2
#define ABNORMAL_SUPPLEMENTARY   4

void abnormal_filter (const char *sam_file, const char *gff_file, const char *db_path);

#endif /* abnormal.h */
