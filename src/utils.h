#ifndef UTILS_H
#define UTILS_H

char * chomp     (char *str);
char * trim      (char *str);
char * path_dir  (const char *path);
char * path_file (const char *path, int rm_ext);

#endif /* utils.h */
