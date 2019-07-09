#ifndef CHR_H
#define CHR_H

#include "hash.h"

typedef Hash ChrStd;

ChrStd     * chr_std_new    (void);
const char * chr_std_lookup (ChrStd *cs, const char *chr);
void         chr_std_free   (ChrStd *cs);

#endif /* chr.h */
