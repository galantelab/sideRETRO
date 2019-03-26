#ifndef TYPES_H
#define TYPES_H

typedef void (*DestroyNotify) (void *data);
typedef int  (*CompareFunc)   (const void *a, const void *b);
typedef int  (*EqualFun)      (const void *a, const void *b);
typedef void (*Func)          (void *data, void *user_data);

#endif /* types.h */
