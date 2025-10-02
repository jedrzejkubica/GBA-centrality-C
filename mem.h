#ifndef _MEM_H_
#define _MEM_H_

#include <stddef.h>

/* allocate memory, or print errMess and die if malloc fails */
void *mallocOrDie(size_t size, char *errMess) ;

#endif
