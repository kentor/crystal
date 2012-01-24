#ifndef LL_H
#define LL_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define null (-1)

#define new_ll(ll, size) \
do { \
   ll = malloc(size * sizeof(int)); \
   memset(ll, null, size * sizeof(int)); \
} while (0) 

#define ll_insert(ll, value, head, counter) \
do { \
   if (ll[value] != null) { printf("Illegal insertion, ll[%d] = %d\n", value, ll[value]); exit(1); } \
   if (counter == 0 && head != null) { perror("Head should be null when inserting first element."); exit(1); } \
   ll[value] = head; \
   head = value; \
   counter++; \
} while (0)

#define ll_remove(ll, value, head, counter) \
do { \
   int i = head; \
   if (i == value) head = ll[value]; \
   else while (ll[i] != value) { \
      i = ll[i]; \
      if (i == null) { printf("wtf\n"); exit(2); } \
   } \
   ll[i] = ll[value]; \
   ll[value] = null; \
   counter--; \
} while (0)

#endif
