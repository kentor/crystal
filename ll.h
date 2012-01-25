#ifndef LL_H
#define LL_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define null -1

typedef struct list {
   int head;
   int size;
} list_t;

#define new_ll(_ll, size) \
do { \
   _ll = malloc(size * sizeof(int)); \
   memset(_ll, -1, size * sizeof(int)); \
} while (0) 

#define ll_insert(ll, value, list) \
do { \
   if (ll[value] != null) { printf("Illegal insertion, ll[%d] = %d\n", value, ll[value]); exit(1); } \
   if (list.size == 0 && list.head != null) { perror("Head should be null when inserting first element."); exit(1); } \
   ll[value] = list.head; \
   list.head = value; \
   list.size++; \
} while (0)

#define ll_remove(ll, value, list) \
do { \
   int i = list.head; \
   if (i == value) list.head = ll[value]; \
   else while (ll[i] != value) i = ll[i]; \
   ll[i] = ll[value]; \
   ll[value] = null; \
   list.size--; \
} while (0)

#endif
