#ifndef __ll_h__
#define __ll_h__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define null -1

typedef struct list {
   int head;
   int size;
} list_t;

#define new_ll(_ll, _size) \
do { \
   _ll = malloc(_size * sizeof(int)); \
   memset(_ll, null, _size * sizeof(int)); \
} while (0) 

#define ll_insert(_ll, _value, _list) \
do { \
   if (_ll[_value] != null) { \
      perror("Illegal insertion error."); \
      abort(); \
   } \
   if (_list.size == 0 && _list.head != null) { \
      perror("Head should be null when inserting first element."); \
      abort(); \
   } \
   _ll[_value] = _list.head; \
   _list.head = _value; \
   _list.size++; \
} while (0)

#define ll_remove(_ll, _value, _list) \
do { \
   int i = _list.head; \
   if (i == _value) _list.head = _ll[_value]; \
   else while (_ll[i] != _value) i = _ll[i]; \
   _ll[i] = _ll[_value]; \
   _ll[_value] = null; \
   _list.size--; \
} while (0)

#endif
