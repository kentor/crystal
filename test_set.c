#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "set.h"

int main()
{
   int a = 1, b = 2, c = 3, d = 4, e = 5, f = 6;

   set *num_tbl = set_new();
   set_insert(num_tbl, &a);
   set_insert(num_tbl, &c);
   set_insert(num_tbl, &e);
   assert(set_include(num_tbl, &a));
   assert(set_include(num_tbl, &c));
   assert(set_include(num_tbl, &e));
   assert(set_size(num_tbl) == 3);

   set_iter iter;
   gpointer key, value;

   set_iter_init(iter, num_tbl);
   while (set_iter_next(iter, key, value)) {
      printf("%d\n", *((int *) value));
   }

   set_insert(num_tbl, &b);
   assert(set_size(num_tbl) == 4);

   set_iter_init(iter, num_tbl);
   while (set_iter_next(iter, key, value)) {
      printf("%d\n", *((int *) value));
   }

   set_remove(num_tbl, &a);
   assert(set_size(num_tbl) == 3);

   set_iter_init(iter, num_tbl);
   while (set_iter_next(iter, key, value)) {
      printf("%d\n", *((int *) value));
   }

   printf("%u %p %p %p %p %p %p\n", num_tbl, &a, &b, &c, &d, &e, &f);
   set_destroy(num_tbl);
   printf("%u %p %p %p %p %p %p\n", num_tbl, &a, &b, &c, &d, &e, &f);

   return 0;
}
