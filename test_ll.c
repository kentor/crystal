#include <stdio.h>
#include <assert.h>
#include "ll.h"

int main()
{
   describe_ll: {
      int *ll;
      list_t list_a = { .head = null, .size = 0 }, list_b = { .head = null, .size = 0};
         
      describe_new_ll: {
         it_should_create_a_new_linked_list_with_null_elements: {
            new_ll(ll, 108000);
            for (int i = 0; i < 108000; i++) {
               assert(ll[i] == null);
            }
            puts(".");
         }
      }

      describe_ll_insert: {
         it_should_make_the_head_of_chain_to_inserted_value_and_increment_correctly: {
            ll_insert(ll, 2, list_a);
            assert(list_a.head == 2);
            assert(list_a.size == 1);
            ll_insert(ll, 4, list_a);
            assert(list_a.head == 4);
            assert(list_a.size == 2);
            ll_insert(ll, 6, list_a);
            assert(list_a.head == 6);
            assert(list_a.size == 3);

            ll_insert(ll, 1, list_b);
            assert(list_b.head == 1);
            assert(list_b.size == 1);
            ll_insert(ll, 3, list_b);
            assert(list_b.head == 3);
            assert(list_b.size == 2);
            puts(".");
         }

         it_should_correctly_link_the_values: {
            assert(ll[6] == 4);
            assert(ll[4] == 2);
            assert(ll[2] == null);
            assert(ll[3] == 1);
            assert(ll[1] == null);
            puts(".");
         }
      }

      test_ll_looping: {
         for (int i = list_a.head; i != null; i = ll[i]) {
            printf("%d ", i);
         }
         printf("\n");
         for (int i = list_b.head; i != null; i = ll[i]) {
            printf("%d ", i);
         }
         printf("\n");
      }

      describe_ll_remove: {
         it_should_decrease_counter_and_set_removed_to_null: {
            ll_remove(ll, 4, list_a);
            assert(list_a.size == 2);
            assert(ll[4] == null);
            puts(".");
         }

         it_should_link_correctly_after_deletion: {
            assert(ll[6] == 2);
            assert(ll[2] == null);
            puts(".");
         }

         it_should_set_new_head_if_deleting_first_element: {
            ll_remove(ll, 3, list_b);
            assert(list_b.head == 1);
            assert(list_b.size == 1);
            puts(".");
         }

         it_should_set_head_null_when_deleting_last_element: {
            ll_remove(ll, 1, list_b);
            assert(list_b.head == null);
            assert(list_b.size == 0);
            puts(".");
         }
      }
   }
   return 0;
}
