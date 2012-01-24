#include <stdio.h>
#include <assert.h>
#include "ll.h"

int main()
{
   describe_ll: {
      int *ll;
      int head_a = null, counter_a = 0, head_b = null, counter_b = 0;
         
      describe_new_ll: {
         it_should_create_a_new_linked_list_with_null_elements: {
            new_ll(ll, 1000);
            for (int i = 0; i < 1000; i++) {
               assert(ll[i] == null);
            }
            puts(".");
         }
      }

      describe_ll_insert: {
         it_should_make_the_head_of_chain_to_inserted_value_and_increment_correctly: {
            ll_insert(ll, 2, head_a, counter_a);
            assert(head_a == 2);
            assert(counter_a == 1);
            ll_insert(ll, 4, head_a, counter_a);
            assert(head_a == 4);
            assert(counter_a == 2);
            ll_insert(ll, 6, head_a, counter_a);
            assert(head_a == 6);
            assert(counter_a == 3);

            ll_insert(ll, 1, head_b, counter_b);
            assert(head_b == 1);
            assert(counter_b == 1);
            ll_insert(ll, 3, head_b, counter_b);
            assert(head_b == 3);
            assert(counter_b == 2);
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
         for (int i = head_a; i != null; i = ll[i]) {
            printf("%d ", i);
         }
         printf("\n");
         for (int i = head_b; i != null; i = ll[i]) {
            printf("%d ", i);
         }
         printf("\n");
      }

      describe_ll_remove: {
         it_should_decrease_counter_and_set_removed_to_null: {
            ll_remove(ll, 4, head_a, counter_a);
            assert(counter_a == 2);
            assert(ll[4] == null);
            puts(".");
         }

         it_should_link_correctly_after_deletion: {
            assert(ll[6] == 2);
            assert(ll[2] == null);
            puts(".");
         }

         it_should_set_new_head_if_deleting_first_element: {
            ll_remove(ll, 3, head_b, counter_b);
            assert(head_b == 1);
            assert(counter_b == 1);
            puts(".");
         }

         it_should_set_head_null_when_deleting_last_element: {
            ll_remove(ll, 1, head_b, counter_b);
            assert(head_b == null);
            assert(counter_b == 0);
            puts(".");
         }
      }
   }
   return 0;
}
