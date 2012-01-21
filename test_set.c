#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "set.h"

int main()
{
   describe_set: {
      set *silver = new_set(1000);

      it_should_start_with_size_of_zero: {
         assert(set_size(silver) == 0);
         puts(".");
      }

      describe_set_insert: {
         it_should_increase_size_when_item_is_inserted: {
            set_insert(silver, 1);
            assert(set_size(silver) == 1);
            set_insert(silver, 2);
            assert(set_size(silver) == 2);
            set_insert(silver, 104);
            assert(set_size(silver) == 3);
            puts(".");
         }

         it_should_do_nothing_when_inserting_oob_item: {
            int size = set_size(silver);
            set_insert(silver, -1);
            assert(size - set_size(silver) == 0);
            set_insert(silver, -3094);
            assert(size - set_size(silver) == 0);
            set_insert(silver, 1001);
            assert(size - set_size(silver) == 0);
            set_insert(silver, 1000);
            assert(size - set_size(silver) == 0);
            puts(".");
         }

         it_should_not_increase_size_when_inserting_existing_item: {
            int size = set_size(silver);
            set_insert(silver, 1);
            assert(size - set_size(silver) == 0);
            set_insert(silver, 2);
            assert(size - set_size(silver) == 0);
            puts(".");
         }
      }

      describe_set_delete: {
         it_should_decrease_size_when_item_is_deleted: {
            int size = set_size(silver);
            set_delete(silver, 1);
            assert(size - set_size(silver) == 1);
            puts(".");
         }

         it_should_do_nothing_when_deleting_non_existing_item: {
            int size = set_size(silver);
            set_delete(silver, 1);
            assert(size - set_size(silver) == 0);
            puts(".");
         }
      }

      describe_set_include:
      {
         it_should_return_boolean: {
            assert(set_include(silver, 2) == true);
            assert(set_include(silver, 104) == true);
            assert(set_include(silver, 1) == false);
            assert(set_include(silver, -1) == false); 
            puts(".");
         }
      }
   }

   describe_set_enum: {
      set *silver = new_set(50);
      set_enum *silver_enum = set_to_enum(silver);

      describe_set_enum_next: {
         it_should_be_false_if_set_is_empty: {
            assert(set_enum_next(silver_enum) == false);
            puts(".");
         }

         set_enum_free(silver_enum);

         it_should_give_next_item_and_increment_index: {
            for (int i = 0; i < 50; i = i+5) {
               set_insert(silver, i);
            }
            silver_enum = set_to_enum(silver);
            assert(set_enum_next(silver_enum) == true);
            assert(silver_enum->value == 0);
            assert(silver_enum->index == 0);
            assert(set_enum_next(silver_enum) == true);
            assert(silver_enum->value == 5);
            assert(silver_enum->index == 1);
            assert(set_enum_next(silver_enum) == true);
            assert(silver_enum->value == 10);
            assert(silver_enum->index == 2);
            while(set_enum_next(silver_enum));
            assert(silver_enum->value == 45);
            assert(silver_enum->index == silver_enum->set->size - 1);
            puts(".");
         }

         set_enum_free(silver_enum);
      }
   }
   return 0;
}
