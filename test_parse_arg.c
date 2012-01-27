#include <stdio.h>
#include <assert.h>
#include "parse_arg.h"

int main()
{
   describe_parse_arg: {
      it_should_set_var_if_strings_match: {
         int x = 0;
         char *string = "string";
         parse_arg("string", "string", x, 1);
         assert(x == 1);
         parse_arg(string, "string", x, 2);
         assert(x == 2);
         parse_arg("string", string, x, 3);
         assert(x == 3);
         parse_arg(string, string, x, 4);
         assert(x == 4);
         puts(".");
      }

      it_should_not_set_var_if_strings_do_not_match: {
         int x = 0;
         char *string = "poop";
         parse_arg("1", "2", x, 1);
         assert(x == 0);
         parse_arg(string, "poopoop", x, 1);
         assert(x == 0);
         puts(".");
      }
   }
   return 0;
}
