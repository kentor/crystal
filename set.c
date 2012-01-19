#include <stdlib.h>
#include <stdbool.h>

typedef struct {
   int max_value;
   int size;
   bool table[];
} set;

typedef struct {
   set *set;
   int index;
   int value;
   int counter;
} set_enum;

set *new_set(int size);
void set_insert(set *s, int value);
void set_delete(set *s, int value);
bool set_include(set *s, int value);
int set_size(set *s);

set_enum *set_to_enum(set *s);
bool set_enum_next(set_enum *s_enum);
void set_enum_free(set_enum *s_enum);

set *new_set(int size)
{
   set *s = malloc(sizeof(set) + size*sizeof(int));

   for (int i = 0; i < size; i++) {
      s->table[i] = false;
   }
   s->max_value = size - 1;
   s->size = 0;

   return s;
}

void set_insert(set *s, int value)
{
   if (value >= 0 && value <= s->max_value && !s->table[value]) {
      s->table[value] = true;
      s->size++;
   }
}

void set_delete(set *s, int value)
{
   if (s->table[value]) {
      s->table[value] = false;
      s->size--;
   }
}

bool set_include(set *s, int value)
{
   if (value < 0 || value > s->max_value) {
      return false;
   }
   return s->table[value];
}

int set_size(set *s)
{
   return s->size;
}

set_enum *set_to_enum(set *s)
{
   set_enum *s_enum = malloc(sizeof(set_enum));
   s_enum->set = s;
   s_enum->index = -1;
   s_enum->value = -1;
   s_enum->counter = 0;

   return s_enum;
}

bool set_enum_next(set_enum *s_enum)
{
   if (s_enum->counter >= s_enum->set->size) {
      return false;
   }

   do {
      s_enum->value++;
   }
   while (!s_enum->set->table[s_enum->value]);

   s_enum->counter++;
   s_enum->index++;

   return true;
}

void set_enum_free(set_enum *s_enum)
{
   free(s_enum);
}
