#include <stdbool.h>

typedef struct {
   int max_value;
   int size;
   bool table[];
} set;

typedef struct {
   const set *set;
   int index;
   int value;
   int counter;
} set_enum;

set *new_set(int size);
void set_insert(set *s, int value);
void set_delete(set *s, int value);
bool set_include(const set *s, int value);
int set_size(const set *s);

set_enum *set_to_enum(const set *s);
bool set_enum_next(set_enum *s_enum);
void set_enum_free(set_enum *s_enum);
