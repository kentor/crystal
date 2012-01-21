#include <string.h>
#ifndef parse_arg
#define parse_arg(flag, match, var, value) do { if (!strcmp(flag, match)) var = (value); } while (0)
#endif
