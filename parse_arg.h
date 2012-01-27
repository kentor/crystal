#ifndef __parse_arg_h__
#define __parse_arg_h__
#include <string.h>

#define parse_arg(flag, match, var, value) do { if (!strcmp(flag, match)) var = (value); } while (0)

#endif
