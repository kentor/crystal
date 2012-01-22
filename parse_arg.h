#ifndef PARSE_ARG_H
#define PARSE_ARG_H
#include <string.h>

#define parse_arg(flag, match, var, value) do { if (!strcmp(flag, match)) var = (value); } while (0)

#endif
