#include <stdarg.h>
#include <string.h>
#include "flags.h"

int flag_match(char *flag, int nflags, ...)
{
    int i;
    va_list ap;

    va_start(ap, nflags);
    for (i = 0; i < nflags; i++)
        if (strcmp(flag, va_arg(ap, char *)) == 0)
        {
            va_end(ap);
            return 1;
        }
    
    va_end(ap);
    return 0;
}
