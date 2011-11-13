#include <stdarg.h>
#include <string.h>

int flag_match(char *flag, int nflags, ...);
void print_progress(int current, int final);

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

void print_progress(int current, int final)
{
    float f;
    int i, bars;
    int max_bars = 50;
    
    f = (float) current / final;

    bars = current / (final / max_bars);

    printf("  ");
    printf("%3d%% ", (int)(100*f));
    printf("[");
    for (i = 0; i < bars; i++)
        printf("=");
    for (i = bars; i < max_bars; i++)
        printf(" ");
    printf("] ");
    printf("%d/%d ", current, final);
    printf("\r");
    fflush(stdout);
}
