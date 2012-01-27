#include <stdio.h>
#include <stdarg.h>
#include <string.h>

void print_progress(int current, int final);

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
