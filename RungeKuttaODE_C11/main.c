#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "prototype_declarations.h"

int k_max, count;
float *t_p, **x_p, dt_save;

int main(int argc, char *argv[])
{
    int i = 0;
    LOOP(i,0,10) printf("%d", i);

    return 0;
}
