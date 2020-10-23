#include <stdio.h>
#include <math.h>
#include "nrutil.h"

int main(int argc, char *argv[])
{
    float x_1 = 2.3;
    float x_2 = -0.1;

    printf("%f", SIGN(x_1, x_2));

    return 0;
}
