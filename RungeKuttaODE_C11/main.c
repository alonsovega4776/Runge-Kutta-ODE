#include <stdio.h>
#include <math.h>
#include "prototype_declarations.h"

struct ODEsolution odeSol;

int main(int argc, char *argv[])
{
   struct ODE_IVP problem;
   ODE_IVP_constructor(&problem,2,0.0,2.0,0.1,0.1,0.2);
   problem.x_0[1] = 0.0;
   problem.x_0[2] = 0.0;

   ODEsolution_constructor(&odeSol,20,0.1,2);

    ODE_driver(problem.x_0,problem.N_var,problem.t_1,problem.t_2,
               problem.TOL,problem.Δ_1,problem.Δ_min,problem.n_good,problem.n_bad,
               f_ODE_IVP,rungeKutta_stepper);

    ODE_IVP_destructor(&problem);
    ODEsolution_destructor(&odeSol);
    return 0;
}

void f_ODE_IVP(float t, float x[], float **xDot, int function_number)
{
        xDot[1][function_number] = 2 * x[1] + 5.0;
        xDot[2][function_number] = x[1] * x[2];
}