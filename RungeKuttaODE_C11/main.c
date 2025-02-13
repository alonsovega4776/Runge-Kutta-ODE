#include <stdio.h>
#include <math.h>
#include "prototype_declarations.h"

struct ODEsolution odeSol;

int main(int argc, char *argv[])
{
   struct ODE_IVP problem;
   ODE_IVP_constructor(&problem,2,0.0,10.0,0.001,0.0001,0.01);
   problem.x_0[1] = 0.0;
   problem.x_0[2] = 0.0;

   ODEsolution_constructor(&odeSol,1000,0.01,2);

    ODE_driver(problem.x_0,problem.N_var,problem.t_1,problem.t_2,
               problem.TOL,problem.Δ_1,problem.Δ_min,
               &problem.n_good,&problem.n_bad,
               f_ODE_IVP,rungeKutta_stepper);

    ODEsolution_printToFile(&odeSol);
    ODEsolution_print(&odeSol);

    ODE_IVP_destructor(&problem);
    ODEsolution_destructor(&odeSol);
    return 0;
}

void f_ODE_IVP(float t, float x[], float **xDot, int function_number)
{
        xDot[1][function_number] = 2*x[1] - 5.0;
        xDot[2][function_number] = x[1]*sin(x[2]);
}