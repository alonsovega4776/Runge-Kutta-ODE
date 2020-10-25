//
// Created by Alonso Vega on 10/25/20.
//

#include "prototype_declarations.h"

void ODEsolution_constructor(struct ODEsolution *odeSol, int maxCount, float precision, int x_dimension)
{
    odeSol->K_max = maxCount;
    odeSol->K = 0;
    odeSol->dt = precision;
    odeSol->x_dim = x_dimension;

    odeSol->t_vector = vector(1,maxCount);
    odeSol->X_matrix = matrix(1,x_dimension,1,maxCount);
}
void ODEsolution_destructor(struct ODEsolution *odeSol)
{
    free_vector(odeSol->t_vector, 1, odeSol->K_max);
    free_matrix(odeSol->X_matrix,1,odeSol->x_dim,1,odeSol->K_max);
}

void ODE_IVP_constructor(struct ODE_IVP *odeIvp, int x_dimension, float inital_t, float final_t,float initialStep, float minStep, float tolerance)
{
    odeIvp->x_0 = vector(1,x_dimension);

    odeIvp->N_var = x_dimension;
    odeIvp->t_1 = inital_t;
    odeIvp->t_2 = final_t;
    odeIvp->Δ_1 = initialStep;
    odeIvp->Δ_min = minStep;
    odeIvp->TOL = tolerance;
}
void ODE_IVP_destructor(struct ODE_IVP *odeIvp)
{
    free_vector(odeIvp->x_0,1,odeIvp->N_var);
}