//
// Created by Alonso Vega on 10/25/20.
//

#include "prototype_declarations.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//----------------------------------------------------ODEsolution-------------------------------------------------------
void ODEsolution_constructor(struct ODEsolution *odeSol, int maxCount, float precision, int x_dimension)
{
    odeSol->K_max = maxCount;
    odeSol->K = 1;
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

void ODEsolution_print(struct ODEsolution *odeSol)
{
    int i_1 = 0, i_2 = 0;

    printf("\n t: ");
    LOOP(i_1, 1, odeSol->K_max) printf("%f, ", odeSol->t_vector[i_1]);

    printf("\n");
    LOOP(i_1,1,odeSol->x_dim)
    {
        printf("x_%d: ", i_1);
        LOOP(i_2, 1, odeSol->K_max)
        {
            printf("%f, ", odeSol->X_matrix[i_1][i_2]);
        }
        printf("\n");
    }
}

void ODEsolution_printToFile(struct ODEsolution *odeSol)
{
    int i_1 = 0, i_2 = 0;

    FILE *f_out;
    f_out = fopen("output.txt", "w");

    fprintf(f_out, "t: ");
    LOOP(i_1, 1, odeSol->K_max) fprintf(f_out, "%f, ", odeSol->t_vector[i_1]);

    fprintf(f_out, "\n");
    LOOP(i_1,1,odeSol->x_dim)
    {
        fprintf(f_out, "x_%d: ", i_1);
        LOOP(i_2, 1, odeSol->K_max)
        {
            fprintf(f_out,"%f, ", odeSol->X_matrix[i_1][i_2]);
        }
        fprintf(f_out,"\n");
    }


    if (fclose(f_out) != 0) nrerror("error in closing files");

}
//----------------------------------------------------ODEsolution-------------------------------------------------------



//------------------------------------------------------ODE_IVP---------------------------------------------------------
void ODE_IVP_constructor(struct ODE_IVP *odeIvp, int x_dimension,
        float inital_t, float final_t,float initialStep, float minStep, float tolerance)
{
    odeIvp->x_0 = vector(1,x_dimension);

    odeIvp->N_var = x_dimension;
    odeIvp->t_1 = inital_t;
    odeIvp->t_2 = final_t;
    odeIvp->Δ_1 = initialStep;
    odeIvp->Δ_min = minStep;
    odeIvp->TOL = tolerance;
    odeIvp->n_bad = 0;
    odeIvp->n_good = 0;
}

void ODE_IVP_destructor(struct ODE_IVP *odeIvp)
{
    free_vector(odeIvp->x_0,1,odeIvp->N_var);
}
//------------------------------------------------------ODE_IVP---------------------------------------------------------