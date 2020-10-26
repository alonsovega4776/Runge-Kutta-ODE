/*
    Alonso Vega
	Embedded Runga-Kutta Method Driver
	6/3/20
	
	Input:	x[1..n]				    variable
			xDot[1..n]			    derivative dx/dt
			f(t,x,xDot)	derivative  dx/dt
			Δ_1					    step size
			Δ_min				    step size
			TOL                     min allowable local error
			[t_1, t_2]              interval of integration
			
	Output:	x[1...n]			    new y
            e[1..n]			        est. of local truncation error
            n_good                  number of good steps taken
            n_bad                   number of bad steps taken

*/

#include <math.h>
#include "nrutil.h"
#include "prototype_declarations.h"
#include <stdio.h>

#define MAX_STEP   10000
#define SCALE_BIAS 1.0e-30
#define A_TOL       0.01
#define R_TOL       0.01

extern struct ODEsolution odeSol;

void ODE_driver(float *x_0, int N_var, float t_1, float t_2, float TOL, float Δ_1, float Δ_min, int *n_good, int *n_bad,
                void (*f)(float, float [], float**, int),
                void (*rungeKutta_stepper)(float [], float**, int, float *, float, float, float [], float *, float *,
                        void (*)(float, float [], float**, int)))
{
	int i = 0, n_step = 0;
    float t_save, t_n, Δ_next, Δ_did, Δ;
	float *x_scale_n, *x_n, **xDot_n;

    x_scale_n = vector(1, N_var);
    x_n       = vector(1, N_var);
    xDot_n    = matrix(1, N_var,1,1);

    t_n = t_1;
    Δ = SIGN(Δ_1, t_2 - t_1);
    *n_good = (*n_bad) = 0;
    odeSol.K = 0;
	LOOP(i, 1, N_var) x_n[i] = x_0[i];

	if (odeSol.K_max > 0) t_save = t_n - odeSol.dt * 2.0; // to save first step
	
	LOOP(n_step, 1, MAX_STEP)
	{
		(*f)(t_n, x_n, xDot_n, 1);                    // note: note not passing a matrix, passing xDot in R^N_vars

		/*------------------------------Using absolute and relative tolerance
        LOOP(i, 1, N_var) x_scale_n[i] = fabs(x_n[i]) + fabs(xDot_n[i][1] * Δ) + SCALE_BIAS;
        //*///--------------------------Using absolute and relative tolerance

        ///*------------------------------Using absolute and relative tolerance
        if(odeSol.K_max)
            LOOP(i,1,N_var) x_scale_n[i] = A_TOL + R_TOL * FMAX(fabs(x_n[i]), fabs(odeSol.X_matrix[i][odeSol.K]));
        else nrerror("NEED K_max_main != 0 IN ODE_driver.c");
        //*///--------------------------Using absolute and relative tolerance

        if (odeSol.K_max > 0 && odeSol.K < odeSol.K_max - 1 && fabs(t_n - t_save) > fabs(odeSol.dt))
		{
            odeSol.t_vector[++odeSol.K] = t_n;
			LOOP(i, 1, N_var) odeSol.X_matrix[i][odeSol.K] = x_n[i];        // store intermediate results
            t_save = t_n;
		}
		if ((t_n + Δ - t_2) * (t_n + Δ - t_1) > 0.0) Δ = t_2 - t_n; // overshoot in [t_1, t_2] protection
		
		(*rungeKutta_stepper)(x_n, xDot_n, N_var, &t_n, Δ, TOL, x_scale_n, &Δ_did, &Δ_next, f);
		if (Δ_did == Δ) ++(*n_good); else ++(*n_bad);
		
		if ((t_n - t_2)*(t_2 - t_1) >= 0.0)                          // done?
		{
			LOOP(i, 1, N_var) x_0[i] = x_n[i];                  // modify IC
			if (odeSol.K_max)
			{
                odeSol.t_vector[++odeSol.K] = t_n;
				LOOP(i, 1, N_var) odeSol.X_matrix[i][odeSol.K] = x_n[i];    // final result
			}

			free_matrix(xDot_n   , 1, N_var,1,1);
			free_vector(x_n      , 1, N_var);
			free_vector(x_scale_n, 1, N_var);
			return;	
		}
		if (fabs(Δ_next) <= Δ_min)
		{
		    printf("step number: %d\n", n_step);
            printf("time: %f\n"       , t_n);
            printf("next step: %f\n"  , Δ_next);
            nrerror("STEP SIZE TOO SMALL IN ODE_driver.c: ");
		}
        Δ = Δ_next;
	}
	nrerror("TOO MANY STEPS IN ODE_driver.c");
}