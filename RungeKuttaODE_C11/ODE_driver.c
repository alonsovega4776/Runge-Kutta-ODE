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

#define MAX_STEP 10000
#define TINY     1.0e-30

extern int k_max, count;
extern float *t_p, **x_p, dt_save;

void ODE_driver(float *x_0, int N_var, float t_1, float t_2, float TOL, float Δ_1, float Δ_min, int *n_good, int *n_bad,
                void (*f)(float, float [], float [], int),
                void (*rungeKutta_stepper)(float [], float [], int, float *, float, float, float [], float *, float *,
                        void (*)(float, float [], float [], int)))
{
	int i = 0, n_step = 0;
    float t_save, t, Δ_next, Δ_did, Δ;
	float *x_scale, *x, *xDot;

    x_scale = vector(1, N_var);
    x       = vector(1, N_var);
    xDot    = vector(1, N_var);

    t = t_1;
    Δ = SIGN(Δ_1, t_2 - t_1);
    *n_good = (*n_bad) = count = 0;
	LOOP(i, 1, N_var) x[i] = x_0[i];
	
	if (k_max > 0) t_save = t - dt_save*2.0;
	
	LOOP(n_step, 1, MAX_STEP)
	{
		(*f)(t, x, xDot, 0);                    // note: note passing a matrix, passing xDot in R^N_vars
		LOOP(i, 1, N_var) x_scale[i] = fabs(x[i]) + fabs(xDot[i] * Δ) + TINY;
		if (k_max > 0 && count < k_max - 1 && fabs(t - t_save) > fabs(dt_save))
		{
            t_p[++count] = t;
			LOOP(i, 1, N_var) x_p[i][count] = x[i];
            t_save = t;
		}
		if ((t + Δ - t_2)*(t + Δ - t_1) > 0.0) Δ = t_2 - t;
		
		(*rungeKutta_stepper)(x, xDot, N_var, &t, Δ, TOL, x_scale, &Δ_did, &Δ_next, f);
		if (Δ_did == Δ) ++(*n_good); else ++(*n_bad);
		
		if ((t - t_2)*(t_2 - t_1) > 0.0)
		{
			LOOP(i, 1, N_var) x_0[i] = x[i];
			if (k_max)
			{
                t_p[++count] = t;
				LOOP(i, 1, N_var) x_p[i][count] = x[i];
			}

			free_vector(xDot   , 1, N_var);
			free_vector(x      , 1, N_var);
			free_vector(x_scale, 1, N_var);
			return;	
		}
		if (fabs(Δ_next) <= Δ_min) nrerror("STEP SIZE TOO SMALL IN ODE_driver.c");
        Δ = Δ_next;
	}
	nrerror("TOO MANY STEPS IN odeint");
}