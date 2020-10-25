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

#define MAX_STEP   10000
#define SCALE_BIAS 1.0e-30
#define A_TOL       0.001
#define R_TOL       0.001

extern int K_max_main, K_main;
extern float *t_main, **X_main, dt_main;

void ODE_driver(float *x_0, int N_var, float t_1, float t_2, float TOL, float Δ_1, float Δ_min, int *n_good, int *n_bad,
                void (*f)(float, float [], float [], int),
                void (*rungeKutta_stepper)(float [], float [], int, float *, float, float, float [], float *, float *,
                        void (*)(float, float [], float [], int)))
{
	int i = 0, n_step = 0;
    float t_save, t_n, Δ_next, Δ_did, Δ;
	float *x_scale_n, *x_n, *xDot_n;

    x_scale_n = vector(1, N_var);
    x_n       = vector(1, N_var);
    xDot_n    = vector(1, N_var);

    t_n = t_1;
    Δ = SIGN(Δ_1, t_2 - t_1);
    *n_good = (*n_bad) = K_main = 0;
	LOOP(i, 1, N_var) x_n[i] = x_0[i];
	
	if (K_max_main > 0) t_save = t_n - dt_main * 2.0; // to save first step
	
	LOOP(n_step, 1, MAX_STEP)
	{
		(*f)(t_n, x_n, xDot_n, 0);                    // note: note not passing a matrix, passing xDot in R^N_vars

		///*------------------------------Using absolute and relative tolerance
        LOOP(i, 1, N_var) x_scale_n[i] = fabs(x_n[i]) + fabs(xDot_n[i] * Δ) + SCALE_BIAS;
        //*///--------------------------Using absolute and relative tolerance

        /*------------------------------Using absolute and relative tolerance
        if(K_max_main) LOOP(i,1,N_var) x_scale_n[i] = A_TOL + R_TOL * FMAX(fabs(x_n[i]), fabs(X_main[i][K_main]));
        else nrerror("NEED K_max_main != 0 IN ODE_driver.c");
        //*///--------------------------Using absolute and relative tolerance

        if (K_max_main > 0 && K_main < K_max_main - 1 && fabs(t_n - t_save) > fabs(dt_main))
		{
            t_main[++K_main] = t_n;
			LOOP(i, 1, N_var) X_main[i][K_main] = x_n[i];        // store intermediate results
            t_save = t_n;
		}
		if ((t_n + Δ - t_2) * (t_n + Δ - t_1) > 0.0) Δ = t_2 - t_n; // overshoot in [t_1, t_2] protection
		
		(*rungeKutta_stepper)(x_n, xDot_n, N_var, &t_n, Δ, TOL, x_scale_n, &Δ_did, &Δ_next, f);
		if (Δ_did == Δ) ++(*n_good); else ++(*n_bad);
		
		if ((t_n - t_2)*(t_2 - t_1) > 0.0)                          // done?
		{
			LOOP(i, 1, N_var) x_0[i] = x_n[i];                  // modify IC
			if (K_max_main)
			{
                t_main[++K_main] = t_n;
				LOOP(i, 1, N_var) X_main[i][K_main] = x_n[i];    // final result
			}

			free_vector(xDot_n   , 1, N_var);
			free_vector(x_n      , 1, N_var);
			free_vector(x_scale_n, 1, N_var);
			return;	
		}
		if (fabs(Δ_next) <= Δ_min) nrerror("STEP SIZE TOO SMALL IN ODE_driver.c");
        Δ = Δ_next;
	}
	nrerror("TOO MANY STEPS IN ODE_driver.c");
}