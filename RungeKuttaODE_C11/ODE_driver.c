/*
    Alonso Vega
	Embedded Runga-Kutta Method Driver
	6/3/20
	
	Input:	x[1..n]				    variable
			xDot[1..n]			    derivative dx/dt
			f(t,x,xDot)	derivative  dx/dt
			Δ					    step size
			
	Output:	x[1...n]			    new y
            e[1..n]			        est. of local truncation error

*/

#include <math.h>
#include "nrutil.h"
#include "prototype_declarations.h"

#define MAX_STEP 10000
#define TINY     1.0e-30

extern int k_max, kount;
extern float *x_p, **y_p, dx_save;

void ODE_driver(float *x_0, int N_var, float t_1, float t_2, float TOL, float Δ_1, float Δ_min, int *n_ok, int *n_bad,
                void (*f)(float, float [], float [], int),
                void (*rungeKutta_stepper)(float [], float [], int, float *, float, float, float [], float *, float *,
                        void (*)(float, float [], float [], int)))
{
	int i = 0, n_step = 0;
    float t_save, t, Δ_next, Δ_did, Δ;
	float *y_scale, *y, *dydx;
	
	y_scale = vector(1, N_var);
	y = vector(1, N_var);
	dydx = vector(1, N_var);

    t = t_1;
    Δ = SIGN(Δ_1, t_2 - t_1);
    *n_ok = (*n_bad) = kount = 0;
	
	LOOP(i, 1, N_var) y[i] = x_0[i];
	
	if (k_max > 0) t_save = t - dx_save * 2.0;
	
	LOOP(n_step, 1, MAX_STEP)
	{
		(*f)(t, y, dydx, 1);
		LOOP(i, 1, N_var) y_scale[i] = fabs(y[i]) + fabs(dydx[i] * Δ) + TINY;
		if (k_max > 0 && kount < k_max - 1 && fabs(t - t_save) > fabs(dx_save))
		{
			x_p[++kount] = t;
			LOOP(i, 1, N_var) y_p[i][kount] = y[i];
            t_save = t;
		}
		if ((t + Δ - t_2) * (t + Δ - t_1) > 0.0) Δ = t_2 - t;
		
		(*rungeKutta_stepper)(y, dydx, N_var, &t, Δ, TOL, y_scale, &Δ_did, &Δ_next, f);
		
		if (Δ_did == Δ) ++(*n_ok); else ++(*n_bad);
		
		if ((t - t_2) * (t_2 - t_1) > 0.0)
		{
			LOOP(i, 1, N_var) x_0[i] = y[i];
			if (k_max)
			{
				x_p[++kount] = t;
				LOOP(i, 1, N_var) y_p[i][kount] = y[i];
			}
			free_vector(dydx, 1, N_var);
			free_vector(y, 1, N_var);
			free_vector(y_scale, 1, N_var);
			return;	
		}
		if (fabs(Δ_next) <= Δ_min) nrerror("STEP SIZE TOO SMALL IN odeint");
        Δ = Δ_next;
	}
	nrerror("TOO MANY STEPS IN odeint");
}