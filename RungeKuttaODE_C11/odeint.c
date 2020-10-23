/*
	Cash-Karp Parameters for Embedded Runga-Kutta Method
	6/3/20
	
	Input:	y[1..n]				variable 
			dydx[1..n]			derivative dy/dx known @ x
			derivs(x,y,dydx)	derivative dy/dx @ x = f(t, y)
			h					step size 
			
	Output:	y[1...n]			new y
			y_err[1..n]			est. of local truncation error 

*/

#include <math.h>
#include "nrutil.h"
#include "prototype_declarations.h"

#define MAX_STEP 10000
#define TINY 1.0e-30
#define LOOP(var, i, f) for (int var = i; var <= f; var++)

extern int k_max, kount;
extern float *x_p, **y_p, dx_save;

void odeint(float y_start[], int n_var, float x_1, float x_2, float ε, float h_1, float h_min, int *nok, int *n_bad, void (*dervis)(float, float [], float []), void (*rkqs)(float [], float [], int, float *, float, float, float [], float *, float *, void (*)(float, float [], float [])))
{
	float x_save, x, h_next, h_did, h;
	float *y_scale, *y, *dydx;
	
	y_scale = vector(1, n_var);
	y = vector(1, n_var);
	dydx = vector(1, n_var);
	
	x = x_1;
	h = SIGN(h_1, x_2 - x_1);
	*nok = (*n_bad) = kount = 0;
	
	LOOP(i, 1, n_var) y[i] = y_start[i];
	
	if (k_max > 0) x_save = x - dx_save*2.0;
	
	LOOP(n_step, 1, MAX_STEP)
	{
		(*dervis)(x, y, dydx);
		LOOP(i, 1, n_var) y_scale[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
		if (k_max > 0 && kount < k_max - 1 && fabs(x - x_save) > fabs(dx_save))
		{
			x_p[++kount] = x;
			LOOP(i, 1, n_var) y_p[i][kount] = y[i];
			x_save = x;
		}
		if ((x + h - x_2)*(x + h - x_1) > 0.0) h = x_2 - x;
		
		(*rkqs)(y, dydx, n_var, &x, h, ε, y_scale, &h_did, &h_next, dervis);
		
		if (h_did == h) ++(*nok); else ++(*n_bad);
		
		if ((x - x_2)*(x_2 - x_1) > 0.0)
		{
			LOOP(i, 1, n_var) y_start[i] = y[i];
			if (k_max)
			{
				x_p[++kount] = x;
				LOOP(i, 1, n_var) y_p[i][kount] = y[i];
			}
			free_vector(dydx, 1, n_var); 
			free_vector(y, 1, n_var); 
			free_vector(y_scale, 1, n_var); 
			return;	
		}
		if (fabs(h_next) <= h_min) nrerror("STEP SIZE TOO SMALL IN odeint");
		h = h_next; 	
	}
	nrerror("TOO MANY STEPS IN odeint");
}