/*
	Fifth-order Runge-Kutta Step Monitoring
	
	6/3/20
	Input:	y[1..n]				variable 
			dydx[1..n]			derivative dy/dx known @ x
			derivs(x,y,dydx)	derivative dy/dx @ x = f(t, y)
			h_try 				attempting stepping size
			ε					accuracy 
			yscal[1..n]			against which error is scaled
			
	Output:	y[1...n]			new y
			x[1..n]				new x
			h_did				accomplished stepsize 
			h_next 				estimated next stepsize
*/


#include <math.h>
#include "nrutil.h"
#include "prototype_declarations.h"

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(float y[], float dydx[], int n, float *x, float h_try, float ε, float y_scal[], float *h_did, float *h_next, void (*derivs)(float, float[], float[]))
{	
	int i;
	float err_max, h, h_temp, x_new, *y_err, *y_temp;
	
	y_err = vector(1, n);
	y_temp = vector(1, n);
	
	h = h_try;												// trial
	for(;;)
	{
		rkck(y, dydx, n, *x, h, y_err, y_temp, derivs);		// take a step 
		
		err_max = 0.0;										// eval. acc.
		for (i = 1;i<=n;i++) err_max = FMAX(err_max, fabs(y_err[i]/y_scal[i]));
		
		err_max /= ε;										// scaled relative to tolerance 
		if (err_max <= 1.0) break;							// step succeeded 
		
		h_temp = SAFETY*h*pow(err_max, PSHRNK);
		
		// truncations error too large
		// reduce stepsize 
		h = (h >= 0.0 ? FMAX(h_temp, 0.1*h) : FMIN(h_temp, 0.1*h));
		
		//no more than a factor of 10
		x_new = (*x) + h;
		if (x_new == *x) nrerror("stepsize underflow in rkqs");	
	}
	
	if (err_max > ERRCON) *h_next = SAFETY*h*pow(err_max, PGROW);
	else *h_next = 5.0*h;									// no more than a factor of 5
	
	*x += (*h_did = h);
	for (i = 1;i<=n;i++) y[i] = y_temp[i];
	
	free_vector(y_temp, 1, n);
	free_vector(y_err, 1, n);
}