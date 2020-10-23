#include "nrutil.h"
#include "prototype_declarations.h"

/*
	5/29/20
	Input:	y[1...n]			variable 
			dydx[1...n]			derivative dy/dx known @ x
			h 					stepping size
			derivs(x,y,dydx)	derivative dy/dx @ x
	Output:	yout[1...n]
*/

void rk4(float y[], float dydx[], int n, float x, float h, float yout[], void (*derivs)(float, float [], float []))
{
	int i;
	float x_half, h_half, h_6, *diff_y_next, *diff_y_next_2, *y_next;
	
	diff_y_next_2 = vector(1, n);		//
	diff_y_next = vector(1, n);			// derivative of y_next 
	y_next = vector(1, n);				// = y_i + k_1/2
	
	h_half = 0.5*h;						// = h/2 (half step)
	h_6 = h/6.0;						// = h/6
	
	x_half = x + h_half;				// = x + h/2
	
	for (i=1;i <= n; i++) y_next[i] = y[i] + h_half*dydx[i];			//1st step
	
	(*derivs)(x_half, y_next, diff_y_next);								//2nd step 
	
	for (i = 1; i <= n; i++) y_next[i] = y[i] + h_half*diff_y_next[i];
	
	(*derivs)(x_half, y_next, diff_y_next_2);							//3rd step 
		
	for (i = 1; i <= n; i++) 
	{
		y_next[i] = y[i] + h*diff_y_next_2[i];
		diff_y_next_2[i] += diff_y_next[i];
	}
	(*derivs)(x+h, y_next, diff_y_next);								//4th step 
	
	for (i = 1; i <= n; i++) 
	{
		yout[i] = y[i] + h_6*(dydx[i] + diff_y_next[i] + 2.0*diff_y_next_2[i]);
	}
	
	free_vector(y_next, 1, n);
	free_vector(diff_y_next, 1, n);
	free_vector(diff_y_next_2, 1, n);	
}

