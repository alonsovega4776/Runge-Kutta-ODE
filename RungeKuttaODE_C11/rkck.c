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

#include "nrutil.h"
#include "prototype_declarations.h"

#define LOOP(var, i, f) for (int var = i; var <= f; var++)

void rkck(float y[], float dydx[], int n, float x, float h, float y_out[], float y_err[], void (*dervis)(float, float[], float[]))
{
	int i;
	static float a_2 = 0.2, a_3 = 0.3, a_4 = 0.6, a_5 = 1.0, a_6 = 0.875; 
	
	static float b_21 = 0.2, \
				 b_31 = 3.0/40.0, b_32 = 9.0/40.0, \
				 b_41 = 0.3, b_42 = -0.9, b_43 = 1.2, \
				 b_51 = -11.0/54.0, b_52 = 2.5, b_53 = -70.0/27.0, b_54=35.0/27.0, \
				 b_61 = 1631.0/55296.0, b_62 = 175.0/512.0, b_63 = 575.0/13824.0, b_64 = 44275.0/110592.0, b_65 = 253.0/4096.0; 
					
	static float c_1 = 37.0/378.0, c_3 = 250.0/621.0, c_4 = 125.0/594.0, c_6 = 512.0/1771.0;
	
	float δ_1 = c_1-2825.0/27648.0, δ_3 = c_3-18575.0/48384.0, δ_4 = c_4-13525.0/55296.0, δ_5 = -277.00/14336.0, δ_6 = c_6-0.25;
	
	float *k_func_2, *k_func_3, *k_func_4, *k_func_5, *k_func_6, *y_temp;
	
	k_func_2=vector(1,n); 
	k_func_3=vector(1,n); 
	k_func_4=vector(1,n); 
	k_func_5=vector(1,n); 
	k_func_6=vector(1,n);
	y_temp = vector(1,n);
	
	LOOP(i, 1, n) y_temp[i] = y[i] + h*(b_21*dydx[i]);	
	(*dervis)(x + a_2*h, y_temp, k_func_2);						
	
	LOOP(i, 1, n) y_temp[i] = y[i] + h*(b_31*dydx[i] + b_32*k_func_2[i]);
	(*dervis)(x + a_3*h, y_temp, k_func_3);
	
	LOOP(i, 1, n) y_temp[i] = y[i] + h*(b_41*dydx[i] + b_42*k_func_2[i] + b_43*k_func_3[i]);
	(*dervis)(x + a_4*h, y_temp, k_func_4);
	
	LOOP(i, 1, n) y_temp[i] = y[i] + h*(b_51*dydx[i] + b_52*k_func_2[i] + b_53*k_func_3[i] + b_54*k_func_4[i]);
	(*dervis)(x + a_5*h, y_temp, k_func_5);
	
	LOOP(i, 1, n) y_temp[i] = y[i] + h*(b_61*dydx[i] + b_62*k_func_2[i] + b_63*k_func_3[i] + b_64*k_func_4[i] + b_65*k_func_5[i]);
	(*dervis)(x + a_6*h, y_temp, k_func_6);
	
	
	LOOP(i, 1, n) y_out[i] = y[i] + h*(c_1*dydx[i] + c_3*k_func_3[i] + c_4*k_func_4[i] + c_6*k_func_6[i]);
	LOOP(i, 1, n) y_err[i] = h*(δ_1*dydx[i] + δ_3*k_func_3[i] + δ_4*k_func_4[i] + δ_5*k_func_5[i] + δ_6*k_func_6[i]);
	
	free_vector(y_temp,1,n); 
	free_vector(k_func_6,1,n); 
	free_vector(k_func_5,1,n); 
	free_vector(k_func_4,1,n); 
	free_vector(k_func_3,1,n); 
	free_vector(k_func_2,1,n);
}