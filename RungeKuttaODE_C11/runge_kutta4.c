/*
    Alonso Vega
	5/29/20
	Input:	x[1...n]			variable
			xDot[1...n]			derivative dx/dt
			Δ 					stepping size
			f(x,t)	            derivative dx/dt
	Output:	x_out[1...n]
*/

#include "nrutil.h"
#include "prototype_declarations.h"

#define N_w  4

void rungeKutta_4(float *x_n, float *xDot_n, int N, float t_n, float Δ, float *x_nPlus1, void (*f)(float, float [], float [], int))
{
	int i = 0, j = 0;
	float t_Δ, Δ_div2, *x_next, **xDot_next, *w, *xDot_vect;

	w = vector(1,N_w);
	LOOP(i,1,N_w) w[i] = 0.0;
	w[1] = 1.0/6.0;
    w[2] = 1.0/1.0;
    w[3] = 1.0/3.0;
    w[4] = 1.0/6.0;

	x_next    = vector(1, N);
	xDot_next = matrix(1, 3,1,N);

	Δ_div2  = Δ/2.0;		// = h/2 (half step)
	t_Δ = t_n + Δ_div2;		// = t + Δ/2

	LOOP(i, 1, N) x_next[i] = x_n[i] + Δ_div2 * xDot_n[i];		   //1st step
	(*f)(t_Δ, x_next, xDot_next, 2);								       //2nd step
    LOOP(i, 1, N) x_next[i] = x_n[i] + Δ_div2 * xDot_next[1][i];
	(*f)(t_Δ, x_next, xDot_next, 3);							           //3rd step
    LOOP(i, 1, N) x_next[i] = x_n[i] + Δ * xDot_next[2][i];
	(*f)(t_n + Δ, x_next, xDot_next, 4);								   //4th step

    xDot_vect = vector(1,N_w);
    LOOP(i, 1, N)
    {
        xDot_vect[1] = xDot_n[i];
        LOOP(j,1,N_w-1) xDot_vect[j+1] = xDot_next[j][i];

        x_nPlus1[i] = x_n[i] + Δ*inner_product(w, xDot_vect, N_w);
    }

	free_vector(x_next   ,1, N);
    free_vector(w        ,1,4);
    free_matrix(xDot_next,1,3,1,N);
}



