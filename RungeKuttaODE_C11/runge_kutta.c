/*
	Alonso Vega
    Embedded Runga-Kutta Method
	6/3/20
	
	Input:	x[1..n]				variable
			xDot[1..n]			derivative dx/dt
			f(t,x,xDot)	        derivative dx/dt
			Δ					step size
			
	Output:	x[1...n]			new x
			e[1..n]			    est. of local truncation error

*/

#include "nrutil.h"
#include "prototype_declarations.h"

#define ORDER_PLUS1 6

void rungeKutta(float *x_n, float *xDot, int N, float t, float Δ, float *x_nPlus1, float *e,
                void (*f)(float, float[], float[], int))
{
	int i_1 = 0, i_2 = 0, i = 0;
	static float *w_Δ, **W, *w, *w_star, *δ;

	w_Δ = vector(1,ORDER_PLUS1);
	w_Δ[1] = 0.0;
    w_Δ[2] = 0.2;
    w_Δ[3] = 0.3;
    w_Δ[4] = 0.6;
    w_Δ[5] = 1.0;
    w_Δ[6] = 0.875;

    W = matrix(1,ORDER_PLUS1,1,ORDER_PLUS1);
    LOOP(i_1,1,ORDER_PLUS1)
        LOOP(i_2,1,ORDER_PLUS1) W[i_1][i_2] = 0.0;
    W[2][1] = 0.2;
    W[3][1] = 0.5; W[3][2] = 0.2;
    W[4][1] = 0.2; W[4][2] = 0.2; W[4][3] = 0.2;
    W[5][1] = 0.2; W[5][2] = 0.2; W[5][3] = 0.2; W[5][4] = 0.2;
    W[6][1] = 0.2; W[6][2] = 0.2; W[6][3] = 0.2; W[6][4] = 0.2; W[6][5] = 0.2;

	w_star = vector(1,ORDER_PLUS1);
    w_star[1] = 0.0;
    w_star[2] = 0.2;
    w_star[3] = 0.3;
    w_star[4] = 0.6;
    w_star[5] = 1.0;
    w_star[6] = 0.875;

    w = vector(1,ORDER_PLUS1);
    w[1] = 0.0;
    w[2] = 0.2;
    w[3] = 0.3;
    w[4] = 0.6;
    w[5] = 1.0;
    w[6] = 0.875;


	float **F;
	F = matrix(1,N,1,ORDER_PLUS1);
    LOOP(i_1,1,ORDER_PLUS1)
        LOOP(i_2,1,N) F[i_1][i_2] = 0.0;
    LOOP(i,1,N) F[i][1] = xDot[i];

    float *x_temp;
    x_temp = vector(1, N);
    LOOP(i_1,2,ORDER_PLUS1)
    {
        LOOP(i_2, 1, N) x_temp[i_2] = x_n[i_2] + Δ*inner_product(W[i_1], F[i_2], i_1 - 1);
        (*f)(t + w_Δ[i_1] * Δ, x_temp, F, i_1);
    }

    LOOP(i, 1, N) x_nPlus1[i] = x_n[i] + Δ*inner_product(w,F[i],ORDER_PLUS1);

    δ = vector(1,ORDER_PLUS1);
    LOOP(i,1,ORDER_PLUS1) δ[i] = w[i] - w_star[i];
    LOOP(i,1, N) e[i] = Δ*inner_product(δ,F[i],ORDER_PLUS1);

    free_vector(w     ,1   ,ORDER_PLUS1);
    free_vector(w_star,1   ,ORDER_PLUS1);
    free_vector(w_Δ   ,1   ,ORDER_PLUS1);
	free_vector(x_temp,1   ,N);
	free_vector(δ     ,1   ,ORDER_PLUS1);
	free_matrix(F     ,1  ,N           ,1,ORDER_PLUS1);
    free_matrix(W     ,1  ,ORDER_PLUS1 ,1, ORDER_PLUS1);
}