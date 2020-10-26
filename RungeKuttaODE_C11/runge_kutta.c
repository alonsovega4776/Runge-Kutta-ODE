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


void rungeKutta(float *x_n, float **xDot, int N, float t, float Δ, float *x_nPlus1, float *e,
                void (*f)(float, float[], float**, int))
{
	int i_1 = 0, i_2 = 0, i = 0;
	static float *w_Δ, **W, *w, *w_star, *δ;

	w_Δ = vector(1,ORDER_PLUS1);
	w_Δ[1] = 0.0;
    w_Δ[2] = 1.0/5.0;
    w_Δ[3] = 3.0/10.0;
    w_Δ[4] = 3.0/5.0;
    w_Δ[5] = 1.0;
    w_Δ[6] = 7.0/8.0;

    W = matrix(1,ORDER_PLUS1,1,ORDER_PLUS1);
    LOOP(i_1,1,ORDER_PLUS1)
        LOOP(i_2,1,ORDER_PLUS1) W[i_1][i_2] = 0.0;
    W[2][1] = 1.0/5.0;
    W[3][1] = 3.0/40.0;       W[3][2] = 9.0/40.0;
    W[4][1] = 3.0/10.0;       W[4][2] = -9.0/10.0;   W[4][3] = 6.0/5.0;
    W[5][1] = -11.0/54.0;     W[5][2] = 5.0/2.0;     W[5][3] = -70.0/27.0;    W[5][4] = 35.0/27.0;
    W[6][1] = 1631.0/55296.0; W[6][2] = 175.0/512.0; W[6][3] = 575.0/13824.0; W[6][4] = 44275.0/110592.0; W[6][5] = 253.0/4096.0;

	w_star = vector(1,ORDER_PLUS1);
    w_star[1] = 2825.0/27648.0;
    w_star[2] = 0.0;
    w_star[3] = 18575.0/48384.0;
    w_star[4] = 13525.0/55396.0;
    w_star[5] = 277.0/14336.0;
    w_star[6] = 1.0/4.0;

    w = vector(1,ORDER_PLUS1);
    w[1] = 37.0/378.0;
    w[2] = 0.0;
    w[3] = 250.0/621.0;
    w[4] = 125.0/594.0;
    w[5] = 0.0;
    w[6] = 512.0/1771.0;


	float **F;
	F = matrix(1,N,1,ORDER_PLUS1);
    LOOP(i_1,1,N)
        LOOP(i_2,1,ORDER_PLUS1) F[i_1][i_2] = 0.0;
    LOOP(i,1,N) F[i][1] = xDot[i][1];

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