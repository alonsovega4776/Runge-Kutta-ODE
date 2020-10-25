/*
	Runge-Kutta Step Monitoring
	
	6/3/20
	Input:	x[1..n]				variable
			xDot[1..n]			derivative dx/dt
			f(t,x,xDot)	        derivative dx/dt
			Δ_try 				attempting stepping size
			ε					accuracy 
			x_scale[1..n]		against which error is scaled
			
	Output:	x[1...n]			new x
			t[1..n]				new t
			Δ_did				accomplished stepsize
			Δ_next 				estimated next stepsize
*/


#include <math.h>
#include "nrutil.h"
#include "prototype_declarations.h"

#define S                   0.9
#define p                   4.0
#define EXP_stepControl     (-1.0/(p+1))
#define EXP_unitStepControl (-1.0/p)
#define ERRCON              pow(5.0/S, 1/EXP_stepControl)
#define FACTOR_MAX          5.0
#define FACTOR_MIN          (1.0/10.0)

void rungeKutta_stepper(float *x, float **xDot, int N, float *t, float Δ_try, float TOL,
                        float *x_scale, float *Δ_did, float *Δ_next,
                        void (*f)(float, float[], float**, int))
{	
	int i = 0;
	float ε, Δ, Δ_temp, Δ_max, t_new, *e, *x_temp, *ε_vect;

    e      = vector(1, N);
    x_temp = vector(1, N);
    ε_vect = vector(1, N);

    Δ = Δ_try;			// trial
	while(1)
	{
        rungeKutta(x, xDot, N, *t, Δ, e, x_temp, f);	            // take a step
		
		ε = 0.0;										            // eval. acc. error
        LOOP(i,1,N) ε_vect[i] = e[i] / x_scale[i];

        ///*---------------------------------------Worst Offender
        LOOP(i,1,N) ε = FMAX(ε, fabs(ε_vect[i]));           // get worst error
        //*/---------------------------------------Worst Offender

        /*---------------------------------------L2 Norm
        ε = inner_product(ε_vect, ε_vect, N);
        ε /= N;
        ε = sqrt(ε);
        //*///-------------------------------------L2 Norm

        ε /= TOL;										            // scaled relative to overall tolerance
        if (ε <= 1.0) break;							            // step succeeded
        //ELSE truncations error too large
        // reduce stepsize_________________________________________________________________________________
		Δ_temp = S*Δ*pow(ε, EXP_unitStepControl);
        Δ_max = FACTOR_MIN*Δ;
		Δ = (Δ >= 0.0 ? FMIN(Δ_temp, Δ_max) : FMAX(Δ_temp, Δ_max)); //no more than a factor of 1/10
		t_new = (*t) + Δ;
        // reduce stepsize_________________________________________________________________________________

		if (t_new == *t) nrerror("stepsize underflow in rungeKutta_stepper.c");
	}
    // increase stepsize_________________________________________________________________________________
	if (ε > ERRCON) *Δ_next = S*Δ*pow(ε, EXP_stepControl);
	else            *Δ_next = FACTOR_MAX*Δ;     					// no more than a factor of 5
	*t += (*Δ_did = Δ);
    // increase stepsize_________________________________________________________________________________

	LOOP(i,1,N) x[i] = x_temp[i];
	
	free_vector(x_temp, 1, N);
	free_vector(e     , 1, N);
	free_vector(ε_vect, 1, N);
}