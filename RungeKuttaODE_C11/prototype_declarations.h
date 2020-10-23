//
// Created by Alonso Vega on 10/23/20.
//

#ifndef RUNGEKUTTAODE_C11_PROTOTYPE_DECLARATIONS_H
#define RUNGEKUTTAODE_C11_PROTOTYPE_DECLARATIONS_H

void rk4(float y[], float dydx[], int n, float x, float h, float yout[], void (*derivs)(float, float [], float []));

void rkck(float y[], float dydx[], int n, float x, float h, float yout[], float yerr[], void (*derivs)(float, float [], float []));

void rkqs(float y[], float dydx[], int n, float *x,float htry, float eps, float yscal[], float *hdid, float *hnext, void (*derivs)(float, float [], float []));

void odeint(float ystart[], int nvar, float x1, float x2,float eps, float h1, float hmin, int *nok, int *nbad,void (*derivs)(float, float [], float []),void (*rkqs)(float [], float [], int, float *, float, float,float [], float *, float *, void (*)(float, float [], float [])));

#endif //RUNGEKUTTAODE_C11_PROTOTYPE_DECLARATIONS_H


