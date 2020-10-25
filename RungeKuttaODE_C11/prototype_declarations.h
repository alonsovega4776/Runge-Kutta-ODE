//
// Created by Alonso Vega on 10/23/20.
//

#ifndef RUNGEKUTTAODE_C11_PROTOTYPE_DECLARATIONS_H
#define RUNGEKUTTAODE_C11_PROTOTYPE_DECLARATIONS_H

void rungeKutta_4(float *x_n, float *xDot_n, int N, float t_n, float Δ, float *x_nPlus1,
                  void (*f)(float, float [], float [], int));
void rungeKutta(float *x_n, float *xDot, int N, float t, float Δ, float *x_nPlus1, float *e,
                void (*f)(float, float [], float [], int));
void rungeKutta_stepper(float *x, float *xDot, int N, float *t, float Δ_try, float TOL, float *x_scale, float *Δ_did, float *Δ_next,
                        void (*f)(float, float [], float [], int));
void ODE_driver(float *x_0, int N_var, float t_1, float t_2, float TOL, float Δ_1, float Δ_min, int *n_ok, int *n_bad,
                void (*f)(float, float [], float [], int),
                void (*rungeKutta_stepper)(float [], float [], int, float *, float, float, float [], float *, float *,
                        void (*)(float, float [], float [], int)));

#endif //RUNGEKUTTAODE_C11_PROTOTYPE_DECLARATIONS_H


