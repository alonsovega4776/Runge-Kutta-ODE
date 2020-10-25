cmake_minimum_required(VERSION 3.17)
project(RungeKuttaODE_C11 C)

set(CMAKE_C_STANDARD 11)

add_executable(RungeKuttaODE_C11 main.c prototype_declarations.h nrutil.h nrutil.c complex.c ODE_driver.c runge_kutta4.c runge_kutta.c runge_kutta_stepper.c ODE.c)