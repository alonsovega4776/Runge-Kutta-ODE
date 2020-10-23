cmake_minimum_required(VERSION 3.17)
project(RungeKuttaODE_C11 C)

set(CMAKE_C_STANDARD 11)

add_executable(RungeKuttaODE_C11 main.c prototype_declarations.h nrutil.h nrutil.c complex.c odeint.c rk4.c rkck.c rkqs.c)