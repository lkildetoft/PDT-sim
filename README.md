# PDT-sim
A Monte-Carlo simulation of the light distribution in photodynamic therapy. Based on MCML by Steven Jaques et. al.

# How to compile:
This program uses OMP for parallelization. 
An example of how to compile using gcc:

gcc Simulation.c main.c Distrs.c -o mc.o -fopenmp -O3

#### NOTE: Repository will be updated with more details as development progresses.
