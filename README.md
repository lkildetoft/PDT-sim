# PDT-sim
A Monte-Carlo simulation of the light distribution in photodynamic therapy (see for example https://www.cancer.gov/about-cancer/treatment/types/photodynamic-therapy). Based on MCML by Steven Jaques et. al (https://www.sciencedirect.com/science/article/pii/016926079501640F).
#### Performed as part of a project in the Medical Optics course at Lund University. Full report about the development can be found at
https://www.overleaf.com/read/vffgvwkdnzwj


# How to compile:
This program uses OMP for parallelization. 

An example of how to compile using gcc:
```
gcc Simulation.c main.c Distrs.c -o mc.o -fopenmp -O3
```
#### NOTE: Repository will be updated with more details as development progresses.
