//
//  Simulation.h
//  PDT-SIM
//
//  Created by Love Kildetoft on 2022-11-23.
//

#ifndef Simulation_h
#define Simulation_h

typedef struct photon {
    float x;
    float y;
    //float z;
    float w;
    int outside;
} photon;

typedef struct medium {
    float mua;
    float mus;
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    //float zmin;
    //float zmax;
} medium;


float randChoice(float *randlist, float *weights, int len);
void genPhaseVals(float *plist, int size, float lim, float g);
void genExpVals(float *glist, int size, float lim, float mu);
void interact(struct photon *p, struct medium *m, float *paths, float *expvals, float *angles, float *phasevals, int len);
void reflect(struct photon *p, struct medium *m, float *paths, float *expvals, float *angles, float *phasevals, int len);
#endif /* Simulation_h */
