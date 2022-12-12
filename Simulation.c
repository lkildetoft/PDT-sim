#include "Distrs.h"
#include "Simulation.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*
Author: Love Kildetoft
*/
float randChoice(float *randlist, float *weights, int len) {
    /*
    Generate weighted random numbers.  
    */
    int low = 0;
    int mid = 0;
    float wsum = 0;
    float cweights[len];
    
    for (int i = 0; i < len; i++) {
        wsum += weights[i];
        cweights[i] = wsum;
    }
    
    float rnd = ((float)rand()/(RAND_MAX))*wsum;
    
    while (low < len) {
        mid = (low + len - 1)/2;
        if (rnd == cweights[mid]) {
            break;
        } else if (rnd > cweights[mid]) {
            low = mid + 1;
        } else {
            len = mid - 1;
        }
    }
    return randlist[mid];
};

void genPhaseVals(float *plist, int size, float lim, float g) {
    /*
    Modify an empty float array to contain scattering phase function
    values. 
    */
    int n = 0;
    for (float i = -1*lim; i < lim; i+=(2*lim/size)) {
        plist[n++] = phaseFunc(i, g);
    }
};

void genExpVals(float *glist, int size, float lim, float mue) {
    /*
    Modify an empty float array to contain exponentially distributed
    path values.
    */
    int n = 0;
    for (float i = 0; i < lim; i+=(lim/size)) {
        glist[n++] = expFunc(i, mue);
    }
};

void interact(struct photon *p, struct medium *m, float *paths, float *expvals, float *angles, float *phasevals, int len) {
    /*
    Simulates the interaction of photons in a medium with appropriate parameters.
    A "photon" is absorbed and scattered at each call.
    */
    float w = p->w;
    
    float mua = m->mua;
    float mus = m->mus;
    float mue = mua + mus;
    
    float theta = randChoice(angles, phasevals, len);

    float l = -1*log(randChoice(paths, expvals, len))/(mue);
    float xs = l*sin(theta);
    float ys = l*cos(theta);

    p->w -= (mua/mue)*w;
    p->x += xs;
    p->y += ys;
};
