//
//  Simulation.c
//  PDT-SIM
//
//  Created by Love Kildetoft on 2022-11-24.
//

#include "Distrs.h"
#include "Simulation.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float randChoice(float *randlist, float *weights, int len) {
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
    int n = 0;
    for (float i = -1*lim; i < lim; i+=(2*lim/size)) {
        plist[n++] = phaseFunc(i, g);
    }
};

void genExpVals(float *glist, int size, float lim, float mue) {
    int n = 0;
    for (float i = 0; i < lim; i+=(lim/size)) {
        glist[n++] = expFunc(i, mue);
    }
};

void interact(struct photon *p, struct medium *m, float *paths, float *expvals, float *angles, float *phasevals, int len) {
    float w = p->w;
    
    float mua = m->mua;
    float mus = m->mus;
    float mue = mua + mus;
    
    float theta = randChoice(angles, phasevals, len);
    //float phi = ((float)rand()/(RAND_MAX))*2*pi;
    float l = -1*log(randChoice(paths, expvals, len))/(mue);
    float xs = l*sin(theta);
    float ys = l*cos(theta);
    //float zs = l*cos(theta);

    p->w -= (mua/mue)*w;
    p->x += xs;
    p->y += ys;
    //p->z += zs;
    
};
