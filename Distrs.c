//
//  Distrs.c
//  PDT-SIM
//
//  Created by Love Kildetoft on 2022-11-23.
//

#include "Distrs.h"
#include <math.h>
#include <stdio.h>

const float pi = 3.141592653589793;

float phaseFunc(float theta, float g) {
    float p = (1/(4*pi))*((1-pow(g, 2))/(pow(1+pow(g, 2)-2*g*cos(theta), 3/2)));
    return p;
};

float expFunc(float x, float mue) {
    float p = mue*exp(-mue*x);
    return p;
};
