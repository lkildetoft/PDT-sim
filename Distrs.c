#include "Distrs.h"
#include <math.h>
#include <stdio.h>

/*
Author: Love Kildetoft
*/

const float pi = 3.141592653589793;

float phaseFunc(float theta, float g) {
    /*
    Henvey-Greenstein scattering phase function. 
    */
    float p = (1/(4*pi))*((1-pow(g, 2))/(pow(1+pow(g, 2)-2*g*cos(theta), 3/2)));
    return p;
};

float expFunc(float x, float mue) {
    /*
    Exponentially decaying distribution (see Beer-Lambert law).
    */
    float p = mue*exp(-mue*x);
    return p;
};
