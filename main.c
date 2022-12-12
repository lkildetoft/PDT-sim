#include <stdio.h>
#include <stdlib.h>
#include "Distrs.h"
#include <time.h>
#include "Simulation.h"
#include <string.h>
#include <math.h>
#include <omp.h>

int main(int argc, const char * argv[]) {
    /* 
    
    PDT-SIM: A Monte-Carlo simulation of light propagation in tissues 
    Based on MCML by Steven Jaques et. al. 
    
    Author: Love Kildetoft 
    
    Usage: 

    niters: sets number of iterations for the simulation. Needs to be set for the simulation to run.
    IMPORTANT: This variable greatly affects output file size as the data is 
    currently saved as .txt files. Tread carefully!

    nvals: number of values to be used for paths/angles. 1000 is usually good so can be left as is
    with good conscience. 

    nphotons: sets number of photons to be used. 
    IMPORTANT: This variable greatly affects output file size as the data is 
    currently saved as .txt files. Tread carefully!

    Each medium is specified as a struct which is given a scattering coefficient, an absorption
    coefficient and their minimum/maximum positions. Tumour and sensitizer is assumed spherical.

    Each photon can be given a start position below. To see the exact usage, see the simulation.h
    header. 

    Outputs four text files containing the x and y positions of scattered and absorbed photons
    at each iteration respectively. This can then be read and plotted using the included Python 
    routines. 

    */
    int niters = 500; 
    int nvals = 1000;
    int nphotons = 1000000;
    srand((unsigned)time(NULL));
    
    FILE *xPosScat = fopen("/Users/lovekildetoft/MedOpts/PDT-SIM/xPosScat.txt", "w");
    FILE *yPosScat = fopen("/Users/lovekildetoft/MedOpts/PDT-SIM/yPosScat.txt", "w");

    FILE *xPosAbs = fopen("/Users/lovekildetoft/MedOpts/PDT-SIM/xPosAbs.txt", "w");
    FILE *yPosAbs = fopen("/Users/lovekildetoft/MedOpts/PDT-SIM/yPosAbs.txt", "w");
    
    photon *photons = calloc(nphotons, sizeof(photon));
    medium *tissue = &(medium){0.56, (9.2-0.9), 0, 3, 0, 4}; 
    medium *tumor = &(medium){1.3, (10-0.9), 0.5, 2.5, 1, 3};
    medium *sensitizer = &(medium){20, (10-0.9), 0.6, 2.4, 1.1, 2.9};
    
    float rTumor = (tumor->xmax - tumor->ymin)/2;
    float rSensi = (sensitizer->xmax - sensitizer->ymin)/2;

    float *paths = calloc(nvals, sizeof(float));
    float *theta = calloc(nvals, sizeof(float));
    float *dermiProb = calloc(nvals, sizeof(float));
    float *tumorProb = calloc(nvals, sizeof(float));
    float *sensiProb = calloc(nvals, sizeof(float));
    float *plist = calloc(nvals, sizeof(float));
    
    genExpVals(dermiProb, nvals, 1, tissue->mua + tissue->mus);
    genExpVals(tumorProb, nvals, 1, tumor->mua + tumor->mus);
    genExpVals(sensiProb, nvals, 1, sensitizer->mua + sensitizer->mus);

    genPhaseVals(plist, nvals, pi, 0.9);
    
    int n = 0;
    for (float i = 0; i < 1; i+=(1/(float)nvals)) {
        paths[n++] = i;
    }

    int k = 0;
    for (float i = -1*pi; i < pi; i+=(2*pi)/nvals) {
        theta[k++] = i;
    }

    for (int i = 0; i < nphotons; i++) {
        photons[i] = (photon){1.5, 1.1, 1, 0};
    }

    printf("Ready to go, with %d photons \n", nphotons);

    float *xPosScatTmp = calloc(nphotons, sizeof(float));
    float *yPosScatTmp = calloc(nphotons, sizeof(float));
    float *xPosAbsTmp = calloc(nphotons, sizeof(float));
    float *yPosAbsTmp = calloc(nphotons, sizeof(float));

    for (int i = 0; i < niters; i++) {
        printf("Running simulation, %f %% done \n", ((float)i/niters)*100);

        #pragma omp parallel for
        for (int q = 0; q < nphotons; q++) {
            photon *p;
            p = &photons[q];

            if ((p->w > 0.2) && (p->outside == 0)) {
                xPosScatTmp[q] = p->x;
                yPosScatTmp[q] = p->y;
                if (pow((p->x - (sensitizer->xmax)/2), 2) + pow((p->y - (sensitizer->ymax)/2), 2) <=  pow(rSensi, 2)) {
                    interact(p, sensitizer, paths, sensiProb, theta, plist, nvals);
                } else if (pow((p->x - (tumor->xmax)/2), 2) + pow((p->y - (tumor->ymax)/2), 2) <=  pow(rTumor, 2)) {
                    interact(p, tumor, paths, tumorProb, theta, plist, nvals);
                } else {
                    interact(p, tissue, paths, dermiProb, theta, plist, nvals);
                }

                if ((p->x < tissue->xmin || p->x > tissue->xmax) || (p->y < tissue->ymin || p->y > tissue->ymax)) {
                    p->outside = 1;
                } else {}
            }
            
            if ((p->w <= 0.2) && (p->outside == 0)) {
                xPosAbsTmp[q] = p->x;
                yPosAbsTmp[q] = p->y;
            }
        }
        for (int i = 0; i < nphotons; i++) {
            fprintf(xPosScat, "%f, ", xPosScatTmp[i]);
            fprintf(yPosScat, "%f, ", yPosScatTmp[i]);
            fprintf(xPosAbs, "%f, ", xPosAbsTmp[i]);
            fprintf(yPosAbs, "%f, ", yPosAbsTmp[i]);
        }
        fprintf(xPosScat, "\n");
        fprintf(yPosScat, "\n");
        fprintf(xPosAbs, "\n");
        fprintf(yPosAbs, "\n");
        
        memset(xPosScatTmp, 0, nphotons*sizeof(float));
        memset(yPosScatTmp, 0, nphotons*sizeof(float));
        memset(xPosAbsTmp, 0, nphotons*sizeof(float));
        memset(yPosAbsTmp, 0, nphotons*sizeof(float));
    }
    printf("Done \n");

    free(photons);
    free(paths);
    free(theta);
    free(dermiProb);
    free(tumorProb);
    free(sensiProb);
    free(plist);
    free(xPosScatTmp);
    free(yPosScatTmp);
    free(xPosAbsTmp);
    free(yPosAbsTmp);

    fclose(xPosScat);
    fclose(yPosScat);
    fclose(xPosAbs);
    fclose(yPosAbs);

    return 0;
}
