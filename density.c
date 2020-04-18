#include <stdlib.h>
#include "datatypes.h"
#include "density.h"
#include "kernel.h"
#include "vector.h"

// Calculates the normalized density for each particle (Equation 4.35 Liu 2003 - Randies and Libersky, 1996; Chen et al., 1999a; 2000)
void normSumDensity(particles_t *parts, int_pairs_t *pairs) {

    int i, j, k;

    double *wi = (double*) malloc(parts->quant*sizeof(double)); // For normalization

    // First calculates the integration of the weight function over the space
    for(k = 0; k < parts->quant; k++) {
        wi[k] = weight(0.0) * parts->particle[k].m / parts->particle[k].rho;
    }
    for(k = 0; k < pairs->quant; k++) {
        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;
        wi[i] += pairs->int_pair[k].w * parts->particle[i].m / parts->particle[i].rho;
        wi[j] += pairs->int_pair[k].w * parts->particle[j].m / parts->particle[j].rho;
    }
    
    // Second, it calculates the integration of density over space
    sumDensity(parts, pairs);

    // Third, calculate the normalized density
    for(k = 0; k < parts->quant; k++) {
        parts->particle[k].rho = parts->particle[k].rho/wi[k];
    }

    free(wi);
}

// Calculates density by simple sum (Equation 4.26 Liu 2003)
void sumDensity(particles_t *parts, int_pairs_t *pairs) {

    int k;
    int i, j; // Particles

    for(k = 0; k < parts->quant; k++) { // Sum the particle itself
        parts->particle[k].rho = weight(0.0) * parts->particle[k].m;
    }

    for(k = 0; k < pairs->quant; k++) { // Sums the other interacting particles

        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;
        
        parts->particle[i].rho += parts->particle[i].m * pairs->int_pair[k].w;
        parts->particle[j].rho += parts->particle[j].m * pairs->int_pair[k].w;
    }
}

// Calculates the rate of variation of density by continuity of the density (Equation 4.31 Liu 2003)
void conDensity(particles_t *parts, int_pairs_t *pairs) {

    int k, beta;
    int i, j; // Particles
    double dv; // Speed difference
    double vcc; // Speed differential

    for(k = 0; k < parts->quant; k++) { // Reset the rate of variation of density in all particles
        parts->particle[k].drhodt = 0.0;
    }

    for(k = 0; k < pairs->quant; k++) { // Calculates the rate of variation of density in all particles

        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        // Calculates the speed differential (vcc)
        vcc = 0.0;
        for(beta = 0; beta < DIM; beta++) {
            dv = parts->particle[i].v[beta] - parts->particle[j].v[beta];
            vcc += dv*pairs->int_pair[k].dwdx[beta];
        }

        parts->particle[i].drhodt += parts->particle[j].m * vcc;
        parts->particle[j].drhodt += parts->particle[i].m * vcc;
    }
}

// Copies the current density to the initial density
void copy2InitDensity(particles_t *parts) {

    int k;

    for(k = 0; k < parts->quant; k++) { // Reset the rate of variation of density in all particles
        parts->particle[k].rho_0 = parts->particle[k].rho;
    }
}