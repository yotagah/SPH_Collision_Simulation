#include "datatypes.h"
#include "force.h"
#include "vector.h"

// Reset the acceleration of all particles
void resetForces(particles_t *parts) {

	int k, alfa;

	for(k = 0; k < parts->quant; k++) {
		for(alfa = 0; alfa < DIM; alfa++) {
			parts->particle[k].a[alfa] = 0;
		}
		parts->particle[k].dedt = 0;
	}
}

// Calculates the acceleration resulting from internal material forces
void internalForce(particles_t *parts, int_pairs_t *pairs) {
	int k, alfa, beta;
    int i, j; // Particles

    double aux_a, aux_e;

    double dv[DIM]; // Speed difference
    double rhoirhoj;

	for(k = 0; k < pairs->quant; k++) { // All interaction pairs

        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        rhoirhoj = parts->particle[i].rho * parts->particle[j].rho;

        subVector(parts->particle[i].v, parts->particle[j].v, dv); // Speed of i minus j calculated in dv

        aux_e = 0.0;
    	for(alfa = 0; alfa < DIM; alfa++) {

			aux_e += (parts->particle[i].p + parts->particle[j].p) * dv[alfa] * pairs->int_pair[k].dwdx[alfa] / rhoirhoj;  // Eq. 8.26 Liu 2003

    		aux_a = 0.0;
    		for(beta = 0; beta < DIM; beta++) {
    			aux_a += (parts->particle[i].sigma[alfa][beta] + parts->particle[j].sigma[alfa][beta]) * pairs->int_pair[k].dwdx[beta] / rhoirhoj; // Eq. 8.22 Liu 2003
    		}

    		parts->particle[i].a[alfa] += parts->particle[j].m * aux_a; // Eq. 8.22 Liu 2003
    		parts->particle[j].a[alfa] -= parts->particle[i].m * aux_a; // Eq. 8.22 Liu 2003
        }

   		parts->particle[i].dedt += 0.5 * parts->particle[j].m * aux_e; // Eq. 8.26 Liu 2003
   		parts->particle[j].dedt += 0.5 * parts->particle[i].m * aux_e; // Eq. 8.26 Liu 2003
    }

    for(k = 0; k < parts->quant; k++) {
    	for(alfa = 0; alfa < DIM; alfa++) {
	   		for(beta = 0; beta < DIM; beta++) {
    			parts->particle[k].dedt += parts->particle[k].tau[alfa][beta] * parts->particle[k].epsilon[alfa][beta] / parts->particle[k].rho; // Eq. 8.27 2003
    		}
    	}
    }
}

// Calculates acceleration due to external forces
void externalForce(particles_t *parts) {

    // Code for external forces not used in my simulation

}