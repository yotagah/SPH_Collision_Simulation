#include <stdlib.h>
#include <math.h>
#include "datatypes.h"
#include "vector.h"
#include "stress.h"

// Calculates the total tension tensor (sigma)
void totalStressTensor(particles_t *parts, int_pairs_t *pairs, double mi) {

	// mi(G): shear module
	// J_0: known von Mieses yield stress

	int k, i, j, alfa, beta, gamma;

	double aux_epsilon;
	double aux_R;
	double dv[DIM];
	double mprhoi, mprhoj;

	// Calculates the rotation and strain rate tensors (Equations 8.24 and 8.25, Liu 2003)

	for(k = 0; k < parts->quant; k++) {
		for(alfa = 0; alfa < DIM; alfa++) { // alfa coordinate
			for(beta = 0; beta < DIM; beta++) { // beta coordinate
				parts->particle[k].epsilon[alfa][beta] = 0;
				parts->particle[k].R[alfa][beta] = 0;
			}
		}
	}

	for(k = 0; k < pairs->quant; k++) { // All iteraction pairs

		i = pairs->int_pair[k].i;
		j = pairs->int_pair[k].j;

		mprhoi = (parts->particle[i].m / parts->particle[i].rho)*0.5;
		mprhoj = (parts->particle[j].m / parts->particle[j].rho)*0.5;

		aux_epsilon = 0.0; // Initialy zerates the sum
		aux_R = 0.0; // Initialy zerates the sum

   		subVector(parts->particle[j].v, parts->particle[i].v, dv); // Speed difference to dv variable

		for(alfa = 0; alfa < DIM; alfa++) { // alfa coordinate
			for(beta = 0; beta < DIM; beta++) { // beta coordinate

	       		// Partial sums
        		aux_epsilon = dv[alfa]*pairs->int_pair[k].dwdx[beta] + dv[beta]*pairs->int_pair[k].dwdx[alfa];
        		aux_R = dv[alfa]*pairs->int_pair[k].dwdx[beta] - dv[beta]*pairs->int_pair[k].dwdx[alfa];

				// Sum in the total
				parts->particle[i].epsilon[alfa][beta] += mprhoj*aux_epsilon;
				parts->particle[j].epsilon[alfa][beta] += mprhoi*aux_epsilon;
				parts->particle[i].R[alfa][beta] += mprhoj*aux_R;
				parts->particle[j].R[alfa][beta] += mprhoi*aux_R;
			}
		}
	}

	// Calculate tau dot, Von Mises limit stress and finally Sigma

	double avarage;

	for(k = 0; k < parts->quant; k++) { // All particles

		// Main diagonal average (Traceless part of epsilon)
		avarage = 0;
		for(alfa = 0; alfa < DIM; alfa++) {
			avarage += parts->particle[k].epsilon[alfa][alfa];
		}
		avarage /= (double) DIM;

		for(alfa = 0; alfa < DIM; alfa++) { // alfa coordinate
			for(beta = 0; beta < DIM; beta++) { // beta coordinate
				
				// Calculates the tensor shear stress rate (tau_dot)

				if(beta == alfa) { // Eq. 8.5 (8.3) Liu 2003
					parts->particle[k].tau_dot[alfa][beta] = mi * (parts->particle[k].epsilon[alfa][beta] - avarage);
				} else {
					parts->particle[k].tau_dot[alfa][beta] = mi * parts->particle[k].epsilon[alfa][beta];
				}
				for(gamma = 0; gamma < DIM; gamma++) { // coordenate gamma Eq. 8.5 Liu 2003
					parts->particle[k].tau_dot[alfa][beta] += parts->particle[k].tau[alfa][gamma]*parts->particle[k].R[beta][gamma]	+ parts->particle[k].tau[beta][gamma]*parts->particle[k].R[alfa][gamma];
				}

				// Eq. 8.2 Liu 2003
				parts->particle[k].sigma[alfa][beta] = parts->particle[k].tau[alfa][beta];
				if(alfa == beta) { // Kroeniger delta
					parts->particle[k].sigma[alfa][beta] -= parts->particle[k].p;
				}
			}
		}
	}
}

// Check that the traction does not exceed the limit and scale it to avoid it (von Mieses, Liu 2003 eq. 8.8)
void plasticYieldModel(particles_t *parts, double J_0) {

	int k, alfa, beta;
	double J;

	for(k = 0; k < parts->quant; k++) { // All particles

		J = 0;
		for(alfa = 0; alfa < DIM; alfa++) { // alfa coordinate
			for(beta = 0; beta < DIM; beta++) { // beta coordinate
				// For the calculation of the von Mises yield stress (Liu 2003 Eq. 8.7)
				J += parts->particle[k].tau[alfa][beta]*parts->particle[k].tau[alfa][beta]; // Eq. 8.7 Liu 2003
			}
		}

		J = sqrt(J);

		if(J > J_0) {
			for(alfa = 0; alfa < DIM; alfa++) { // alfa coordinate
				for(beta = 0; beta < DIM; beta++) { // beta coordinate
					parts->particle[k].tau[alfa][beta] *= sqrt(J_0/(3*J*J)); // Eq. 8.8 Liu 2003
				}
			}
		}
	}
}