#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "datatypes.h"
#include "particles.h"
#include "integration.h"
#include "stress.h"

double norm; // Normalization of the weight function
double dt; // Integration interval

int main() {
	int i, j, alfa, beta;

	double t;
	double h_size; // Size of the interaction area between particles

	particles_t particles; // Partícles
	int_pairs_t int_pairs; // Interaction pairs

	particles = newParticles();
	int_pairs = newIntPairs();

	srand48(time(NULL)); // No simulation will be the same

	int c; // Counter for selection of "frames to print"

	h_size = 0.00038*1.2; // Define the interaction size

	// The normalization constant depends on the number of dimensions of the problem and the weight function choosed
	if(DIM == 2)
		norm = 15.0/(7.0*M_PI*h_size*h_size);
	else if(DIM == 3)
		norm = 3.0/(2.0*M_PI*h_size*h_size*h_size);

	// Create the material that will collide
	for(j = 0; j < 67; j++) {
		for(i = 0; i < 20; i++) {
			addParticle(&particles, i*0.00038 - 0.00361, j*0.00038 + 0.00019, 0.0, -221.0, 1.13354e-3, 7850, h_size);
		}
	}

	// Create the wall wher the material will collide (virtual particles are not integrated)
	for(j = 0; j < 5; j++) {
		for(i = 0; i < 76; i++) {
			addParticle(&particles, i*0.00038 - 0.01425, j*0.00038 - 0.00171, 0.0, 0.0, 1.13354e-3, 7850, h_size);
			particles.particle[particles.quant-1].virt = 1;
		}
	}

	dt = 0.00000001; // Define the integration step
	t = 0.0; // Initial time

	// Do the first step to calculate forces
	singleStep(&particles, &int_pairs);

	// Calculate initial values for begin with LeapFrog
	for(i = 0; i < particles.quant; i++) {
		if(!particles.particle[i].virt) {
			particles.particle[i].rho += dt*particles.particle[i].drhodt/2.0;
			particles.particle[i].e += dt*particles.particle[i].dedt/2.0;
			for(alfa = 0; alfa < DIM; alfa++) {
				particles.particle[i].v[alfa] += dt*particles.particle[i].a[alfa]/2.0 + particles.particle[i].av[alfa];
				particles.particle[i].r[alfa] += dt*particles.particle[i].v[alfa];
				for(beta = 0; beta < DIM; beta++) {
					particles.particle[i].tau[alfa][beta] += dt*particles.particle[i].tau_dot[alfa][beta]/2.0;
				}
			}
		}
	}

	c = 0; // Counter for print
	for(t = 0.0 + dt; t < 0.0001; t += dt) {

		for(i = 0; i < particles.quant; i++) {
			if(!particles.particle[i].virt) {
				particles.particle[i].rho_prev = particles.particle[i].rho;
				particles.particle[i].e_prev = particles.particle[i].e;
				particles.particle[i].rho += dt*particles.particle[i].drhodt/2.0;
				particles.particle[i].e += dt*particles.particle[i].dedt/2.0;
				for(alfa = 0; alfa < DIM; alfa++) {
					particles.particle[i].v_prev[alfa] = particles.particle[i].v[alfa];
					particles.particle[i].v[alfa] += dt*particles.particle[i].a[alfa]/2.0;
					for(beta = 0; beta < DIM; beta++) {
						particles.particle[i].tau_prev[alfa][beta] = particles.particle[i].tau[alfa][beta];
						particles.particle[i].tau[alfa][beta] += dt*particles.particle[i].tau_dot[alfa][beta]/2.0;
					}
				}
			}
		}

		singleStep(&particles, &int_pairs);

		// LeapFrog Integration. Liu 2003 (pag 209)
		for(i = 0; i < particles.quant; i++) {
			if(!particles.particle[i].virt) {
				particles.particle[i].rho = particles.particle[i].rho_prev + dt*particles.particle[i].drhodt;
				particles.particle[i].e = particles.particle[i].e_prev + dt*particles.particle[i].dedt;
				for(alfa = 0; alfa < DIM; alfa++) {
					particles.particle[i].v[alfa] = particles.particle[i].v_prev[alfa] + dt*particles.particle[i].a[alfa] + particles.particle[i].av[alfa];
					particles.particle[i].r[alfa] += dt*particles.particle[i].v[alfa];
					for(beta = 0; beta < DIM; beta++) {
						particles.particle[i].tau[alfa][beta] = particles.particle[i].tau_prev[alfa][beta] + dt*particles.particle[i].tau_dot[alfa][beta];
					}
				}
			}
		}

		// The next line enable or disable the plasticity of the material
		plasticYieldModel(&particles, 5.0e8); // von Mieses - Armco Iron

		// Print the progress to Standard Error every 100 steps
		if(c%100 == 0) fprintf(stderr, "%0.1f %%\n", (t/0.0001)*100);

		// Print the data to Standard Output every 100 steps
		for(i = 0; i < particles.quant; i++) {
			if(!particles.particle[i].virt && c%100 == 0) { // Print if not a virtual particle
				printf("%g %g\n", particles.particle[i].r[0], particles.particle[i].r[1]);
			}
		}
		c++;
	}

	return 0;
}
