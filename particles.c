#include <stdlib.h>
#include "datatypes.h"
#include "particles.h"

// Presets for particles
particles_t newParticles() {
	particles_t new;

	new.quant = 0;
	new.alocated = ALOC_STEP;
	new.particle = (particle_t*) malloc(sizeof(particle_t)*ALOC_STEP);

	return new;
}

// Add new particle
void addParticle(particles_t *parts, double x, double y, double vx, double vy, double m, double rho_0, double h_size) {
//void addParticle(particles_t *parts, double x, double y, double z, double vx, double vy, double vz, double m, double rho_0, double h_size) {

	int alfa, beta;

	particle_t new;

	new.r[0] = x;
	new.r[1] = y;
	//new.r[2] = z;

	new.v[0] = vx;
	new.v[1] = vy;
	//new.v[2] = vz;

	new.m = m;
	new.rho = rho_0;
	new.rho_0 = rho_0;
	new.e = 0;
	new.dedt = 0;

	new.h = h_size;

	new.virt = 0;

	for(alfa = 0; alfa < DIM; alfa++) {
		for(beta = 0; beta < DIM; beta++) {
			new.tau[alfa][beta] = 0;
		}
		new.av[alfa] = 0.0;
	}

	parts->quant += 1;

    if(parts->quant > parts->alocated) { // If there is no memory allocated for that amount
    	alocMoreParticles(parts); // Alocate more
    }

	parts->particle[parts->quant-1] = new;
}

// Allocates more memory for particles
void alocMoreParticles(particles_t *parts) {
	parts->alocated += ALOC_STEP;
	parts->particle = (particle_t*) realloc(parts->particle, sizeof(particle_t)*parts->alocated);
}

// Presets for interaction pairs
int_pairs_t newIntPairs() {
	int_pairs_t new;

	new.quant = 0;
	new.alocated = ALOC_STEP;
	new.int_pair = (int_pair_t*) malloc(sizeof(int_pair_t)*ALOC_STEP);

	return new;
}

// Allocates more memory for interaction pairs
void alocMorePairs(int_pairs_t *pairs) {
	pairs->alocated += ALOC_STEP;
	pairs->int_pair = (int_pair_t*) realloc(pairs->int_pair, sizeof(int_pair_t)*pairs->alocated);
}