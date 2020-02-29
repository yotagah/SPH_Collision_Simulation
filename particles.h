#ifndef __PARTICLES__
#define __PARTICLES__

#define ALOC_STEP 100 // Allocates memory to each quantity of particles / pairs

// Presets for particles
particles_t newParticles();

// Add new particle
void addParticle(particles_t *parts, double x, double y, double vx, double vy, double m, double rho_0, double h_size);

// Allocates more memory for particles
void alocMoreParticles(particles_t *parts);

// Presets for interaction pairs
int_pairs_t newIntPairs();

// Allocates more memory for interaction pairs
void alocMorePairs(int_pairs_t *pairs);

#endif
