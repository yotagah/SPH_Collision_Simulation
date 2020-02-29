#ifndef __STRESS__
#define __STRESS__

// Calculates the total tension tensor (sigma)
void totalStressTensor(particles_t *parts, int_pairs_t *pairs, double mi);

// Check that the traction does not exceed the limit and scale it to avoid it (von Mieses, Liu 2003 eq. 8.8)
void plasticYieldModel(particles_t *parts, double J_0);

#endif