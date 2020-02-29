#ifndef __DENSITY__
#define __DENSITY__

// Calculates the normalized density for each particle (Equation 4.35 Liu 2003 - Randies and Libersky, 1996; Chen et al., 1999a; 2000)
void normSumDensity(particles_t *parts, int_pairs_t *pairs);

// Calculates density by simple sum (Equation 4.26 Liu 2003)
void sumDensity(particles_t *parts, int_pairs_t *pairs);

// Calculates the rate of variation of density by continuity of the density (Equation 4.31 Liu 2003)
void conDensity(particles_t *parts, int_pairs_t *pairs);

// Copies the current density to the initial density
void copy2InitDensity(particles_t *parts);

#endif