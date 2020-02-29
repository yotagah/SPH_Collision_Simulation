#ifndef __PRESSURE__
#define __PRESSURE__

// Calculates the pressure and speed of sound for the particles (Equation 8.9 and 8.17 Liu, Mie-Gruneisen equation for solids)
void solidsPressureAndSoundSpeed(particles_t *parts, double gamma, double sound_speed, double slope,  double mi);

#endif