#include "datatypes.h"
#include "integration.h"
#include "vector.h"
#include "particles.h"
#include "kernel.h"
#include "density.h"
#include "stress.h"
#include "pressureandsound.h"
#include "force.h"
#include "viscosity.h"
#include "temperature.h"

#include <stdio.h>

// Calculate the values for a time step
void singleStep(particles_t *particles, int_pairs_t *int_pairs) {

    // Interaction pairs
    directFind(particles, int_pairs);

    // Equations
    conDensity(particles, int_pairs);
    solidsPressureAndSoundSpeed(particles, 1.81, 3630.0, 1.8, 8.0e10); // gamma, sound_speed, slope, shear_modulus - Armco Iron
    totalStressTensor(particles, int_pairs, 8.0e10); // shear modulus - Armco Iron

    // Forces
    zerateForces(particles);
    artificialViscosity(particles, int_pairs);
    //artificialHeat(particles, int_pairs);
    internalForce(particles, int_pairs);
    //externalForce(particles);

    hUpgrade(particles);
    //avarageVelocity(particles, int_pairs);
}

// Find the particles that will make the interaction (that are within the SPH radius) and calculate the weight function for each pair
void directFind(particles_t *particles, int_pairs_t *int_pairs) {

    int i, j;
    int scale_k; // Scale factor for the SPH radius
    double dist[DIM]; // Distance vector
    double m_dist; // Distance module
    double mh; // Average particle interaction distance
    
    scale_k = 2; // Scale factor for: cubic spline kernel by W4 - Spline (Monaghan 1985)

    int_pairs->quant = 0; // Zerate number of interaction pairs

    for(i = 0; i < particles->quant; i++) { // All particles
        for(j = i+1; j < particles->quant; j++) { // All pairs
            subVector(particles->particle[i].r, particles->particle[j].r, dist); // Find distance vector
            m_dist = modVector(dist); // Find distance module

            mh = (particles->particle[i].h + particles->particle[j].h) / 2.0;

            if(m_dist < scale_k*mh) { // IF is within the SPH radius
                if(int_pairs->quant+1 > int_pairs->alocated) // If there is no memory alocated for the pairs
                    alocMorePairs(int_pairs); // Alocate more memory

                int_pairs->int_pair[int_pairs->quant].i = i; // First particle of the interaction pair
                int_pairs->int_pair[int_pairs->quant].j = j; // Second particle of the interaction pair

                kernel(m_dist, dist, &(int_pairs->int_pair[int_pairs->quant]), mh); // Calculates weight functoin and derivative

                int_pairs->quant += 1; // Increase number of pairs
            }
        }
    }
}

// Calculates average speed to fix speed and prevent penetration (Monaghan, 1992)
void avarageVelocity(particles_t *parts, int_pairs_t *pairs) {

    int i, j, k, alfa;
    double aux, dv, epsilon;

    // epsilon -- a small constants chosen by experence, may lead to instability.
    // for example, for the 1 dimensional shock tube problem, the E <= 0.3
    epsilon = 0.1;

    for(k = 0; k < parts->quant; k++) {
        for(alfa = 0; alfa < DIM; alfa++) {
            parts->particle[k].av[alfa] = 0.0;
        }
    }

    for(k = 0; k < pairs->quant; k++) {
        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        aux = pairs->int_pair[k].w / (parts->particle[i].rho + parts->particle[j].rho);

        for(alfa = 0; alfa < DIM; alfa++) {
            dv = parts->particle[i].v[alfa] - parts->particle[j].v[alfa];
            parts->particle[i].av[alfa] -= parts->particle[j].m * dv * aux;
            parts->particle[j].av[alfa] += parts->particle[i].m * dv * aux;
        }
    }

    epsilon *= 2.0;
    for(k = 0; k < parts->quant; k++) {
        for(alfa = 0; alfa < DIM; alfa++) {
            parts->particle[k].av[alfa] = epsilon * parts->particle[k].av[alfa];
        }
    }
}

// Update the h size of smooth according to the article of Liu 2005 Eq. 11
void hUpgrade(particles_t *parts) {

    int k;
    double aux;

    for(k = 0; k < parts->quant; k++) {
        aux = -(parts->particle[k].h / (DIM * parts->particle[k].rho)) * parts->particle[k].drhodt;
        parts->particle[k].h += aux*dt;
        if(parts->particle[k].h <= 0.0)
            parts->particle[k].h -= aux*dt;
    }
}