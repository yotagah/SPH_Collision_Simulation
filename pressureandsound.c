#include <math.h>
#include "datatypes.h"
#include "pressureandsound.h"

// Calculates the pressure and speed of sound for the particles (Equation 8.9 and 8.17 Liu, Mie-Gruneisen equation for solids)
void solidsPressureAndSoundSpeed(particles_t *parts, double gamma, double sound_speed, double slope,  double mi) {

	// gamma: Grunaisen parameter
	// mi(G): shear module

	int k;

	double a_0, b_0, c_0; // Hugoniot curve constants
	double p_H; // pressure curve Hugoniot
	double eta; // rho/rho_0 Eq. 8.10 Liu 2003
	//double dpdrho; // Eq. 8.17 Liu 2003

	for(k = 0; k < parts->quant; k++) { // All particles
		a_0 = parts->particle[k].rho_0 * sound_speed * sound_speed; // Eq. 8.13 Liu 2003
		b_0 = a_0 * (1.0 + 2.0*(slope - 1.0)); // Eq. 8.14 Liu 2003
		c_0 = a_0 * (2.0*(slope - 1.0) + 3.0*(slope - 1.0)*(slope - 1.0)); // Eq. 8.15 Liu 2003

		eta = (parts->particle[k].rho / parts->particle[k].rho_0) - 1.0; // Eq. 8.10 Liu 2003

		if(eta > 0.0) { // Eq. 8.11 Liu 2003
			p_H = a_0*eta + b_0*eta*eta + c_0*eta*eta*eta;
			//dpdrho = gamma*parts->particle[k].e + (a_0 + eta*(2.0*b_0 - a_0*gamma) + eta*eta*(3.0*c_0 - 3.0*b_0*gamma/2.0) + 2*eta*eta*eta*c_0*gamma); // dp/drho (derivative of Eq. 8.9 Liu 2003)
		} else {
			p_H = a_0*eta;
			//dpdrho = gamma*parts->particle[k].e + (a_0 - eta*a_0*gamma); // dp/drho (derivative of Eq. 8.9 Liu 2003)
		}

		parts->particle[k].p = (1 - gamma*eta/2.0)*p_H + gamma*parts->particle[k].rho*parts->particle[k].e; // Eq. 8.9 Liu 2003

		parts->particle[k].c = sqrt(4*mi/(3*parts->particle[k].rho_0)); // Eq. 8.17 Liu 2003 (removed: + dpdrho)
	}

}