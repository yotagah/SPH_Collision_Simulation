#ifndef __DATATYPES__
#define __DATATYPES__

#define DIM 2 // Number of dimensions of the problem

extern double norm; // Weight function normalization
extern double dt; // Integration interval

// Particle
typedef struct particle_t {
	double a[DIM];				// Acceleration		(m/s^2)
	double v[DIM];				// Speed			(m/s)
	double r[DIM];				// Position			(m)

	double v_prev[DIM];			// Previous speed (for LeapFrog integration)
	double av[DIM];				// Average speed (for correction)

	double m; 					// Mass				(kg)
	double p;					// Pressure			(N/m^3)
	double rho;					// Density			(kg/m^3)
	double rho_0;				// Initial density	(kg/m^3)
	double drhodt;				// Density var.		(kg/s*m^3)
	double e;					// Energy			(J)
	double dedt;				// Energy var.		(J/s)
	double c;					// Sound speed 		(m/s)

	double rho_prev;				// Previous density
	double e_prev;				// Previous energy

	double h;					// Smoothing size

	double epsilon[DIM][DIM];   // Strain rate tensioner
	double tau[DIM][DIM];		// Shear stress tensioner
	double tau_dot[DIM][DIM];	// Shear stress rate tensor
	double sigma[DIM][DIM];		// Total tension tensioner

	double tau_prev[DIM][DIM];	// Previous shear stress tensioner
	double R[DIM][DIM];

	short virt;					// Virtual particle flag
} particle_t;

// Particles
typedef struct particles_t {
	int quant; 				// Particles quantity
	particle_t *particle; 	// Particles pointer
	int alocated;			// Alocated space quantity
} particles_t;

// Interaction pair
typedef struct int_pair_t {
	double w;				// Weight function (0.0 to 1.0)
	double dwdx[2];			// Derivative weight function (-1.0 to 1.0)
	int i;					// First particle of the pair
	int j;					// Second particle of the pair
} int_pair_t;

// Interation pairs
typedef struct int_pairs_t {
	int quant; 				// Pairs quantity
	int_pair_t *int_pair;	// Pairs pointer
	int alocated;			// Alocated space quantity
} int_pairs_t;

#endif