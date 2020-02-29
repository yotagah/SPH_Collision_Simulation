#include <math.h>
#include "vector.h"
#include "datatypes.h"

// subtract two vectors
void subVector(double a[DIM], double b[DIM], double result[DIM]) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    //result[2] = a[2] - b[2];
}

// Calculates the modulus of a vector
double modVector(double v[DIM]) {
    return sqrt(v[0]*v[0] + v[1]*v[1]);
    //return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}