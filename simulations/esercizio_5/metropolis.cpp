#include "metropolis.h"

double groundstate(double r[3]){ // phi 100
	double rho = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	return exp(-2*rho) / M_PI; // p
}

double excitedstate(double r[3]){ // phi 210
	double rho = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	double theta = atan (r[1] / r[0]);
	return rho * exp(-rho) * cos(theta) / (32*M_PI);
}

