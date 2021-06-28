#ifndef __SIMANN__
#define __SIMANN__

#include "random.h"
#include <armadillo>

// RNG
int seed[4];
Random rnd;

// Constants
const double pi=3.1415927;

// Simulation parameters
int n_cities, n_steps, n_blocks;
bool mode, report, reload;
int attempted, accepted;
int improved;
// Simulation Data
arma::mat cities;
arma::urowvec currentpath;
arma::urowvec newpath;
arma::urowvec bestpath;
double temp, rate;

// Simulation
void Initialize();
void Move();
void Report(int);
void Record(int);
void Anneal(int);

// Mutations
arma::urowvec Mutate(arma::urowvec &, double);
arma::urowvec Swap(arma::urowvec &, int, int);
arma::urowvec Invert(arma::urowvec &, int, int);

// Metropolis
double Weight(arma::urowvec &);
double Cost(arma::urowvec &);
double Norm(int, int);

// Cities
void CitiesOnCircumference();
void CitiesInSquare();

// Other
void Welcome();
void ReadInput();

#endif
