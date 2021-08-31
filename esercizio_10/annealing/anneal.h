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

// Simulation Data
arma::mat cities;
arma::urowvec currentpath, oldpath, bestpath;
double starting_temp, temp, rate;

// Simulation
void Initialize();
void Move();
void Report(int);
void Record(int);
void Anneal(int);

// Mutations
void Mutate(arma::urowvec &, double);
void Swap(arma::urowvec &, int, int);
void Invert(arma::urowvec &, int, int);

// Metropolis
double Cost(arma::urowvec &);
double Norm(int, int);

// Cities
void CitiesOnCircumference();
void CitiesInSquare();

// Other
void Welcome();
void ReadInput();

#endif
