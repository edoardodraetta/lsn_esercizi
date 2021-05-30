#ifndef __GENALGO__
#define __GENALGO__

#include "random.h"
#include <armadillo>

// Top level
void Initialize(std::string filename = "input.dat");
void Reset();
void Rank();
void Save();
void Select();
// void Crossover();
void Mutate();
void Report();

// Mutation
void Swap(int, int);
void Shift(int, int, int, int);
void Permute(int, int, int);
void Invert(int, int, int);

// Other
void Welcome();

void ReadInput(std::string);
void Print(std::string);
double Cost(int);
int Pbc(int);
void CitiesOnCircumference();


// Parameters
int seed[4];
Random rnd;

int n_cities, n_replicas, n_dim, n_generations;

// Data
arma::umat population;
arma::mat cities;
arma::vec fitness;
arma::uvec selection;

// Output
arma::vec best_path;
arma::vec avg_distance;

// Constants
const double pi=3.1415927;

#endif
