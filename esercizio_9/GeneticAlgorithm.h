#ifndef __GENALGO__
#define __GENALGO__

#include "random.h"
#include <armadillo>


// RNG
int seed[4];
Random rnd;

// Params
int n_cities, n_replicas, n_generations;
bool mode, report, reload;

// Data
arma::mat cities;
arma::field<arma::urowvec> population;
arma::field<arma::urowvec> oldpopulation;

// Init
void Initialize();
void ReadInput();
void AllRandom();
void Reset();

// Fitness
arma::rowvec fitness;
arma::uvec selection;
void Rank();
double Cost(int);
double Norm(int, int);

// Reproduction
arma::rowvec pdf;
arma::rowvec cdf;
void Select();
void Crossover(int, int, int, int, int);

// Mutation
double p_mut; // unit of mutation chance
int mutations;
void Mutate();
void Swap(int, int, int);
void Shift(int, int, int, int);
void Permute(int, int, int, int);
void Invert(int, int, int);
void SwapAB(int, int, int);

// Messages and Output
arma::urowvec best_path;
double blk_avg;
double least_cost;
void Welcome();
void Report(int);
void Save(int);
void Check(int, std::string);
void PrintPopulation(int);

// Other
void CitiesOnCircumference();
void CitiesInSquare();
bool IsIn(const arma::urowvec &, int);

// Constants
const double pi=3.1415927;

#endif
