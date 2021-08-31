#ifndef __PARSEARCH__
#define __PARSEARCH__

#include "random.h"
#include <armadillo>

// RNG
int seed[4];
Random rnd;

// Parallel Search
int n_migr;
int mpi_size, mpi_rank;
arma::urowvec glob_best_path;
double glob_least_cost;
void Migrate();
void Exchange();

// Output files
std::string bestpathfile, leastcostfile, avgcostfile;

// Params
int n_cities, n_replicas, n_generations;
bool mode, report, reload;

// Data
arma::mat cities;
arma::field<arma::urowvec> population;
arma::field<arma::urowvec> oldpopulation;

// Init
void Initialize();
void InitRNG();
void ReadInput();
void AllRandom();
void Reset();

// Messages and Output
arma::urowvec best_path;
double blk_avg;
double least_cost;
void Welcome();
void Report(int);
void Save(int);

arma::rowvec pdf;
arma::rowvec cdf;

arma::rowvec fitness;
arma::uvec selection;
void RankFitness();
double Cost(int);
double Norm(int, int);

// Genetics
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


// Other
void MakeCities(int);
bool IsIn(const arma::urowvec &, int);

// Constants
const double pi=3.1415927;

#endif
