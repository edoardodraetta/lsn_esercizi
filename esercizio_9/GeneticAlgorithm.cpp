#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include "GeneticAlgorithm.h"

using namespace arma;
using namespace std;

int main(){

  Initialize();
  for (int igen = 1; igen <= n_generations; ++igen) {
    Reset(); // Resets observables
    Rank(); // Evaluate fitness
    Save(); // Measure and print observables
    Select(); // Select and cross parents
    Mutate(); // Mutate population
    Report(); // Print to terminal
  }
  return 0;
}

// Top level

void Initialize(string filename){
  ReadInput(filename);

	// Seed random number generator
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  // City coordinates
  cities = zeros<mat>(n_cities, n_dim);
  CitiesOnCircumference();

  // Matrix of fitness
  fitness = zeros<vec>(n_replicas);

  // Randomize starting population
  population = zeros<umat>(n_replicas, n_cities);

  for (int i = 0; i < n_replicas; ++i){
    for (int j = 0; j < n_cities; ++j){
      population(i, j) = j;
    }
  }
  int n_init = 1000;
  for (int i = 0; i < n_init; ++i){
    for (int irep = 0; irep < n_replicas; ++irep){
      // Choose a random city (not the first!)
      int rand_city = (int) rnd.Rannyu(1, n_cities-1);
      Swap(irep, rand_city);
    }
  }

  population.brief_print("Starting population : ");
  cities.brief_print("City coordinates : ");
  population.save("population0.dat", raw_ascii);
  cities.save("cities.dat", raw_ascii);
}

void Reset(){}

void Rank(){
  for (int irep = 0; irep < n_replicas; ++irep){
    fitness(irep) = Cost(irep);
  }
  selection = sort_index(fitness, "ascend");

  fitness.save("fitness.dat", raw_ascii);
  selection.save("selection.dat", raw_ascii);
}

void Select(){
  umat old = population;
  for (int irep = 0; irep < n_replicas; ++irep){
    for (int icity = 0; icity < n_cities; ++icity){
      population(irep, icity) = old(selection(irep), icity);
    }
  }
}

void Mutate(){
  int r = (int)rnd.Rannyu(0, 100);

  if (r <= 5) {
    int n = (int)rnd.Rannyu(1, n_cities-2);
    int m = (int)rnd.Rannyu(1, n_cities-2-n);
    int l = (int)rnd.Rannyu(1, n_cities-2-n-m);
    cout << population.row(irep) << endl;
    cout << l << ", " << m << ", " << n << endl;
    Shift(irep, l, m, n);
    cout << population.row(irep) << endl;
  }
}

void Report(){}

// Mutation

void Swap(int i, int j){ // Swaps two cities

  // Swap it with the following one
  int u = population(i, j);
  population(i, j) = population(i, j+1);
  population(i, j+1) = u;
}

void Shift(int irep, int l, int m, int n){ // shift n cities by m

  umat path = population.row(irep);
  umat take = zeros<umat>(m);

  cout << path << endl << endl;

  for (int i = 0; i < m; ++i){
    take(i) = path(l);
    path.shed_col(l);
  }

  if (l+n>n_cities-l) path.insert_cols(1, take.t() );
  else path.insert_cols(l+n, take.t() );

  for (int i = 0; i < n_cities; ++i){
    population(irep, i) = path(i);
  }

  cout << path << endl << endl;
}

void Permute(int irep, int l, int m){
  // m < N/2
  umat path = population.row(irep);
  umat take_a = zeros<umat>(m);
  umat take_b = zeros<umat>(m);

  cout << path << endl << endl;

  for (int i = 0; i < m; ++i){
    take_a(i) = path(l);
    path.shed_col(l);
  }
  for (int i = 0; i < m; ++i){
    take_b(i) = path(l);
    path.shed_col(l);
  }

  path.insert_cols(l, take_b.t() );
  path.insert_cols(l+m, take_a.t() );

  cout << path << endl << endl;
  for (int i = 0; i < n_cities; ++i){
    population(irep, i) = path(i);
  }
}

void Invert(int irep, int l, int m){
  umat path = population.row(irep);
  umat take = zeros<umat>(m);

  cout << path << endl << endl;

  for (int i = 0; i < m; ++i){
    take(i) = path(l);
    path.shed_col(l);
  }

  cout << path << endl << endl;
  path.insert_cols(l, fliplr(take.t()));
  cout << path << endl << endl;

  for (int i = 0; i < n_cities; ++i){
    population(irep, i) = path(i);
  }
}

// Messages and file saving

void ReadInput(std::string filename){
  ifstream input;
  input.open(filename);
  input >> n_cities;
  input >> n_replicas;
  input >> n_dim;
  input.close();
}

void Welcome(){}

void Print(string filename){

  population.save(filename, raw_ascii);
}

// Other

double Cost(int i){ // Calculate the cost of a path
  double cost = 0;
  int a, b;
  for (int icity = 0; icity < n_cities; ++icity){
    a = population(i, icity);
    b = population(i, Pbc(icity+1));
    cost += abs(cities(a,0)+cities(a,1)-cities(b,0)-cities(b,1));
  }
  return cost;
}

int Pbc(int icity){ // periodic boundary conditions
  if (icity==n_cities) icity = 0;
  return icity;
}

void CitiesOnCircumference(){
  double x, y, theta, R;
  R = 10;
  for (int i = 0; i < n_cities; ++i){
    theta = rnd.Rannyu(0, 2*pi);
    x = R * cos(theta);
    y = R * sin(theta);

    cities(i, 0) = x;
    cities(i, 1) = y;
  }
}
