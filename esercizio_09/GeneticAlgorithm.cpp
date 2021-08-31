#include "GeneticAlgorithm.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace arma;
using namespace std;


int main() {
  Initialize();
  for (int igen = 1; igen <= n_generations; ++igen) {
    // generation 1 is already randomized
    Reset();    // Reset observables
    if (igen != 1){
      Select();
      Mutate();
    }
    Rank();     // Evaluate fitness
    Save(igen); // Save results to file
    if (report == 1 && igen % 5 == 0)
      Report(igen); // Print to terminal
  }
  cities.save("cities.reload"); // To enable restarting
  return 0;
}

/*
Seeds RNG, generates city coordinates, initializes data,
randomizes initial population,
*/
void Initialize() {
  ReadInput(); // read input.dat
  // Seed random number generator
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();

  // Generate city coordinates
  cities = zeros<mat>(n_cities, 2); // X, Y ; first = last
  if (reload == 1) {
    cities.load("cities.reload");
  } else {
    if (mode == 0)
      CitiesOnCircumference();
    else
      CitiesInSquare();
  }

  cities.save("cities.dat", raw_ascii);

  // Matrix of fitness
  best_path = zeros<urowvec>(n_cities - 1); // Imply city 0

  fitness = zeros<rowvec>(n_replicas); // Cost of each replica
  pdf = zeros<rowvec>(n_replicas); // e^(-b*Cost)
  cdf = zeros<rowvec>(n_replicas + 1); // CDF for rolling

  // Matrix of paths
  population = field<urowvec>(n_replicas);
  oldpopulation = field<urowvec>(n_replicas);

  // populate and then shuffle
  int ic = 0;
  urowvec path(n_cities - 1);
  path.imbue([&]() { // 1 ... n_city-1
    ++ic;
    return ic;
  });
  for (int irep = 0; irep < n_replicas; ++irep) {
    population(irep) = shuffle(path); // randomize elements
  }


  if (report == 1) Welcome(); // Welcome message to terminal

  // Clear files
  ofstream of;
  of.open("./best_path.dat");
  of.close();
  of.open("./least_cost.dat");
  of.close();
  of.open("./avg_cost.dat");
  of.close();
}

void ReadInput() { // Read input.dat
  ifstream input;
  input.open("./input.dat");
  input >> n_cities;
  input >> n_replicas;
  input >> n_generations;
  input >> p_mut;
  input >> mode;
  input >> report;
  input >> reload;
  input.close();
}

void Reset() { // Reset observables
  best_path.zeros();
  blk_avg = 0.0;
  mutations = 0;
}

// ============== Reproduction ===========================

/*
Keeps 4 fittest members, then selects two parents (with
probability proportional to exp(-b*Cost)) and then crosses
them to populate the next generation with Crossover()
*/
void Select() { // Select parents
  int a, b, splice, n_kids;
  double r;
  oldpopulation = population;
  n_kids = 0;
  // Keep the best parents
  population.at(0) = oldpopulation.at(selection(0));
  population.at(1) = oldpopulation.at(selection(1));
  population.at(2) = oldpopulation.at(selection(2));
  population.at(3) = oldpopulation.at(selection(3));
  n_kids += 4;
  // Fill new population with
  while (n_kids < n_replicas - 1) {
    // Find two parents
    r = rnd.Rannyu();
    for (int i = 0; i < n_replicas; ++i) {
      if (r < cdf(i)) {
        a = selection(i);
        break;
      }
    }
    r = rnd.Rannyu();
    for (int j = 0; j < n_replicas; ++j) {
      if (r < cdf(j)) {
        b = selection(j);
        break;
      }
    }
    // Make sure they're different
    while (a == b) {
      r = rnd.Rannyu();
      for (int j = 0; j < n_replicas; ++j) {
        if (r < cdf(j)) {
          b = selection(j);
          break;
        }
      }
    }
    // populate
    splice = (int)rnd.Rannyu(2, n_cities - 2);
    Crossover(a, b, n_kids, n_kids + 1, splice);
    n_kids += 2;
  }
}

/*
Given indices of mom and dad, and location at which to
splice their "DNA", creates son and daughter in new
population
*/
void Crossover(int m, int d, int s, int dg, int splice) {

  // son and daughter
  (population.at(s)).zeros();
  (population.at(dg)).zeros();
  // ! mom = (oldpopulation.at(m));
  // ! dad = (oldpopulation.at(d));

  // Each child takes a bit of a parent
  for (int i = 0; i < splice; ++i) {
    (population.at(s)).at(i) = (oldpopulation.at(m)).at(i);
    (population.at(dg)).at(i) = (oldpopulation.at(d)).at(i);
  }

  // And is afterwards completed (without repeats) and according to
  // crossover rules
  int j = 0;
  for (int i = 0; i < n_cities - 1; ++i) {
    if (!IsIn(population.at(s), (oldpopulation.at(d)).at(i))) {
      (population.at(s)).at(splice + j) = (oldpopulation.at(d)).at(i);
      ++j;
    }
  }
  j = 0;
  for (int i = 0; i < n_cities - 1; ++i) {
    if (!IsIn(population.at(dg), (oldpopulation.at(m)).at(i))) {
      (population.at(dg)).at(splice + j) = (oldpopulation.at(m)).at(i);
      ++j;
    }
  }
}

// ============== Fitness ==============================

/*
Get:
  fitness of each path -> fitness(irep)
  indices of best paths -> selection(irep)
  (e.g. best path is at selection(0))
  normalized exp(-b*Cost) of the paths -> pdf, b = 1/00
  cumulative distribution of the paths -> cdf
*/
void Rank() {
  double accumulator;

  for (int irep = 0; irep < n_replicas; ++irep) {
    fitness(irep) = Cost(irep);
  }
  // Get indices of top paths
  selection = sort_index(fitness, "ascend");
  best_path = population.at(selection(0));
  least_cost = Cost(selection(0));
  // Calculate probabilities weighted by fitness
  pdf = exp(-fitness / 100);
  pdf /= accu(pdf);
  pdf = sort(pdf, "descend");
  // Create a cdf to cycle through
  accumulator = 0;
  for (int i = 0; i < n_replicas; ++i) {
    accumulator += pdf(i);
    cdf(i + 1) = accumulator; // keep zero at i=0
  }
  // cdf.save("./cdf.dat", raw_ascii);
}

// Calculates cost of a replica using Norm()
double Cost(int ir) {
  double cost;
  int a, b;
  cost = 0;
  // first and last
  a = 0;
  b = (population.at(ir)).at(0);
  cost += Norm(a, b);
  a = (population.at(ir)).at(n_cities - 2);
  b = 0;
  cost += Norm(a, b);
  for (int ic = 0; ic < n_cities - 2; ++ic) {
    a = (population.at(ir)).at(ic);
    b = (population.at(ir)).at(ic + 1);
    cost += Norm(a, b);
  }
  return cost;
}

// L2 Norm
double Norm(int a, int b) {
  double sum;
  sum = (cities(a, 0) - cities(b, 0)) * (cities(a, 0) - cities(b, 0));
  sum += (cities(a, 1) - cities(b, 1)) * (cities(a, 1) - cities(b, 1));
  return abs(sum);
}

// ================= Mutation =============================

/*
Attempts to mutate each path, each mutation has p_mut of
happening.

! Never touches the home city

Mutations:
- Swap(replica, o, u) : Swap city o with city u
- Invert(rep, o, m) : Reverse m cities at o
- Shift(rep, o, m, n) : Shift n cities by m at o
- Permute(rep, o, u, m) : Swap m cities at o with m cities at u
*/
void Mutate() {
  double r;
  int o, u, m, n;
  // string msg; // See Check()
  for (int irep = 1; irep < n_replicas; ++irep) {
    r = rnd.Rannyu(); // U[0,1)
    if (r < 1 * p_mut) { // Swap two cities
      o = (int)rnd.Rannyu(0, n_cities - 1);
      u = (int)rnd.Rannyu(0, n_cities - 1);
      while (o == u) {
        u = (int)rnd.Rannyu(0, n_cities - 1);
      }
      Swap(irep, o, u);
      ++mutations;
    } else if (r < 2 * p_mut) { // Invert contiguous cities
      m = (int)rnd.Rannyu(3, n_cities - 1);
      o = (int)rnd.Rannyu(0, n_cities - 1 - m);
      Invert(irep, o, m);
      ++mutations;
    } else if (r < 3 * p_mut) { // Shuffle the path
      population.at(irep) = shuffle(population.at(irep));
      ++mutations;
    } else if (r < 4 * p_mut) {    // Shift a set of cities up
      n = (int)rnd.Rannyu(1, n_cities - 3);         // number of cities to shift
      m = (int)rnd.Rannyu(1, n_cities - 3 - n);     // size of shift
      o = (int)rnd.Rannyu(0, n_cities - 3 - n - m); // starting index
      Shift(irep, o, m, n);
      ++mutations;
    } else if (r < 5 * p_mut) {      // Swap two sets of cities
      m = (int)rnd.Rannyu(2, (n_cities - 1) / 2 - 1); // size of set
      o = (int)rnd.Rannyu(0, (n_cities - 1) / 2 - m); // index of first set
      u = (int)rnd.Rannyu((n_cities - 1) / 2, n_cities - 1 - m); // index of second set
      Permute(irep, o, u, m);
      ++mutations;
    }
  }
}

void AllRandom() {
  for (int ir = 0; ir < n_replicas; ++ir) {
    population.at(ir) = shuffle(population.at(ir));
  }
}

// Swap city a with city b
void Swap(int ir, int ic_a, int ic_b) {
  int u = (population.at(ir)).at(ic_a);
  (population.at(ir)).at(ic_a) = (population.at(ir)).at(ic_b);
  (population.at(ir)).at(ic_b) = u;
}

// Invert (reverse) m cities at ic
void Invert(int ir, int ic, int m) {
  urowvec copy(m);
  for (int i = 0; i < m; ++i) {
    copy(i) = (population.at(ir)).at(i + ic);
  }
  copy = reverse(copy);
  for (int i = 0; i < m; ++i) {
    (population.at(ir)).at(i + ic) = copy(i);
  }
}

// Shift n cities by m at ic
void Shift(int ir, int ic, int m, int n) {
  urowvec path = population.at(ir);
  urowvec take = zeros<urowvec>(n);
  // take n cities
  for (int i = 0; i < n; ++i) {
    take(i) = path(ic);
    path.shed_col(ic);
  }
  // place them at l+m
  path.insert_cols(ic + m, take);
  for (int i = 0; i < n_cities - 1; ++i) {
    (population.at(ir)).at(i) = path(i);
  }
}

// Swap m cities at ic_a with m cities at ic_b
void Permute(int ir, int ic_a, int ic_b, int m) {
  urowvec path = population.at(ir);
  urowvec take_a = zeros<urowvec>(m);
  urowvec take_b = zeros<urowvec>(m);
  // take the first set of m cities at ic_a
  for (int i = 0; i < m; ++i) {
    take_a(i) = path(ic_a);
    path.shed_col(ic_a);
  }
  // take the second set at ic_b
  for (int i = 0; i < m; ++i) {
    take_b(i) = path(ic_b - m);
    path.shed_col(ic_b - m);
  }
  // Switch the two
  path.insert_cols(ic_a, take_b);
  path.insert_cols(ic_b, take_a);
  for (int i = 0; i < n_cities - 1; ++i) {
    (population.at(ir)).at(i) = path(i);
  }
}

// ============== Input and Output =======================

void Save(int ig) {
  ofstream of;
  // print current best path to file
  of.open("./best_path.dat", ios::app);
  best_path.save(of, raw_ascii);
  of.close();
  of.open("./least_cost.dat", ios::app);
  of << setprecision(8) << least_cost << endl;
  of.close();
  // Average of best paths
  fitness = sort(fitness);
  for (int ir = 0; ir < n_replicas / 2; ++ir) {
    blk_avg += fitness(ir);
  }
  blk_avg /= (long double)n_replicas / 2;
  of.open("./avg_cost.dat", ios::app);
  of << ig << " " << blk_avg << endl;
  of.close();
}

void Report(int ig) {
  cout << "Generation : " << ig << endl;
  // best_path.print("Best solution :");
  cout << "Avg cost : " << blk_avg << endl;
  cout << "Least cost : " << setprecision(8) << least_cost << endl;
  cout << "Mutations: " << mutations << endl;
  cout << "=========================" << endl;
}

void Welcome() {
  cout << endl;
  cout << "Traveling Salesman with Genetic Algorithm" << endl;
  cout << n_cities << " cities placed ";
  if (mode == 0)
    cout << "randomly on a circumference.";
  else
    cout << "randomly inside a square";
  cout << endl;
  cout << n_replicas << " replicas" << endl;
  cout << n_generations << " generations" << endl;
  cout << "Chance of mutations is " << 5 * p_mut * 100 << " percent" << endl;
  cout << endl << endl;
}

// ============== Other Functions =======================

void CitiesOnCircumference() {
  double x, y, theta, R;
  R = 10;
  for (int i = 0; i < n_cities; ++i) {
    theta = rnd.Rannyu(0, 2 * pi);
    x = R * cos(theta);
    y = R * sin(theta);
    cities(i, 0) = x;
    cities(i, 1) = y;
  }
}

void CitiesInSquare() {
  double a = 10;
  for (int i = 0; i < n_cities; ++i) {
    cities(i, 0) = rnd.Rannyu(-a, a);
    cities(i, 1) = rnd.Rannyu(-a, a);
  }
}

bool IsIn(const urowvec &vec, int city) {
  // Returns true
  for (int ic = 0; ic < n_cities - 1; ++ic) {
    if (city == vec(ic)) {
      return true;
    }
  }
  return false;
}

void Check(int ir, string msg) {
  ofstream check;
  check.open("check.out", ios::app);
  check << ir << endl;
  check << msg << endl << endl;
  check << population.at(ir) << endl;
}

void PrintPopulation(int igen) { // print the current sorted population
  int n_best = 5;
  ofstream pop;
  pop.open("./evolution.dat", ios::app);
  pop << "Generation : " << igen << endl;
  for (int ir = 0; ir < n_best; ++ir) {
    (population.at(selection(ir))).save(pop, raw_ascii);
  }
  pop << endl;
}
