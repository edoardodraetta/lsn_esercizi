#include "anneal.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace arma;
using namespace std;

int main() {
  Initialize();
  for (int iblock = 1; iblock <= n_blocks; ++iblock) {
    // Reset();
    for (int istep = 1; istep <= n_steps; ++istep) {
      Move();
    }
    Report(iblock);  // to terminal
    Record(iblock);  // to file
    Anneal(iblock);  //
  }
  return 0;
}

void Initialize() {
  // Read input file
  ReadInput();

  // Generate city coordinates
  cities = zeros<mat>(n_cities + 1, 2);
  if (reload == 1) {
    cities.load("cities.reload");
  } else {
    if (mode == 0)
      CitiesOnCircumference();
    else
      CitiesInSquare();
  }

  cities.save("./cities.dat", raw_ascii);

  // Random starting path
  currentpath = zeros<urowvec>(n_cities - 1);
  int ic = 0;
  currentpath.imbue([&]() {  // 1 ... n_city-1
    ++ic;
    return ic;
  });
  currentpath = shuffle(currentpath);
  bestpath = currentpath;

  // output files reset
  ofstream reset;
  reset.open("./best_paths.dat");
  reset.close();

  improved = 0;
  // welcome message
  Welcome();
}

void Move() {  // Move with the Metropolis algorithm
  double alpha;

  // 1. Generate a new path with a symmetric transition

  newpath = Mutate(currentpath, rnd.Rannyu());

  // 2. Accept or reject move

  attempted++;
  if (Cost(newpath) < Cost(currentpath)) {
    currentpath = newpath;
    accepted++;
  } else {
    alpha = Weight(newpath) / Weight(currentpath);
    if (alpha > rnd.Rannyu()) {
      currentpath = newpath;
      accepted++;
    }
  }

  // 3. Keep the best path
  if (Cost(newpath) < Cost(bestpath)) {
    bestpath = newpath;
    improved++;
  }
}

void Report(int iblk) {  // Print to terminal
  cout << "Block : " << iblk << endl;
  if (n_cities < 12)
    bestpath.raw_print("Current Best Path : ");
  else
    bestpath.brief_print("Current Best Path : ");
  cout << "Fitness : " << Cost(bestpath) << endl;
  cout << "Acceptance rate : ";
  cout << (double)accepted / (double)attempted << endl;
  cout << "temp : " << temp << endl;
  cout << "Times improved : " << improved << endl;
  cout << "======================================" << endl;
}

void Record(int iblk) {  // Save to file
  ofstream file;
  file.open("./best_paths.dat", ios::app);
  bestpath.save(file, raw_ascii);
  file.close();
  if (iblk == n_blocks) cities.save("./cities.reload");
}

void Anneal(int iblk) { temp *= exp(-rate * iblk); }

// ====== Mutations =====================

urowvec Mutate(urowvec& path, double r) {
  int o, u, m, n;
  if (r < 0.33) {  // Swap two random cities
    o = (int)rnd.Rannyu(0, n_cities - 1);
    u = (int)rnd.Rannyu(0, n_cities - 1);
    while (o == u) {
      u = (int)rnd.Rannyu(0, n_cities - 1);
    }
    return Swap(path, o, u);
  } else if (r < .66) {  // Invert m cities at l
    m = (int)rnd.Rannyu(3, n_cities - 1);
    o = (int)rnd.Rannyu(0, n_cities - 1 - m);
    return Invert(path, o, m);
  } else if (r < 1) {  // Shuffle
    return shuffle(path);
  }
}

urowvec Swap(urowvec& path, int a, int b) {  // Swap two cities in the path
  int u = path(a);
  path(a) = path(b);
  path(b) = u;
  return path;
}

urowvec Invert(urowvec& path, int l, int m) {
  urowvec copy(m);
  for (int i = 0; i < m; ++i) {  // take m from l
    copy(i) = path(i + l);
  }
  copy = reverse(copy);
  for (int i = 0; i < m; ++i) {  // replace in reverse order
    path.at(i + l) = copy(i);
  }
  return path;
}

urowvec Permute(urowvec& path, int l, int k, int m) {
  urowvec take_a = zeros<urowvec>(m);
  urowvec take_b = zeros<urowvec>(m);

  // take the first set of m cities at l
  for (int i = 0; i < m; ++i) {
    take_a(i) = path(l);
    path.shed_col(l);
  }
  // take the second set at k
  for (int i = 0; i < m; ++i) {
    take_b(i) = path(k - m);
    path.shed_col(k - m);
  }
  // Switch the two
  path.insert_cols(l, take_b);
  path.insert_cols(k, take_a);

  return path;
}

// ======= Norm, Cost, Weight ===============

double Weight(urowvec& path) {
  // Returns boltzmann weight of a path
  double cost = Cost(path);
  return exp(-cost / temp);
}

double Cost(urowvec& path) {
  // Given path, calculates its current cost
  double cost;
  int a, b;
  cost = 0;

  // first and last
  a = 0;
  b = path(0);
  cost += Norm(a, b);

  a = path(n_cities - 2);
  b = 0;
  cost += Norm(a, b);

  // the rest
  for (int ic = 0; ic < n_cities - 2; ++ic) {
    a = path(ic);
    b = path(ic + 1);
    cost += Norm(a, b);
  }
  return cost;
}

double Norm(int a, int b) {
  // Given index of two cities, calculates distance
  double sum;
  sum = (cities(a, 0) - cities(b, 0)) * (cities(a, 0) - cities(b, 0));
  sum += (cities(a, 1) - cities(b, 1)) * (cities(a, 1) - cities(b, 1));
  return abs(sum);
}

// ======= I/O ============

void ReadInput() {  // Read from input.dat
  ifstream input;
  input.open("./input.dat");
  input >> n_cities;
  input >> n_steps;
  input >> n_blocks;
  input >> temp;
  input >> rate;
  input >> mode;
  input >> report;
  input >> reload;
  input.close();

  // Seed random number generator
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();
  input.open("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();
}

void Welcome() {  // Welcome message
  cout << endl;
  cout << "Traveling Salesman with Metropolis Algorithm" << endl;
  cout << n_cities << " cities placed ";
  if (mode == 0)
    cout << "randomly on a circumference.";
  else
    cout << "randomly inside a square";
  cout << endl;
  cout << n_steps << " steps" << endl;
  cout << n_blocks << " blocks" << endl;
  cout << "Temperature : " << temp << endl;
  cout << "Cooling rate : " << rate << endl;
  cout << endl << endl;
}

// ======= Cities ===========

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
  cities(n_cities, 0) = cities(0, 0);
  cities(n_cities, 1) = cities(0, 1);
}

void CitiesInSquare() {
  double a = 10;
  for (int i = 0; i < n_cities; ++i) {
    cities(i, 0) = rnd.Rannyu(-a, a);
    cities(i, 1) = rnd.Rannyu(-a, a);
  }
  cities(n_cities, 0) = cities(0, 0);
  cities(n_cities, 1) = cities(0, 1);
}
