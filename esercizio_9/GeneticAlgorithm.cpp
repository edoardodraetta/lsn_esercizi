#include <string>
#include <cmath>
#include <iostream>
#include "GeneticAlgorithm.h"
#include <ctime>

using namespace arma;
using namespace std;

int main(){
  clock_t start, stop;
  double erank = 0;
  double esave = 0;
  double eselect = 0;
  double emutate = 0;
  double erep = 0;

  Initialize();
  AllRandom(); // Randomize population

  for (int igen = 1; igen <= n_generations; ++igen){

    Reset(); // Reset observables

    start = clock();
    Rank(); // Evaluate fitness
    stop = clock();
    erank += (double) (stop - start) / CLOCKS_PER_SEC;

    start = clock();
    Save(igen); // Save results to file
    stop = clock();
    esave += (double) (stop - start) / CLOCKS_PER_SEC;

     // Create new generation
    if (igen != n_generations){
      start = clock();
      Select();
      stop = clock();
      eselect += (double) (stop - start) / CLOCKS_PER_SEC;
      start = clock();
      Mutate();
      stop = clock();
      emutate += (double) (stop - start) / CLOCKS_PER_SEC;
    }

    start = clock();
    if (report==1) Report(igen); // Print to terminal
    stop = clock();
    erep = (double) (stop - start) / CLOCKS_PER_SEC;
  }

  printf("\nTime to rank: %.5f\n", erank);
  printf("\nTime to save: %.5f\n", esave);
  printf("\nTime to select: %.5f\n", eselect);
  printf("\nTime to mutate: %.5f\n", emutate);
  printf("\nTime to report: %.5f\n", erep);
  printf("\nTime to Search: %.5f\n", timer);
  printf("\nTime to cross: %.5f\n", timer2);
  printf("\nTime to breed: %.5f\n", timer3);
  printf("\nTime to complete: %.5f\n", timer4);
  return 0;
}

// ============== Initialize ===========================

void Initialize(){

  ReadInput();

	// Seed random number generator
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  mutations = 0;
  best_path = zeros<urowvec>(n_cities-1);

  // Generate city coordinates
  cities = zeros<mat>(n_cities+1, 2);
  if (mode==0) CitiesOnCircumference();
  else CitiesInSquare();

  // Matrix of fitness
  fitness = zeros<rowvec>(n_replicas);
  pdf = zeros<rowvec>(n_replicas);
  cdf = zeros<rowvec>(n_replicas+1);

  // Matrix of paths
  urowvec path(n_cities-1);
  population = field<urowvec>(n_replicas);
  oldpopulation = field<urowvec>(n_replicas);

  int ic = 0;
  path.imbue([&]() {++ic; return ic;} ); // 1 ...
  for (int irep = 0; irep < n_replicas; ++irep){
    population(irep) = shuffle(path); // randomize elements
  }

  // Unit of mutation probability
  p_mut = 0.05;
  cities.save("cities.dat", raw_ascii);
  if (report==1) Welcome();

  ofstream of;
  of.open("./best_path.dat");
  of.close();
  of.open("./least_cost.dat");
  of.close();
  of.open("./avg_cost.dat");
  of.close();

  timer = 0;
  timer2 = 0;
  timer3 = 0;
  timer4 = 0;
}

void ReadInput(){

  ifstream input;
  input.open("./input.dat");
  input >> n_cities;
  input >> n_replicas;
  input >> n_generations;
  input >> mode;
  input >> report;
  input.close();
}

void Reset(){
  best_path.zeros();
  blk_avg = 0.0;
  mutations = 0;
}

// ============== Reproduction =========================

void Select(){
  int a, b, splice, n_kids;
  double r;
  oldpopulation = population;

  clock_t start, stop;
  n_kids = 0;
  // Keep best!
  population.at(0) = oldpopulation.at(selection(0));
  n_kids+=1;
  while (n_kids < n_replicas-2){

    // Find two parents
    r = rnd.Rannyu();

    start = clock();
    for (int i = 0; i < n_replicas; ++i){
      if (r < cdf(i)) {
        a = selection(i);
        break;
      }
    }

    r = rnd.Rannyu();
    for (int j = 0; j < n_replicas; ++j){
      if (r < cdf(j)){
        b = selection(j);
        break;
      }
    }

    // Make sure they're different
    while (a==b) {
      r = rnd.Rannyu();
      for (int j = 0; j < n_replicas; ++j){
        if (r < cdf(j)){
          b = selection(j);
          break;
        }
      }
    }

    stop = clock();
    timer += (double) (stop - start) / CLOCKS_PER_SEC;

    // populate
    start = clock();
    splice = (int) rnd.Rannyu(2, n_cities-2);
    Crossover(a, b, n_kids, n_kids+1, splice);
    n_kids+=2;
    stop = clock();
    timer2 += (double) (stop - start) / CLOCKS_PER_SEC;
  }
}

void Crossover(int m, int d, int s, int dg, int splice){

  clock_t start, stop;

  // son and daughter
  (population.at(s)).zeros();
  (population.at(dg)).zeros();

  // mom and dad
  // mom = (oldpopulation.at(m));
  // dad = (oldpopulation.at(d));

  start = clock();
  // Each child takes a bit of a parent
  for (int i = 0; i < splice; ++i){
    (population.at(s)).at(i) = (oldpopulation.at(m)).at(i);
    (population.at(dg)).at(i) = (oldpopulation.at(d)).at(i);
  }
  stop = clock();
  timer3 += (double) (stop - start) / CLOCKS_PER_SEC;

  start = clock();
  // And is afterwards completed
  int j = 0;
  for (int i = 0; i < n_cities-1; ++i) {
    if (!IsIn(population.at(s), (oldpopulation.at(d)).at(i))){
      (population.at(s)).at(splice+j) = (oldpopulation.at(d)).at(i);
      ++j;
    }
  }
  j = 0;
  for (int i = 0; i < n_cities-1; ++i) {
    if (!IsIn(population.at(dg), (oldpopulation.at(m)).at(i))){
      (population.at(dg)).at(splice+j) = (oldpopulation.at(m)).at(i);
      ++j;
    }
  }
  stop = clock();
  timer4 += (double) (stop - start) / CLOCKS_PER_SEC;

}

// ============== Fitness ==============================

void Rank(){
  double accumulator;

  // Calculate fitness of each path
  for (int irep = 0; irep < n_replicas; ++irep){
    fitness(irep) = Cost(irep);
  }

  // Get indices of top paths
  selection = sort_index(fitness, "ascend");

  best_path = population.at(selection(0));
  least_cost = Cost(selection(0));

  // Calculate probabilities weighted by fitness

  pdf = exp(- fitness / 100);
  pdf /= accu(pdf);
  pdf = sort(pdf, "descend");

  accumulator = 0;
  for (int i = 0; i < n_replicas; ++i){
    accumulator += pdf(i);
    cdf(i+1) = accumulator; // keep zero at i=0
  }

  cdf.save("./cdf.dat", raw_ascii);
}

double Cost(int ir){
  double cost;
  int a, b;

  cost = 0;

  // first and last
  a = 0;
  b = (population.at(ir)).at(0);
  cost += Norm(a, b);

  a = (population.at(ir)).at(n_cities-2);
  b = 0;
  cost += Norm(a, b);

  for (int ic = 0; ic < n_cities-2; ++ic){
    a = (population.at(ir)).at(ic);
    b = (population.at(ir)).at(ic+1);
    cost += Norm(a, b);
  }

  return cost;
}

double Norm(int a, int b){
  double sum;
  sum=cities(a,0)+cities(a,1)-cities(b,0)-cities(b,1);
  return abs(sum);
}

// ============== Mutation =============================

void Mutate(){
  double r;
  int o, u, m;
  // int n;
  string msg;

  for (int irep = 1; irep < n_replicas; ++irep){
    r = rnd.Rannyu(); // U[0,1)

    if (r > 0 && r < 1*p_mut){ // Swap two cities

      o = (int) rnd.Rannyu(1,n_cities-1);
      u = (int) rnd.Rannyu(1,n_cities-1);
      while (o==u) {
        u = (int) rnd.Rannyu(1,n_cities-1);
      }
      Swap(irep, o, u);
      ++mutations;

    } else if (r > 1*p_mut && r < 2*p_mut){

      m = (int) rnd.Rannyu(3,n_cities-1);
      o = (int) rnd.Rannyu(0, n_cities-1-m);
      Invert(irep, o, m);
      ++mutations;

    } else if (r > 2*p_mut && r < 5*p_mut){
      population.at(irep) = shuffle(population.at(irep));
      ++mutations;
    }
  }
}

void AllRandom(){
  for (int ir = 0; ir < n_replicas; ++ir){
    population.at(ir) = shuffle(population.at(ir));
  }
}

void Swap(int ir, int ic_a, int ic_b){
  // Swap two cities in the path
  int u = ( population.at(ir) ).at(ic_a);
  (population.at(ir) ).at(ic_a) = (population.at(ir) ).at(ic_b);
  (population.at(ir) ).at(ic_b) = u;
}

void Invert(int ir, int l, int m){
  urowvec copy(m);
  for (int i = 0; i < m; ++i){
    copy(i) = (population.at(ir)).at(i+l);
  }
  copy = reverse(copy);
  for (int i = 0; i < m; ++i){
    (population.at(ir)).at(i+l) = copy(i);
  }
}

// NOT FUNCTIONAL
//void Shift(int irep, int l, int m, int n){ // shift n cities by m
//
//  // umat path = population.row(irep);
//  // umat take = zeros<umat>(m);
//
//  // for (int i = 0; i < m; ++i){
//  //   take(i) = path(l);
//  //   path.shed_col(l);
//  // }
//
//  // if (l+n>n_cities-l) path.insert_cols(1, take.t() );
//  // else path.insert_cols(l+n, take.t() );
//
//  // for (int i = 0; i < n_cities; ++i){
//  //   population(irep, i) = path(i);
//  // }
//}
//
//void Permute(int irep, int l, int m){
//  // // m < N/2
//  // umat path = population.row(irep);
//  // umat take_a = zeros<umat>(m);
//  // umat take_b = zeros<umat>(m);
//
//  // for (int i = 0; i < m; ++i){
//  //   take_a(i) = path(l);
//  //   path.shed_col(l);
//  // }
//  // for (int i = 0; i < m; ++i){
//  //   take_b(i) = path(l);
//  //   path.shed_col(l);
//  // }
//
//  // path.insert_cols(l, take_b.t() );
//  // path.insert_cols(l+m, take_a.t() );
//
//  // for (int i = 0; i < n_cities; ++i){
//  //   population(irep, i) = path(i);
//  // }
//}

// ============== Input and Output =======================

void Save(int ig){
  ofstream of;

  // print current best path to file

  of.open("./best_path.dat", ios::app);
  of << 0; // Only necessary until i rewrite the plotting func
  best_path.save(of, raw_ascii);
  of.close();

  of.open("./least_cost.dat", ios::app);
  of << setprecision(8) << least_cost << endl;
  of.close();

  // Average of best paths
  fitness = sort(fitness);
  for (int ir = 0; ir < n_replicas/2; ++ir) {
    blk_avg += fitness(ir);
  }
  blk_avg /= (long double) n_replicas/2;

  of.open("./avg_cost.dat", ios::app);
  of << ig << " " << blk_avg << endl;
  of.close();
}

void Report(int ig){
  cout << "Generation : " << ig << endl;
  best_path.raw_print("Best solution :");
  cout << "Avg cost : " << blk_avg << endl;
  cout << "Least cost : "<< setprecision(8) << least_cost << endl;
  cout << "Mutations: " << mutations << endl;
  cout << "=========================" << endl << endl;
}

void Welcome(){

  cout << endl;

  cout << "Traveling Salesman with Genetic Algorithm" << endl;

  cout << n_cities << " cities placed ";
  if (mode==0) cout << "randomly on a circumference.";
  else cout << "randomly inside a square";
  cout << endl;

  cout << n_replicas << " replicas" << endl;
  cout << n_generations << " generations" << endl;

  cout << "Chance of mutations is " << 5*p_mut*100 << " percent" << endl;

  cout << endl << endl;
}

// ============== Other Functions =======================

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
  cities(n_cities, 0) = cities(0, 0);
  cities(n_cities, 1) = cities(0, 1);
}

void CitiesInSquare(){
  double a = 10;
  for (int i = 0; i < n_cities; ++i){
    cities(i, 0) = rnd.Rannyu(-a, a);
    cities(i, 1) = rnd.Rannyu(-a, a);
  }
  cities(n_cities, 0) = cities(0, 0);
  cities(n_cities, 1) = cities(0, 1);
}

bool IsIn(const urowvec & vec, int city){
  // Returns true
  for (int ic = 0; ic < n_cities-1; ++ic){
    if (city == vec(ic)) {
      return true;
    }
  }
  return false;
}

void Check(int ir, string msg){
  ofstream check;
  check.open("check.out", ios::app);
  check << ir << endl;
  check << msg << endl << endl;
  check << population.at(ir) << endl;
}
