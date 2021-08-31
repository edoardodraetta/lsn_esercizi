#include "ParallelSearch.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>

using namespace arma;
using namespace std;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  Initialize();

  for (int igen = 1; igen <= n_generations; ++igen) {

    if (igen % n_migr == 0){
      Migrate();
      RankFitness();
      if (mpi_rank == 0) Report(igen); // Print to terminal
    }

    Reset(); // Reset observables

    if (igen != 1) { // Create new generation
      Select();
      Mutate();
    }

    RankFitness(); // Evaluate fitness --> New name?
    Save(igen); // Save results to file
  }

  if (mpi_rank == 0 ){
    cities.save("cities.reload");
    cities.save("cities.dat", raw_ascii);
  }

  string seedout = "seed_" + to_string(mpi_rank) + ".out";
  rnd.SaveSeed(seedout);
  MPI_Finalize();
  return 0;
}

void Migrate(){

  // 1. Track the global least cost on the root node
  double send_cost = least_cost;
  double* recv_cost = new double[mpi_size];

  MPI_Gather(&send_cost, 1, MPI_DOUBLE, recv_cost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int best_node = 0;
  if (mpi_rank == 0 ){
    for (int i = 0; i < mpi_size; ++i) {
      cout << recv_cost[i] << " ";
      if (recv_cost[i] < glob_least_cost) {
        glob_least_cost = recv_cost[i];
        best_node = i; // keep track of where the best path is
      }
    }
    cout << endl;
  }

  // cue the best node to send its best path
  int msg_size = n_cities - 1;
  int* cue_send = new int[mpi_size];
  int cue_recv = 0;

  if (mpi_rank == 0) {
    for (int i = 0; i < mpi_size; ++i) {
      if (i == best_node) cue_send[i] = 1; // fittest node gets "1"
      else cue_send[i] = 0; // others get "0"
    }
  }

  // scatter the cues
  MPI_Scatter(&cue_send[0], 1, MPI_INT, &cue_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // better behavior with std::vec
  vector<int> send_path = conv_to<vector<int> >::from(best_path);
  vector<int> get_path = conv_to<vector<int> >::from(glob_best_path);

  if (cue_recv == 1) {
    MPI_Send(&send_path[0], msg_size, MPI_INT, 0, 6, MPI_COMM_WORLD);
  }
  if (mpi_rank == 0 && best_node != 0){
    MPI_Recv(&get_path[0], msg_size, MPI_INT, best_node, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    glob_best_path = conv_to<urowvec>::from(get_path);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  Exchange(); // Exchange best paths between nodes.
  delete[] recv_cost;
  delete[] cue_send;
}

void Exchange() {
  // Each of the four nodes needs to exchange part of its population.

  int msg_size = n_cities - 1;
  vector<int> get_path(msg_size);
  vector<int> send_path;
  urowvec good_path;
  int n_send = 4; // paths to exchange

  for (int i = 0; i < n_send; ++i) {

    good_path = population.at(selection(i)); // path to send
    send_path = conv_to<vector<int> >::from(good_path);

    if (mpi_rank == 0){
      MPI_Send(&send_path[0], msg_size, MPI_INT, 1, 1, MPI_COMM_WORLD);
      MPI_Recv(&get_path[0], msg_size, MPI_INT, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (mpi_rank == 1){

      MPI_Recv(&get_path[0], msg_size, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&send_path[0], msg_size, MPI_INT, 0, 2, MPI_COMM_WORLD);

    } else if (mpi_rank == 2){

      MPI_Send(&send_path[0], msg_size, MPI_INT, 3, 3, MPI_COMM_WORLD);
      MPI_Recv(&get_path[0], msg_size, MPI_INT, 3, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    else if (mpi_rank == 3){

      MPI_Recv(&get_path[0], msg_size, MPI_INT, 2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&send_path[0], msg_size, MPI_INT, 2, 4, MPI_COMM_WORLD);

    }

    population.at(selection(n_replicas-1-i)) = conv_to<urowvec>::from(get_path); // copy to worst slots of population
  }

}

void Report(int ig) {
  cout << "Generation : " << ig << endl;
  // cout << "(This is node " << mpi_rank << " reporting.)" << endl;
  cout << "The best path so far has a cost of " << glob_least_cost << endl;
  // best_path.print("Best solution :");
  cout << "=========================" << endl;
}

void InitRNG() {
  int p1, p2;

  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  input.close();

  // unique seed for each node
  for (int i = 0; i < 4; ++i) seed[i] *= mpi_rank;

  rnd.SetRandom(seed, p1, p2);
  arma_rng::set_seed(0); // seed armadillo
}

void Initialize() {

  ReadInput(); // read input.dat

  InitRNG(); // Seed random number generator

  // Generate city coordinates
  cities = zeros<mat>(n_cities, 2);
  if (reload == 1) {
    cout << "Reload" << endl;
    cities.load("cities.reload");
  } else if (mpi_rank == 0){
    MakeCities(mode);
  }

  // Matrix of fitness
  best_path = zeros<urowvec>(n_cities - 1);
  fitness = zeros<rowvec>(n_replicas);

  // For choosing parents
  pdf = zeros<rowvec>(n_replicas);
  cdf = zeros<rowvec>(n_replicas + 1);

  // tracking
  glob_least_cost = 100000.0;
  glob_best_path = zeros<urowvec>(n_cities - 1);

  // Matrix of paths
  urowvec path(n_cities - 1);
  population = field<urowvec>(n_replicas);
  oldpopulation = field<urowvec>(n_replicas);

  // populate and then shuffle
  int ic = 0;
  path.imbue([&]() { // 1 ... n_city-1
    ++ic;
    return ic;
  });
  for (int irep = 0; irep < n_replicas; ++irep) {
    population(irep) = shuffle(path); // randomize elements
  }

  mutations = 0; // number of mutations

  // Output files

  bestpathfile = "./best_path_" + to_string(mpi_rank) + ".dat";
  leastcostfile = "./least_cost_" + to_string(mpi_rank) + ".dat";
  avgcostfile = "./avg_cost_" + to_string(mpi_rank) + ".dat";

  ofstream of;
  of.open(bestpathfile);
  of.close();
  of.open(leastcostfile);
  of.close();
  of.open(avgcostfile);
  of.close();

  if (mpi_rank == 0) {
    of.open("glob_least_cost.dat");
    of.close();
    of.open("glob_best_path.dat");
    of.close();
    Welcome();
  }
}

void Reset() { // Reset observables
  best_path.zeros();
  blk_avg = 0.0;
  mutations = 0;
}

void RankFitness() {
  double accumulate;

  // Calculate fitness of each path
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
  accumulate = 0;
  for (int i = 0; i < n_replicas; ++i) {
    accumulate += pdf(i);
    cdf(i + 1) = accumulate; // keep zero at i=0
  }
}

void Select() { // Select parents
  int a, b, splice, n_kids;
  double r;

  oldpopulation = population;
  // Keep the best parents
  for (int i = 0; i < 4; ++i){
    population.at(i) = oldpopulation.at(selection(i));
  }

  n_kids = 4;

  // Fill new population with children
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
    splice = (int)rnd.Rannyu(2, n_cities - 2); // splicing location
    Crossover(a, b, n_kids, n_kids + 1, splice);
    n_kids += 2;
  }
}

void Mutate() {
  double r;
  int o, u, m, n;
  string msg;

  for (int irep = 1; irep < n_replicas; ++irep) {
    r = rnd.Rannyu();             // U[0,1)

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

// Save results to file
void Save(int ig) {
  ofstream Bp, Lc, Ac, Glc, Gbp;

  // print current best path to file
  Bp.open(bestpathfile, ios::app);
  best_path.save(Bp, raw_ascii);
  Bp.close();

  // ... and its cost
  Lc.open(leastcostfile, ios::app);
  Lc << least_cost << endl;
  Lc.close();

  // global best path and least cost
  if (mpi_rank == 0 ){
    if (ig == 1) {
      glob_least_cost = least_cost;
      glob_best_path = best_path;
    }
    Glc.open("glob_least_cost.dat", ios::app);
    Glc << glob_least_cost << endl;
    Glc.close();

    Gbp.open("glob_best_path.dat", ios::app);
    glob_best_path.save(Gbp, raw_ascii);
    Gbp.close();
  }

  // Average of top half paths paths
  fitness = sort(fitness);
  for (int ir = 0; ir < n_replicas / 2; ++ir) {
    blk_avg += fitness(ir);
  }
  blk_avg /= (long double)n_replicas / 2;

  Ac.open(avgcostfile, ios::app);
  Ac << " " << blk_avg << endl;
  Ac.close();
}

// Crossover operator for generating new generation
void Crossover(int m, int d, int s, int dg, int splice) {
  // son and daughter
  (population.at(s)).zeros();
  (population.at(dg)).zeros();
  // mom and dad
  // mom = (oldpopulation.at(m));
  // dad = (oldpopulation.at(d));
  // Each child takes a bit of a parent
  for (int i = 0; i < splice; ++i) {
    (population.at(s)).at(i) = (oldpopulation.at(m)).at(i);
    (population.at(dg)).at(i) = (oldpopulation.at(d)).at(i);
  }
  // And is afterwards completed
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

// Cost function for fitness
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

// L2 norm
double Norm(int a, int b) {
  double sum;
  sum = (cities(a, 0) - cities(b, 0)) * (cities(a, 0) - cities(b, 0));
  sum += (cities(a, 1) - cities(b, 1)) * (cities(a, 1) - cities(b, 1));
  return abs(sum);
}

// === MUTATION OPERATORS ===

// Swap two cities in the path
void Swap(int ir, int ic_a, int ic_b) {
  int u = (population.at(ir)).at(ic_a);
  (population.at(ir)).at(ic_a) = (population.at(ir)).at(ic_b);
  (population.at(ir)).at(ic_b) = u;
}

// Invert m cities at l
void Invert(int ir, int l, int m) {
  urowvec copy(m);
  for (int i = 0; i < m; ++i) {
    copy(i) = (population.at(ir)).at(i + l);
  }
  // copy = reverse(copy);
  for (int i = 0; i < m; ++i) {
    (population.at(ir)).at(i + l) = copy(m-i-1);
  }
}

// Shift n cities by m
void Shift(int ir, int l, int m, int n) {
  urowvec path = population.at(ir);
  urowvec take = zeros<urowvec>(n);
  // take n cities
  for (int i = 0; i < n; ++i) {
    take(i) = path(l);
    path.shed_col(l);
  }
  // place them at l+m
  path.insert_cols(l + m, take);
  for (int i = 0; i < n_cities - 1; ++i) {
    (population.at(ir)).at(i) = path(i);
  }
}

// Swap two sets of cities
void Permute(int ir, int l, int k, int m) {
  urowvec path = population.at(ir);
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
  for (int i = 0; i < n_cities - 1; ++i) {
    (population.at(ir)).at(i) = path(i);
  }
}

// ===========================

void ReadInput() { // Read input.dat
  ifstream input;
  input.open("./input.dat");
  input >> n_cities;
  input >> n_replicas;
  input >> n_generations;
  input >> p_mut;
  input >> mode;
  input >> reload;
  input >> n_migr;
  input.close();
}

void Welcome() { // Welcome msg
  cout << endl;
  cout << "Traveling Salesman with Genetic Algorithm with Parallel Search" << endl;
  cout << "(This is node " << mpi_rank << " of " << mpi_size << " speaking.)" << endl;
  cout << n_cities << " cities placed ";
  if (mode == 0) cout << "randomly on a circumference.";
  else cout << "randomly inside a square";
  cout << endl;
  if (reload == 1) cout << "Configuration reloaded from file." << endl;
  cout << n_replicas << " replicas" << endl;
  cout << n_generations << " generations" << endl;
  cout << "Chance of mutation is " << 5 * p_mut * 100 << " percent." << endl;
  cout << endl << endl;
}

void MakeCities(int MODE) {
  // Make cities on root node
  if (mpi_rank == 0) {
    if (MODE == 0) { // on a circle
      double x, y, theta, R;
      R = 10;
      for (int i = 0; i < n_cities; ++i) {
        theta = rnd.Rannyu(0, 2 * pi);
        x = R * cos(theta);
        y = R * sin(theta);
        cities(i, 0) = x;
        cities(i, 1) = y;
      }
    } else if (MODE == 1) { // in a square
      double a = 10;
      for (int i = 0; i < n_cities; ++i) {
        cities(i, 0) = rnd.Rannyu(-a, a);
        cities(i, 1) = rnd.Rannyu(-a, a);
      }
    }
  }

  // Broadcast city coordinates to all nodes
  MPI_Bcast(&cities[0], n_cities*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Returns true if city is in vec
bool IsIn(const urowvec &vec, int city) {
  for (int ic = 0; ic < n_cities - 1; ++ic) {
    if (city == vec(ic)) {
      return true;
    }
  }
  return false;
}
