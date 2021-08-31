/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "statistics.h"
#include "RandomWalk.h"

using namespace std;


int main (int argc, char *argv[]){

   // Random number generator

   Random rnd;
   int seed[4];
   int p1, p2;

   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

/*===========================================================================*/

   // 2.1.1 1D Integral via Monte Carlo

   int M = 1000000; // Rolls
   int N = 100;    // Blocks
   int L = M/N;    // Rolls per block
   double r;

   vector<double> ave(N,0);
   vector<double> av2(N,0);

   // Integrate

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
         ave[i] += (M_PI/2) * cos(M_PI * 0.5 * rnd.Rannyu()); // Integrand
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   string datafile = "basic_sampling.dat";
   blocked_stats(ave, av2, N, datafile);

   // 2.1.2 Importance Sampling with 1-x distribution

   // rnd.SetRandom(seed,p1,p2);

   M = 1000000; // Rolls
   N = 100;    // Blocks
   L = M/N;    // Rolls per block

   fill(ave.begin(), ave.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   // Integrate

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
         r = rnd.ImportanceSample();
         ave[i] +=-(M_PI/4)*cos((M_PI/2)*r) / (r - 1); // Integrand
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   datafile = "importance_sampling.dat";
   blocked_stats(ave, av2, N, datafile);

   // 2.1.3? Unimportance Sampling with x distribution

   // rnd.SetRandom(seed,p1,p2);

   M = 1000000; // Rolls
   N = 100;    // Blocks
   L = M/N;    // Rolls per block

   fill(ave.begin(), ave.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   // Integrate

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
         r = rnd.UnImportanceSample();
         ave[i] += - (M_PI/4) * cos((M_PI/2) * r) / (r);
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   datafile = "unimportance_sampling.dat";
   blocked_stats(ave, av2, N, datafile);

/*===========================================================================*/

   // 2.2.1 3D Random Walk on a Cubic Lattice

   int n_walks = 10000; // number of random walks per experiment
   int n_steps = 100;   // max number of steps

   fill(ave.begin(), ave.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   vector<int> lattice(3); // current position

   for (int i = 0; i < n_steps; i++){ // 100 random walks of steps 1, 2, ... , 100

      for (int j = 0; j < n_walks; j++) { // take many walks

         lattice[0] = 0; // start at origin
         lattice[1] = 0;
         lattice[2] = 0;

         // walk
         Discrete_Random_Walk(lattice, rnd, i+1);

         // calculate distance from origin
         ave[i] += Distance_Formula_Lattice(lattice);
      }

      ave[i] /= n_walks;
      av2[i] = ave[i] * ave[i];
   }

   // Statistics

   datafile = "rw_lattice.dat";
   blocked_stats(ave, av2, n_steps, datafile);

   // 2.2.2

   // rnd.SetRandom(seed,p1,p2);

   n_walks = 10000; // number of random walks per experiment
   n_steps = 100;   // max number of steps

   vector<double> position(3); // current position

   fill(ave.begin(), ave.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   for (int i = 0; i < n_steps; i++){ // 100 random walks of steps 1, 1, 2, ... 100

      for (int j = 0; j < n_walks; j++) { // take many walks

         position[0] = 0;   // start at origin
         position[1] = 0;
         position[2] = 0;

         Continuous_Random_Walk(position, rnd, i+1);
         ave[i] += Distance_Formula(position); // calculate distance from origin
      }

      ave[i] /= n_walks;
      av2[i] = ave[i] * ave[i];
   }

   datafile = "rw_continuum.dat";
   blocked_stats(ave, av2, n_steps, datafile);

   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
