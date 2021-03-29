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
#include "functions.h"
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

   cout << "Primes: " << p1 << " " << p2 << endl;
   cout << "Seed: ";
   for (int i=0; i<4; i++){
      cout << seed[i] << " ";
   }
   cout << endl;

   // 2.1.1 1D Integral via Monte Carlo

   int M = 100000;
   int N = 100;
   int L = M/N;
   vector<double> ave(N,0);
   vector<double> av2(N,0);

   // Integrate

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
         ave[i] += (M_PI/2) * cos(M_PI * (.5) * rnd.Rannyu());
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   string datafile = "./data/stats_2.1.1.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.2 Importance Sampling

   // 2.2.1 3D Random Walk on a Cubic Lattice

   // test

   for (int i = 0; i < 10; i++){

   }




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
