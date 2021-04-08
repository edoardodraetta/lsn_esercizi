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
#include <algorithm>

#include "random.h"
#include "statistics.h"
#include "GeometricBrownianMotion.h"

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

   // Params

   float S = 100;
   float S_t;
   float T = 1;
   float K = 100;
   float r = 0.1; // mu = r
   float sigma = 0.25;

   // Black-Scholes Solution to Option Pricing

   cout << "== Black-Scholes ==" << endl;
   cout << "Call Price: " << call_option(S,T,K,sigma,r) << endl;
   cout << "Put Price: " << put_option(S,T,K,sigma,r) << endl;

   // 3.1.1 - Direct Sampling of GBM

   int M = 100000;
   int N = 100;

   vector<double> ave(N,0);
   vector<double> av2(N,0);

   for (int i = 0; i < N; i++){
      S_t = 100;
      for (int j = 0; j < M; j++){
         S_t = GBM(S,T,r,sigma,rnd);

         ave[i] += exp(-r*T)* max(0,K-S_t);
      }

      ave[i] /= M; // average price

      av2[i] = ave[i]*ave[i];
   }

   // Cumulative blocked statistics

   string datafile = "../data/stats_3.1.1.dat";
   blocked_stats(ave, av2, N, datafile);

   // 3.1.2 - Path Sampling of GBM

   M = 10000;
   int timesteps = 100;
   float t = T / timesteps;

   fill(ave.begin(), ave.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   for (int i = 0; i < N; i++){

      for (int j = 0; j < M; j++){
         S_t = 100;
         for (int k = 0; k < timesteps; k++){   
            S_t = GBM(S_t,t,r,sigma,rnd);
         }
         ave[i] += exp(-r*T)* max(0, S_t-K); // max doesn't work this way
      }
      ave[i] /= M;
      av2[i] = ave[i]*ave[i];
   }

   datafile = "../data/stats_3.1.2.dat";
   blocked_stats(ave, av2, N, datafile);

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
