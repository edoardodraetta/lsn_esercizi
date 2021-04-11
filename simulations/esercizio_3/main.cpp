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

   double S = 100;
   double S_t;
   double T = 1;
   double K = 100;
   double r = 0.1; // mu = r
   double sigma = 0.25;

   // Black-Scholes Solution to Option Pricing

   cout << "== Black-Scholes ==" << endl;
   cout << "Call Price: " << call_option(S,T,K,sigma,r) << endl;
   cout << "Put Price: " << put_option(S,T,K,sigma,r) << endl;

   // 3.1.1 - Direct Sampling of GBM

   int M = 10000;
   int N = 100;

   double diff;

   vector<double> call_ave(N,0);
   vector<double> call_av2(N,0);
   vector<double> put_ave(N,0);
   vector<double> put_av2(N,0);

   for (int i = 0; i < N; i++){
      S_t = 100;
      for (int j = 0; j < M; j++){
         S_t = GBM(S,T,r,sigma,rnd);

         diff = S_t-K;
         call_ave[i] += exp(-r*T)*max(0., diff); // Call Price

         diff = K-S_t;
         put_ave[i] += exp(-r*T)*max(0., diff); // Put Price
      }

      call_ave[i] /= M;
      call_av2[i] = call_ave[i]*call_ave[i];

      put_ave[i] /= M;
      put_av2[i] = put_ave[i]*put_ave[i];
   }

   cout << "== Direct Sampling of S(T) ==" << endl;
   cout << "Call Price: " << call_ave.back() << endl;
   cout << "Put Price: " << put_ave.back() << endl;

   // Cumulative blocked statistics

   string datafile = "../../data/stats_3.1.1_call.dat";
   blocked_stats(call_ave, call_av2, N, datafile);

   datafile = "../../data/stats_3.1.1_put.dat";
   blocked_stats(put_ave, put_av2, N, datafile);

   // 3.1.2 - Path Sampling of GBM

   rnd.SetRandom(seed,p1,p2);

   M = 10000;
   int timesteps = 100;
   float t = T / timesteps;

   fill(call_ave.begin(), call_ave.end(), 0);
   fill(call_av2.begin(), call_av2.end(), 0);
   fill(put_ave.begin(), put_ave.end(), 0);
   fill(put_av2.begin(), put_av2.end(), 0);

   for (int i = 0; i < N; i++){ // Blocks

      for (int j = 0; j < M; j++){ // Sampling the future price

         S_t = 100;
         for (int k = 0; k < timesteps; k++){ // Path simulation
            S_t = GBM(S_t, t, r, sigma, rnd);
         }

         diff = S_t-K;
         call_ave[i] += exp(-r*T)*max(0., diff); // Call Price

         diff = K-S_t;
         put_ave[i] += exp(-r*T)*max(0., diff); // Put Price

      }
      call_ave[i] /= M;
      call_av2[i] = call_ave[i]*call_ave[i];

      put_ave[i] /= M;
      put_av2[i] = put_ave[i]*put_ave[i];
   }

   cout << "== Path Sampling of S(T) ==" << endl;
   cout << "Call Price: " << call_ave.back() << endl;
   cout << "Put Price: " << put_ave.back()<< endl;

   // Statistics

   datafile = "../../data/stats_3.1.2_call.dat";
   blocked_stats(call_ave, call_av2, N, datafile);

   datafile = "../../data/stats_3.1.2_put.dat";
   blocked_stats(put_ave, put_av2, N, datafile);


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
