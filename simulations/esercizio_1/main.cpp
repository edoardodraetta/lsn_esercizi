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
            rnd.SetRandom(seed,p1,p2); // Initialize
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

   // 1.1.1 First Integral

   int M = 100000; // Rolls
   int N = 100;    // Blocks
   int L = M/N;    // Rolls per Block
   vector<double> ave(N,0);
   vector<double> av2(N,0);

   // Integrate

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
         ave[i] += rnd.Rannyu(); // Integrand
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   string datafile = "../../data/stats_1.1.1.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.2 Second Integral

   // rnd.SetRandom(seed,p1,p2);       // Reinitialize
   fill(ave.begin(), av2.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   // Integrate

   double r;
   for (int i=0; i < N; i++) {
      for (int j=0; j < L; j++) {
         r = rnd.Rannyu();
         ave[i] += (r-0.5) * (r-0.5); // Integrand
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Statistics

   datafile = "../../data/stats_1.1.2.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.3 Chi Squared Test

   // rnd.SetRandom(seed,p1,p2);

   int n = 10000;
   M = 100;
   int E = n/M;
   // double r;
   int count;
   double a,b;   // Lower and upper bounds
   double chi_squared;

   ofstream statsfile;
   statsfile.open("../../data/stats_1.1.3.dat");

   for (int k = 0; k < M; k++ ){       // iterations of chi^2 test
      chi_squared = 0;
      for (int i = 0 ; i < M; i++){    // sub-intervals

         a = (double) i / M;
         b = (double) (i+1) / M;
         count = 0;

         for (int j = 0; j < n; j++){  // count hits
            r = rnd.Rannyu();
            if ( (a <= r ) && ( r < b ) ) count += 1; // hits
         }

         chi_squared += ((count - E) * (count - E));
      }

      chi_squared /= E;
      statsfile << chi_squared << endl;
   }

   statsfile.close();

   // 1.2 Testing the Central Limit Theorem

   // rnd.SetRandom(seed,p1,p2);

   int n_sums = 10000;     // Number of sums to compute
   vector<int> sums = {1,2,10,100};

   // Uniform Distribtution

   statsfile.open("../../data/stats_1.2.1.dat");
   for (int i = 0; i < n_sums; i++){
      for (int j = 0; j < sums.size(); j++){
         statsfile << rnd.sum_uniform(sums[j],0,1) << " ";
      }
      statsfile << endl;
   }
   statsfile.close();

   // Exponential Distribtution

   statsfile.open("../../data/stats_1.2.2.dat");
   for (int i = 0; i < n_sums; i++){
      for (int j = 0; j < sums.size(); j++){
         statsfile << rnd.sum_exponential(sums[j],1) << " ";
      }
      statsfile << endl;
   }
   statsfile.close();

   // Lorentzian Distribtution

   statsfile.open("../../data/stats_1.2.3.dat");
   for (int i = 0; i < n_sums; i++){
      for (int j = 0; j < sums.size(); j++){
         statsfile << rnd.sum_lorentzian(sums[j],0,1) << " ";
      }
      statsfile << endl;
   }
   statsfile.close();

   // 1.3 Buffon Experiment

   // rnd.SetRandom(seed,p1,p2); // Reinitialize rng

   M = 10000000;     // throws
   N = 100;          // blocks
   int T = M/N;      // throws per block
   float d = 10;     // distance between gridlines
   float length = 5; // needle length

   double t,y1,y2; // random numbers
   int hits;       // intersections in each block

   vector<double> pi_estimate(N,0);
   vector<double> pi2_estimate(N,0);

   for (int i = 0; i < N; i++){
      hits = 0;
      for (int j = 0; j < T; j++){

         t = rnd.RanTheta();        // needle rotation
         y1 = d * rnd.Rannyu();     // place one end of the needle
         y2 = y1 + length * sin(t); // calculate the other end

         if ((y1<0 || y1>10) || (y2<0 || y2>10)) { // check for hits
            hits += 1;
         }
      }
      pi_estimate[i] = make_pi(hits, T, d, length);
      pi2_estimate[i] = pi_estimate[i] * pi_estimate[i];
   }

   // Statistics

   datafile = "../../data/stats_1.3.dat";
   blocked_stats(pi_estimate, pi2_estimate, N, datafile);

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
