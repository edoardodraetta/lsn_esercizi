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

   // 1.1.1 First Integral

   int M = 100000;
   int N = 100;
   int L = M/N;
   vector<double> ave(N,0);
   vector<double> av2(N,0);

   // Integrate

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
         ave[i] += rnd.Rannyu();
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   string datafile = "./data/stats_1.1.1.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.2 Second Integral

   rnd.SetRandom(seed,p1,p2);
   fill(ave.begin(), av2.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   // Integrate

   double r;
   for (int i=0; i < N; i++) {
      for (int j=0; j < L; j++) {
         r = rnd.Rannyu();
         ave[i] += (r-0.5) * (r-0.5);
      }
      ave[i] /= L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Statistics

   datafile = "./data/stats_1.1.2.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.3 Chi Squared Test

   M = 1000000; // throws
   N = 100;     // blocks
   L = M/N;     // throws per block
   int n = L/N; // Expectation value for chi2 test
   int count;   // Experimental value for chi2 test
   float a,b;   // Lower and upper bounds

   float chi2_sum = 0;
   vector<double> chi2_prog(N,0);

   ofstream statsfile;
   statsfile.open("./data/stats_1.1.3.dat");

   for (int i=0; i<N; i++){
      count = 0;
      float r; // ok ?
      for (int j=0; j<L; j++){

         r = rnd.Rannyu();
         a = (float)i / N;
         b = (float)(i+1) / N;

         if ( (a <= r ) && ( r < b ) ){
            count +=1;
         }
      }
      chi2_sum  += ((count - n) * (count - n)) /n;
      chi2_prog[i] = chi2_sum;
      statsfile << count << " " << chi2_sum << " " << chi2_prog[i] << endl;
   }

   statsfile.close();

   // 1.2 The Central Limit Theorem

   rnd.SetRandom(seed,p1,p2);

   int n_sums = 10000; // We sum m random numbers n_sums times
   vector<int> m_numbers = {1,2,10,100};

   // Uniform
   statsfile.open("./data/stats_1.2.1.dat");
   for (int i = 0; i < n_sums; i++){
      for (int j = 0; j < m_numbers.size(); j++){
         statsfile << rnd.sum_uniform(m_numbers[j],0,1) << " ";
      }
      statsfile << endl;
   }
   statsfile.close();

   // Exponential
   statsfile.open("./data/stats_1.2.2.dat");
   for (int i = 0; i < n_sums; i++){
      for (int j = 0; j < m_numbers.size(); j++){
         statsfile << rnd.sum_exponential(m_numbers[j],1) << " ";
      }
      statsfile << endl;
   }
   statsfile.close();

   // Lorentzian
   statsfile.open("./data/stats_1.2.3.dat");
   for (int i = 0; i < n_sums; i++){
      for (int j = 0; j < m_numbers.size(); j++){
         statsfile << rnd.sum_lorentzian(m_numbers[j],0,1) << " ";
      }
      statsfile << endl;
   }
   statsfile.close();

   // 1.3 Buffon Experiment

   rnd.SetRandom(seed,p1,p2);

   M = 10000000; // throws
   N = 100;      // blocks
   int T = M/N;  // throws per block
   float d = 10;     // distance between gridlines
   float length = 5; // needle length

   double t,y1,y2; // random numbers
   int hits;       // intersections in each block

   vector<double> pi_estimate(N,0);
   vector<double> pi2_estimate(N,0);

   // Run the experiment
   for (int i = 0; i < N; i++){
      hits = 0;
      for (int j = 0; j < T; j++){
         t = rnd.RanTheta();
         y1 = d * rnd.Rannyu();
         y2 = y1 + length * sin(t);
         if ((y1<0 || y1>10) || (y2<0 || y2>10)) {
            hits += 1;
         }
      }
      pi_estimate[i] = make_pi(hits, T, d, length);
      pi2_estimate[i] = pi_estimate[i] * pi_estimate[i];
   }

   // Statistics
   datafile = "./data/stats_1.3.dat";
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
