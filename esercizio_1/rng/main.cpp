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

using namespace std;

float make_pi(int hits, int n_throws, float d, float length){
   if (hits==0) return 0;
   else {
      return (2*length*n_throws) / (hits * d);
   }
}

float error(vector<float> AV, vector<float> AV2, int n){
   if (n == 0){
      return 0;
   } else {
      return sqrt( (AV2[n] - AV[n]*AV[n]) / n);
   }
}

void blocked_stats(vector<float>  AV, vector<float>  AV2, int N, string filename){
   ofstream outfile;
   outfile.open(filename);

   vector<float> sum(N,0);
   vector<float> sum2(N,0);
   vector<float> err(N,0);

   for (int i=0; i<N; i++){
      for (int j=0; j<(i+1); j++){
         sum[i] += AV[j]; // cumulative sum of averages
         sum2[i] += AV2[j]; // sum of square averages
      }
      sum[i] /= (i+1); // cumulative average
      sum2[i] /= (i+1); // cumulative square average
      err[i] = error(sum, sum2, i); // uncertainty
      outfile << sum[i] << " " << sum2[i] << " " << err[i] << endl;
   }
   outfile.close();
}

int main (int argc, char *argv[]){


   Random rnd;
   int seed[4];
   int p1, p2;

   // Read Primes and Seeds.in

   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   cout << "Primes: " << p1 << " " << p2 << endl;

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

   cout << "Seed: ";
   for (int i=0; i<4; i++){
      cout << seed[i] << " ";
   }
   cout << endl;

   // 1.1.1

   int M = 100000;
   int N = 100;
   int L = M/N;

   vector<float> ave(N,0);
   vector<float> av2(N,0);
   vector<float> sum_prog(N,0);
   vector<float> su2_prog(N,0);
   vector<float> err_prog(N,0);

   float sum, r;
   string datafile;

   // Integrate

   for (int i = 0; i < N; i++) {
      sum = 0;
      for (int j = 0; j < L; j++) {
         sum += rnd.Rannyu();
      }
      ave[i] = sum/L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Cumulative blocked statistics

   datafile = "./data/stats_1.1.1.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.2

   rnd.SetRandom(seed,p1,p2);

   fill(ave.begin(), av2.end(), 0);
   fill(av2.begin(), av2.end(), 0);

   // Integrate

   for (int i=0; i < N; i++) {
      sum = 0;
      for (int j=0; j < L; j++) {
         r = rnd.Rannyu();
         sum += (r-0.5) * (r-0.5);
      }
      ave[i] = sum/L;
      av2[i] = (ave[i]*ave[i]);
   }

   // Statistics

   datafile = "./data/stats_1.1.2.dat";
   blocked_stats(ave, av2, N, datafile);

   // 1.1.3 Chi Squared Test

   // M = 1000000;
   // L = M/N;
   // float chi2_sum = 0;
   // float chi2_prog[N] = {0};
   // int n = L/N;
   // int count;


   // statsfile.open("./data/stats_1.1.3.dat");

   // for (int i=0; i<N; i++){
   //    count = 0;
   //    float r; // ok ?
   //    for (int j=0; j<L; j++){

   //       r = rnd.Rannyu();
   //       a = (float)i / (float)N;
   //       b = (float)(i+1) / (float)N;

   //       if ( (a <= r ) && ( r < b ) ){
   //          count +=1;
   //       }
   //    }
   //    // cout << count << endl;
   //    chi2_sum  += ((count - n) * (count - n)) /n;
   //    chi2_prog[i] = chi2_sum;
   //    statsfile << count << " " << chi2_sum << " " << chi2_prog[i] << endl;
   // }

   // statsfile.close();

   // // 1.2 The Central Limit Theorem

   rnd.SetRandom(seed,p1,p2);

   int n_sums = 10000;

   int M1 = 1;
   int M2 = 2;
   int M3 = 10;
   int M4 = 100;

   ofstream statsfile;
   statsfile.open("./data/stats_1.2.1.dat");

   for (int i = 0; i < n_sums; i++){
      statsfile << rnd.sum_uniform(M1,0,1) << " ";
      statsfile << rnd.sum_uniform(M2,0,1) << " ";
      statsfile << rnd.sum_uniform(M3,0,1) << " ";
      statsfile << rnd.sum_uniform(M4,0,1);
      statsfile << endl;
   }

   statsfile.close();

   statsfile.open("./data/stats_1.2.2.dat");
   for (int i = 0; i < n_sums; i++){
      statsfile << rnd.sum_exponential(M1,1) << " ";
      statsfile << rnd.sum_exponential(M2,1) << " ";
      statsfile << rnd.sum_exponential(M3,1) << " ";
      statsfile << rnd.sum_exponential(M4,1) << endl;

   }
   statsfile.close();

   statsfile.open("./data/stats_1.2.3.dat");
   for (int i = 0; i < n_sums; i++){
      statsfile << rnd.sum_lorentzian(M1,0,1) << " ";
      statsfile << rnd.sum_lorentzian(M2,0,1) << " ";
      statsfile << rnd.sum_lorentzian(M3,0,1) << " ";
      statsfile << rnd.sum_lorentzian(M4,0,1) << endl;

   }
   statsfile.close();

   // // 1.3 Buffon Experiment

   // rnd.SetRandom(seed,p1,p2);

   // // Test

   // int n_throws = 1000000;  // Needles
   // float d = 10;   // Grid spacing
   // float length = d/2;   // Needle length
   // int hits = 0;

   // float y1,y2,t;

   // for (int i = 0; i < n_throws; i++){
   //    t = rnd.RanTheta();
   //    y1 = d * rnd.Rannyu();
   //    y2 = y1 + length * sin(t);
   //    if ((y1<0 || y1>10) || (y2<0 || y2>10)) {
   //       hits += 1;
   //    }
   // }

   // cout << make_pi(hits,n_throws, d, length) << endl;

   // // Blocked Experiment

   // rnd.SetRandom(seed,p1,p2);
   // n_throws = 10000000;
   // // N = 100
   // int T = n_throws/N;

   // float pi_estimate[N] = {0};
   // float pi2_estimate[N] = {0};

   // memset(sum_prog, 0, sizeof(sum_prog));
   // memset(su2_prog, 0, sizeof(su2_prog));
   // memset(err_prog, 0, sizeof(err_prog));


   // statsfile.open("./data/stats_1.3_pi.dat");

   // for (int i = 0; i < N; i++){
   //    hits = 0;
   //    for (int j = 0; j < T; j++){
   //       t = rnd.RanTheta();
   //       y1 = d * rnd.Rannyu();
   //       y2 = y1 + length * sin(t);
   //       if ((y1<0 || y1>10) || (y2<0 || y2>10)) {
   //          hits += 1;
   //       }
   //    }


   //    pi_estimate[i] = make_pi(hits, T, d, length);
   //    pi2_estimate[i] = pi_estimate[i] * pi_estimate[i];
   //    statsfile << hits << " " << pi_estimate[i] << " " << pi2_estimate[i] << endl;
   // }

   // statsfile.close();

   // // Statistics

   // statsfile.open("./data/stats_1.3.dat");

   // for (int i=0; i<N; i++){
   //    for (int j=0; j<(i+1); j++){
   //       sum_prog[i] += pi_estimate[j]; // cumulative sum of averages
   //       su2_prog[i] += pi2_estimate[j]; // sum of square averages
   //    }
   //    sum_prog[i] /= (i+1); // cumulative average
   //    su2_prog[i] /= (i+1); // cumulative square average
   //    err_prog[i] = error(sum_prog, su2_prog, i); // uncertainty

   //    statsfile << sum_prog[i] << " " << su2_prog[i] << " " << err_prog[i] << endl;
   // }

   // statsfile.close();

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
