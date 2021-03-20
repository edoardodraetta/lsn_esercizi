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
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;

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

   // Uniform:

   // pass to a function?
   string rng_file = "./random.out";
   ofstream outFile;
   outFile.open(rng_file);

   int n_random = 1000000;

   // importante fare così, cioè non tenerseli in memoria...
   for (int i=0; i<n_random; i++){
      outFile << rnd.Rannyu() << endl;
   }
   outFile.close();

   cout << "Uniformly distributed random numbers are in '" << rng_file << "'." << endl;

   // Lorenzian:

   rng_file = "./random_lor.out";
   outFile.open(rng_file);

   for (int i=0; i<n_random; i++){
      outFile << rnd.Lorentzian(0,1) << endl;
   }
   outFile.close();
   cout << "Lorentzian random numbers are in '" << rng_file << "'." << endl;

   // Exponential:

   rng_file = "./random_exp.out";
   outFile.open(rng_file);

   for (int i=0; i<n_random; i++){
      outFile << rnd.Exponential(1) << endl;
   }
   outFile.close();
   cout << "Exponential random numbers are in '" << rng_file << "'." << endl;

   // Theta:

   rng_file = "./random_theta.out";
   outFile.open(rng_file);

   for (int i=0; i<n_random; i++){
      outFile << rnd.RanTheta() << endl;
   }
   outFile.close();
   cout << "Random 2D angles are in '" << rng_file << "'." << endl;

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
