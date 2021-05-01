#include <vector>
#include <iostream>
#include <fstream>
#include "random.h"
#include "metropolis.h"

using namespace std;

void printpos(double vec[3]){
	cout << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
	cout << endl;
}

int main(){

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

 /*===========================================================================*/

	int N = 10000;
	double r0[3] = {1.,0.,0.};
	double rold[3];
	double r[3];
	vector<vector<double>> d(10);
	for (int i=0; i<3; i++) r[i] = r0[i];

	double alpha;
	int count = 0;

	ofstream printpos("./wavewalk.dat");

	for (int i = 0; i < N; i++){
		// 1. keep old position
		for (int j=0; j<3; j++) rold[j] = r[j];
		// 2. try a new position
		rnd.RanUniform3d(.7, r);
		// 3. calculate acceptance rate
		alpha = groundstate(r) / groundstate(rold);
		// 4. keep or discard
		if (rnd.Rannyu() > alpha) {
			for (int j=0; j<3; j++) r[j] = rold[j];
		} else {
			count++;
		}
		// 5. print to file
		for (int j = 0; j < 3; j++){
			printpos << r[j] << " ";
		}
		printpos << endl;
	}
	printpos.close();
	cout << "acceptance rate:" << (double)count / N << endl;

	count = 0;
	printpos.open("./wavewalk2.dat");
	for (int i = 0; i < N; i++){
		// 1. keep old position
		for (int j=0; j<3; j++) rold[j] = r[j];
		// 2. try a new position
		rnd.RanUniform3d(1., r);
		// 3. calculate acceptance rate
		alpha = excitedstate(r) / excitedstate(rold);
		// 4. keep or discard
		if (rnd.Rannyu() > alpha) {
			for (int j=0; j<3; j++) r[j] = rold[j];
		} else {
			count++;
		}
		// 5. print to file
		for (int j = 0; j < 3; j++){
			printpos << r[j] << " ";
		}
		printpos << endl;
	}

	cout << "acceptance rate:" << (double)count / N << endl;
	return 0;
}


