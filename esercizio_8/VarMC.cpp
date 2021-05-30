#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "random.h"
#include "VarMC.h"

int main(){

	Initialize(); // read from input.dat
	for (int iblock = 1; iblock <= m_blocks; ++iblock){
		Reset(iblock);

		for (int istep = 1; istep <= nsteps; ++istep){

			Move(); // Try a move

			if(istep%imeasure==0) Accumulate(iblock); // measure

			if(istep%iprint==0) PrintPos(r);
		}

		Average(iblock); // compute average for current block
		if(iblock%10==0) Report(iblock);
	}

	return 0;
}

void Initialize(){

	// Initialize RNG
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

  // Input File
  ifstream ReadInput("./input.dat");
  ReadInput >> state;
  ReadInput >> mode;
  ReadInput >> dr;
  ReadInput >> nsteps;
  ReadInput >> imeasure;
  ReadInput >> iprint;
  ReadInput.close();

  // Initialize params
  if (state==0) {
    r0[0] = 1.2;
    r0[1] = 0;
    r0[2] = 0;
  } else {
    r0[0] = 0;
    r0[1] = 0;
    r0[2] = 5;
  }

  for (int i=0; i<3; ++i) r[i] = r0[i];

  Welcome(); // print summary to terminal
}

void Welcome(){

	cout << "Sampling Wavefunctions of the Hydrogen Atom with the Metropolis Algorithm" << endl;

  if(state==0) {
  	cout << "Sampling ground state wavefunction.";
  } else {
  	cout << "Sampling excited state wavefunction.";
  }
  cout << endl;

  cout << "Walking " << nsteps << " steps in each block." << endl;
  cout << "Size of step : " << dr << endl;
  cout << "Recording every " << iprint << " steps" << endl;
  cout << "Measuring every " << imeasure << " steps" << endl;

  if(mode==0) {
  	cout << "Moves attempted with uniform probability.";
  } else {
  	cout << "Moves attempted with Gaussian probability.";
  }
  cout << endl;

  cout << "Starting position : ";
  cout << "a_0 * (" << r0[0] << ", " << r0[1] << ", " << r0[2] << ")";

  cout << endl << endl;
}

void Reset(int iblock){
	accepted = 0;
	attempted = 0;
	ofstream pos;
	pos.open("output_positions.dat");
	pos.close();
	ave_pos[iblock-1] = 0;
	av2_pos[iblock-1] = 0;
}

void Report(int iblock){
	cout << "Block " << iblock << ", "  << "Acceptance Rate: ";
	cout << (double)accepted/attempted << endl;
	cout << "---------------------------------" << endl;
}

void Move(){ // attempt move with Metropolis Algorithm

	for (int j=0; j<3; j++) rold[j] = r[j];

  if (mode==0) rnd.RanUniform3d(dr, r); // Uniform Transition Matrix
  if (mode==1) rnd.RanUniform3d(dr, r); // Gaussian Transition Matrix

	if (state==0) alpha = GroundState(r) / GroundState(rold); // Acceptance Rate
	if (state==1) alpha = ExcitedState(r) / ExcitedState(rold);

	attempted++;

	if (rnd.Rannyu() > alpha) { // keep old pos
		for (int j=0; j<3; j++) r[j] = rold[j];
	} else { // keep new position
		accepted++;
	}
}

void Accumulate(int iblock){ // collect positions for averaging
	// cout << DistanceFormula(r) << endl;
	ave_pos[iblock-1] += DistanceFormula(r);
}

void Average(int iblock){ // compute average and print to file
	ave_pos[iblock-1] /= (nsteps/imeasure);
	av2_pos[iblock-1] = ave_pos[iblock-1] * ave_pos[iblock-1];
	if(iblock==m_blocks){
		string datafile = "./stats.dat";
		BlockedStats(datafile, ave_pos, av2_pos, m_blocks);
	}
}

double DistanceFormula(double r[3]){

	return sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
}

void PrintPos(double vec[3]){
	ofstream pos;
	pos.open("output_positions.dat",ios::app);
	pos << vec[0] << " " << vec[1] << " " << vec[2] << endl;
	pos.close();
}

void BlockedStats(string filename, double AV[], double AV2[], int N){
  ofstream outFile;
  outFile.open(filename);

  vector<double> sum(N,0);
  vector<double> sum2(N,0);
  vector<double> err(N,0);

  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum[i] += AV[j]; // cumulative sum of averages
      sum2[i] += AV2[j]; // sum of square averages
    }
    sum[i] /= (i+1); // cumulative average
    sum2[i] /= (i+1); // cumulative square average
    err[i] = Error(sum, sum2, i); // uncertainty
    outFile << sum[i] << " " << sum2[i] << " " << err[i] << endl;
  }
  outFile.close();
}

double Error(vector<double> AV, vector<double> AV2, int n){
  if (n == 0){
    return 0;
  } else {
    return sqrt( (AV2[n] - AV[n] * AV[n]) / n);
  }
}

double TrialState(double mu, double sigma, double x){
	double psi = exp(-((x-mu)*(x-mu))/(2*sigma));
	psi += exp(-((x+mu)*(x+mu))/(2*sigma));
	return psi * psi; // returns P
}
