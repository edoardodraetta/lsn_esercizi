#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "metropolis.h"
#include "random.h" // RNG

using namespace std;

int main() {
  Initialize();  // read from input.dat
  for (int iblock = 1; iblock <= nblocks; ++iblock) {
    Reset(iblock);
    for (int istep = 1; istep <= nsteps; ++istep) {
      Move();  // Try a move

      // Measure (for avg position)
      if (istep % imeasure == 0) Accumulate(iblock);  // measure

      // Print (for plotting of pdf)
      if (istep % iprint == 0) PrintPos(r);
    }
    Average(iblock);  // compute average for current block
    if (iblock % 5 == 0) Report(iblock);
  }
  return 0;
}

// ========== Startup ===============================================

void Initialize() {
  // Initialize RNG
  ifstream Primes("Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);  // Initialize
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;

  // Input File
  ifstream ReadInput("./input.dat");
  ReadInput >> state;
  ReadInput >> mode;
  ReadInput >> dr;
  ReadInput >> nsteps;
  ReadInput >> nblocks;
  ReadInput >> imeasure;
  ReadInput >> iprint;
  ReadInput.close();

  // Starting positions
  if (state == 0) {  // ground state
    r0[0] = 1.2;
    r0[1] = 0;
    r0[2] = 0;
  } else {  // excited state
    r0[0] = 0;
    r0[1] = 0;
    r0[2] = 5;
  }

  for (int i = 0; i < 3; ++i) r[i] = r0[i];

  Welcome();  // print summary to terminal
}

void Welcome() {
  cout << "Sampling Wavefunctions of the Hydrogen Atom with the Metropolis "
          "Algorithm"
       << endl;

  if (state == 0) {
    cout << "Sampling ground state wavefunction.";
  } else {
    cout << "Sampling excited state wavefunction.";
  }
  cout << endl;

  cout << "Walking " << nsteps << " steps in each block." << endl;
  cout << "For " << nblocks << " blocks." << endl;
  cout << "Size of step : " << dr << endl;
  cout << "Recording every " << iprint << " steps" << endl;
  cout << "Measuring every " << imeasure << " steps" << endl;

  if (mode == 0) {
    cout << "Moves attempted with uniform probability.";
  } else {
    cout << "Moves attempted with Gaussian probability.";
  }
  cout << endl;

  cout << "Starting position : ";
  cout << "a_0 * (" << r0[0] << ", " << r0[1] << ", " << r0[2] << ")";

  cout << endl << endl;
}

void Reset(int iblock) {  // Resets observables w/in each block

  if (iblock == 1) {
    glob_ave = 0.0;
    glob_av2 = 0.0;
  }

  accepted = 0;
  attempted = 0;
  ave_pos = 0.0;
  av2_pos = 0.0;
  blk_norm = 0.0;
}

// ========== Algorithm ==========================================

void Move() {  // Metropolis Algorithm

  // 1. Keep old position
  for (int j = 0; j < 3; j++) rold[j] = r[j];

  // 2. Generate a new position
  if (mode == 0) rnd.RanUniform3d(dr, r);  // Uniform Transition Matrix
  if (mode == 1) {                         // Gaussian Transition Matrix
    for (int j = 0; j < 3; j++) r[j] = rnd.Gauss(r[j], dr);
  }

  // 3. Calculate the acceptance probability

  attempted++;

  if (state == 0) alpha = GroundState(r) / GroundState(rold);
  if (state == 1) alpha = ExcitedState(r) / ExcitedState(rold);

  // 4. Accept or reject move
  if (rnd.Rannyu() > alpha) {
    for (int j = 0; j < 3; j++) r[j] = rold[j];
  } else {  // keep new position
    accepted++;
  }
}

double GroundState(double r[3]) {  // phi 100
  double rho = DistanceFormula(r);
  return exp(-2 * rho) / M_PI;  // probability
}

double ExcitedState(double r[3]) {  // phi 210
  double p;
  double rho = DistanceFormula(r);
  double theta = acos(r[2] / rho);
  p = rho * rho * exp(-rho);
  p *= cos(theta) * cos(theta);
  p /= (32 * M_PI);
  return p;  // probability
}

// ========== Measurement ==========================================

void Accumulate(int iblock) {  // collect positions for averaging
  ave_pos += DistanceFormula(r);
  blk_norm += 1.0;
}

double DistanceFormula(double r[3]) {
  return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}

void PrintPos(double vec[3]) {
  ofstream pos;
  pos.open("positions.dat", ios::app);
  pos << vec[0] << " " << vec[1] << " " << vec[2] << endl;
  pos.close();
}

void Average(int iblock) {  // compute average and print to file
  ofstream pos;

  // Compute ave and ave2 within block
  ave_pos /= (blk_norm);
  av2_pos = (ave_pos * ave_pos);

  // Compute cumulative average
  glob_ave += ave_pos;
  glob_av2 += av2_pos;
  err = Error(glob_ave, glob_av2, iblock);

  // Print to file
  pos.open("./average_position.dat", ios::app);
  pos << iblock << " " << glob_ave / (double)iblock << " " << err << endl;
  pos.close();
}

double Error(double sum, double sum2, int iblk) {
  if (iblk == 1)
    return 0.0;
  else {
    return sqrt(
        (sum2 / (double)iblk - (sum / (double)iblk) * (sum / (double)iblk)) /
        (double)(iblk - 1));
  }
}

void Report(int iblock) {
  cout << "Block " << iblock << ", "
       << "Acceptance Rate: ";
  cout << (double)accepted / attempted << endl;
  cout << "Estimate : " << glob_ave / (double)iblock << endl;
  cout << "---------------------------------" << endl;
}

