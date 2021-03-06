/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Monte_Carlo_ISING_1D.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>

using namespace std;

int main() {
  Input();                                 // Inizialization
  for (int iblk = 1; iblk <= nblk; ++iblk) // Simulation
  {
    Reset(iblk); // Reset block averages
    for (int istep = 1; istep <= nstep; ++istep) {
      Move(metro);
      Measure();
      Accumulate(); // Update block averages
    }
    Averages(iblk); // Print results for current block
  }
  ConfFinal(); // Write final configuration

  return 0;
}

void Input(void) // Initialize
{
  ifstream ReadInput;

  // Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();

  // Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0 / temp;
  ReadInput >> nspin;
  ReadInput >> J;
  ReadInput >> h;
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> restart;
  ReadInput.close();

  // Prepare arrays for measurements
  iu = 0; // Energy
  ic = 1; // Heat capacity
  im = 2; // Magnetization
  ix = 3; // Magnetic susceptibility

  n_props = 4; // Number of observables

  // initial configuration
  Restart();

  // Evaluate energy etc. of the initial configuration
  Measure();
  Welcome();
}

void Welcome() // Welcome message
{
  cout << "Classical 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  cout << "Temperature = " << temp << endl;
  cout << "Number of spins = " << nspin << endl;
  cout << "Exchange interaction = " << J << endl;
  cout << "External field = " << h << endl << endl;

  if (metro == 1)
    cout << "The program perform Metropolis moves" << endl;
  else
    cout << "The program perform Gibbs moves" << endl;

  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;

  if (restart == 1)
    cout << "Restarting from previous configuration." << endl;
  else
    cout << "Generating infinite temperature configuration." << endl;

  // Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu] / (double)nspin << endl;
  cout << "Initial specific heat = " << walker[ic] << endl;
  cout << "Initial susceptibility = " << walker[ix] << endl;
  cout << "Initial magnetization = " << walker[im] << endl;
  // cout <<  ... << endl;
  cout << endl;
}

void Restart() // Set initial configuration
{
  if (restart == 0) {
    for (int i = 0; i < nspin; ++i) {
      if (rnd.Rannyu() >= 0.5)
        s[i] = 1;
      else
        s[i] = -1;
    }
  } else {
    ifstream confin;
    confin.open("config.final");
    for (int i = 0; i < nspin; ++i) {
      confin >> s[i];
    }
    confin.close();
  }
}

void PrintState(int a) // debugging
{
  cout << "Current state: " << endl;
  for (int i = 0; i < nspin; ++i) {
    cout << s[i];
    if (i == a)
      cout << " <--";
    cout << endl;
  }
}

void Move(int metro) {
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for (int i = 0; i < nspin; ++i) {
    o = (int)(rnd.Rannyu() * nspin); // Select randomly a particle
    if (metro == 1)                  // Metropolis
    {
      ++attempted;
      if (s[o] == 1)
        sm = -1;
      else
        sm = 1;
      energy_old = Boltzmann(s[o], o);
      energy_new = Boltzmann(sm, o);
      if (energy_new < energy_old) {
        s[o] = sm;
        ++accepted;
      } else {
        p = exp(-beta * (energy_new - energy_old));
        if (p >= rnd.Rannyu()) {
          s[o] = sm;
          ++accepted;
        }
      }
    } else if (metro == 0) // Gibbs sampling
    {
      energy_up = Boltzmann(1, i);
      energy_down = Boltzmann(-1, i);
      p = 1 / (1 + exp(-beta * (energy_up - energy_down)));
      if (p >= rnd.Rannyu())
        s[i] = -1;
      else
        s[i] = +1;
    }
  }
}

double Boltzmann(int sm, int ip) {
  double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
  return ene;
}

void Measure() {
  int bin;
  double u = 0.0, c = 0.0, m = 0.0, x = 0.0;

  // cycle over spins
  for (int i = 0; i < nspin; ++i) { // collect energy and magnetiz.
    u += -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
    m += s[i];
  }
  x = m * m; // square magnetization
  c = u * u; // square energy
  walker[iu] = u;
  walker[ic] = c;
  walker[im] = m;
  walker[ix] = x;
}

void Reset(int iblk) // Reset block averages
{
  if (iblk == 1) {
    for (int i = 0; i < n_props; ++i) {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for (int i = 0; i < n_props; ++i) {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void) // Update block averages
{
  for (int i = 0; i < n_props; ++i) {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) // Print results for current block
{
  ofstream Ene, Heat, Mag, Chi;
  const int wd = 14; // delimeter

  cout << "Block number " << iblk << endl;

  cout << "Acceptance rate ";
  if (metro == 1)
    cout << accepted / attempted << endl;
  else
    cout << "100 percent (Gibbs)" << endl;

  // Internal energy
  Ene.open("output.ene.0", ios::app);
  stima_u = blk_av[iu] / blk_norm / (double)nspin; // Energy

  glob_av[iu] += stima_u;
  glob_av2[iu] += stima_u * stima_u;
  err_u = Error(glob_av[iu], glob_av2[iu], iblk);
  Ene << iblk << setw(wd) << stima_u << setw(wd);
  Ene << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
  Ene.close();

  // Specific Heat
  Heat.open("output.heat.0", ios::app);
  stima_c = blk_av[ic] / blk_norm;
  stima_c -= blk_av[iu] * blk_av[iu] / blk_norm / blk_norm;
  stima_c *= beta * beta;
  stima_c /= (double)nspin;

  glob_av[ic] += stima_c;
  glob_av2[ic] += stima_c * stima_c;
  err_c = Error(glob_av[ic], glob_av2[ic], iblk);
  Heat << iblk << setw(wd) << stima_c << setw(wd);
  Heat << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
  Heat.close();

  // Magnetization
  Mag.open("output.mag.0", ios::app);
  stima_m = blk_av[im] / blk_norm / (double)nspin;

  glob_av[im] += stima_m;
  glob_av2[im] += stima_m * stima_m;
  err_m = Error(glob_av[im], glob_av2[im], iblk);
  Mag << iblk << setw(wd) << stima_m << setw(wd);
  Mag << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
  Mag.close();

  // Susceptibility
  Chi.open("output.chi.0", ios::app);
  stima_x = beta * blk_av[ix] / blk_norm / (double)nspin;

  glob_av[ix] += stima_x;
  glob_av2[ix] += stima_x * stima_x;
  err_x = Error(glob_av[ix], glob_av2[ix], iblk);
  Chi << iblk << setw(wd) << stima_x << setw(wd);
  Chi << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
  Chi.close();

  cout << "T : " << temp << endl;
  cout << "E/N : " << stima_u << endl;
  cout << "C/N : " << stima_c << endl;
  cout << "M/N : " << stima_m << endl;
  cout << "X/N : " << stima_x << endl;
  cout << "----------------------------" << endl;
}

void ConfFinal(void) {
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i = 0; i < nspin; ++i) {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i) // Algorithm for periodic boundary conditions
{
  if (i >= nspin)
    i = i - nspin;
  else if (i < 0)
    i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk) {
  if (iblk == 1)
    return 0.0;
  else
    return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) /
                (double)(iblk - 1));
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
