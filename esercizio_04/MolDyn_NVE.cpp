/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>      // std::setprecision
#include <string>
#include <vector>
#include "MolDyn_NVE.h"

using namespace std;

int main() {
  Input(); // Inizialization
  int nconf = 1;
  for (int iblock = 1; iblock <= nblocks; ++iblock) {
    Reset(iblock);
    for (int istep = 1; istep <= nstep; ++istep) {
      Move(); // Move particles with Verlet algorithm
      if (istep % iprint == 0) { // Write actual configuration in XYZ format:
        // ConfXYZ(nconf); // Commented to avoid "filesystem full"!
        nconf += 1;
      }

      Measure();    // properties
      Accumulate(); // accumulate to compute averages

      if ((istep == nstep - 1) && (iblock == nblocks))
        ConfPenult();
    }
    Averages(iblock); // compute and output
  }
  ConfFinal(); // Write final configuration to restart
  return 0;
}


void Input() { // Prepare all stuff for the simulation
  ifstream ReadInput, ReadConf;
  // double ep, ek, pr, et, vir;

  seed = 1;    // Set seed for random numbers
  srand(seed); // Initialize random number generator

  ReadInput.open("input.dat"); // Read input

  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblocks;
  ReadInput >> iprint;
  ReadInput >> restart;
  ReadInput >> rescale;
  ReadInput >> rescaletemp;

  vol = (double)npart / rho;
  box = pow(vol, 1.0 / 3.0);

  ReadInput.close();
  Welcome(); // Welcome message

  // Prepare array for measurements
  iv = 0;      // Potential energy
  ik = 1;      // Kinetic energy
  ie = 2;      // Total energy
  it = 3;      // Temperature
  iw = 4;      // Pressure (virial)
  n_props = 5; // Number of observables

  // measurement of g(r): da 0 a L/2
  igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box / 2.0) / (double)nbins;

  // Restart or rescale?
  if (restart == 1) {
    Restart();
    if (rescale == 1) {
      Rescale();
    }
  } else {
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i = 0; i < npart; ++i) {
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
    Prepare();
  }
}

void Reset(int iblock) { // zero averages
  if (iblock == 1) {
    for (int i = 0; i < n_props; ++i) {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }
  for (int i = 0; i < n_props; ++i) {
    blk_av[i] = 0;
  }
  blk_norm = 0;
}

void Move(void) { // Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for (int i = 0; i < npart; ++i) { // Force acting on particle i
    fx[i] = Force(i, 0);
    fy[i] = Force(i, 1);
    fz[i] = Force(i, 2);
  }

  for (int i = 0; i < npart; ++i) { // Verlet integration scheme

    xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
    ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
    znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

    vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
    vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
    vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

void Measure() { // Properties measurement
  double bin;
  double v, w, t;
  double vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("output_epot.dat", ios::app);
  Ekin.open("output_ekin.dat", ios::app);
  Temp.open("output_temp.dat", ios::app);
  Press.open("output_press.dat", ios::app);
  Etot.open("output_etot.dat", ios::app);

  v = 0.0; // reset observables
  w = 0.0;
  t = 0.0;

  // reset the hystogram of g(r)
  for (int k = igofr; k < igofr + nbins; ++k)
    walker[k] = 0.0;

  // cycle over pairs of particles
  for (int i = 0; i < npart - 1; ++i) {
    for (int j = i + 1; j < npart; ++j) {

      dx = Pbc(xold[i] - xold[j]); // here I use old configurations [old = r(t)]
      dy = Pbc(yold[i] - yold[j]); // to be compatible with EKin which uses v(t)
      dz = Pbc(zold[i] - zold[j]); // => EPot should be computed with r(t)

      dr = dx * dx + dy * dy + dz * dz;
      dr = sqrt(dr);

      // update of the histogram of g(r)
      for (int ibin = 0; ibin < nbins; ++ibin) {
        bin = bin_size * (ibin + 1.0);
        if (dr < bin) {
          walker[igofr + ibin] += 2;
          break;
        }
      }

      if (dr < rcut) {
        vij = 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
        wij = 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6);
        v += vij; // Potential energy
        w += wij; // virial
      }
    }
  }

  // Kinetic energy
  for (int i = 0; i < npart; ++i)
    t += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

  // Update walkers
  walker[iv] = 4.0 * v;
  walker[ik] = 0.5 * t;
  walker[ie] = (0.5 * t + 4.0 * v);
  walker[it] = t / 3.0;
  walker[iw] = 48.0 * w / 3.0;

  // Instantaneous Measurements
  stima_pot = 4.0 * v / (double)npart;                // Potential energy per particle
  stima_kin = 0.5 * t / (double)npart;                // Kinetic energy per particle
  stima_temp = (1.0 / 3.0) * t / (double)npart; // Temperature
  stima_press = rho * stima_temp + walker[iw] / vol;  // Pressure
  stima_etot = (0.5 * t + 4.0 * v) / (double)npart;   // Total energy per particle

  Epot << stima_pot << endl;
  Ekin << stima_kin << endl;
  Temp << stima_temp << endl;
  Press << stima_press << endl;
  Etot << stima_etot << endl;

  current_temp = stima_temp;

  Epot.close();
  Ekin.close();
  Temp.close();
  Press.close();
  Etot.close();

  return;
}

void Accumulate(){
  for (int i = 0; i < n_props; ++i) {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblock) { // compute averages
  ofstream Pot, Kin, Etot, Temp, Pres, Gofr, Gave;
  double r, dr, gdir, dvol;
  int pr = 9;

  cout << "Block number : " << iblock << endl;

  Pot.open("./output.ave_epot.dat", ios::app);
  Kin.open("./output.ave_ekin.dat", ios::app);
  Etot.open("./output.ave_etot.dat", ios::app);
  Temp.open("./output.ave_temp.dat", ios::app);
  Pres.open("./output.ave_pres.dat", ios::app);
  Gofr.open("./output.ave_gofr.dat", ios::app);

  // Potential Energy
  stima_pot = blk_av[iv] / blk_norm / (double)npart;
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot * stima_pot;
  err_pot = Error(glob_av[iv], glob_av2[iv], iblock);
  // Kinetic Energy
  stima_kin = blk_av[ik] / blk_norm / (double)npart;
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin * stima_kin;
  err_kin = Error(glob_av[ik], glob_av2[ik], iblock);
  // Total Energy
  stima_etot = blk_av[ie] / blk_norm / (double)npart;;
  glob_av[ie] += stima_etot;
  glob_av2[ie] += stima_etot * stima_etot;
  err_etot = Error(glob_av[ie], glob_av2[ie], iblock);
  // Temperature
  stima_temp = blk_av[it] / blk_norm / (double)npart;;
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp * stima_temp;
  err_temp = Error(glob_av[it], glob_av2[it], iblock);
  // Pressure
  stima_press = rho * stima_temp + (blk_av[iw] / blk_norm) / vol;
  glob_av[iw] += stima_press;
  glob_av2[iw] += stima_press * stima_press;
  err_press = Error(glob_av[iw], glob_av2[iw], iblock);

  // Output for the above
  Pot << iblock << " " << setprecision(pr) << glob_av[iv]/(double)iblock;
  Pot << " " << setprecision(pr) << err_pot << endl;
  Kin << iblock << " " << setprecision(pr) << glob_av[ik]/(double)iblock;
  Kin << " " << setprecision(pr) << err_kin << endl;
  Etot << iblock << " " << setprecision(pr) << glob_av[ie]/(double)iblock;
  Etot << " " << setprecision(pr) << err_etot << endl;
  Temp << iblock << " " << setprecision(pr) << glob_av[it]/(double)iblock;
  Temp << " " << setprecision(pr) << err_temp  << endl;
  Pres << iblock << " " << setprecision(pr) << glob_av[iw]/(double)iblock;
  Pres << " " << setprecision(pr) << err_press << endl;

  // g(r)
  dr = bin_size;
  for (int ibin = 0; ibin < nbins; ++ibin) {

    // Delta Vol = (4 pi / 3) [(r+dr)^3 - r^3 ]
    r = dr * ibin;
    dvol = 4 * M_PI * (pow(r + dr, 3) - pow(r, 3)) / 3;

    gdir = (blk_av[igofr + ibin] / blk_norm) / (double)npart / rho / dvol;
    glob_av[igofr + ibin] += gdir;
    glob_av2[igofr + ibin] += gdir * gdir;
    err_gdir = Error(glob_av[igofr + ibin], glob_av2[igofr + ibin], iblock);

    Gofr << iblock << " " << ibin + 1 << " " << r << " ";
    Gofr << gdir << " " << glob_av[igofr + ibin] / (double)iblock << " "
         << err_gdir << " " << endl;

    if (iblock == nblocks){ // Final g(r)
      Gave.open("./output.gave.dat", ios::app);
      Gave << ibin << " " << r << " ";
      Gave << glob_av[igofr + ibin] / (double)iblock << " ";
      Gave << err_gdir << " " << endl;
      Gave.close();
    }
  }

  Pot.close();
  Kin.close();
  Etot.close();
  Pres.close();
  Gofr.close();

  cout << "T(t): " << current_temp << endl;
  cout << "<T>: " << stima_temp << endl;
  cout << "<KE> : " << stima_kin << endl;
  cout << "<PE> : " << stima_pot << endl;
  cout << "<E> : " << stima_etot << endl;
  cout << "<P> : " << stima_press << endl;
  cout << "----------------------------" << endl;
}

void ConfPenult(void) { // Write final configuration
  ofstream WriteConf;

  WriteConf.open("old.penult");
  for (int i = 0; i < npart; ++i) {
    WriteConf << x[i] / box << "   " << y[i] / box << "   ";
    WriteConf << z[i] / box << endl;
  }
  WriteConf.close();
}

void ConfFinal(void) { // Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file old.final " << endl << endl;

  WriteConf.open("old.final");
  for (int i = 0; i < npart; ++i) {
    WriteConf << x[i] / box << "   " << y[i] / box << "   ";
    WriteConf << z[i] / box << endl;
  }
  WriteConf.close();

  return;
}

/*=========================EQUILIBRAZIONE=====================================*/

void Restart(void) {
  cout << "Restarting from previous configuration. " << endl << endl;
  ifstream ReadConf;
  ReadConf.open("old.final"); // final config
  for (int i = 0; i < npart; ++i) {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box; // LJ reduced units
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  ReadConf.open("old.penult"); // penultimate config
  for (int i = 0; i < npart; ++i) {
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box; // LJ reduced units
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();
}

void Rescale(void) {
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
  double sumv2 = 0, fs;
  cout << "Rescale temperature for quick equilibration." << endl << endl;

  for (int i = 0; i < npart; ++i) { // Force acting on particle i
    fx[i] = Force(i, 0);
    fy[i] = Force(i, 1);
    fz[i] = Force(i, 2);
  }

  for (int i = 0; i < npart; ++i) { // Verlet integration scheme for r(t+dt)

    // 1. Compute r(t+dt) with Verlet algorithm
    xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
    ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
    znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

    // 2. calculate velocity at v(t + dt/2)
    vx[i] = Pbc(xnew - x[i]) / (delta);
    vy[i] = Pbc(ynew - y[i]) / (delta);
    vz[i] = Pbc(znew - z[i]) / (delta);

    // Calculate average square velocity
    sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

    // Keep
    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }

  sumv2 /= (double)npart;
  if (rescaletemp==0){
    fs = sqrt(3 * temp / sumv2); // T2/T1
  } else {
    fs = sqrt(temp / rescaletemp); // T2/T1
  }

  cout << "Desired temp: " << temp << endl;
  cout << "Computed temp: " << sumv2 / 3. << endl;
  cout << "Rescaling factor: " << fs << endl << endl;

  for (int i = 0; i < npart; ++i) {
    // 3. Rescale velocities according to desired temperature
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
    // 4. Estimate an OLD spatial config, r(t-dt)
    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
}

void Prepare(void) {
  cout << "Prepare random velocities with c.o.m. velocity equal to zero";
  cout << endl << endl;

  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < npart; ++i) {
    vx[i] = rand() / double(RAND_MAX) - 0.5;
    vy[i] = rand() / double(RAND_MAX) - 0.5;
    vz[i] = rand() / double(RAND_MAX) - 0.5;
    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }

  for (int idim = 0; idim < 3; ++idim)
    sumv[idim] /= (double)npart;
  double sumv2 = 0.0, fs;
  for (int i = 0; i < npart; ++i) {
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
  }
  sumv2 /= (double)npart;

  fs = sqrt(3 * temp / sumv2); // fs = velocity scale factor
  for (int i = 0; i < npart; ++i) {
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
  return;
}

/*===========================================================================*/

void Welcome() { // prints welcome message
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl;
  cout << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]";
  cout << endl;
  cout << endl;

  cout << "The program uses Lennard-Jones units " << endl;

  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << endl;
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "The program integrates Newton equations with the Verlet method ";
  cout << endl;

  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblocks << endl;
  cout << "Restart = " << restart << endl;
  cout << "Rescale = " << rescale << endl;
  cout << endl;
}

double Force(int ip, int idir) { // Compute forces as -Grad_ip V(r)
  double f = 0.0;
  double dvec[3], dr;

  for (int i = 0; i < npart; ++i) {
    if (i != ip) {
      dvec[0] = Pbc(x[ip] - x[i]); // distance ip-i in pbc
      dvec[1] = Pbc(y[ip] - y[i]);
      dvec[2] = Pbc(z[ip] - z[i]);

      dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
      dr = sqrt(dr);

      if (dr < rcut) { // -Grad_ip V(r)
        f += dvec[idir] *
             (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8));
      }
    }
  }

  return f;
}

void ConfXYZ(int nconf) { // Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i = 0; i < npart; ++i) {
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   ";
    WriteXYZ << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Error(double AV, double AV2, int N) {
  if (N == 1) {
    return 0;
  } else {
    return sqrt((AV2 / (double)N - (AV / (double)N) * (AV / (double)N)) /
                (double)(N - 1));
  }
}

double Pbc(double r) {
  // Algorithm for periodic boundary conditions with side L=box
  return r - box * rint(r / box);
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
