#include <iostream>
#include <fstream>
#include <string>

#include "statistics.h"
#include "MolDyn.h"

using namespace std;

MolDyn::MolDyn(string input){

	ifstream ReadInput;
	srand(seed); // Initialize random number generator

	// === PARAMS ===
	inputfile = input;
	ReadInput.open(inputfile); // Read input

	ReadInput >> temp;
	ReadInput >> npart;
	ReadInput >> rho;
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> restart; // new
	ReadInput >> quick_rescale; // new

	ReadInput.close();

	// === MESSAGE ===
	Welcome();

	// === INITIALIZE ===
	Initialize();
}

MolDyn::~MolDyn(){}

void MolDyn::Welcome(){
	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	cout << "Simulation parameters from " << inputfile << "" << endl;
	cout << "Number of particles = " << npart << endl;
	cout << "Density of particles = " << rho << endl;
	cout << "Volume of the simulation box = " << vol << endl;
	cout << "Edge of the simulation box = " << box << endl;

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
}

/*=======================Simulation Methods==================================*/

void MolDyn::Simulate(){

	cout << "=== Simulating ===" << endl << endl;

	int nconf = 1;
	int iblock = 0;
	int L = nstep / nblocks;

	for(int istep=1; istep <= nstep; ++istep){
		Move(); // Verlet
		if(istep%iprint == 0) { // Report progress
			cout << "Number of time-steps: " << istep << endl;
		}

		// Measure observables and write current config.xyz
	  if(istep%10 == 0){
      // Measure();
      // ConfXYZ(nconf);
      nconf += 1;
	  }

	  // Measurement for blocked Statistics
	  Measure();
		ave_etot[iblock] += stima_etot;
		ave_epot[iblock] += stima_pot;
		ave_ekin[iblock] += stima_kin;
		ave_temp[iblock] += stima_temp;
		if(istep%L ==0) iblock += 1;

	  // Print penultimate configuration
	  if(istep == nstep-1) {
	  	ConfOut("old.penult");
	  }
	  if(istep == nstep) {
	  	ConfOut("old.final");
		}
	}

	// Compute averages and square averages
	for (int i = 0; i < nblocks; i++){
		ave_etot[i] /= L;
		av2_etot[i] = ave_etot[i] * ave_etot[i];

		ave_epot[i] /= L;
		av2_epot[i] = ave_epot[i] * ave_epot[i];

		ave_ekin[i] /= L;
		av2_ekin[i] = ave_ekin[i] * ave_ekin[i];

		ave_temp[i] /= L;
		av2_temp[i] = ave_temp[i] * ave_temp[i];
	}
}

void MolDyn::PrintStats(string datadir){
	string datafile;

 	datafile = datadir + "/ave_etot.out";
 	blocked_stats(ave_etot, av2_etot, nblocks, datafile);

	datafile = datadir + "/ave_epot.out";
 	blocked_stats(ave_epot, av2_epot, nblocks, datafile);

 	datafile = datadir + "/ave_ekin.out";
 	blocked_stats(ave_ekin, av2_ekin, nblocks, datafile);

 	datafile = datadir + "/ave_temp.out";
 	blocked_stats(ave_temp, av2_temp, nblocks, datafile);
}

/*=======================Initialization Methods==============================*/

void MolDyn::Initialize(){
	// restart = 1 --> restart from previous configuration
	// rescale = 1 --> rescale velocities to match temperature
	if (restart==1){
		Restart();
		if (quick_rescale==1){
			Rescale();
		}
	} else{
		PrepareVelocities();
	}
}

void MolDyn::Restart(){
	cout << "=== Restarting from previous configuration === " << endl << endl;
	ifstream ReadConf;
	ReadConf.open("old.final"); // final config
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box; // LJ reduced units
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();
	ReadConf.open("old.penult"); // penultimate config
	for (int i=0; i<npart; ++i){
	 ReadConf >> xold[i] >> yold[i] >> zold[i];
	 xold[i] = xold[i] * box; // LJ reduced units
	 yold[i] = yold[i] * box;
	 zold[i] = zold[i] * box;
  }
  ReadConf.close();

}

void MolDyn::Rescale(){
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
	double sumv2 = 0, fs;
	cout << "=== Rescaling temperature for quick equilibration === " << endl << endl;

	for(int i=0; i<npart; ++i){ // Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){ // Verlet integration scheme for r(t+dt)

		// 1. Compute r(t+dt) with Verlet algorithm
		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		// 2. calculate velocity at v(t + dt/2)
		vx[i] = Pbc(xnew - x[i])/(delta);
		vy[i] = Pbc(ynew - y[i])/(delta);
		vz[i] = Pbc(znew - z[i])/(delta);

		// Calculate average square velocity
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];

		// Keep
		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}

	sumv2 /= (double) npart;
	fs = sqrt(3*temp / sumv2); // T2/T1

	cout << "RESCALE: Desired temp: " << temp << endl;
	cout << "RESCALE: Computed temp: " << sumv2/3. << endl;
	cout << "RESCALE: Rescaling factor: " << fs << endl;

	for (int i = 0; i < npart; ++i){

		// 3. Rescale velocities according to desired temperature
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;

		// 4. Estimate an OLD spatial config, r(t-dt)
		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}

	sumv2 = 0;
	for (int i = 0; i < npart; i++){
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}
	sumv2 /= (double) npart;
}

void MolDyn::PrepareVelocities(){
	cout << "=== Starting from single configuration file  === " << endl << endl;
	ifstream ReadConf;

  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box; // LJ reduced units
	  y[i] = y[i] * box;
	  z[i] = z[i] * box;
  }
  ReadConf.close();

	cout << "Prepare random velocities with center of mass velocity equal to zero ";
	cout << endl << endl;
	double sumv[3] = {0.0, 0.0, 0.0};
	for (int i=0; i<npart; ++i){
		vx[i] = rand()/double(RAND_MAX) - 0.5;
		vy[i] = rand()/double(RAND_MAX) - 0.5;
		vz[i] = rand()/double(RAND_MAX) - 0.5;
		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}
	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	double sumv2 = 0.0, fs;
	for (int i=0; i<npart; ++i){
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;
  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor T2 / T1
  for (int i=0; i<npart; ++i){
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;
		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}
}

/*==========================Physics Methods==================================*/

void MolDyn::Move(){ // Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){ // Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){ // Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}
void MolDyn::Measure(){ //Properties measurement
	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	ofstream Epot, Ekin, Etot, Temp;

	Epot.open("output_epot.dat",ios::app);
	Ekin.open("output_ekin.dat",ios::app);
	Temp.open("output_temp.dat",ios::app);
	Etot.open("output_etot.dat",ios::app);

	v = 0.0; //reset observables
	t = 0.0;

	// cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
	for (int j=i+1; j<npart; ++j){
		 dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
		 dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
		 dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

		 dr = dx*dx + dy*dy + dz*dz;
		 dr = sqrt(dr);

		 if(dr < rcut){
			vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

			  // Potential energy
			  v += vij;
			}
	}
	}

	// Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

	stima_pot = v/(double)npart; // Potential energy per particle
	stima_kin = t/(double)npart; // Kinetic energy per particle
	stima_temp = (2.0 / 3.0) * t/(double)npart; // Temperature
	stima_etot = (t+v)/(double)npart; // Total energy per particle

	Epot << stima_pot  << endl;
	Ekin << stima_kin  << endl;
	Temp << stima_temp << endl;
	Etot << stima_etot << endl;

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
}

double MolDyn::Force(int ip, int idir){ // Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
		if(i != ip){
		  dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
		  dvec[1] = Pbc( y[ip] - y[i] );
		  dvec[2] = Pbc( z[ip] - z[i] );

		  dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
		  dr = sqrt(dr);

		  if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
		  }
		}
  }
  return f;
}

double MolDyn::Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}

/*=======================Configuration Methods===============================*/

void MolDyn::ConfOut(string filename){ //Write final configuration
  ofstream WriteConf;

  // cout << "Print configuration to file " << filename << endl;
  WriteConf.open(filename);

  for (int i=0; i<npart; ++i){
	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void MolDyn::ConfFinal(){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void MolDyn::ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
	WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}
