#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <vector>

#include "MolDyn.h" 		// original program -> MolDyn class
#include "statistics.h"

using namespace std;

int main() {

	MolDyn MD;
	MD.Input();

	int nconf = 1;
	int ns = MD.nstep;
	int ip = MD.iprint;

	const int N = MD.nblocks;
	int iblock = 0;
	int L = ns/N; // measured steps per block

	vector<double> ave_epot(N,0), ave_ekin(N,0), ave_etot(N,0), ave_temp(N,0);
	vector<double> av2_epot(N,0), av2_ekin(N,0), av2_etot(N,0), av2_temp(N,0);

	cout << ns << " time steps, " << N << " blocks." << endl;
	cout << "Measure at every step, " <<  ns << " measurements." << endl;
	cout << ns/(N) << " measurements per block." << endl << endl;

	int count = 0;
	for(int istep=1; istep <= ns; ++istep){

		MD.Move(); // Verlet

		// Report progress
		if(istep%ip == 0) {
			cout << "Number of time-steps: " << istep << endl;
		}

		// Measure observables and write current config.xyz
	  if(istep%10 == 0){
      // MD.Measure();
      // ConfXYZ(nconf);
      nconf += 1;
	  }

	  // Print penultimate configuration
	  if(istep == ns-1) {
	  	MD.ConfOut("old.0"); // "config.penult"
		}

  	// Measurement for blocked Statistics

		MD.Measure();
		ave_etot[iblock] += MD.stima_etot;
		ave_epot[iblock] += MD.stima_pot;
		ave_ekin[iblock] += MD.stima_kin;
		ave_temp[iblock] += MD.stima_temp;
		count += 1;
		if(istep%L ==0) iblock += 1;
	}

	// Compute averages and square averages
	for (int i = 0; i < N; i++){
		ave_etot[i] /= L;
		av2_etot[i] = ave_etot[i] * ave_etot[i];

		ave_epot[i] /= L;
		av2_epot[i] = ave_epot[i] * ave_epot[i];

		ave_ekin[i] /= L;
		av2_ekin[i] = ave_ekin[i] * ave_ekin[i];

		ave_temp[i] /= L;
		av2_temp[i] = ave_temp[i] * ave_temp[i];
	}

 	// Blocked Statistics
 	string datafile;

 	datafile  = "../../data/ave_etot.out";
 	blocked_stats(ave_etot, av2_etot, N, datafile);

	datafile = "../../data/ave_epot.out";
 	blocked_stats(ave_epot, av2_epot, N, datafile);

 	datafile = "../../data/ave_ekin.out";
 	blocked_stats(ave_ekin, av2_ekin, N, datafile);

 	datafile = "../../data/ave_temp.out";
 	blocked_stats(ave_temp, av2_temp, N, datafile);

 	// End
 	MD.ConfFinal();
 	cout << "Printed penultimate configuration to file old.0" << endl;

}
