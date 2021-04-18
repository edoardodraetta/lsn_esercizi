#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <vector>
#include "MolDyn.h" // original program -> MolDyn class
#include "statistics.h"

using namespace std;

int main() {

	MolDyn MD;
	MD.Input();

	int nconf = 1;
	int ns = MD.nstep;
	int ip = MD.iprint;

	for(int istep=1; istep <= ns; ++istep){
		MD.Move(); // Verlet

		if(istep%ip == 0) {
			cout << "Number of time-steps: " << istep << endl;
		}

	   if(istep%10 == 0){
	      MD.Measure(); // Properties measurement
	      // ConfXYZ(nconf); // Write configuration in XYZ format
	      nconf += 1;
	   }

	   if(istep == ns-1) MD.ConfOut("old.0");
	}

	cout << endl;
 	MD.ConfFinal();
 	cout << "Also, printed final-1 configuration to file old.0" << endl;



}
