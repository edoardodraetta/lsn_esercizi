#include <iostream>
#include "MolDyn.h"

using namespace std;

int main(){

	// MolDyn MD("./input/Ar.solid");
	// MD.Simulate();
	// MD.PrintStats("../../data/solidphase");

	MolDyn MD("./input/Ar.solid");
	MD.Simulate();
	MD.PrintStats("../../data/4.3.1");
	return 0;

}
