#include <iostream>
#include "MolDyn.h"

using namespace std;

int main(){

	MolDyn MD;
	MD.Simulate();
	MD.PrintStats();
	return 0;

}
