#include "metropolis.h"

using namespace std;

int main(){

	Initialize(); // read from input.dat

	for (int iblock = 1; iblock <= nblocks; ++iblock){
		for (int istep = 1; istep < nsteps; ++istep){
			Move(); // Try a move
			Accumulate(iblock); // measure
		}
		Average(iblock); // compute average for current block
	}
}
