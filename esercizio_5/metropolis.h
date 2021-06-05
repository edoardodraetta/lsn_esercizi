
#include "random.h"

// Positions
double r0[3];
double rold[3], r[3];

// Metropolis
double alpha; // acceptance probability
double dr; // size of step
int nsteps, nblocks;
int accepted, attempted;
int iprint, imeasure;

// Params
bool state; // ground or excited state
bool mode;

// Observables
double ave_pos, av2_pos;
double glob_ave, glob_av2, err;
double blk_norm;

// RNG
Random rnd;
int seed[4];
int p1, p2;

void Initialize();
void Welcome();
void Move();
void Accumulate(int);
void Average(int);
void Reset(int);
void Report(int);
void PrintPos(double []);

double GroundState(double []);
double ExcitedState(double []);
double DistanceFormula(double []);
double Error(double, double, int);
