#include "random.h"

double r0[3];
double rold[3], r[3]; // positions
double alpha; // acceptance probability
double dr; // size of step

const int m_blocks = 100;
int nsteps, nblocks;
int accepted, attempted;
double ave_pos[m_blocks], av2_pos[m_blocks];
int iprint, imeasure;
bool state; // ground or excited state
bool mode;

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

double TrialState(double, double, double);
double DistanceFormula(double []);

void BlockedStats(std::string, double [], double [], int);
double Error(std::vector<double>, std::vector<double>, int);
