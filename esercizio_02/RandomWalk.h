#ifndef __RandomWalk__
#define __RandomWalk__

using namespace std;

void Discrete_Random_Step(vector<int> &, Random &);
void Discrete_Random_Walk(vector<int> &, Random &, int);
void Continuous_Random_Step(vector<double> &, Random &);
void Continuous_Random_Walk(vector<double> &, Random &, int);
double Distance_Formula_Lattice(vector<int>);
double Distance_Formula(vector<double>);

#endif
