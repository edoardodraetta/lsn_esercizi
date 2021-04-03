#ifndef __functions__
#define __functions__

using namespace std;

void blocked_stats(vector<double>, vector<double>, int, string);
float error(vector<double>, vector<double>, int);
void Discrete_Random_Step(vector<int> &, Random &);
void Discrete_Random_Walk(vector<int> &, Random &, int);
void Continuous_Random_Step(vector<double> &, Random &);
void Continuous_Random_Walk(vector<double> &, Random &, int);
double Distance_Formula_Lattice(vector<int>);
double Distance_Formula(vector<double>);

#endif
