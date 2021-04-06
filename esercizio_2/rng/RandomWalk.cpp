#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "random.h"
#include "RandomWalk.h"

void Discrete_Random_Walk(vector<int> & pos, Random & rnd, int N){

   for (int i = 0; i < N; i++){
      Discrete_Random_Step(pos, rnd);
   }
   
}

void Discrete_Random_Step(vector<int> & pos, Random & rnd){

  double r;
  r = rnd.Rannyu(); // r ~ U[0,1)

   if (r < 1./6){
      pos[0] += 1;

   } else if (1./6 <= r && r < 2./6) {
      pos[0] -= 1;

   } else if (2./6 <= r && r < 3./6) {
      pos[1] += 1;

   } else if (3./6 <= r && r < 4./6){
      pos[1] -= 1;

   } else if (4./6 <= r && r < 5./6){
      pos[2] += 1;

   } else {
      pos[2] -= 1;
   }
}

void Continuous_Random_Walk(vector<double> & pos, Random & rnd, int N){

   for (int i = 0; i < N; i++){
      Continuous_Random_Step(pos, rnd);
   }

}

void Continuous_Random_Step(vector<double> & pos, Random & rnd){

   double theta = rnd.RanTheta3d();
   double phi = rnd.Rannyu(0,2*M_PI);

   pos[0] += sin(theta) * cos(phi);
   pos[1] += sin(theta) * cos(phi);
   pos[2] += cos(theta);

}

double Distance_Formula_Lattice(vector<int> pos){
   float R = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
   return R;
}

double Distance_Formula(vector<double> pos){
   float R = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
   return R;
}
