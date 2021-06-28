b#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "random.h"
#include "functions.h"

void Discrete_Random_Walk(vector<int> & pos, Random & rnd, int N){
   for (int i = 0; i < N; i++){
      Discrete_Random_Step(pos, rnd);
   }
}

void Discrete_Random_Step(vector<int> & pos, Random & rnd){

  double r;
  r = rnd.Rannyu(); // r ~ U[0,1)

   if (r < 1./6) pos[0] += 1;
   else if (r < 2./6) pos[0] -= 1;
   else if (r < 3./6) pos[1] += 1;
   else if (r < 4./6) pos[1] -= 1;
   else if (r < 5./6) pos[2] += 1;
   else pos[2] -= 1;
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

void blocked_stats(vector<double>  AV, vector<double>  AV2, int N, string filename){
   ofstream outfile;
   outfile.open(filename);

   vector<double> sum(N,0);
   vector<double> sum2(N,0);
   vector<double> err(N,0);

   for (int i=0; i<N; i++){
      for (int j=0; j<(i+1); j++){
         sum[i] += AV[j]; // cumulative sum of averages
         sum2[i] += AV2[j]; // sum of square averages
      }
      sum[i] /= (i+1); // cumulative average
      sum2[i] /= (i+1); // cumulative square average
      err[i] = error(sum, sum2, i); // uncertainty
      outfile << sum[i] << " " << sum2[i] << " " << err[i] << endl;
   }
   outfile.close();
}

float error(vector<double> AV, vector<double> AV2, int n){
   if (n == 0){
      return 0;
   } else {
      return sqrt( (AV2[n] - AV[n]*AV[n]) / n);
   }
}

