#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "statistics.h"

void blocked_stats(double AV[], double AV2[], int N, string filename){
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

double error(vector<double> AV, vector<double> AV2, int n){
   if (n == 0){
      return 0;
   } else {
      return sqrt( (AV2[n] - AV[n]*AV[n]) / n);
   }
}

