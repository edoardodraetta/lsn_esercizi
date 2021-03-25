

float error(float* AV, float* AV2, int n){
   if (n == 0){
      return 0;
   } else {
      return sqrt( (AV2[n] - AV[n]*AV[n]) / n);
   }
}

