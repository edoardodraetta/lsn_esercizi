#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"


float GBM(float S, float T, float mu, float sigma , Random & rnd){
	return S*exp( (mu-0.5 * sigma*sigma)*T + sigma*rnd.Gauss(0,1)*sqrt(T));
}

// Black-Scholes:

float N(float x){
	return 0.5 * (1 + erf( x/sqrt(2) ) );
}

float call_option(float S, float T, float K, float sigma, float r){

	float d1 = ( log(S/K) + (r + 0.5*sigma*sigma*T) ) / (sigma*sqrt(T));
	float d2 = d1 - sigma*sqrt(T);
	return S*N(d1) - K*exp(-r*T)*N(d2);
}

float put_option(float S, float T, float K, float sigma, float r){
	float d1 = ( log(S/K) + (r + 0.5*sigma*sigma*T) ) / (sigma*sqrt(T));
	float d2 = d1 - sigma*sqrt(T);
	return S*(N(d1)-1) - K*exp(-r*T)*(N(d2)-1);
}
