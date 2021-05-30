#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"


double GBM(double S, double T, double mu, double sigma , Random & rnd){
	return S*exp( (mu-0.5 * sigma*sigma)*T + sigma*rnd.Gauss(0,1)*sqrt(T));
}

// Black-Scholes:

double N(double x){
	return 0.5 * (1 + erf( x/sqrt(2) ) );
}

double call_option(double S, double T, double K, double sigma, double r){

	double d1 = ( log(S/K) + (r + 0.5*sigma*sigma*T) ) / (sigma*sqrt(T));
	double d2 = d1 - sigma*sqrt(T);
	return S*N(d1) - K*exp(-r*T)*N(d2);
}

double put_option(double S, double T, double K, double sigma, double r){
	double d1 = ( log(S/K) + (r + 0.5*sigma*sigma*T) ) / (sigma*sqrt(T));
	double d2 = d1 - sigma*sqrt(T);
	return S*(N(d1)-1) - K*exp(-r*T)*(N(d2)-1);
}
