#ifndef __GeometricBrownianMotion__
#define __GeometricBrownianMotion__


double GBM(double, double, double, double, Random &);

// Black-Scholes:
double call_option(double, double, double, double, double);
double put_option(double, double, double, double, double);
double N(double);

#endif 