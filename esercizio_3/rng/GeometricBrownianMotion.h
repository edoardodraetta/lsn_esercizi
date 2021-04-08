#ifndef __GeometricBrownianMotion__
#define __GeometricBrownianMotion__


float GBM(float, float, float, float, Random &);

// Black-Scholes:
float call_option(float, float, float, float, float);
float put_option(float, float, float, float, float);
float N(float);

#endif 