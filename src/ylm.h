#ifndef ylm_h
#define ylm_h ylm_h
#include<iostream>
#include<cmath>
#include<complex>
#include<factorial.h>
#define complex complex<double>

using namespace std;

complex ylm(long l, long m, double theta, double phi);
complex ylm2(long l, long m, double theta, double phi);

#endif // ylm_h





