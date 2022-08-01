#ifndef blackscholes_cuh
#define blackshcoles_cuh
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuda.h>
//
//
#define sqrt2 1.41421356237
//
double cdfNormal(double x)
{
	return 0.5*erfc(-x/sqrt2);
}
//
//
class BlackScholes{
protected:
	double S;
	double X;
	double tau;
	double r;
	double sigma;
public:
	__host__ __device__ BlackScholes(double S,double X,double tau,double r,double sigma)
	{
		this->S = S;
		this->X = X;
		this->tau = tau;
		this->r = r;
		this->sigma = sigma;
	}
	
	__host__ __device__ double callprice()
	{
		double d1 = (log(this->S/this->X) + (this->r + 0.5*this->sigma*this->sigma)*this->tau)/(this->sigma*sqrt(this->tau));
		double d2 = d1 - this->sigma*sqrt(this->tau);
		return this->S*cdfNormal(d1) - exp(-this->r*this->tau)*this->X*cdfNormal(d2);
	}
};



#endif
