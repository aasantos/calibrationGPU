#ifndef integrator_cuh
#define integrator_cuh
//
//
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <math.h>
#include <stdlib.h>
#include "legendrepolynomial.cuh"
//
//
class Function{
public:
	__host__ __device__ Function()
	{
	}

	__host__ __device__ double value(double x)
	{
		return x*x;
	}

};

template <class Obj>
class Integrator{
protected:
	//
	Obj *obj;
	//
	double a;
        double b;
        int n;
        double width;
        double mean;
        //
        LegendrePolynomial *lp;
        double *mWeight;
        double *mRoot;
	//
	//
public:
	__host__ __device__ Integrator(Obj *obj,double a,double b,int n)
	{
	    this->obj = obj;
	    this->a = a;
	    this->b = b;
	    this->n = n;
	    this->width = 0.5*(b-a);
	    this->mean = 0.5*(a+b);
	    this->lp = new LegendrePolynomial(a,b,n);
	    this->mWeight = this->lp->getWeight();
	    this->mRoot = this->lp->getRoot();
	}

	__host__ __device__ ~Integrator()
	{
		if(lp) delete lp;
	}

	__host__ __device__ void setObj(Obj *obj)
	{
		if(obj) delete obj;
		this->obj = obj;
	}

	__host__ __device__ double value()
	{
	double gaussLegendre = 0;
        for(int step = 1; step <= n; step++) {
            gaussLegendre += mWeight[step]*this->obj->value(this->width * mRoot[step] + this->mean);
        }
        //
        return  gaussLegendre * this->width;
	}

};


#endif
