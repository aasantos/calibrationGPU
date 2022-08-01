#ifndef legendrepolynomial_cuh
#define legendrepolynomial_cuh

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <math.h>
#include <stdlib.h>


class LegendrePolynomial{
protected:
    const double EPSILON = 1e-15;

    struct Result {
        double value;
        double derivative;

        __host__ __device__ Result() : value(0), derivative(0) {}
        __host__ __device__ Result(double val, double deriv) : value(val), derivative(deriv) {}
    };

    double mLowerBound;
    double mUpperBound;
    size_t mNumberOfIterations;
    double *mWeight;
    double *mRoot;

public:
    //
    __host__ __device__ LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations)
    {
        this->mLowerBound = lowerBound;
        this->mUpperBound = upperBound;
        this->mNumberOfIterations = numberOfIterations;
        this->mWeight = new double[numberOfIterations + 1];
        this->mRoot = new double[numberOfIterations + 1];
        //
	calculateWeightAndRoot();
    }

    __host__ __device__ ~LegendrePolynomial()
    {
	    if(mWeight) delete[] mWeight;
	    if(mRoot)   delete[] mRoot;
    }

    __host__ __device__ double* getWeight() {
        return this->mWeight;
    }

    __host__ __device__ double* getRoot(){
        return this->mRoot;
    }

        __host__ __device__ void calculateWeightAndRoot() {
        for(int step = 0; step <= mNumberOfIterations; step++) {
            double root = cos(M_PI * (step-0.25)/(mNumberOfIterations+0.5));
            Result result = calculatePolynomialValueAndDerivative(root);

            double newtonRaphsonRatio;
            do {
                newtonRaphsonRatio = result.value/result.derivative;
                root -= newtonRaphsonRatio;
                result = calculatePolynomialValueAndDerivative(root);
            } while (fabs(newtonRaphsonRatio) > EPSILON);

            this->mRoot[step] = root;
            this->mWeight[step] = 2.0/((1-root*root)*result.derivative*result.derivative);
        }
    }

    __host__ __device__ Result calculatePolynomialValueAndDerivative(double x) {

        Result result(x, 0);
        double value_minus_1 = 1;
        const double f = 1/(x*x-1);
        for(int step = 2; step <= mNumberOfIterations; step++) {
            const double value = ((2*step-1)*x*result.value-(step-1)*value_minus_1)/step;
            result.derivative = step*f*(x*value - result.value);

            value_minus_1 = result.value;
            result.value = value;
        }
        return result;
    }


};

#endif
