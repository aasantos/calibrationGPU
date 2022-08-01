#ifndef heston_cuh
#define heston_cuh


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "complex.cuh"
#include "legendrepolynomial.cuh"
#include "integrator.cuh"


class CFH{
protected:
    //
    //
    double S;
    double X;
    double tau;
    double r;
    double kappa;
    double theta;
    double sigma;
    double rho;
    double v0;
    //
    //
public:
    __host__ __device__ CFH(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
    {
        //
        //
        this->S = S;
        this->X = X;
        this->tau = tau;
        this->r = r;
        this->kappa = kappa;
        this->theta = theta;
        this->sigma = sigma;
        this->rho = rho;
        this->v0 = v0;
        //
        //
    }
    //
    //
    __host__ __device__ double value(double t)
    {
        double a = this->kappa*this->theta;
        double u = -0.5;
        double b = this->kappa;
        Complex<double> c1 = Complex<double>(-b,this->rho*this->sigma*t).power(2.0);
        Complex<double> c2 = Complex<double>(this->sigma*this->sigma*t*t,-this->sigma*this->sigma*2.0*u*t);
        Complex<double> dd = (c1 + c2).squareRoot();
        Complex<double> c3 = Complex<double>(b,-this->rho*this->sigma*t);
        Complex<double> gg = (c3 + dd)/(c3 - dd);
        Complex<double> cc = gg.inverse();
        Complex<double> c4 = (-this->tau*dd).exponential();
        Complex<double> DD = ((c3 - dd)/this->sigma/this->sigma)*((1.0 - c4)/(1.0 - cc*c4));
        Complex<double> GG = (1.0 - cc*c4)/(1.0 - cc);
        Complex<double> CC = Complex<double>(0.0,this->r*this->tau*t) + (a/this->sigma/this->sigma)*(this->tau*(c3 - dd) - 2.0*GG.logarithm());
        Complex<double> ff = (CC + this->v0*DD + Complex<double>(0.0,t*log(this->S))).exponential();
        Complex<double> num = Complex<double>(0.0,-t*log(this->X)).exponential()*ff;
        Complex<double> den = Complex<double>(t*t,t);
        return (num/den).real();
    }
    //
    //
};



class Heston{
protected:
    //
    //
    double S;
    double X;
    double tau;
    double r;
    double kappa;
    double theta;
    double sigma;
    double rho;
    double v0;
    CFH *cfh;
    Integrator<CFH> *ii;
    
public:
    //
    //
    __host__ __device__ Heston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
    {
        //
        this->S = S;
        this->X = X;
        this->tau = tau;
        this->r = r;
        this->kappa = kappa;
        this->theta = theta;
        this->sigma = sigma;
        this->rho = rho;
        this->v0 = v0;
        //
        //
        this->cfh = new CFH(S,X,tau,r,kappa,theta,sigma,rho,v0);
        this->ii = new Integrator<CFH>(cfh,0.000001,50.0,100);
    }
    //
    //
    __host__ __device__ ~Heston()
    {
        if(cfh) delete cfh;
        if(ii) delete ii;
    }
    //
    //
    __host__ __device__ double callprice()
    {
        double prob = this->ii->value();
        double p = 0.5 + prob/M_PI;
        return this->S - exp(-this->r*this->tau)*this->X*p;
    }
    //
    //
    __host__ __device__ double putprice()
    {
        double prob = this->ii->value();
        return  exp(-this->r*this->tau)*this->X*(0.5 - prob/M_PI);
        
    }
};




#endif
