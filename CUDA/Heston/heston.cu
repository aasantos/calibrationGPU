#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "complex.cuh"
#include "constvectors1000.cuh"
//
__device__ const double a = 0.0;
__device__ const double b = 500.0;
__device__ const int n = 1000;
//
//
__device__ double valueheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double t)
{
	double a = kappa*theta;
        double u = -0.5;
        double b = kappa;
        Complex<double> c1 = Complex<double>(-b,rho*sigma*t).power(2.0);
        Complex<double> c2 = Complex<double>(sigma*sigma*t*t,-sigma*sigma*2.0*u*t);
        Complex<double> dd = (c1 + c2).squareRoot();
        Complex<double> c3 = Complex<double>(b,-rho*sigma*t);
        Complex<double> gg = (c3 + dd)/(c3 - dd);
        Complex<double> cc = gg.inverse();
        Complex<double> c4 = (-tau*dd).exponential();
        Complex<double> DD = ((c3 - dd)/sigma/sigma)*((1.0 - c4)/(1.0 - cc*c4));
        Complex<double> GG = (1.0 - cc*c4)/(1.0 - cc);
        Complex<double> CC = Complex<double>(0.0,r*tau*t) + (a/sigma/sigma)*(tau*(c3 - dd) - 2.0*GG.logarithm());
        Complex<double> ff = (CC + v0*DD + Complex<double>(0.0,t*log(S))).exponential();
        Complex<double> num = Complex<double>(0.0,-t*log(X)).exponential()*ff;
        Complex<double> den = Complex<double>(t*t,t);
        return (num/den).real();

}
//
//
__device__ double integrateheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
{
    double width = 0.5*(b-a);
    double mean = 0.5*(a+b);
    double gaussLegendre = 0;
    for(int step = 1; step <= n; step++) {
        gaussLegendre += weight[step]*valueheston(S,X,tau,r,kappa,theta,sigma,rho,v0,width*root[step] + mean);
    }
    return gaussLegendre*width;
}
//
//
class Heston{
        protected:
            //
            //
            double S;
            double X;
            double tau;
            double r;
            //
            //
	    double kappa;
	    double theta;
            double sigma;
            double rho;
            double v0;
            //
	    //
        public:
        __device__ Heston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
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
        __device__ double callprice()
        {
                double prob = integrateheston(this->S,this->X,this->tau,this->r,this->kappa,this->theta,this->sigma,this->rho,this->v0);
                return this->S - exp(-this->r*this->tau)*this->X*(0.5 + prob/M_PI);
        }
        //
        //
        __device__ double putprice()
        {
                double prob = integrateheston(this->S,this->X,this->tau,this->r,this->kappa,this->theta,this->sigma,this->rho,this->v0);
                return  exp(-this->r*this->tau)*this->X*(0.5 - prob/M_PI);

        }
        //
        //
};
//
__device__ double callheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
{
        Heston *heston = new Heston(S,X,tau,r,kappa,theta,sigma,rho,v0);
        double price = heston->callprice();
        delete heston;
        return price;
}

__device__ double putheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
{
        Heston *heston = new Heston(S,X,tau,r,kappa,theta,sigma,rho,v0);
        double price = heston->putprice();
        delete heston;
        return price;
}


__global__ void kernelheston(double *price,double *stockprice,double *strike,double *tau,int *type,int n,double r,double kappa,double theta,double sigma,double rho,double v0)
{
        int idx = blockDim.x*blockIdx.x + threadIdx.x;
        if(idx < n){
                if(type[idx] == 0){
                    price[idx] = callheston(stockprice[idx],strike[idx],tau[idx],r,kappa,theta,sigma,rho,v0);
                }else{
                    price[idx] = putheston(stockprice[idx],strike[idx],tau[idx],r,kappa,theta,sigma,rho,v0);
                }
        }
}
