#ifndef heston_hpp
#define heston_hpp

#include <iostream>
#include <vector>
#include <stdio.h>
#include "complex.hpp"
#include "const_vector.hpp"
//

using namespace std;

const double a = 0.0;
const double b = 500.0;
const int n = 1000;


double valueheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double t)
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
double integrateheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
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
        Heston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
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
        double callprice()
        {
                double prob = integrateheston(this->S,this->X,this->tau,this->r,this->kappa,this->theta,this->sigma,this->rho,this->v0);
                return this->S - exp(-this->r*this->tau)*this->X*(0.5 + prob/M_PI);
        }
        //
        //
        double putprice()
        {
                double prob = integrateheston(this->S,this->X,this->tau,this->r,this->kappa,this->theta,this->sigma,this->rho,this->v0);
                return  exp(-this->r*this->tau)*this->X*(0.5 - prob/M_PI);

        }
        //
        //
};
//
double callheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
{
        Heston *heston = new Heston(S,X,tau,r,kappa,theta,sigma,rho,v0);
        double price = heston->callprice();
        delete heston;
        return price;
}
//
//
double putheston(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0)
{
        Heston *heston = new Heston(S,X,tau,r,kappa,theta,sigma,rho,v0);
        double price = heston->putprice();
        delete heston;
        return price;
}
//
double fitfunction(vector<double> param,vector<double> stockprice,vector<double> strike,
                   vector<double> tau,vector<int> typed,vector<double> price)
{
    size_t n = price.size();
    //
    double kappa = param[0];
    double theta = param[1];
    double sigma = param[2];
    double rho = param[3];
    double v0 = param[4];
    //
    double fit = 0.0;
    for(int i=0;i<n;i++){
        if(typed[i] == 0){
            double pp = callheston(stockprice[i], strike[i], tau[i], 0.02,
                                  kappa, theta, sigma, rho, v0);
            fit += fabs(pp - price[i])/price[i];
        }else{
            double pp = putheston(stockprice[i], strike[i], tau[i], 0.02,
                                  kappa, theta, sigma, rho, v0);
            fit += fabs(pp - price[i])/price[i];
        }
    }
    return fit;
}
//
#endif
