//
//  bates.hpp
//  OptionPricingHestonBates
//
//  Created by Antonio Santos on 01/08/2022.
//  Copyright © 2022 Antonio Santos. All rights reserved.
//

#ifndef bates_hpp
#define bates_hpp

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
//
class Bates{
protected:
    //
    double S;
    double X;
    double tau;
    //
    //
    double kappa;
    double theta;
    double sigma;
    double rho;
   	double v0;
    double lambda;
    double muj;
    double sigmaj;
    //
    double r;
    //
    //
public:
    //
    //
    Bates(double S,double X,double tau,double r,
          double kappa,double theta,double sigma,double rho,double v0,
          double lambda,double muj,double sigmaj)
    {
        this->S = S;
        this->X = X;
        this->tau = tau;
        this->kappa = kappa;
        this->theta = theta;
        this->sigma = sigma;
        this->rho = rho;
        this->v0 = v0;
        this->lambda = lambda;
        this->muj = muj;
        this->sigmaj = sigmaj;
        this->r = r;
    }
    //
    //
    //
    double callprice()
    {
        double prob = integratebates(this->S,this->X,this->tau,this->r,
                                     this->kappa,this->theta,this->sigma,this->rho,this->v0,
                                     this->lambda,this->muj,this->sigmaj);
        return this->S - exp(-this->r*this->tau)*this->X*(0.5 + prob/M_PI);
    }
    //
    //
    //
    double putprice()
    {
        double prob = integratebates(this->S,this->X,this->tau,this->r,
                                     this->kappa,this->theta,this->sigma,this->rho,this->v0,
                                     this->lambda,this->muj,this->sigmaj);
        return  exp(-this->r*this->tau)*this->X*(0.5 - prob/M_PI);
    }
    //
};
//
//
double valuebates(double S,double X,double tau,double r,
                  double kappa,double theta,double sigma,double rho,double v0,
                  double lambda,double muj,double sigmaj,double phi)
{
    double a = kappa*theta;
    double u = -0.5;
    double b = kappa;
    //
    Complex<double> c1 = Complex<double>(-b,rho*sigma*phi);
    Complex<double> c2 = Complex<double>(sigma*sigma*phi*phi,-2.0*u*sigma*sigma*phi);
    Complex<double> dd = (c1*c1 + c2).squareRoot();
    Complex<double> c3 = Complex<double>(b,-rho*sigma*phi) + dd;
    Complex<double> c4 = Complex<double>(b,-rho*sigma*phi) - dd;
    Complex<double> gg = c3/c4;
    Complex<double> cc = gg.inverse();
    Complex<double> c5 = (dd*(-tau)).exponential();
    Complex<double> DD = c4/sigma/sigma * ((1.0 - c5)/(1.0 - cc*c5));
    Complex<double> GG = (1.0 - cc*c5)/(1.0 - cc);
    Complex<double> CC = Complex<double>(0.0,(r - lambda*(exp(muj + 0.5*sigmaj*sigmaj) -1.0))*phi*tau ) + (a/sigma/sigma)*(tau*c4 - 2.0*GG.logarithm());
    Complex<double> ff = (CC + v0*DD + Complex<double>(0.0,phi*log(S))).exponential();
    Complex<double> cf = (tau*lambda*(Complex<double>(-0.5*sigmaj*sigmaj*phi*phi,muj*phi).exponential() - 1.0)).exponential();
    Complex<double> FF = Complex<double>(0.0,-phi*log(X)).exponential()*ff*cf/Complex<double>(phi*phi,phi);
    double re = FF.real();
    return re;
}
//
//
//
double integratebates(double S,double X,double tau,double r,
                      double kappa,double theta,double sigma,double rho,double v0,
                      double lambda,double muj,double sigmaj)
{
    double width = 0.5*(b-a);
    double mean = 0.5*(a+b);
    double gaussLegendre = 0;
    for(int step = 1; step <= n; step++) {
        gaussLegendre += weight[step]*valuebates(S,X,tau,r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj,width*root[step] + mean);
    }
    return gaussLegendre*width;
}


//
double callbates(double S,double X,double tau,double r,
                 double kappa,double theta,double sigma,double rho,double v0,
                 double lambda,double muj,double sigmaj)
{
    Bates *bates = new Bates(S,X,tau,r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
    double price = bates->callprice();
    delete bates;
    return price;
}
//
//
//
double putbates(double S,double X,double tau,double r,
                double kappa,double theta,double sigma,double rho,double v0,
                double lambda,double muj,double sigmaj)
{
    Bates *bates = new Bates(S,X,tau,r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
    double price = bates->putprice();
    delete bates;
    return price;
}

double fitfunction(vector<double> param,vector<double> stockprice,vector<double> strike,
                   vector<double> tau,vector<int> typed,vector<double> price)
{
    size_t n = price.size();
    double kappa = param[0];
    double theta = param[1];
    double sigma = param[2];
    double rho = param[3];
    double v0 = param[4];
    double lambda = param[5];
    double muj = param[6];
    double sigmaj = param[7];
    double fit = 0.0;
    for(int i=0;i<n;i++){
        if(typed[i] == 0){
            double pp = callbates(stockprice[i], strike[i], tau[i], 0.02,
                                  kappa, theta, sigma, rho, v0, lambda, muj, sigmaj);
            fit += fabs(pp - price[i])/price[i];
        }else{
            double pp = putbates(stockprice[i], strike[i], tau[i], 0.02,
                                  kappa, theta, sigma, rho, v0, lambda, muj, sigmaj);
            fit += fabs(pp - price[i])/price[i];
        }
    }
    return fit;
}
//
//
#endif /* bates_hpp */
