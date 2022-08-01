#include "complex.cuh"
#include "constvectors1000.cuh"
//
__device__ const double a = 0.0;
__device__ const double b = 500.0;
__device__ const int n = 1000;
//
//
__device__ double valuebates(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double lambda,double muj,double sigmaj,double phi)
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
__device__ double integratebates(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double lambda,double muj,double sigmaj)
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
//
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
	__device__ Bates(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double lambda,double muj,double sigmaj)
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
	__device__ double callprice()
	{
	double prob = integratebates(this->S,this->X,this->tau,this->r,this->kappa,this->theta,this->sigma,this->rho,this->v0,this->lambda,this->muj,this->sigmaj);
	return this->S - exp(-this->r*this->tau)*this->X*(0.5 + prob/M_PI);
	}
	//
	//
	//
	__device__ double putprice()
	{
	  double prob = integratebates(this->S,this->X,this->tau,this->r,this->kappa,this->theta,this->sigma,this->rho,this->v0,this->lambda,this->muj,this->sigmaj);
	  return  exp(-this->r*this->tau)*this->X*(0.5 - prob/M_PI);
	}
	//
};
//
//
//
__device__ double callbates(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double lambda,double muj,double sigmaj)
{
        Bates *bates = new Bates(S,X,tau,r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
        double price = bates->callprice();
        delete bates;
        return price;
}
//
//
//
__device__ double putbates(double S,double X,double tau,double r,double kappa,double theta,double sigma,double rho,double v0,double lambda,double muj,double sigmaj)
{
        Bates *bates = new Bates(S,X,tau,r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
        double price = bates->putprice();
        delete bates;
        return price;
}
//
//
//
__global__ void kernelbates(double *price,double *stockprice,double *strike,double *tau,int *type,int n,double r,double kappa,double theta,double sigma,double rho,double v0,double lambda,double muj,double sigmaj)
{
        int idx = blockDim.x*blockIdx.x + threadIdx.x;
        if(idx < n){
                if(type[idx] == 0){
                    price[idx] = callbates(stockprice[idx],strike[idx],tau[idx],r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
                }else{
                    price[idx] = putbates(stockprice[idx],strike[idx],tau[idx],r,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
                }
        }
}
//
//

