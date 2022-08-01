#ifndef complex_cuh
#define complex_cuh
//
//
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <math.h>
#include <stdlib.h>
//
//
template <typename T>
class Complex{
protected:
    //
    //
    T x; //real part
    T y; //imaginary part
    T r;
    T theta;
    //
    //
public:
    //
    __host__ __device__ Complex(){
    }
    //
    __host__ __device__ Complex(T x,T y)
    {
        //
        //
        this->x = x;
        this->y = y;
        this->r = sqrt(this->x*this->x + this->y*this->y);
        this->theta = atan2(this->y, this->x);
        //
        //
    }

    //
    //
    __host__ __device__ T real()
    {
        return this->x;
    }
    //
    //
    __host__ __device__ T imaginary()
    {
        return this->y;
    }
    //
    //
    __host__ __device__ T modulus()
    {
        return this->r;
    }
    //
    //
    __host__ __device__ T arg()
    {
        return this->theta;
    }
    //
    //
    __host__ __device__ Complex complexReal()
    {
        return Complex(this->x,0.0);
    }

     //
    __host__ __device__ Complex complexImaginary()
    {
        return Complex(0.0,this->y);
    }
    //
    //
    __host__ __device__ Complex conjugate()
    {
        return Complex(this->x,-1.0*this->y);
    }
    //
    //
    __host__ __device__ Complex exponential()
    {
        return Complex(exp(this->x)*cos(this->y),exp(this->x)*sin(this->y));
    }
    //
    __host__ __device__ Complex logarithm()
    {
        return Complex(log(this->r),this->theta);
    }
    //
    __host__ __device__ Complex add(Complex z)
    {
        T xc = this->x + z.real();
        T yc = this->y + z.imaginary();
        return Complex(xc,yc);
    }
    //
    //
    __host__ __device__ Complex subtract(Complex z)
    {
        T xc = this->x - z.real();
        T yc = this->y - z.imaginary();
        return Complex(xc,yc);
    }
    //
    //
    __host__ __device__ Complex multiply(Complex z)
    {
        T zx = this->x*z.real() - this->y*z.imaginary();
        T zy = this->x*z.imaginary() + z.real()*this->y;
        return Complex(zx,zy);
    }
    //
    //
    __host__ __device__ Complex divide(Complex z)
    {
        T den = z.real()*z.real() + z.imaginary()*z.imaginary();
        T xc = this->x*z.real() + this->y*z.imaginary();
        T yc = this->y*z.real() - this->x*z.imaginary();
        return Complex(xc/den,yc/den);
    }
    //
    //   Arithmetic operations
    __host__ __device__ Complex operator+ (const Complex& other) const
    {
        return Complex(this->x + other.x, this->y + other.y);
    }
    //
    //
    __host__ __device__ Complex operator- (const Complex& other) const
    {
        return Complex(this->x - other.x, this->y - other.y);
    }
    //
    //
    __host__ __device__ Complex operator* (const Complex& other) const
    {
        return Complex(this->x * other.x - this->y * other.y,
            this->x * other.y + this->y * other.x);
    }
    //
    //
    __host__ __device__ Complex operator/ (const Complex& other) const
    {
        const double denominator = other.x * other.x + other.y * other.y;
        return Complex((this->x * other.x + this->y * other.y) / denominator,
            (this->y * other.x - this->x * other.y) / denominator);
    }
    //
    //
    //
    __host__ __device__ Complex operator+ (const double val) const
    {
        return Complex(this->x + val, this->y);
    }
    //
    //
    __host__ __device__ Complex operator- (const double val) const
    {
        return Complex(this->x - val, this->y);
    }
    //
    //
    __host__ __device__ Complex operator* (const double val) const
    {
        return Complex(this->x * val, this->y * val);
    }
    //
    //
    __host__ __device__ Complex operator/ (const double val) const
    {
        return Complex(this->x / val, this->y / val);
    }
    //
    //
    __host__ __device__ Complex& operator+= (const double val)
    {
        this->x += val;
        return *this;
    }
    //
    //
    __host__ __device__ Complex& operator-= (const double val)
    {
        this->x -= val;
        return *this;
    }
    //
    //
    //
    __host__ __device__ Complex& operator*= (const double val)
    {
        this->x *= val;
        this->y *= val;
        return *this;
    }
    //
    //
    __host__ __device__ Complex& operator/= (const double val)
    {
        this->x /= val;
        this->y /= val;
        return *this;
    }
    //
    //
    __host__ __device__ friend Complex operator+ (const double left, const Complex& right)
    {
        return Complex(left + right.x, right.y);
    }
    //
    //
    __host__ __device__ friend Complex operator- (const double left, const Complex& right)
    {
        return Complex(left - right.x, -right.y);
    }
    //
    //
    __host__ __device__ friend Complex operator* (const double left, const Complex& right)
    {
        return Complex(left * right.x, left * right.y);
    }
    //
    //
    __host__ __device__ friend Complex operator/ (const double left, const Complex& right)
    {
        const double denominator = right.x * right.x + right.y * right.y;
        return Complex(left * right.x / denominator,
            -left * right.y / denominator);
    }
    //
    //
    __host__ __device__ Complex power(T z)
    {
        T rt = pow(this->r,z);
        T tt = z*this->theta;
        return Complex(rt*cos(tt),rt*sin(tt));
    }
    //
    //
    __host__ __device__ Complex squareRoot()
    {
        return Complex(this->x,this->y).power(0.5);
    }
    //
    //
    __host__ __device__ Complex inverse()
    {
        T rr = this->x*this->x + this->y*this->y;
        return Complex(this->x/rr,-this->y/rr);
    }
    //
    //
    __host__ __device__ void print()
    {
        printf("(%.6f,%.6f)\n",this->x,this->y);
    }
    //
    //
    __host__ __device__ void printPolar()
    {
        printf("(%.6f,%.6f)\n",this->r,this->theta);
    }
    //
    //
};


#endif
