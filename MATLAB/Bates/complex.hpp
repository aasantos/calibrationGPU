//
//  complex.hpp
//  Jul2022
//
//  Created by Antonio Santos on 15/07/2022.
//  Copyright © 2022 Antonio Santos. All rights reserved.
//

#ifndef complex_hpp
#define complex_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <math.h>
#include <stdlib.h>
//
//
//    Complex object
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
    Complex(){
        
    }
    //
    Complex(T x,T y)
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
    T real()
    {
        return this->x;
    }
    //
    //
    T imaginary()
    {
        return this->y;
    }
    //
    //
    T modulus()
    {
        return this->r;
    }
    //
    //
    T arg()
    {
        return this->theta;
    }
    //
    //
    Complex complexReal()
    {
        return Complex(this->x,0.0);
    }
    //
    //
    Complex complexImaginary()
    {
        return Complex(0.0,this->y);
    }
    //
    //
    Complex conjugate()
    {
        return Complex(this->x,-1.0*this->y);
    }
    //
    //
    Complex exponential()
    {
        return Complex(exp(this->x)*cos(this->y),exp(this->x)*sin(this->y));
    }
    //
    Complex logarithm()
    {
        return Complex(log(this->r),this->theta);
    }
    //
    Complex add(Complex z)
    {
        T xc = this->x + z.real();
        T yc = this->y + z.imaginary();
        return Complex(xc,yc);
    }
    //
    //
    Complex subtract(Complex z)
    {
        T xc = this->x - z.real();
        T yc = this->y - z.imaginary();
        return Complex(xc,yc);
    }
    //
    //
    Complex multiply(Complex z)
    {
        T zx = this->x*z.real() - this->y*z.imaginary();
        T zy = this->x*z.imaginary() + z.real()*this->y;
        return Complex(zx,zy);
    }
    //
    //
    Complex divide(Complex z)
    {
        T den = z.real()*z.real() + z.imaginary()*z.imaginary();
        T xc = this->x*z.real() + this->y*z.imaginary();
        T yc = this->y*z.real() - this->x*z.imaginary();
        return Complex(xc/den,yc/den);
    }
    //
    //   Arithmetic operations
    Complex operator+ (const Complex& other) const
    {
        return Complex(this->x + other.x, this->y + other.y);
    }
    //
    //
    Complex operator- (const Complex& other) const
    {
        return Complex(this->x - other.x, this->y - other.y);
    }
    //
    //
    Complex operator* (const Complex& other) const
    {
        return Complex(this->x * other.x - this->y * other.y,
                       this->x * other.y + this->y * other.x);
    }
    //
    //
    Complex operator/ (const Complex& other) const
    {
        const double denominator = other.x * other.x + other.y * other.y;
        return Complex((this->x * other.x + this->y * other.y) / denominator,
                       (this->y * other.x - this->x * other.y) / denominator);
    }
    //
    //
    Complex operator+ (const double val) const
    {
        return Complex(this->x + val, this->y);
    }
    
    Complex operator- (const double val) const
    {
        return Complex(this->x - val, this->y);
    }
    
    Complex operator* (const double val) const
    {
        return Complex(this->x * val, this->y * val);
    }
    
    Complex operator/ (const double val) const
    {
        return Complex(this->x / val, this->y / val);
    }
    
    Complex& operator+= (const double val)
    {
        this->x += val;
        return *this;
    }
    
    Complex& operator-= (const double val)
    {
        this->x -= val;
        return *this;
    }
    //
    //
    Complex& operator*= (const double val)
    {
        this->x *= val;
        this->y *= val;
        return *this;
    }
    //
    //
    Complex& operator/= (const double val)
    {
        this->x /= val;
        this->y /= val;
        return *this;
    }
    //
    friend Complex operator+ (const double left, const Complex& right)
    {
        return Complex(left + right.x, right.y);
    }
    
    friend Complex operator- (const double left, const Complex& right)
    {
        return Complex(left - right.x, -right.y);
    }
    
    friend Complex operator* (const double left, const Complex& right)
    {
        return Complex(left * right.x, left * right.y);
    }
    
    friend Complex operator/ (const double left, const Complex& right)
    {
        const double denominator = right.x * right.x + right.y * right.y;
        return Complex(left * right.x / denominator,
                       -left * right.y / denominator);
    }
    //
    Complex power(T z)
    {
        T rt = pow(this->r,z);
        T tt = z*this->theta;
        return Complex(rt*cos(tt),rt*sin(tt));
    }
    //
    //
    Complex squareRoot()
    {
        return Complex(this->x,this->y).power(0.5);
    }
    //
    //
    Complex inverse()
    {
        T rr = this->x*this->x + this->y*this->y;
        return Complex(this->x/rr,-this->y/rr);
    }
    //
    //
};
//

#endif /* complex_hpp */
