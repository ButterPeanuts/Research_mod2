#ifndef _INTEGRAL_H_
#define _INTEGRAL_H_
#include<vector>
#include<functional>
using namespace std;
double Semi_Para(const double x);
//double Trapezoidal(const double a,const double b,double(&Function)(double),const int n);
vector<double> Trapezoidal(const double a, const double b, std::function<double(double)> Function, const int n, vector<double> Memo);
//double Simpson(const double a,const double b,double(&Function)(double),const int n);
double Romberg(const double a, const double b, const int k, const int n,std::function<double(double)> Function);
#endif
