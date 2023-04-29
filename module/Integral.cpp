#include"Integral.h"
#include<cmath>
#include<iostream>
#include<vector>
#include<functional>
using namespace std;
double Semi_Para(const double x){
	return x * sqrt(x * x + 1);
}
/*
double Trapezoidal(const double a,const double b,double (&Function)(double),const int n){
	if(n < 0)return 0;
	if(n == 0)return (b - a) * (Function(a) + Function(b)) / 2;
	double Tndec;
	//1段階荒い積分値を求める
	Tndec = Trapezoidal(a,b,Function,n - 1);
	//Tndecを1段階上の精度にすべく計算
	double d = 0;
	for(int i = 0;i < (int)pow(2.0,(double)n - 1);i++){
		d += Function(a + (1 + 2 * i) * (b - a) / pow(2.0,(double)n));
	}
	return (Tndec / 2 + d * (b - a) / pow(2.0,(double)n));
}
*/
//メモ化対応版
vector<double> Trapezoidal(const double a,const double b,std::function<double(double)> Function,const int n,vector<double> Memo){
	if(n < 0)return Memo;
	if (n == 0) {
		double T0 = (b - a) * (Function(a) + Function(b)) / 2;
		Memo.push_back(T0);
		return Memo;
	}
	Memo = Trapezoidal(a,b,Function,n - 1,Memo);
	double d = 0;
	for(int i = 0;i < (int)pow(2.0,(double)n - 1);i++){
		d += Function(a + (1.0l + 2.0l * (double)i) * (b - a) / pow(2.0,(double)n));
	}
	double Tn = (Memo[n-1] / 2 + d * (b - a) / pow(2.0,(double)n));
	Memo.push_back(Tn);
	return Memo;
}
/*
double Simpson(const double a,const double b,double (&Function)(double),const int n){
	return (4 * Trapezoidal(a,b,Function,n) - Trapezoidal(a,b,Function,n - 1)) / 3;
}
*/
double Romberg(const double a, const double b, const int k, const int n, std::function<double(double)> Function) {
	//n < kの場合ロンバーグ積分不能
	if (n < k)return 0;
	vector<vector<double>> Rnk;
	vector<double> Emp;
	Rnk.push_back(Emp);
	Rnk[0] = Trapezoidal(a, b, Function, n, Rnk[0]);
	//関数でのメモ化がめんd
	//よってここから手続き型風
	for (int ik = 1; ik <= k; ik++) {
		vector<double> Empn(n + 1);
		for (int in = ik + n - k; in <= n; in++) {
			Empn[in] = (pow(4, ik) * Rnk[ik - 1][in] - Rnk[ik - 1][in - 1]) / (pow(4, ik) - 1);
		}
		Rnk.push_back(Empn);
	}
	return Rnk[k][n];
}
