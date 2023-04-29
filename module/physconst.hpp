/*!
 * @file physconst.hpp
 * @brief 物理定数などの計算部分
 * @author ButterPeanuts
 * @date 2023-04-26
*/
#pragma once
#include<random>
#include<cmath>
#include<iostream>
/*!
 * @brief 物理定数, 他ヘルパー
 * @details 物理定数やヘルパー関数を収めるクラス
*/
class physconst {
	private:
		physconst();
	public:
		/*! ディラック定数 \f$ \hbar = \frac{h}{2\pi} [J\cdot s]\f$*/
		static inline const double dirac = 1.054571817e-34;
		/*! 円周率 */
		static inline const double pi = 3.14159265358979323846;
		/*! ボルツマン定数 \f$ k_B [J/K] \f$*/
		static inline const double boltzmann = 1.380649e-23;
		/*! 一様乱数発生器 */
		static std::mt19937_64 mtrand;
		/*!
		 * @brief ノイマンの棄却法による乱数発生器
		 * @param (*f)(double) 乱数分布, 確率変数を引数とし, 確率密度を返す
		 * @param xi 確率変数の下限
		 * @param xs 確率変数の上限
		 * @param fm 確率密度の上限
		*/
		static double vonNeumann_rejection(double (*f)(double),double xi, double xs, double fm);
};
