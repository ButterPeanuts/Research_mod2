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
#include<utility>
#include<functional>
/*!
 * @brief 物理定数, 他ヘルパー
 * @details 物理定数やヘルパー関数を収めるクラス
*/
class physconst {
	private:
		physconst();
	public:
		/*! ディラック定数 \f$ \hbar = \frac{h}{2\pi} [J\cdot s]\f$*/
		constexpr static inline double dirac = 1.054571817e-34;
		/*! ボルツマン定数 \f$ k_B [J/K] \f$*/
		constexpr static inline double boltzmann = 1.380649e-23;
		/*! 一様乱数発生器 */
		static std::mt19937_64 mtrand;
		/*!
		 * @brief ノイマンの棄却, 関数を呼び出すたびに棄却か採用かを返す
		 * @param f 乱数分布, 確率変数を引数とし, 確率密度を返す関数
		 * @param xdist 確率変数の範囲
		 * @param fdist 確率密度の範囲
		 * @param engine 乱数器
		 * @return <結果, 採用ならその時の確率変数>のタプル
		*/
		static std::pair<bool, double> vonNeumann_rejection(std::function<double(double)>& f, std::uniform_real_distribution<>&& xdist, std::uniform_real_distribution<>& fdist, std::mt19937_64& engine);
		
		/*!
		 * @brief 非推奨
		 * @details Bose-Einstein統計関数
		 * @param ang_freq 各周波数\f$\omega\f$
		 * @param temp 温度\f$T\f$
		*/
		static double bedist(double ang_freq, double temp);
		
		/*!
		 * @brief Bose-Einstein統計関数
		 * @details Bose-Einstein統計関数の定数倍を返す
		 * @param ang_freq 各周波数\f$\omega\f$
		 * @param temp 温度\f$T\f$
		 * @param left_const 統計関数の計算結果にかける定数
		*/
		static double bedist2(double ang_freq, double temp, double left_const);
		
		/*!
		 * @brief メッシュ座標(index)から規格化座標を返す
		 * @details 0からndivの値を0から1に直す
		 * @param std::tuple<int, int, int> メッシュ座標(index)
		 * @param int 分割数(ndiv)
		 * @return std::tuple<double, double, double> 規格化座標
		*/
		static std::tuple<double, double, double> indextostd(std::tuple<int, int, int> const &, int const);
		
		/*!
		 * @brief エウクレイデス距離(ユークリッド距離)を求める
		 * @param std::tuple<double, double, double> 座標
		 * @return エウクレイデス距離
		*/
		static double eukleideia_metrike(std::tuple<double, double, double> const &);
};
