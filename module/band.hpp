/*!
 * @file band.hpp
 * @brief 粒子が存在する? バンドについての情報
 * @author ButterPeanuts
 * @data 2023-05-10
*/
#pragma once
#include "modeenum.hpp"
#include<random>

/*!
 * @brief 粒子が存在する? バンドについての情報のクラス
 * @details 粒子が存在する? 従う? バンドについての情報のクラス
 * バンドが関連する物質定数, 定関数の類を表としてまとめる
*/
class band{
	public:
	/*!
	 * @brief 状態密度(Density of State)を取得する関数
	 * @param omega 角周波数\f$\omega\f$
	*/
	virtual double dos_getter(double omega) = 0;
	
	/*!
	 * @brief dos_getterのomegaとして指定しうる範囲を取得する関数
	 * @return dos_getterのomegaとして指定しうる範囲のdistribution
	*/
	virtual std::uniform_real_distribution<> dos_omega_distribution_getter() = 0;
	
	/*!
	 * @brief 状態密度の範囲を取得する関数
	 * @return 状態密度の範囲のdistribution
	*/
	virtual std::uniform_real_distribution<> dos_distribution_getter() = 0;
	
	/*!
	 * @brief (角周波数-)MC粒子密度の範囲を取得する関数
	 * @param double 温度
	 * @return MC粒子密度の範囲のdistribution
	*/
	virtual std::uniform_real_distribution<> domcp_distribution_getter(double) = 0;
	
	/*!
	 * @brief dos_getterのomegaとして指定しうる範囲の下端を取得する関数
	 * @return dos_getterのomegaとして指定しうる範囲の下端(double)
	*/
	virtual double dos_leftedge() = 0;
	
	/*!
	 * @brief dos_getterのomegaとして指定しうる範囲の上端を取得する関数
	 * @return dos_getterのomegaとして指定しうる範囲の上端(double)
	*/
	virtual double dos_rightedge() = 0;
	
	/*!
	 * @brief 群速度(Group velocity)を取得する関数
	 * @param omega 角周波数\f$\omega\f$
	*/
	virtual double gvelocity_getter(double omega) = 0;
	
	/*!
	 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)を取得する関数
	 * @details 詳細は久木田(2013)の式2.13を参照
	 * @param omega 角周波数\f$\omega\f$
	*/
	//virtual double domcp_getter(double omega) = 0;
	
	/*!
	 * @brief domcp_getterのomegaとして指定しうる範囲を取得する関数
	 * @return domcp_getterのomegaとして指定しうる範囲のdistribution
	*/
	//virtual std::uniform_real_distribution<> domcp_omega_distribution_getter() = 0;
	
	/*!
	 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)の範囲を取得する関数
	 * @details かなり純粋なgetter これがないと棄却法にO(n)かかる
	 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)の範囲を表すdistribution
	*/
	//virtual std::uniform_real_distribution<> domcp_distribution_getter() = 0;
	
	/*!
	 * @brief 縦波か横波かを取得する関数
	*/
	virtual wave_direction directions_getter() = 0;
	
	/*!
	 * @brief 音響波か光学波か取得する関数
	*/
	virtual wave_mode mode_getter() = 0;
	
	/*!
	 * @brief ウムクラップ散乱の散乱係数a
	*/
	virtual double a() = 0;
	
	/*!
	 * @brief ウムクラップ散乱の散乱係数b
	*/
	virtual double b() = 0;
	
	/*!
	 * @brief ウムクラップ散乱の散乱係数\f$\chi\f$
	*/
	virtual double chi() = 0;
	
	/*!
	 * @brief ウムクラップ散乱の散乱係数\f$\xi\f$
	*/
	virtual double xi() = 0;
	
	/*!
	 * @brief 欠陥散乱の散乱係数c
	*/
	virtual double c() = 0;
	
	/*!
	 * @brief 境界散乱の散乱係数f
	*/
	virtual double f() = 0;
};
