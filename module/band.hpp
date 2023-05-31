/*!
 * @file band.hpp
 * @brief 粒子が存在する? バンドについての情報
 * @author ButterPeanuts
 * @data 2023-05-10
*/
#pragma once
#include "modeenum.hpp"

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
	 * @brief dos_getterのomegaとして指定しうる最小値を取得する関数
	*/
	virtual double dos_omega_min_getter() = 0;
	
	/*!
	 * @brief dos_getterのomegaとして指定しうる最大値値を取得する関数
	*/
	virtual double dos_omega_max_getter() = 0;
	
	/*!
	 * @brief 状態密度の最大値を取得する関数
	*/
	virtual double dos_max_getter() = 0;
	
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
	virtual double domcp_getter(double omega) = 0;
	
	/*!
	 * @brief domcp_getterのomegaとして指定しうる最小値を取得する関数
	*/
	virtual double domcp_omega_min_getter() = 0;
	
	/*!
	 * @brief domcp_getterのomegaとして指定しうる最大値を取得する関数
	*/
	virtual double domcp_omega_max_getter() = 0;
	
	/*!
	 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)の最大値を取得する関数
	 * @details かなり純粋なgetter これがないと棄却法にO(n)かかる
	*/
	virtual double domcp_max_getter() = 0;
	
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
