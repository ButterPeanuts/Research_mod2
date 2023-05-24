/*!
 * @file band_obj.hpp
 * @brief 粒子が存在する? バンドについての情報
 * @author ButterPeanuts
 * @data 2023-05-10
*/
#pragma once
#include <string>
#include "modeenum.hpp"
#include "band.hpp"
#include <curve.hpp>
#include <logger.hpp>

/*!
 * @brief 粒子が存在する? バンドについての情報のクラス
 * @details 粒子が存在する? 従う? バンドについての情報のクラス
 * バンドが関連する物質定数, 定関数の類を表としてまとめる
*/
class band_obj : public band{
	private:
		curve dos, gvelocity, domcp;
		double domcp_max, dos_max;
		wave_direction directions;
		wave_mode mode;
		mc_sim::logger& logger;
	public:
		band_obj(mc_sim::logger& newlogger, curve dos, curve gvelocity, curve domcp, wave_direction direction, wave_mode mode);
		/*!
		 * @brief 状態密度(Density of State)を取得する関数
		 * @param omega 角周波数\f$\omega\f$
		*/
		double dos_getter(double omega) override;
		
		/*!
		 * @brief dos_getterのomegaとして指定しうる最小値を取得する関数
		*/
		double dos_omega_min_getter() override;
		
		/*!
		 * @brief dos_getterのomegaとして指定しうる最大値値を取得する関数
		*/
		double dos_omega_max_getter() override;
		
		/*!
		 * @brief 状態密度の最大値を取得する関数
		*/
		double dos_max_getter() override;
		
		/*!
		 * @brief 群速度(Group velocity)を取得する関数
		 * @param omega 角周波数\f$\omega\f$
		*/
		double gvelocity_getter(double omega) override;
		
		/*!
		 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)を取得する関数
		 * @details 詳細は久木田(2013)の式2.13を参照
		 * @param omega 角周波数\f$\omega\f$
		*/
		double domcp_getter(double omega) override;
		
		/*!
		 * @brief domcp_getterのomegaとして指定しうる最小値を取得する関数
		*/
		virtual double domcp_omega_min_getter() override;
		
		/*!
		 * @brief domcp_getterのomegaとして指定しうる最大値を取得する関数
		*/
		virtual double domcp_omega_max_getter() override;
		
		/*!
		 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)の最大値を取得する関数
		 * @details かなり純粋なgetter これがないと棄却法にO(n)かかる
		*/
		double domcp_max_getter() override;
		
		/*!
		 * @brief 縦波か横波かを取得する関数
		*/
		wave_direction directions_getter() override;
		
		/*!
		 * @brief 音響波か光学波か取得する関数
		*/
		wave_mode mode_getter() override;
};
