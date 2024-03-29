/*!
 * @file band_obj.hpp
 * @brief 粒子が存在する? バンドについての情報
 * @author ButterPeanuts
 * @data 2023-05-10
*/
#pragma once
#include <string>
#include<random>
#include "modeenum.hpp"
#include "band.hpp"
#include <curve.hpp>
#include <logger.hpp>
#include <scatconst.hpp>

/*!
 * @brief 粒子が存在する? バンドについての情報のクラス
 * @details 粒子が存在する? 従う? バンドについての情報のクラス
 * バンドが関連する物質定数, 定関数の類を表としてまとめる
*/
class band_obj : public band{
	private:
		curve& dos, gvelocity, domcpmax;
		double dos_max;
		wave_direction directions;
		wave_mode mode;
		scatconst& bands_scatconst;
		std::shared_ptr<mc_sim::logger> logger;
	public:
		band_obj(std::shared_ptr<mc_sim::logger>& newlogger, curve& dos, curve& gvelocity, curve& domcpmax, wave_direction direction, wave_mode mode, scatconst& bands_scatconst);
		/*!
		 * @brief 状態密度(Density of State)を取得する関数
		 * @param omega 角周波数\f$\omega\f$
		*/
		double dos_getter(double omega) override;
		
		/*!
		 * @brief dos_getterのomegaとして指定しうる範囲を取得する関数
		 * @return dos_getterのomegaとして指定しうる範囲のdistribution
		*/
		std::uniform_real_distribution<> dos_omega_distribution_getter() override;
		
		/*!
		 * @brief 状態密度の範囲を取得する関数
		 * @return 状態密度の範囲のdistribution
		*/
		std::uniform_real_distribution<> dos_distribution_getter() override;
		
		/*!
		 * @brief (角周波数-)MC粒子密度の範囲を取得する関数
		 * @param double 温度
		 * @return MC粒子密度の範囲のdistribution
		*/
		std::uniform_real_distribution<> domcp_distribution_getter(double) override;
		
		/*!
		 * @brief 群速度(Group velocity)を取得する関数
		 * @param omega 角周波数\f$\omega\f$
		*/
		double gvelocity_getter(double omega) override;
		
		/*!
		 * @brief dos_getterのomegaとして指定しうる範囲の下端を取得する関数
		 * @return dos_getterのomegaとして指定しうる範囲の下端(double)
		*/
		double dos_leftedge() override;
		
		/*!
		 * @brief dos_getterのomegaとして指定しうる範囲の上端を取得する関数
		 * @return dos_getterのomegaとして指定しうる範囲の上端(double)
		*/
		double dos_rightedge() override;
		
		/*!
		 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)を取得する関数
		 * @details 詳細は久木田(2013)の式2.13を参照
		 * @param omega 角周波数\f$\omega\f$
		*/
		//double domcp_getter(double omega) override;
		
		/*!
		 * @brief domcp_getterのomegaとして指定しうる範囲を取得する関数
		 * @return domcp_getterのomegaとして指定しうる範囲のdistribution
		*/
		//std::uniform_real_distribution<> domcp_omega_distribution_getter() override;
		
		/*!
		 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)の範囲を取得する関数
		 * @details かなり純粋なgetter これがないと棄却法にO(n)かかる
		 * @brief モンテカルロ粒子密度(Density of MonteCalro Particle)の範囲を表すdistribution
		*/
		//std::uniform_real_distribution<> domcp_distribution_getter() override;
		
		/*!
		 * @brief 縦波か横波かを取得する関数
		*/
		wave_direction directions_getter() override;
		
		/*!
		 * @brief 音響波か光学波か取得する関数
		*/
		wave_mode mode_getter() override;
		
		/*!
		 * @brief ウムクラップ散乱の散乱係数a
		*/
		double a() override;
		
		/*!
		 * @brief ウムクラップ散乱の散乱係数b
		*/
		double b() override;
		
		/*!
		 * @brief ウムクラップ散乱の散乱係数\f$\chi\f$
		*/
		double chi() override;
		
		/*!
		 * @brief ウムクラップ散乱の散乱係数\f$\xi\f$
		*/
		double xi() override;
		
		/*!
		 * @brief 欠陥散乱の散乱係数c
		*/
		double c() override;
		
		/*!
		 * @brief 境界散乱の散乱係数f
		*/
		double f() override;
		
		/*!
		 * @brief フォノン-フォノン散乱確率
		 * @param double 角周波数\f$\omega\f$
		 * @param double 温度\f$T\f$
		 * @param double 経過時間\f$\Delta t\f$
		*/
		double scatprob_u(double, double, double) override;
		/*!
		 * @brief フォノン-欠陥散乱確率
		 * @param double 角周波数\f$\omega\f$
		 * @param double 経過時間\f$\Delta t\f$
		*/
		double scatprob_d(double, double) override;
		/*!
		 * @brief フォノン-境界散乱確率
		 * @param double 群速度\f$\bar(v)\f$
		 * @param double 特徴長さ\f$L\f$
		 * @param double 経過時間\f$\Delta t\f$
		*/
		double scatprob_b(double, double, double) override;
		
		/*!
		 * @brief 全緩和時間の最小を求める
		 * @param double 最大温度\f$T_{max}\f$
		 * @param double 特徴長さ\f$L\f$
		*/
		double mintau(double, double) override;
};
