/*!
 * @file mcparticles.hpp
 * @brief モンテカルロ粒子の定義に関わる部分
 * @author ButterPeanuts
 * @date 2023-04-12
*/
#pragma once
#include<vector>
#include<band.hpp>
#include<logger.hpp>
#include<memory>
/**
 * @brief モンテカルロ粒子のクラス
 * @details モンテカルロ粒子 略してMCParticle
 * 久木田(2014)(https://doi.org/10.18910/34463)のp.21で語られる概念
 * 簡単にするならばエネルギーEnergyを持つ粒子複数個の集まり
 * 一つのMC粒子について持つエネルギー(つまりEnergyと粒子の数の積)は一定
*/
namespace mc_sim{
	class mc_particle{
		private :
			/*! 空間の次元 何のために一般化しているのかはわかりません */
			static inline const int dimension = 3;
			
			/*! 角周波数\f$\omega(=E / \hbar)\f$ */
			double angular_frequency = 0;
			
			/*! 速度ベクトル方向 */
			std::vector<double> velocity_pointing;
			
			/*! 電荷q */
			double charge = 0;
			
			/*! 現在のバンド 参照する分散関係など */
			std::shared_ptr<band> band_current;
			
			/*! 粒子が取りうるすべてのバンド情報 */
			const std::vector<std::shared_ptr<band>>& banddata;
			
			std::shared_ptr<mc_sim::logger> logger;
			
			/*!
			 * @brief 粒子に弾性散乱を起こす
			 * @details 弾性散乱(速さが変化せず速度が変化するもの)を起こす
			 * 中身は速度方向ベクトルの新造
			 * Scatterで弾性散乱が発生したらこれが呼び出される
			*/
			void elastic_scattering();
			
			/*!
			 * @brief 粒子に非弾性散乱を起こす
			*/
			void inelastic_scattering(double temperature, double dt);
			
			/*!
			 * @brief 粒子の角周波数を再配置する
			 * @param double 温度
			 * @param bool Kirchhoff則適用の有無
			 * @param double Kirchhoff則における経過時間
			*/
			void angfreq_replace(double, bool, double = 1);
			
			/*!
			 * @brief MC粒子分布最大値が最大となるバンドのdomcp_distribution
			 * @param double 温度
			*/
			std::uniform_real_distribution<double> max_dd(double);
		
		public :
			/*! 変位 */
			std::vector<double> position;
			
			/*!
			 * @brief コンストラクタ
			 * @param newlogger ロガー
			 * @param Energy 粒子(MC粒子にあらず)のエネルギー \f$\omega\hbar\f$
			 * @param Temperature 不明 cppを解析されたい
			 * @param bandinj バンドデータの注入
			*/
			mc_particle(const std::shared_ptr<mc_sim::logger>& newlogger, double temperature, const std::vector<std::shared_ptr<band>>& bandinj);
			
			/*!
			 * @brief 時間dtの後の状態にMC粒子を遷移させる
			 * @param dt 格子時間dt
			*/
			void nextstep(double dt);
			
			/*!
			 * @brief 境界散乱Bを起こす
			 * @param max_x シミュレーション空間の大きさ?
			 * @param max_y シミュレーション空間の大きさ?
			 * @param max_z シミュレーション空間の大きさ?
			*/
			void boundaryscatter_b(double max_x, double max_y, double max_z);
			
			/*!
			 * @brief 散乱の判定を行い, 適切な散乱を発生させる
			 * @detail 棄却法によってどの散乱が発生するか決定
			 * その後然るべきメソッドを使用して散乱による変化を粒子に引き起こす
			 * @param Temperature 温度? 要検証 -> おそらくmcparticlesがいるところの温度
			 * @param dt MCParticles::Nextstep(double dt)などで使った時間間隔
			 * @param min_structure 確率P_bに関わるもの Simulationから呼び出せということか
			*/
			void scatter(double Temperature,double dt,double min_structure);
	};
}
