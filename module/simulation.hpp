/*!
 * @file simulation.hpp
 * @brief シミュレーター自体のクラス
 * @author ButterPeanuts
 * @date 2023-08-23
*/
#pragma once
#include<functional>
#include<string>
#include<memory>
#include<mutex>

#include"mcparticle.hpp"
#include<curve.hpp>
#include<logger.hpp>
/**
 * @brief シミュレーター自体のクラス
*/
class simulation {
	private:
		/*! ロガー */
		std::shared_ptr<mc_sim::logger> logger;
		
		/*! 粒子が状態として取りうるバンド */
		std::vector<std::shared_ptr<band>>& banddata;
		
		/*! 温度から内部エネルギー密度を求めるcurve */
		curve& internal_energy;
		
		/*! 内部エネルギー密度から温度を求める関数 比熱ではない */
		curve& heat_cap;
		
		/*!
		 * @brief 粒子の連続座標をmeshの離散三次元座標に直す
		 * @details 粒子の座標からどのmeshにいるかを探す
		 * @param const std::vector<double>& 3次元連続座標
		 * @return std::array<int, 3> 3次元離散座標
		*/
		std::array<int, 3> square(const std::vector<double>&);
		
		/*! 空間メッシュ分割数 */
		std::vector<int> spacemesh;
		
		/*! シミュレーション空間の大きさ */
		std::vector<double> max_r;
		
		/*! メッシュの大きさ */
		std::vector<double> dr;
		
		/*! 各メッシュごとのMC粒子数 */
		std::vector<int> mcp_freqdist;
		
		/*!
		 * @brief 各メッシュにおける粒子の頻度分布を更新する
		 * @details 頻度分布を計算し, mcp_freqdistを更新する
		*/
		void freqdist_construct();
		
		/*!
		 * @brief meshの離散三次元座標を離散1次元座標に直す
		 * @detail temperatureやmcp_freqdistは1次元コンテナなのでそのindexへ直すために
		 * @param const std::array<int, 3>& 3次元座標
		 * @return int 1次元座標
		*/
		int tempcoor_to_fdlinear(const std::array<int, 3>&);
		
		/*!
		 * @brief 各メッシュにおける粒子の温度分布を更新する
		 * @details 温度分布を計算し, temperatureを更新する
		*/
		bool temperature_construct();
	public:
		/*!
		 * @brief コンストラクタ
		 * @param int mc粒子の数
		 * @param std::vector<double>& max_r, シミュレーション空間の大きさ
		 * @param std::vector<int>& シミュレーション空間の分割数
		 * @param double シミュレーション空間の温度が均一になったときの温度
		 * @param curve& 温度から内部エネルギー密度を求めるcurve
		 * @param curve& 内部エネルギー密度から温度を求めるcurve
		 * @param std::shared_ptr<mc_sim::logger>& ロガー
		 * @param std::vector<std::shared_ptr<band>>& band一覧
		*/
		simulation(int, std::vector<double>&, std::vector<int>&, double, curve&, curve&, std::shared_ptr<mc_sim::logger>&, std::vector<std::shared_ptr<band>>&);
		
		/*! シミュレーション空間の粒子が持つエネルギー */
		double U;
		
		/*! シミュレーション空間の体積 */
		double volume;
		
		/*! mc粒子の持つエネルギー */
		double energy_mcparticles;
		
		/*! 粒子一覧 */
		std::vector<mc_sim::mc_particle> mc_particles;
		
		/*! 各メッシュの温度 */
		std::vector<double> Temperature;
		
		/*!
		 * @brief 粒子の移動, 散乱を行う
		 * @param dt 格子時間dt
		*/
		void Particle_move(double dt);
		
		/*!
		 * @brief 粒子の変位をファイルに出力する
		 * @param filename ファイル名
		*/
		void Particle_Disp_output(std::string filename);
		
		/*!
		 * @brief 指定された関数でmc粒子の変位を初期化する
		 * @param const std::function<std::vector<double>()>& 変位を出力する関数
		*/
		void particle_posinit(const std::function<std::vector<double>()>&);
};
