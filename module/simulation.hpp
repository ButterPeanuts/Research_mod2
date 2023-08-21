#pragma once
#include<functional>
#include<string>
#include<memory>
#include<mutex>

#include"mcparticle.hpp"
#include<curve.hpp>
#include<logger.hpp>
class simulation {
	private:
		std::shared_ptr<mc_sim::logger> logger;
		//いる?
		//bandって粒子の種類によって選択肢が決まるもんだし
		std::vector<std::shared_ptr<band>>& banddata;
		//std::vector<double> Internal_energy;
		curve& internal_energy;
		//比熱
		curve& heat_cap;
		
		//どのメッシュにいるのか
		std::array<int, 3> square(const std::vector<double>&);
		//空間メッシュ分割数
		std::vector<int> spacemesh;
		//x方向に進むのを想定することに
		//y方向の薄膜?
		std::vector<double> max_r;
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
		 * @details temperatureは3次元コンテナ, mcp_freqdistは1次元コンテナ
		 * この2つの座標系を橋渡しする
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
		//physconstへ移転
		/* static double Total_energy2(double Temperature); */
		simulation(int, std::vector<double>&, std::vector<int>&, double, curve&, curve&, std::shared_ptr<mc_sim::logger>&, std::vector<std::shared_ptr<band>>&);
		double U, volume, energy_mcparticles;
		//粒子一覧
		std::vector<mc_sim::mc_particle> mc_particles;
		
		
		std::vector<double> Temperature;
		void Particle_move(double dt);
		void Particle_Disp_output(std::string filename);
};
