#pragma once
#include<functional>
#include<string>
#include<memory>

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
		std::vector<int> square(std::vector<double>);
		//空間メッシュ分割数
		std::vector<int> spacemesh;
		//x方向に進むのを想定することに
		//y方向の薄膜?
		std::vector<double> max_r;
		std::vector<double> dr;
	public:
		//physconstへ移転
		/* static double Total_energy2(double Temperature); */
		simulation(int, std::vector<double>&, std::vector<int>&, double, curve&, curve&, std::shared_ptr<mc_sim::logger>&, std::vector<std::shared_ptr<band>>&);
		double U, volume, energy_mcparticles;
		//粒子一覧
		std::vector<mc_sim::mc_particle> mc_particles;
		
		
		std::vector<std::vector<std::vector<double>>> Temperature;
		void Particle_move(double dt);
		bool Temperature_construct();
		void Particle_Disp_output(std::string filename);
};
