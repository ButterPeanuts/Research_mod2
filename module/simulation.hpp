#pragma once
#include<functional>
#include<string>
#include"mcparticles.hpp"
#include<curve.hpp>
#include<logger.hpp>
class simulation {
	private:
		mc_sim::logger& logger;
		std::vector<band>& banddata;
		//std::vector<double> Internal_energy;
		curve internal_energy;
		//比熱
		curve heat_cap;
		
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
		simulation(int, std::vector<double>, std::vector<int>, double, curve, curve, mc_sim::logger&, std::vector<band>&);
		double U, volume, energy_mcparticles;
		//粒子一覧
		std::vector<mc_particles::MCParticles> MCParticles;

		//状態密度
		static std::vector<std::vector<std::vector<std::vector<double>>>> dispersion;
		static std::vector<std::vector<double>> DOS_TA;
		static std::vector<std::vector<double>> DOS_LA;


		std::vector<std::vector<std::vector<double>>> Temperature;
		void Particle_move(double dt);
		bool Temperature_construct();
		void Particle_Disp_output(std::string filename);
};
