#pragma once
#include<functional>
#include<string>
#include"mcparticles.hpp"
#include<curve.hpp>
class simulation {
	private:
	public:
		static double Total_energy2(double Temperature);
		simulation(int NumberofMCparticles,double max_x,double max_y,double max_z,std::vector<int> spacemesh, double, curve);
		double U, volume,Energy_MCParticles;
		//x方向に進むのを想定することに
		//y方向の薄膜?
		double max_x, max_y, max_z;
		//空間メッシュ分割数
		std::vector<int> spacemesh;
		//粒子一覧
		std::vector<mc_particles::MCParticles> MCParticles;

		//状態密度
		static std::vector<std::vector<std::vector<std::vector<double>>>> dispersion;
		static std::vector<std::vector<double>> DOS_TA;
		static std::vector<std::vector<double>> DOS_LA;

		//比熱
		std::vector<double> Heat_cap;
		//std::vector<double> Internal_energy;
		curve internal_energy;

		std::vector<std::vector<std::vector<double>>> Temperature;
		void Particle_move(double dt);
		bool Temperature_construct();
		void Particle_Disp_output(std::string filename);
};
