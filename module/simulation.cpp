#include<functional>
#include<string>
#include<fstream>
#include<algorithm>
#include<omp.h>

#include"simulation.hpp"
#include"physconst.hpp"
#include"massconst.hpp"
#include"Integral.h"
#include"mcparticles.hpp"
simulation::simulation(int numof_mcp, double max_x,double max_y,double max_z,std::vector<int> spacemesh, double tempof_device, curve internal_energy, curve heat_cap, mc_sim::logger& logger, std::vector<band>& band_inj) :  logger(logger), banddata(band_inj), internal_energy(internal_energy), heat_cap(heat_cap) {
	//体積
	this->volume = max_x * max_y * max_z;
	
	
	//physconstへ移管
	/* Internal_energy[0] = (0); */
	/* #pragma omp parallel for */
	/* 	for (int Temperature = 1; Temperature <= massconst::heatcaps_Tempmax ; Temperature++) { */
	/* 		double InE = 0; */
	/* 		InE += simulation::Total_energy2(Temperature); */
	/* 		Internal_energy[Temperature] = InE; */
	/* 	} */
	/* #pragma omp barrier */
	/* std::cout << "Internal_energy table was calculated." << std::endl; */

	//全体のエネルギー
	this->U = this->internal_energy.itpl_getter(tempof_device) * this->volume;

	this->energy_mcparticles = this->U / numof_mcp;

	this->max_x = max_x;
	this->max_y = max_y;
	this->max_z = max_z;
	
	this->MCParticles = std::vector<mc_particles::MCParticles>();
	for (int i = 0; i < numof_mcp; i++) {
		mc_sim::logger& newlogger = this->logger.copy_samesink("mcp" + std::to_string(i));
		this->MCParticles.push_back(mc_particles::MCParticles(newlogger, static_cast<double>(tempof_device), this->banddata));
	}
	this->Temperature = std::vector<std::vector<std::vector<double>>>(spacemesh[0], std::vector < std::vector<double>>(spacemesh[1], std::vector<double>(spacemesh[2], 0.0)));
	this->spacemesh = spacemesh;
	
	this->logger.info("Simulation is constructed");
}

//massconstへ移転
/* double simulation::Total_energy2(double Temperature){ */
/* 	return Romberg(0, Temperature, 10, 10, [Temperature](double T) { */
/* 		if (T <= 0)return massconst::Si_heatcap[0]; */
/* 		if (T >= massconst::heatcaps_Tempmax)return *(massconst::Si_heatcap.end() - 1); */
/* 		double standard = T - std::floor(T); */
/* 		return (1 - standard) * massconst::Si_heatcap[(int)std::floor(T)] */
/* 			+ (standard) * massconst::Si_heatcap[(int)std::floor(T) + 1]; */ 
/* 	}); */
/* }; */

void simulation::Particle_Disp_output(std::string filename) {
	std::ofstream files(filename + ".txt");
	if (!files) {
		std::cout << "保存に失敗しました" << std::endl;
		return;
	}
	for (auto i = MCParticles.begin(); i < MCParticles.end(); i++) {
		files << i->position[0] << ", "<<i->position[1] <<", "<<i->position[2] <<std::endl;
	}
	files.close();
}

bool simulation::Temperature_construct() {
	//各meshのEnergy密度
	std::vector<std::vector<std::vector<double>>> MeshEnergy= std::vector<std::vector<std::vector<double>>>(spacemesh[0], std::vector < std::vector<double>>(spacemesh[1], std::vector<double>((int)spacemesh[2], 0)));
	int index_x,index_y,index_z;
	//MeshEnergyを計算
	std::vector<double> dr = std::vector<double>(3, 0);
	dr[0] = this->max_x / spacemesh[0];
	dr[1] = this->max_y / spacemesh[1];
	dr[2] = this->max_z / spacemesh[2];
	double drpro = dr[0] * dr[1] * dr[2];
	
	#pragma omp parallel for
	for (int inum = 0; inum < this->MCParticles.size(); inum++) {
		auto i = this->MCParticles[inum];
		index_x = static_cast<int>(std::clamp(i.position[0] / dr[0], 0.0, static_cast<double>(spacemesh[0]) - 0.5));
		index_y = static_cast<int>(std::clamp(i.position[1] / dr[1], 0.0, static_cast<double>(spacemesh[1]) - 0.5));
		index_z = static_cast<int>(std::clamp(i.position[2] / dr[2], 0.0, static_cast<double>(spacemesh[2]) - 0.5));
		MeshEnergy[index_x][index_y][index_z] += this->energy_mcparticles / drpro;
	}
	#pragma omp barrier
	
	//dTは多分Tの差分
	//Tの変動が少ないと終わるはず
	double dT = 0;
	for(int i = 0 ; i < spacemesh[0];i++){
		for(int j = 0 ; j < spacemesh[1];j++){
			#pragma omp parallel for
			for(int k = 0 ; k < spacemesh[2];k++){
				double old_Temperature = Temperature[i][j][k];
				Temperature[i][j][k] = this->heat_cap.itpl_getter(MeshEnergy[i][j][k]);
				dT += fabs(old_Temperature - Temperature[i][j][k]);
			}
			#pragma omp barrier
		}
	}
	/* std::cout << dT << std::endl; */
	this->logger.info("dT is " + std::to_string(dT));
	/*if (emesh, 3) / 2)return true;
	else*/ return false;
}

void simulation::Particle_move(double dt) {
	#pragma omp parallel for
	//for (auto i = this->MCParticles.begin(); i < this->MCParticles.end(); i++) {
	for(int j = 0;j < static_cast<int>(MCParticles.size());j++){
		auto i = &(MCParticles[j]);
		(*i).Nextstep(dt);
		(*i).Boundary_Scatter_B(this->max_x, this->max_y, this->max_z);
		int index_x = (int)std::max(0,std::min((int)floor(((*i).position)[0] / this->max_x * this->spacemesh[0]),this->spacemesh[0] - 1));
		int index_y = (int)std::max(0,std::min((int)floor(((*i).position)[1] / this->max_y * this->spacemesh[1]),this->spacemesh[1] - 1));
		int index_z = (int)std::max(0,std::min((int)floor(((*i).position)[2] / this->max_z * this->spacemesh[2]),this->spacemesh[2] - 1));
		(*i).Scatter(Temperature[index_x][index_y][index_z], dt, std::min({ max_x,max_y,max_z }));
	}
	#pragma omp barrier
}
