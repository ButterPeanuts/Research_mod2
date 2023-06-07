#include<functional>
#include<string>
#include<fstream>
#include<algorithm>
#include<numeric>
#include<omp.h>
#include<future>

#include"simulation.hpp"
#include"physconst.hpp"
#include"massconst.hpp"
#include"Integral.h"
#include"mcparticles.hpp"
simulation::simulation(int numof_mcp, std::vector<double> max_r, std::vector<int> spacemesh, double tempof_device, curve internal_energy, curve heat_cap, mc_sim::logger& logger, std::vector<band>& band_inj) :  logger(logger), banddata(band_inj), internal_energy(internal_energy), heat_cap(heat_cap) {
	//デバイス大きさ
	this->max_r = max_r;
	//体積
	this->volume = std::accumulate(this->max_r.begin(), this->max_r.end(), 1.0, std::multiplies<>());
	//メッシュ細かさ
	this->spacemesh = spacemesh;
	
	//メッシュ粗さ[m]
	this->dr = std::vector<double>();
	for (int i = 0; i < std::min(max_r.size(), spacemesh.size()); i++){
		this->dr.push_back(this->max_r[i] / static_cast<double>(this->spacemesh[i]));
	}
	
	
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
	
	this->logger.debug("We will construct mcparticles.");
	
	//MC粒子の初期化開始
	this->MCParticles = std::vector<mc_particles::MCParticles>();
	#pragma omp parallel for
	for (int i = 0; i < numof_mcp; i++) {
		mc_sim::logger& newlogger = this->logger.copy_samesink("mcp" + std::to_string(i));
		auto mcp_part = mc_particles::MCParticles(newlogger, static_cast<double>(tempof_device), this->banddata);
		//本来criticalはあまりスピード的に優越しない
		//ただ今回はmcparticlesのコンストラクトのほうが支配的な時間をかけるという仮定の元やってみる
		#pragma omp critical(mcparticle_pushback)
		{
			this->MCParticles.push_back(mcp_part);
			this->logger.debug("mcp" + std::to_string(i) + " is constructed.");
		}
	}
	#pragma omp barrier
	this->logger.debug("We have just constructed mcparticles.");
	
	this->Temperature = std::vector<std::vector<std::vector<double>>>(spacemesh[0], std::vector < std::vector<double>>(spacemesh[1], std::vector<double>(spacemesh[2], 0.0)));
	
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
	double drpro = std::accumulate(this->dr.begin(), this->dr.end(), 1.0, std::multiplies<>());
	
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
	this->logger.info("Temperature table is updated. dT is " + std::to_string(dT));
	/*if (emesh, 3) / 2)return true;
	else*/ return false;
}

void simulation::Particle_move(double dt) {
	#pragma omp parallel for
	for(int j = 0;j < static_cast<int>(MCParticles.size());j++){
		mc_particles::MCParticles& i = MCParticles[j];
		i.nextstep(dt);
		i.boundaryscatter_b(this->max_r[0], this->max_r[1], this->max_r[2]);
		std::vector<int> index = this->square(i.position);
		i.scatter(Temperature[index[0]][index[1]][index[2]], dt, *std::min_element(this->max_r.begin(), this->max_r.end()));
	}
	#pragma omp barrier
}

std::vector<int> simulation::square(std::vector<double> position){
	std::vector<int> res = std::vector<int>();
	for (int i = 0; i < this->dr.size(); i++){
		res.push_back(static_cast<int>(std::clamp(position[i] / this->dr[i], 0.0, static_cast<double>(spacemesh[i]) - 0.5)));
	}
	return res;
}
