#include<functional>
#include<string>
#include<fstream>
#include<algorithm>
#include<numeric>
#include<omp.h>
#include<future>
#include<iostream>
#include<mutex>
#include<array>

#include"simulation.hpp"
#include"mcparticle.hpp"
simulation::simulation(int numof_mcp, std::vector<double>& max_r, std::vector<int>& spacemesh, double tempof_device, curve& internal_energy, curve& heat_cap, std::shared_ptr<mc_sim::logger>& logger, std::vector<std::shared_ptr<band>>& band_inj) :  logger(logger), banddata(band_inj), internal_energy(internal_energy), heat_cap(heat_cap) {
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
	
	this->logger->debug("We will construct mcparticles.");
	
	//MC粒子の初期化開始
	this->mc_particles = std::vector<mc_sim::mc_particle>();
	#pragma omp parallel for
	for (int i = 0; i < numof_mcp; i++) {
		std::shared_ptr<mc_sim::logger> newlogger = this->logger->copy_samesink("mcp" + std::to_string(i));
		auto mcp_part = mc_sim::mc_particle(newlogger, static_cast<double>(tempof_device), this->banddata);
		//本来criticalはあまりスピード的に優越しない
		//ただ今回はmcparticlesのコンストラクトのほうが支配的な時間をかけるという仮定の元やってみる
		#pragma omp critical(mcparticle_pushback)
		{
			this->mc_particles.push_back(mcp_part);
			this->logger->debug("mcp" + std::to_string(i) + " is constructed.");
		}
	}
	#pragma omp barrier
	this->logger->debug("We have just constructed mcparticles.");
	
	this->Temperature = std::vector<std::vector<std::vector<double>>>(spacemesh[0], std::vector < std::vector<double>>(spacemesh[1], std::vector<double>(spacemesh[2], 0.0)));
	this->mcp_freqdist = std::vector<int>(std::accumulate(this->spacemesh.begin(), this->spacemesh.end(), 1, std::multiplies<>()), 0);
	std::vector<std::mutex> temp_mutex(std::accumulate(this->spacemesh.begin(), this->spacemesh.end(), 1, std::multiplies<>()));
	this->mcp_freqdist_mutex.swap(temp_mutex);
	this->freqdist_construct();
	this->temperature_construct();
	
	this->logger->info("Simulation is constructed");
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
	for (auto i = this->mc_particles.begin(); i < this->mc_particles.end(); i++) {
		files << i->position[0] << ", "<<i->position[1] <<", "<<i->position[2] <<std::endl;
	}
	files.close();
}

void simulation::freqdist_construct(){
	std::vector<std::future<void>> futures;
	for (auto& i: this->mc_particles){
		futures.push_back(std::async(std::launch::async, [&i, this](){
			auto coor = this->square(i.position);
			auto index = this->tempcoor_to_fdlinear({coor[0], coor[1], coor[2]});
			{
				std::lock_guard<std::mutex> lock(this->mcp_freqdist_mutex[index]);
				this->mcp_freqdist[index]++;
			}
		}));
	}
	for (auto& i: futures)i.get();
}

bool simulation::temperature_construct() {
	this->logger->debug("We will construct MeshEnergys.");
	
	//dTは多分Tの差分
	//Tの変動が少ないと終わるはず
	double dT = 0;
	//meshの体積の逆数
	const double drpro_inv = 1.0 / std::accumulate(this->dr.begin(), this->dr.end(), 1.0, std::multiplies<>());
	
	std::vector<std::future<double>> futures;
	for(int i = 0 ; i < spacemesh[0];i++){
		futures.push_back(std::async(std::launch::async, [i, this, &drpro_inv](){
			double dtpart = 0;
			for(int j = 0 ; j < spacemesh[1];j++){
				for(int k = 0 ; k < spacemesh[2];k++){
					auto index = this->tempcoor_to_fdlinear({i, j, k});
					double old_Temperature = Temperature[i][j][k];
					Temperature[i][j][k] = this->heat_cap.itpl_getter(static_cast<double>(this->mcp_freqdist[index]) * this->energy_mcparticles * drpro_inv);
					dtpart += fabs(old_Temperature - Temperature[i][j][k]);
				}
			}
			return dtpart;
		}));
	}
	for (auto& i: futures)dT += i.get();
	/* std::cout << dT << std::endl; */
	this->logger->info("Temperature table is updated. dT is " + std::to_string(dT));
	/*if (emesh, 3) / 2)return true;
	else*/ return false;
}

void simulation::Particle_move(double dt) {
	/* #pragma omp parallel for */
	/* for(int j = 0;j < static_cast<int>(MCParticles.size());j++){ */
	/* 	mc_particles::MCParticles& i = MCParticles[j]; */
	/* 	i.nextstep(dt); */
	/* 	i.boundaryscatter_b(this->max_r[0], this->max_r[1], this->max_r[2]); */
	/* 	std::vector<int> index = this->square(i.position); */
	/* 	i.scatter(Temperature[index[0]][index[1]][index[2]], dt, *std::min_element(this->max_r.begin(), this->max_r.end())); */
	/* } */
	/* #pragma omp barrier */
	std::vector<std::future<void>> futures;
	for (auto i: this->mc_particles){
		futures.push_back(std::async(std::launch::async, [dt, this, &i](){
			auto beforecoor = this->square(i.position);
			auto beforeindex = this->tempcoor_to_fdlinear({beforecoor[0], beforecoor[1], beforecoor[2]});
			
			i.nextstep(dt);
			i.boundaryscatter_b(this->max_r[0], this->max_r[1], this->max_r[2]);
			std::vector<int> index = this->square(i.position);
			i.scatter(Temperature[index[0]][index[1]][index[2]], dt, *std::min_element(this->max_r.begin(), this->max_r.end()));
			
			auto aftercoor = this->square(i.position);
			auto afterindex = this->tempcoor_to_fdlinear({aftercoor[0], aftercoor[1], aftercoor[2]});
			{
				std::lock_guard<std::mutex> lock(this->mcp_freqdist_mutex[beforeindex]);
				this->mcp_freqdist[beforeindex]--;
			}
			{
				std::lock_guard<std::mutex> lock(this->mcp_freqdist_mutex[afterindex]);
				this->mcp_freqdist[afterindex]++;
			}
		}));
	}
	for (auto& i: futures){
		i.get();
	}
	temperature_construct();
}

std::vector<int> simulation::square(std::vector<double> position){
	std::vector<int> res = std::vector<int>();
	for (int i = 0; i < this->dr.size(); i++){
		res.push_back(static_cast<int>(std::clamp(position[i] / this->dr[i], 0.0, static_cast<double>(spacemesh[i]) - 0.5)));
	}
	return res;
}

int simulation::tempcoor_to_fdlinear(const std::array<int, 3>& meshcoor){
	return meshcoor[2] + meshcoor[1] * this->spacemesh[2] + meshcoor[0] * this->spacemesh[2] * this->spacemesh[1];
}
