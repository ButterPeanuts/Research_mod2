#include<iostream>
#include<vector>
#include<fstream>
#include<chrono>
#include"module/physconst.hpp"
#include"module/mcparticles.hpp"
#include"module/Simulation.hpp"
#include"module/Integral.h"
#include"module/massconst.hpp"

int main(){
	massconst::Ndiv = 48;
	//massconst::Si_dispersion_table_construct();
	//massconst::dispersion_table_output("Si_dispersion_48t", massconst::Si_dispersion);
	massconst::dispersion_table_input("Si_dispersion_48t", massconst::Si_dispersion);
	//massconst::Si_DOS_table_construct();
	//massconst::DOS_table_output("Si_DOS_TA_48t", "Si_DOS_LA_48t", massconst::Si_DOS_TA, massconst::Si_DOS_LA);
	massconst::DOS_table_input("Si_DOS_TA_48t", "Si_DOS_LA_48t", massconst::Si_DOS_TA, massconst::Si_DOS_LA);
	
	//massconst::Si_heatcap_table_construct();
	//massconst::heatcap_table_output("Si_heatcap_48t", massconst::Si_heatcap);
	massconst::heatcap_table_input("Si_heatcap_48t", massconst::Si_heatcap);
	
	Simulation mainSimulation(1000, 50.00, 50.00,50.00, {1,1,1});
	//mainSimulation.Temperature_construct();
	
	//最小緩和時間和を求める
	//シミュレーションステップはこれで決めるらしい
	double maxPu1 = massconst::Si_scatter_ALA *
		std::pow((*(massconst::Si_DOS_LA.end() - 1))[0], massconst::Si_scatter_chiLA) *
		std::pow(massconst::heatcaps_Tempmax, massconst::Si_scatter_xiLA) *
		std::exp(-1 * massconst::Si_scatter_BLA / massconst::heatcaps_Tempmax);
	double maxPu2 = massconst::Si_scatter_ATA *
		std::pow((*(massconst::Si_DOS_TA.end() - 1))[0], massconst::Si_scatter_chiTA) *
		std::pow(massconst::heatcaps_Tempmax, massconst::Si_scatter_xiTA) *
		std::exp(-1 * massconst::Si_scatter_BTA / massconst::heatcaps_Tempmax);
	double minTu = 1 / (std::min({ maxPu1,maxPu2 }));
	
	double maxomega = std::max((*(massconst::Si_DOS_LA.end() - 1))[0], (*(massconst::Si_DOS_TA.end() - 1))[0]);
	double minTd = 1 / (massconst::Si_scatter_C * pow(maxomega, 4));
	
	//double minTb = 8.05e-6 * 0.55 / 8480;
	
	double minT = minTu + minTd /*+ minTb*/;
	//今からやるのは半導体の端まで到達しないシミュ
	//じゃあPb -> 0でよくねってなった
	double time_sim = 0;
	std::cout << "Simulation is started" << std::endl;
	std::vector<std::vector<std::vector<double>>> Temperature = { {{300}} };
	mainSimulation.Temperature = Temperature;
	
	std::vector<double> output_thr = { 1e-9,1e-8,1e-7,1e-6,1.6e-6,2.5e-6,4e-6,6.3e-6,1e-5,1.6e-5,2.5e-5,4e-5,6.3e-5,1e-4};
	int output_size = output_thr.size();
	int outputed = 0;
	auto start = std::chrono::system_clock::now();
	for (;;) {
		mainSimulation.Particle_move(minT);
		time_sim += minT;
		//if (mainSimulation.Temperature_construct())break;
		std::cout << time_sim << std::endl;
		//if (time_sim > 2e-6)break;
		if (time_sim >= output_thr[0]) {
			mainSimulation.Particle_Disp_output("Displacement" + std::to_string(outputed));
			output_thr.erase(output_thr.begin());
			outputed++;
			if (outputed == output_size)break;
		}
		
	}
	auto end = std::chrono::system_clock::now();
	auto time = end - start;
	std::cout << std::chrono::duration_cast<std::chrono::seconds>(time).count();
	mainSimulation.Particle_Disp_output("Displacement");
}
