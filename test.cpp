#include<iostream>
#include<vector>
#include<fstream>
#include<chrono>
#include<logger_obj.hpp>
#include<logger.hpp>
#include<band.hpp>
#include<band_obj.hpp>
#include<curve.hpp>
#include<scatconst.hpp>
#include<modeenum.hpp>
#include<brillouin_zone_funcobj.hpp>
#include<brillouin_zone.hpp>
#include<simulation.hpp>
using namespace mc_sim;

int main(){
	logger_obj logobj = logger_obj("main", "log/simulator.log");
	
	
	//laバンドの情報
	std::shared_ptr<logger> lados_logger(logobj.copy_samesink("LA_dos"));
	std::shared_ptr<logger> ladomcpmax_logger(logobj.copy_samesink("LA_domcpmax"));
	std::shared_ptr<logger> lagv_logger(logobj.copy_samesink("LA_gv"));
	std::shared_ptr<logger> lascat_logger(logobj.copy_samesink("LA_scat"));
	curve dos_la = curve(lados_logger, "data/Si_DOS_LA_100.txt");
	curve domcpmax_la = curve(ladomcpmax_logger, "data/Si_DOMCPMAX_LA_100.txt");
	curve gv_la = curve(lagv_logger, "data/Si_gvelocity_LA.txt");
	scatconst scatconst_la = scatconst(lascat_logger, "data/Si_scatconst_LA.txt");
	
	std::shared_ptr<logger> logger_la = logobj.copy_samesink("band_LA");
	auto band_la = std::make_shared<band_obj>(band_obj(logger_la, dos_la, gv_la, domcpmax_la, wave_direction::longitudinal, wave_mode::acoustic, scatconst_la));
	
	
	//taバンドの情報
	std::shared_ptr<logger> tados_logger(logobj.copy_samesink("TA_dos"));
	std::shared_ptr<logger> tadomcpmax_logger(logobj.copy_samesink("TA_domcpmax"));
	std::shared_ptr<logger> tagv_logger(logobj.copy_samesink("TA_gv"));
	std::shared_ptr<logger> tascat_logger(logobj.copy_samesink("TA_scat"));
	curve dos_ta = curve(tados_logger, "data/Si_DOS_TA_100.txt");
	curve domcpmax_ta = curve(tadomcpmax_logger, "data/Si_DOMCPMAX_TA_100.txt");
	curve gv_ta = curve(tagv_logger, "data/Si_gvelocity_TA.txt");
	scatconst scatconst_ta = scatconst(tascat_logger, "data/Si_scatconst_TA.txt");
	
	std::shared_ptr<logger> logger_ta = logobj.copy_samesink("band_TA");
	auto band_ta = std::make_shared<band_obj>(band_obj(logger_ta, dos_ta, gv_ta, domcpmax_ta, wave_direction::transverse, wave_mode::acoustic, scatconst_ta));
	
	
	//全バンドの情報
	std::vector<std::shared_ptr<band>> bandlist;
	bandlist.push_back(band_la);
	bandlist.push_back(band_ta);
	bandlist.push_back(band_ta);
	
	
	//内部エネルギー関数とその逆関数
	std::shared_ptr<logger> internalenergy_logger(logobj.copy_samesink("internalenergy"));
	std::shared_ptr<logger> tempcurve_logger(logobj.copy_samesink("tempcurve"));
	curve internalenergy = curve(internalenergy_logger, "data/Si_internal_100.txt");
	curve tempcurve = curve(tempcurve_logger, "data/Si_temperature_100.txt");
	
	
	//その他パラメータ
	std::shared_ptr<logger> simulator_logger(logobj.copy_samesink("simulator"));
	std::vector<double> max_r = {1, 1, 1};
	std::vector<int> spacemesh = {20, 20, 20};
	double tempof_device = 300;
	
	
	//シミュレーター
	simulation mainsimulation(1000, max_r, spacemesh, tempof_device, internalenergy, tempcurve, simulator_logger, bandlist);
	/* mainsimulation.Temperature_construct(); */
	
	
	//最小緩和時間和を求める
	double mintau = std::numeric_limits<double>::infinity();
	for (auto& i: bandlist){
		double part = i->mintau(internalenergy.right_edge(), *std::min_element(max_r.begin(), max_r.end()));
		if (part < mintau){
			mintau = part;
		}
	}
	
	
	double time_sim = 0;
	std::cout << "Simulation is started" << std::endl;
	
	//std::vector<double> output_thr = { 1e-9,1e-8,1e-7,1e-6,1.6e-6,2.5e-6,4e-6,6.3e-6,1e-5,1.6e-5,2.5e-5,4e-5,6.3e-5,1e-4};
	//auto start = std::chrono::system_clock::now();
	double time_end = 1e-13;
	for (;;) {
		mainsimulation.Particle_move(mintau);
		//Particle_moveに吸収された動作
		/* mainsimulation.Temperature_construct(); */
		time_sim += mintau;
		std::cout << time_sim << std::endl;
		if (time_sim > time_end)break;
	}
	mainsimulation.Particle_Disp_output("Displacement");
}
