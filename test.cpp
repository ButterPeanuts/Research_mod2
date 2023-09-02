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
	scatconst scatconst_la = scatconst(lascat_logger, "data/Si_scatconst_LA_inf.txt");
	
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
	scatconst scatconst_ta = scatconst(tascat_logger, "data/Si_scatconst_TA_inf.txt");
	
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
	std::vector<double> max_r = {7.16e-2, 7.16e-2, 7.16e-2};
	std::vector<int> spacemesh = {1, 1, 1};
	double tempof_device = 300;
	
	
	//シミュレーター
	simulation mainsimulation(1000, max_r, spacemesh, tempof_device, internalenergy, tempcurve, simulator_logger, bandlist);
	/* mainsimulation.Temperature_construct(); */
	mainsimulation.particle_posinit([max_r](){
		std::vector a{max_r[0] / 2, max_r[1] / 2, max_r[2] / 2};
		return a;
	});
	
	
	//最小緩和時間和を求める
	double mintau = std::numeric_limits<double>::infinity();
	for (auto& i: bandlist){
		/* double part = i->mintau(tempcurve.max(), *std::min_element(max_r.begin(), max_r.end())); */
		double part = i->mintau(tempof_device, *std::min_element(max_r.begin(), max_r.end()));
		if (part < mintau){
			mintau = part;
		}
	}
	
	
	double time_sim = 0;
	
	int step = 7 * 40;
	double start = -13;
	double end = -6;
	std::vector<double> output_thr = {};
	for (int i = 0; i <= step; i++){
		output_thr.push_back(std::pow(10, std::lerp(start, end, static_cast<double>(i) / static_cast<double>(step))));
	}
	mainsimulation.Particle_Disp_output("disp/Dispacementinf");
	logobj.info("Simulation is started");
	logobj.info("mintau is " + std::to_string(std::log10(mintau)));
	logobj.debug("Debugtest");
	/* auto start = std::chrono::system_clock::now(); */
	try{
	for (;;) {
		mainsimulation.Particle_move(mintau);
		time_sim += mintau;
		if (time_sim > output_thr[0]){
			mainsimulation.Particle_Disp_output("disp/Dispacement" + std::to_string(std::log10(time_sim)));
			logobj.info("It is " + std::to_string(std::log10(time_sim)));
			output_thr.erase(output_thr.begin());
			if (output_thr.empty()) break;
		}
	}
	}
	catch(...){
		std::cout << "Wahts happen?" << std::endl;
	}
}
