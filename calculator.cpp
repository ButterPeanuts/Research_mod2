#include<massconst.hpp>
#include<logger_obj.hpp>
#include<logger.hpp>
#include<band.hpp>
#include<band_obj.hpp>
#include<curve.hpp>
#include<scatconst.hpp>
#include<modeenum.hpp>
#include<brillouin_zone_funcobj.hpp>
#include<brillouin_zone.hpp>

#include<iostream>
using namespace mc_sim;

void heatcap_calculator(logger_obj& logobj, std::vector<std::shared_ptr<band>>& bandlist);
void internalenergy_calculator(logger_obj& logobj, std::vector<std::shared_ptr<band>>& bandlist);
void dos_domcpmax_calculator(logger_obj& logobj);

int main(){
	logger_obj logobj = logger_obj("calculator", "log/calculator.log");
	
	//dos, domcpmaxを計算する
	dos_domcpmax_calculator(logobj);
	
	std::shared_ptr<logger> lados_logger(logobj.copy_samesink("LA_dos_logger"));
	std::shared_ptr<logger> ladomcpmax_logger(logobj.copy_samesink("LA_domcpmax_logger"));
	std::shared_ptr<logger> lagv_logger(logobj.copy_samesink("LA_gv_logger"));
	std::shared_ptr<logger> lascat_logger(logobj.copy_samesink("LA_scat_logger"));
	curve dos_la = curve(lados_logger, "data/Si_DOS_LA_100.txt");
	curve domcpmax_la = curve(ladomcpmax_logger, "data/Si_DOMCPMAX_LA_100.txt");
	curve gv_la = curve(lagv_logger, "data/Si_gvelocity_LA.txt");
	scatconst scatconst_la = scatconst(lascat_logger, "data/Si_scatconst_LA.txt");
	std::shared_ptr<logger> logger_la = logobj.copy_samesink("band_LA");
	auto band_la = std::make_shared<band_obj>(band_obj(logger_la, dos_la, gv_la, domcpmax_la, wave_direction::longitudinal, wave_mode::acoustic, scatconst_la));
	
	std::shared_ptr<logger> tados_logger(logobj.copy_samesink("TA_dos_logger"));
	std::shared_ptr<logger> tadomcpmax_logger(logobj.copy_samesink("TA_domcpmax_logger"));
	std::shared_ptr<logger> tagv_logger(logobj.copy_samesink("TA_gv_logger"));
	std::shared_ptr<logger> tascat_logger(logobj.copy_samesink("TA_scat_logger"));
	curve dos_ta = curve(tados_logger, "data/Si_DOS_TA_100.txt");
	curve domcpmax_ta = curve(tadomcpmax_logger, "data/Si_DOMCPMAX_TA_100.txt");
	curve gv_ta = curve(tagv_logger, "data/Si_gvelocity_TA.txt");
	scatconst scatconst_ta = scatconst(tascat_logger, "data/Si_scatconst_TA.txt");
	std::shared_ptr<logger> logger_ta = logobj.copy_samesink("band_TA");
	auto band_ta = std::make_shared<band_obj>(band_obj(logger_ta, dos_ta, gv_ta, domcpmax_ta, wave_direction::transverse, wave_mode::acoustic, scatconst_ta));
	
	std::vector<std::shared_ptr<band>> bandlist;
	bandlist.push_back(band_la);
	bandlist.push_back(band_ta);
	bandlist.push_back(band_ta);
	heatcap_calculator(logobj, bandlist);
	internalenergy_calculator(logobj, bandlist);
}

void heatcap_calculator(logger_obj& logobj, std::vector<std::shared_ptr<band>>& bandlist){
	std::shared_ptr<logger> curveslogger(logobj.copy_samesink("heatcap_curve"));
	curve heatcap(massconst::heatcap_curve_construct(bandlist, curveslogger));
	heatcap.file_output("data/Si_heatcap_100.txt");
}

void internalenergy_calculator(logger_obj& logobj, std::vector<std::shared_ptr<band>>& bandlist){
	std::shared_ptr<logger> curves_logger(logobj.copy_samesink("tempcurve"));
	auto res = (massconst::internal_energy_construct(bandlist, curves_logger));
	res.first.file_output("data/Si_internal_100.txt");
	res.second.file_output("data/Si_temperature_100.txt");
}

void dos_domcpmax_calculator(logger_obj& logobj){
	logobj.info("We start calculating dos and domcpmax.");
	std::shared_ptr<logger> curveslogger(logobj.copy_samesink("curves"));
	int ndiv = 100;
	
	brillouin_zone_funcobj bz_la(ndiv, massconst::si_angfreq_100_la, wave_direction::longitudinal, wave_mode::acoustic);
	brillouin_zone_funcobj bz_ta(ndiv, massconst::si_angfreq_100_ta, wave_direction::transverse, wave_mode::acoustic);
	
	auto lares = massconst::dos_domcpmax_tetrahedron(bz_la, curveslogger);
	auto tares = massconst::dos_domcpmax_tetrahedron(bz_ta, curveslogger);
	
	lares.first.file_output("data/Si_DOS_LA_100.txt");
	tares.first.file_output("data/Si_DOS_TA_100.txt");
	lares.second.file_output("data/Si_DOMCPMAX_LA_100.txt");
	tares.second.file_output("data/Si_DOMCPMAX_TA_100.txt");
	logobj.info("Dos and Domcpmax is calculated!");
}
