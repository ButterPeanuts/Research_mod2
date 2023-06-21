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

int main(){
	std::cout << "That is sign" << std::endl;
	logger_obj logobj = logger_obj("calculator", "log/calculator.log");
	std::shared_ptr<logger> curveslogger(logobj.copy_samesink("curves"));
	
	logobj.info("We start calculating dos.");
	std::cout << "That is sign" << std::endl;
	int ndiv = 48;
	brillouin_zone_funcobj bz_la(ndiv, massconst::si_angfreq_100_la, wave_direction::longitudinal, wave_mode::acoustic);
	brillouin_zone_funcobj bz_ta(ndiv, massconst::si_angfreq_100_ta, wave_direction::transverse, wave_mode::acoustic);
	curve dos_la = massconst::doscurve_tetrahedron(bz_la, curveslogger);
	curve dos_ta = massconst::doscurve_tetrahedron(bz_ta, curveslogger);
	dos_la.file_output("data/Si_DOS_LA_48neo.txt");
	dos_ta.file_output("data/Si_DOS_TA_48neo.txt");
	std::cout << dos_la.itpl_getter(1e+10) << std::endl;
	logobj.info("Dos is calculated!");
	std::cout << "That is sign" << std::endl;
	
	//curve dos_la = curve(curveslogger, "data/Si_DOS_LA_48t.txt");
	curve gv_la = curve(curveslogger, "data/Si_gvelocity_LA.txt");
	scatconst scatconst_la = scatconst(curveslogger, "data/Si_scatconst_LA.txt");
	std::shared_ptr<logger> logger_la = logobj.copy_samesink("band_LA");
	auto band_la = std::make_shared<band_obj>(band_obj(logger_la, dos_la, gv_la, wave_direction::longitudinal, wave_mode::acoustic, scatconst_la));
	
	//curve dos_ta = curve(curveslogger, "data/Si_DOS_TA_48t.txt");
	curve gv_ta = curve(curveslogger, "data/Si_gvelocity_TA.txt");
	scatconst scatconst_ta = scatconst(curveslogger, "data/Si_scatconst_TA.txt");
	std::shared_ptr<logger> logger_ta = logobj.copy_samesink("band_TA");
	auto band_ta = std::make_shared<band_obj>(band_obj(logger_ta, dos_la, gv_la, wave_direction::transverse, wave_mode::acoustic, scatconst_la));
	
	std::vector<std::shared_ptr<band>> bandlist;
	bandlist.push_back(band_la);
	bandlist.push_back(band_ta);
	bandlist.push_back(band_ta);
	curve heatcap(massconst::heatcap_curve_construct(bandlist, curveslogger));
	heatcap.file_output("data/Si_heatcap_48neo.txt");
}
