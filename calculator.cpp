#include<massconst.hpp>
#include<logger_obj.hpp>
#include<logger.hpp>
#include<band.hpp>
#include<band_obj.hpp>
#include<curve.hpp>
#include<scatconst.hpp>
#include<modeenum.hpp>
using namespace mc_sim;

int main(){
	logger_obj logobj = logger_obj("calculator", "log/calculator.log");
	
	logger& curveslogger(logobj.copy_samesink("curves"));
	curve dos_la = curve(curveslogger, "data/Si_DOS_LA_48t.txt");
	curve gv_la = curve(curveslogger, "data/Si_gvelocity_LA.txt");
	scatconst scatconst_la = scatconst(curveslogger, "data/Si_scatconst_LA.txt");
	auto band_la = std::make_shared<band>(band_obj(logobj.copy_samesink("band"), dos_la, gv_la, wave_direction::longitudinal, wave_mode::acoustic, scatconst_la));
	
	curve dos_ta = curve(curveslogger, "data/Si_DOS_TA_48t.txt");
	curve gv_ta = curve(curveslogger, "data/Si_gvelocity_TA.txt");
	scatconst scatconst_ta = scatconst(curveslogger, "data/Si_scatconst_TA.txt");
	auto band_ta = std::make_shared<band>(band_obj(logobj.copy_samesink("band"), dos_la, gv_la, wave_direction::transverse, wave_mode::acoustic, scatconst_la));
	
	std::vector<std::shared_ptr<band>> bandlist;
	bandlist.push_back(band_la);
	bandlist.push_back(band_ta);
	bandlist.push_back(band_ta);
	curve heatcap(massconst::heatcap_curve_construct(bandlist, logobj));
	heatcap.file_output("data/Si_heatcap_48neo.txt");
}
