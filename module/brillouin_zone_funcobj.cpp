#include <brillouin_zone_funcobj.hpp>
#include <physconst.hpp>

mc_sim::brillouin_zone_funcobj::brillouin_zone_funcobj(int const & ndiv, std::function<double(std::tuple<double, double, double> const &)> const & disparsion_calc, wave_direction const & direction, wave_mode const & mode) : ndiv(ndiv), direction(direction), mode(mode), disparsion_calc(disparsion_calc){}


int mc_sim::brillouin_zone_funcobj::ndiv_getter(){
	return ndiv;
}

double mc_sim::brillouin_zone_funcobj::angfreq_index(std::tuple<int, int, int> const & kcoor_std){
	return this->disparsion_calc(physconst::indextostd(kcoor_std, this->ndiv));
}

wave_direction mc_sim::brillouin_zone_funcobj::directions_getter(){
	return this->direction;
}

wave_mode mc_sim::brillouin_zone_funcobj::mode_getter(){
	return this->mode;
}
