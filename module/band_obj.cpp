#include "modeenum.hpp"
#include "band_obj.hpp"
#include <logger.hpp>

band_obj::band_obj(mc_sim::logger& newlogger, curve dos, curve gvelocity, curve domcp, wave_direction direction, wave_mode mode) : dos(dos), gvelocity(gvelocity), domcp(domcp), logger(newlogger){
	this->directions = direction;
	this->mode = mode;
	this->domcp_max = this->domcp.max();
	this->dos_max = this->dos.max();
}

double band_obj::dos_getter(double omega){
	return this->dos.itpl_getter(omega);
}

double band_obj::dos_omega_min_getter(){
	return this->dos.left_edge();
}

double band_obj::dos_omega_max_getter(){
	return this->dos.right_edge();
}

double band_obj::dos_max_getter(){
	return this->dos_max;
}

double band_obj::gvelocity_getter(double omega){
	return this->gvelocity.itpl_getter(omega);
}

double band_obj::domcp_getter(double omega){
	return this->domcp.itpl_getter(omega);
}

double band_obj::domcp_omega_min_getter(){
	return this->domcp.left_edge();
}

double band_obj::domcp_omega_max_getter(){
	return this->domcp.right_edge();
}

double band_obj::domcp_max_getter(){
	return this->domcp_max;
}

wave_direction band_obj::directions_getter(){
	return this->directions;
}

wave_mode band_obj::mode_getter(){
	return this->mode;
}
